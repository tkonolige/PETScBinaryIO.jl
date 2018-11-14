module PETScBinaryIO

export writePETSc, readPETSc

using SparseArrays

classids = Dict("Vec"=>1211214, "Mat"=>1211216)
ids_to_class = Dict(zip(values(classids), keys(classids)))

# PETSc IO for binary matrix IO. Format documentation here:
# http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatLoad.html
#  int    MAT_FILE_CLASSID
#  int    number of rows
#  int    number of columns
#  int    total number of nonzeros
#  int    *number nonzeros in each row
#  int    *column indices of all nonzeros (starting index is zero)
#  PetscScalar *values of all nonzeros

"""
### writePETSc(filename, mat :: SparseMatrixCSC; int_type = Int32, scalar_type = Float64)

Write a sparse matrix to `filename` in a format PETSc can understand. `int_type`
and `scalar_type` need to be set to integer and scalar types that PETSc was
configured with.
"""
function writePETSc(filename, mat :: SparseMatrixCSC; int_type = Int32, scalar_type = Float64)
    m = SparseMatrixCSC(transpose(mat)) # transpose for fast row operations
    open(filename, "w") do io
        rows, cols = size(mat)
        write(io, hton(int_type(classids["Mat"]))) # MAT_FILE_CLASSID
        write(io, hton(int_type(rows))) # number of rows
        write(io, hton(int_type(cols))) # number of columns
        write(io, hton(int_type(nnz(mat)))) # number of nonzeros

        # write row lengths
        for i = 1:rows
            write(io, hton(int_type(length(nzrange(m, i)))))
        end

        # write column indices
        cols = rowvals(m)
        for i = 1:rows
            for j in nzrange(m, i)
                write(io, hton(int_type(cols[j] - 1))) # PETSc uses 0-indexed arrays
            end
        end

        # write nonzero values
        vals = nonzeros(m)
        for i = 1:rows
            for j in nzrange(m, i)
                write(io, hton(scalar_type(vals[j])))
            end
        end
    end
end

"""
### writePETSc(filename, mat :: SparseMatrixCSC; int_type = Int32, scalar_type = Float64)

Write a vector to `filename` in a format PETSc can understand. `int_type`
and `scalar_type` need to be set to integer and scalar types that PETSc was
configured with.
"""
function writePETSc(filename, vec :: Vector; int_type = Int32, scalar_type = Float64)
    open(filename, "w") do io
        write(io, hton(int_type(classids["Vec"]))) # MAT_FILE_CLASSID
        write(io, hton(int_type(length(vec)))) # number of rows
        write(io, hton.(scalar_type.(vec)))
    end
end

function read_prefix_vec(io, int_type, scalar_type)
    len = ntoh(read(io, int_type))
    read_vec(io, scalar_type, len)
end

function read_vec(io, ty, sz)
    ntoh.(read!(io, Array{ty}(undef, sz)))
end

function read_mat(io, int_type, scalar_type)
    rows = ntoh(read(io, int_type))
    cols = ntoh(read(io, int_type))
    nnz = ntoh(read(io, int_type))

    row_ptr = Array{int_type}(undef, rows+1)
    row_ptr[1] = 1

    # read row lengths
    row_ptr[2:end] = read_vec(io, int_type, rows)
    cumsum!(row_ptr, row_ptr)

    # write column indices
    colvals = read_vec(io, int_type, nnz) .+ int_type(1)

    # write nonzero values
    vals = read_vec(io, scalar_type, nnz)

    mat = SparseMatrixCSC(cols, rows, row_ptr, colvals, vals)
    SparseMatrixCSC(transpose(mat))
end

function read_single(io, int_type, scalar_type)
    class_id = ntoh(read(io, int_type))
    if !in(class_id, keys(ids_to_class))
        throw("Invalid PETSc binary file $class_id")
    end
    if ids_to_class[class_id] == "Vec"
        read_prefix_vec(io, int_type, scalar_type)
    elseif ids_to_class[class_id] == "Mat"
        read_mat(io, int_type, scalar_type)
    else
        error("Invalid class id $class_id")
    end
end

"""
### readPETSc(filename; int_type = Int32, scalar_type = Float64) :: Union{SparseMatrixCSC, Vector}

Read a sparse matrix in PETSc's binary format from `filename`. `int_type` and
`scalar_type` must be set to the integer and scalar types that PETSc was
configured with.
"""
function readPETSc(filename; int_type = Int32, scalar_type = Float64) :: Union{SparseMatrixCSC, Vector{SparseMatrixCSC}, Vector, Vector{Vector}}
    open(filename) do io
        items = []
        while !eof(io)
            push!(items, read_single(io, int_type, scalar_type))
        end
        items
    end
end

end

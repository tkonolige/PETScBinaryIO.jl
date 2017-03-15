module PETScBinaryIO

export writePETSc, readPETSc

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

Write a sparse matrix to `filename` in a format PETSc can understand. `int_typ`
and `scalar_type` need to be set to integer and scalar types that PETSc was
configured with.
"""
function writePETSc(filename, mat :: SparseMatrixCSC; int_type = Int32, scalar_type = Float64)
    m = transpose(mat) # transpose for fast row operations
    open(filename, "w") do io
        rows, cols = size(mat)
        write(io, hton(int_type(1211216))) # MAT_FILE_CLASSID
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
### readPETSc(filename; int_type = Int32, scalar_type = Float64) :: SparseMatrixCSC

Read a sparse matrix in PETSc's binary format from `filename`. `int_type` and
`scalar_type` must be set to the integer and scalar types that PETSc was
configured with.
"""
function readPETSc(filename; int_type = Int32, scalar_type = Float64) :: SparseMatrixCSC
    open(filename) do io
        class_id = ntoh(read(io, int_type))
        if class_id != 1211216
            throw("Invalid PETSc binary file")
        end
        rows = ntoh(read(io, int_type))
        cols = ntoh(read(io, int_type))
        nnz = ntoh(read(io, int_type))

        row_ptr = Array{int_type}(rows+1)
        row_ptr[1] = 1

        # read row lengths
        row_ptr[2:end] = map(ntoh, read(io, int_type, rows))
        cumsum!(row_ptr, row_ptr)

        # write column indices
        colvals = map(ntoh, read(io, int_type, nnz)) .+ 1

        # write nonzero values
        vals = map(ntoh, read(io, scalar_type, nnz))

        mat = SparseMatrixCSC(cols, rows, row_ptr, colvals, vals)
        transpose(mat)
    end
end

end

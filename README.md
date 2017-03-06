# PETScBinaryIO.jl

A Julia package for reading and writing sparse matrices in a format PETSc understands.

## Exported Functions

```julia
writePETSc(filename, mat :: SparseMatrixCSC; int_type = Int32, scalar_type = Float64)
```

Write a sparse matrix to `filename` in a format PETSc can understand. `int_typ`
and `scalar_type` need to be set to integer and scalar types that PETSc was
configured with.

```julia
readPETSc(filename; int_type = Int32, scalar_type = Float64) :: SparseMatrixCSC
```

Read a sparse matrix in PETSc's binary format from `filename`. `int_type` and
`scalar_type` must be set to the integer and scalar types that PETSc was
configured with.

using PETScBinaryIO
using Test
using SparseArrays

# write your own tests here
@testset "Idempotent IO" begin
    @test begin
        try
            A = sprand(1000,1000,0.1)
            writePETSc("tmp.petsc", A)
            AA = readPETSc("tmp.petsc")
            AA == A
        finally
            rm("tmp.petsc")
        end
    end

    @test begin
        try
            A = rand(1000)
            writePETSc("vec.petsc", A)
            AA = readPETSc("vec.petsc")
            AA == A
        finally
            rm("vec.petsc")
        end
    end
end

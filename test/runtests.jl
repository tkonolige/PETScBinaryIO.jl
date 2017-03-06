using PETScBinaryIO
using Base.Test

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
end

module QuantumTest

using QuantumBayesian
using Test
using LinearAlgebra
using SparseArrays

###
# Fixtures for tests
###

# Capture STDOUT
rd, wr = redirect_stdout()

# Define standard quantum spaces
q = qubit();
f = QFactor(3, "Foo");
o = osc(5);
s = ⊗(q,q,f)
t = (q ⊗ (o ⊗ f)) ⊗ q

include("testQubit.jl")

include("testOscillator.jl")

include("testOperations.jl")

include("testQuantum.jl")

include("testEvolution.jl")

end

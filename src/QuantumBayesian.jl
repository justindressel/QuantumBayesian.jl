__precompile__()

"""
    module QuantumBayesian

A Julia package for efficient simulation of quantum state evolution.
"""
module QuantumBayesian

#################################
# Package dependences and imports
#

# NOTE: there is no simple way to convert v0.5 to v0.6 of julia
#       with the Compat module, since the syntax for inner
#       constructors has fundamentally changed. The present code
#       works only with v0.6+

#using Distributions  # Minimize depedencies
using Interpolations
import Base.product

# Imported solely for method overloading purposes
# import Base.call
import Base.show
import Base.showarray
import Base.length
import Base.size
import Base.getindex
import Base.setindex!
import Base.ndims
import Base.indices
import Base.print_matrix
import Base.IteratorsMD
import Base.LinearIndices
import Base.map
import Statistics.mean
import Statistics.median
import Statistics.std
###

#################################
# Abstract types and type aliases
abstract type QObj <: Any end

# Note implementation is light-weight, using mostly type aliases
# for the features in Base that do all the real work
const Time = Float64
const QComp = ComplexF64
const QName = AbstractString
const QOp = AbstractArray
const QKet = AbstractVector
const QOps = Dict{AbstractString, QOp}
###

#######################
# Include functionality
include("Quantum.jl")
include("QuantumOscillator.jl")
include("QuantumBit.jl")
include("QuantumOperations.jl")
include("QuantumEvolution.jl")
###

#########
# Exports

# QuantumBayesian
export Time, QObj, QComp, QInd, QName, QOp, QKet, QOps
# Quantum
export QFactor, QSpace, QView
export size, length, show, showarray, CartesianIndices, LinearIndices, getindex, setindex!
export name, factors, unview, subview
export superket, unsuperket, superopl, superopr
export ⊗, lift, ptrace, bra
export osc, qubit
export groundvec, ground, projector, transition, coherentvec, coherent
export comm, acomm, ⊖, ⊕, sand, diss, inn
export scomm, sacomm, ssand, sdiss
export expect, expectvec, weakvalue, weakvaluevec
# QuantumEvolution
export Trajectory, Ensemble
export map, mean, median, std
export size, length, ndims, print_matrix
export ham, ham_rk4, lind, lind_rk4
export sham, slind
export meas, trajectory, ensemble, trajectoryToEnsemble
###

end # module

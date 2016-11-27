__precompile__()
module QuantumBayesian

#################################
# Package dependences and imports
#
using Distributions
import Base.product

# Imported solely for method overloading purposes
import Base.call
import Base.show
import Base.showarray
import Base.length
import Base.size
import Base.getindex
import Base.sub2ind
import Base.ind2sub
import Base.dot
###

#################################
# Abstract types and type aliases
abstract QObj <: Any

# Note implementation is light-weight, using mostly type aliases
# for the features in Base that do all the real work

# Define efficient (sparse) representations for Quantum 
# objects that keep track of tensor product factors properly
# for the purposes of partial traces 
typealias QComp Complex128
typealias QInd Int
typealias QName AbstractString
typealias QOp{T,I} SparseMatrixCSC{T, I}
typealias QKet{T,I} SparseVector{T, I}
typealias QOps{T,I} Dict{AbstractString, QOp{T,I}}
###

#######################
# Include functionality
include("Quantum.jl")
include("QuantumEvolution.jl")
###

#########
# Exports

# QuantumBayesian
export QObj, QComp, QInd, QName, QOp, QKet, QOps
# Quantum
export QFactor, QSpace, QView
export size, length, name, show, showarray, sub2ind, ind2sub, getindex, setindex!
export factors, unview, subview
export ⊗, lift, ptrace, dot, ⋅, bra
export osc, qubit
export groundvec, ground, projector, transition, coherentvec, coherent
export comm, acomm, ⊖, ⊕, sand, diss, inn
export expect, expectvec, weakvalue, weakvaluevec 
# QuantumEvolution
export ham, lind, lind_runge, trajectory
###

end # module

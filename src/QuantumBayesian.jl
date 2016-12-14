__precompile__()

"""
    module QuantumBayesian

A Julia package for efficient simulation of quantum state evolution.
"""
module QuantumBayesian

#################################
# Package dependences and imports
#

#using Distributions  # Minimize depedencies 
using Interpolations
import Base.product

# Imported solely for method overloading purposes
import Base.call
import Base.show
import Base.showarray
import Base.length
import Base.size
import Base.getindex
import Base.setindex!
import Base.ndims
import Base.indices
import Base.print_matrix
import Base.sub2ind
import Base.ind2sub
import Base.dot
import Base.map
import Base.mean
import Base.median
import Base.std
###

#################################
# Abstract types and type aliases
abstract QObj <: Any

# Note implementation is light-weight, using mostly type aliases
# for the features in Base that do all the real work

# Define efficient (sparse) representations for Quantum 
# objects that keep track of tensor product factors properly
# for the purposes of partial traces 
typealias Time Float64
typealias QComp Complex128
typealias QInd Int
typealias QName AbstractString
typealias QOp SparseMatrixCSC
typealias QKet SparseVector
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
export Time, QObj, QComp, QInd, QName, QOp, QKet, QOps
# Quantum
export QFactor, QSpace, QView
export size, length, show, showarray, sub2ind, ind2sub, getindex, setindex!
export name, factors, unview, subview
export superket, unsuperket, superopl, superopr
export ⊗, lift, ptrace, dot, ⋅, bra
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
export meas, trajectory, ensemble
###

end # module

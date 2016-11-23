__precompile__

module QuantumBayesian

include("Quantum.jl")
include("QuantumEvolution.jl")

# Quantum
export QObj, QComp, QInd, QName, QOp, QKet, QOps
export QFactor, QSpace, QView
export size, length, name, show, showarray, sub2ind, ind2sub, getindex, setindex!
export factors, unview, subview
export ⊗, lift, ptrace, dot, ⋅, bra
export osc, qubit
export groundvec, ground, projector, transition, coherentvec, coherent
export comm, acomm, ⊖, ⊕, diss, inn, expect, weakvalue, weakvaluevec 
# QuantumEvolution
export ham, lind, lind_runge, trajectory

end # module

### QuantumOperations.jl
# Convenience functionality for handling common operations
###

"""
    bra(a::QKet)

Construct a dual-vector (one-form) from a ket vector `a`.

Such a one-form has type: (b::QKet) -> dot(a,b)
"""
bra(a::QKet) = (b::QKet) -> dot(a,b)

"""
    comm(a, b) = a ⊖ b

Commutator of operators `a` and `b`.

### Returns:
  - Anti-Hermitian operator: a * b' - b * a'
"""
comm(a, b) = a * b' - b * a'
⊖ = comm

"""
    scomm(a)

Superoperator for commutator with operator `a`.
Assumes Hermitian superket.

### Returns:
  - Superoperator: scomm(a) * superket(b) == superket(a * b - b * a')
"""
scomm(a) = superopl(a) - superopr(a')

"""
    acomm(a, b) = a ⊕ b

Anti-commutator of operators `a` and `b`.

### Returns:
  - Hermitian operator: a * b' + b * a'
"""
acomm(a, b) = a * b' + b * a'
⊕ = acomm

"""
    sacomm(a)

Superoperator for anticommutator with operator `a`.
Assumes Hermitian superket.

### Returns:
  - Superoperator: scomm(a) * superket(b) == superket(a * b + b * a')
"""
sacomm(a) = superopl(a) + superopr(a')

"""
    sand(a, b)

Sandwich `b` operator with `a`.

### Returns:
  - Operator: a * b * a'
"""
sand(a, b) = a * b * a'

"""
    ssand(a)

Superoperator for sandwich with operator `a`.

### Returns:
  - Superoperator: ssand(a) * superket(b) == superket(a * b * a')
"""
ssand(a) = superopl(a) * superopr(a')

"""
    diss(a)

Dissipation function for `a` action.

### Returns:
  - Function: ρ -> sand(a, ρ) - acomm(at*a, ρ)/2
"""
diss(a) = ρ -> sand(a, ρ) - acomm(a'*a, ρ)/2

"""
    sdiss(a)

Dissipation superoperator for `a` action.

### Returns:
  - ssand(a) - sacomm(at*a)/2
"""
sdiss(a) = ssand(a) - sacomm(a'*a)/2

"""
    inn(a)

Innovation function for `a` action.

### Returns:
  - Function: ρ -> acomm(a, ρ) - tr(acomm(a, ρ))*ρ
"""
inn(a)  = ρ -> let ac = acomm(a, ρ);  ac - tr(ac)*ρ end

"""
    expect(ρ, op)

Expectation value of `op` in state `ρ`.

### Returns:
  - tr(ρ * op) / tr(ρ)
"""
expect(ρ, op) = tr(ρ * op) / tr(ρ)

"""
    expectvec(ψ, op)

Expectation value of `op` in state vector `ψ`.

### Returns:
  - dot(ψ, op * ψ) / dot(ψ, ψ)
"""
expectvec(ψ, op) = dot(ψ, op * ψ) / dot(ψ, ψ)

"""
    weakvalue(ρi, ρf, op)

Generalized weak value of `op` with initial state `ρi`
and final state `ρf`.

### Returns:
  - tr(ρf * op * ρi) / tr(ρf * ρi)
"""
weakvalue(ρi, ρf, op) = tr(ρf * op * ρi) / tr(ρf * ρi)

"""
    weakvaluevec(ψi, ψf, op)

Weak value of `op` with initial state `ψi` and final state `ψf`.

### Returns:
  - dot(ψf, op * ψi) / dot(ψf, ψi)
"""
weakvaluevec(ψi, ψf, op) = dot(ψf, op * ψi) / dot(ψf, ψi)

### QuantumOscillator.jl
# Convenience functionality for handling finite-dimensional systems,
# specifically related to the (truncated) harmonic oscillator
###

"""
    osc(n::Int[, name="Osc(n)"::QName])

Create harmonic oscillator in Fock basis with `n` levels.

### Default ops for QFactor:
  - 'i' : identity operator
  - 'n' : number operator
  - 'd' : lowering operator
  - 'u' : raising operator
  - 'x' : in-phase quadrature
  - 'y' : out-of-phase quadrature

"""
function osc(levels::Int, name=""::QName)
    if name == ""
        name = "Osc($(levels))"
    end
    s = QFactor(levels, name)
    s.ops['d'] = sparse([x == y - 1 ? sqrt(QComp(x)) : QComp(0) for x=1:levels, y=1:levels])
    s.ops['u'] = s('d')'
    s.ops['n'] = let l=1:levels; sparse(l,l,map(QComp, 0:(levels - 1))) end
    s.ops['x'] = s('d') + s('u')
    s.ops['y'] = (s('d') - s('u')) .* im
    s
end

"""
    coherentvec(o::QObj, α::Number)::QKet

Create a coherent state ket vector with amplitude `α`.
"""
# Generate coherent state α in space o
function coherentvec(o::QObj, α::Number)::QKet
  m = length(o)
  nbar = abs2(α)
  @assert (nbar + 3*sqrt(nbar) <= m) "Mean n of $nbar too large for max n of $m."
  e = exp(-nbar/2)
  cv = map(k -> QComp(e*α^k/sqrt(gamma(k+1))), 0:(m - 1))
  sparsevec(cv / norm(cv))
end

"""
    coherent(o::QObj, α::QComp)::QKet

Create a coherent state projection operator with amplitude `α`.
"""
coherent(o::QObj, α) = projector(coherentvec(o, α))

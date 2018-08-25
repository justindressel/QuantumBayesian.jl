### QuantumBit.jl
# Convenience functionality for working with qubit-based systems
###

"""
    qubit([name="Qubit"::QName])

Create qubit in computational basis.

**Note:** Rows are permuted from usual matrix representation.
This is so the 0 state can represent the ground state of energy
when constructing a qubit Hamiltonian.

**Example:**
```julia
q = qubit()
q('z') * groundvec(q) == sparse([QComp(-1); 0])  # True
H = (ħ*ωq/2)*q('z') # Correct energy structure
```

### Default ops for QFactor:
  - 'i' : identity operator
  - 'd' : lowering operator (``σ_-``)
  - 'u' : raising operator  (``σ_+``)
  - 'x' : Pauli x operator  (``σ_x``)
  - 'y' : Pauli y operator  (``σ_y``)
  - 'z' : Pauli z operator  (``σ_z``)

"""
function qubit(name="Qubit")
    q = osc(2, name)
    # Replace the number operator with a rescaling to get Pauli z
    merge!(q.ops, Dict('z'=> 2 .* q('n') .- q('i')))
    delete!(q.ops,'n')
    q
end

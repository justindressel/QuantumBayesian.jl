### Quantum.jl
# Core functionality for setting up quantum Hilbert spaces
# of operators, handling tensor products, and computing
# partial traces
###


###
# Quantum Hilbert Space Factor
###
"""
    QFactor(dim::Int, name::QName[, ops::QOps])

Single Hilbert space factor for a quantum space.

### Fields:
  - dim : Hilbert space dimension
  - name : String naming the factor for clarity
  - ops  : Dict of named operators belonging to space

### Bases:
  - Indexing as an array with [] indexes from 1
  - Indexing as a function with () indexes from 0
"""
struct QFactor <: QObj
    dim  :: Int
    name :: QName
    ops  :: QOps
    # Input sanity check
    function QFactor(d::Int, n::QName, ops::QOps)
        (d > 0) || error("Dimension must be positive.")
        new(d,n,ops)
    end
end
# Simplified constructor that auto-creates identity operator
QFactor(dim::Int, name::QName) =
    QFactor(dim, name,
            QOps('i' => sparse(QComp(1.0)SparseArrays.I, dim, dim))
            )

# Helper functions
size(s :: QFactor) = s.dim
length(s :: QFactor) = s.dim
name(s :: QFactor) = s.name

# Overload function call of string to get operator
(q::QFactor)(arg::QOpName) = q.ops[arg]

###
# Quantum Tensor Product of Hilbert Spaces
###
"""
    QSpace(dim::Int, name::QName)
    QSpace(f::QFactor)
    QSpace(fs::Vector{QFactor}, ops::QOps)

Tensor product containing quantum Hilbert space factors.

### Fields:
  - factors : Vector of Hilbert space factors
  - ops     : Dict of named operators belonging to space

### Bases:
  - Indexing as an array with [] indexes from 1
  - Indexing as a function with () indexes from 0
"""
struct QSpace <: QObj
    factors :: Vector{QFactor}
    ops     :: QsOps
end
# Convenient constructors
QSpace(f::QFactor) = QSpace(vec([f]),
                            QsOps([(string(k),v) for (k,v) in f.ops]))
QSpace(dim::Int, name::QName) = QSpace(QFactor(dim, name))

# Helper functions
factors(s :: QSpace) = s.factors
@inline size(s :: QSpace) = tuple(map(size, factors(s))...)
length(s :: QSpace) = prod(size(s))
name(s :: QSpace) = join(map(name, factors(s)), " ⊗ ")

# Overload function call
# Memoize creation of tensor product of composite operators for efficiency
(q::QSpace)(arg::AbstractString) = begin
    if haskey(q.ops, arg)
        q.ops[arg]
    elseif length(arg) == length(factors(q))
        q.ops[arg] = ⊗([q.factors[i](c) for (i,c) in enumerate(arg)]...)
    else
        error("QSpace does not contain named operator")
    end
end


###
# Pretty printing
###
function show(io::IO, q::QObj)
    println(io, "$(typeof(q)): ", name(q))
    println(io, "Dims  : ", size(q))
    println(io, "Ops   : '", join(keys(q.ops), "', '"), "'")
end

###
# Subviews
###
"""
    subview(s::QObj, op::QOp)

Return a reshaped view on quantum operator that makes
subsystem indices transparent.

`t[i,j]` is equivalent to ``T_{i;j} = a_{i} ⊗ b_{j}``.

`t[i,k,j,l]` is equivalent to ``T_{ik;jl} = a_{ij} ⊗ b_{kl}``.

**Example:**
```julia
a = ["a11" "a12" ; "a21" "a22"];
b = ["b11" "b12" ; "b21" "b22"];
c = ["c11" "c12" "c13" ; "c21" "c22" "c23"; "c31" "c32" "c33"];
abc = kron(a,b,c);
sys = qubit() ⊗ qubit() ⊗ osc(3);
abcv = subview(sys, abc);
"a21b12c31" == abcv[2,1,3,1,2,1] # True
```
"""
@inline function subview(s :: QSpace, op :: AbstractArray)
    # Compute correct permutation of dimensions
    # for the reshaped view of the array in memory
    dim = size(s)
    d = length(dim)
    rdim = reverse(dim)
    l = length(size(op))
    if l==1
        # Note: for a dimension list (2, 3, 4)
        #       the reshape should be: reshape(o, 4, 3, 2)
        #       with permutation: [3, 2, 1]
        PermutedDimsArray(reshape(op, rdim...), reverse(1:d))
    elseif l==2
        # Note: for a dimension list (2, 3, 4)
        #       the reshape should be: reshape(o, 4, 3, 2, 4, 3, 2)
        #       with permutation: [3, 2, 1, 6, 5 ,4]
        perm = vcat(reverse(1:d), reverse((d+1):(2*d)))
        PermutedDimsArray(reshape(op, rdim..., rdim...), perm)
    else
        error("subview assumes an initial 1D or 2D array")
    end
end
subview(s::QFactor, op::AbstractArray) = subview(QSpace(s), op)

"""
    superket(op::QOp)

Convert a square matrix operator into a ket vector in the superoperator space.
"""
superket(op::QOp) = reshape(op, length(op), 1)

"""
    unsuperket(op::AbstractArray)

Convert a superket back into a square matrix operator.
"""
function unsuperket(op::QOp)
    d = Int(sqrt(length(op)))
    d^2 == length(op) || error("Not a superket")
    reshape(op, d, d)
end

"""
    superopl(op::QOp)

Convert `op` into a superoperator that performs left multiplication
by `op`. `superopl(op)*superket(v) == superket(op * v)`
"""
superopl(op::QOp) = sparse(QComp(1.0)SparseArrays.I, size(op)) ⊗ op

"""
    superopr(op::QOp)

Convert `op` into a superoperator that performs right multiplication
by `op`. `superopr(op)*superket(v) == superket(v * op)`
"""
superopr(op::QOp) = transpose(op) ⊗ sparse(QComp(1.0)SparseArrays.I, size(op))

###################################################
# Workhorse functions
###

###
# Tensor Product
###
"""
    ⊗(q::AbstractArray, s::AbstractArray)
    ⊗(q::QObj, s::QObj)

Tensor product.

  - equivalent to `kron` for arrays
  - constructs tensor-product QSpaces from QFactors

"""
@inline ⊗(q::AbstractArray...) = kron(q...)
@inline function ⊗(q::QSpace...)
    factors = vcat(map(s->s.factors, q)...)
    QSpace(factors, QsOps())
end
@inline ⊗(q::QFactor...) = ⊗(map(QSpace,q)...)
⊗(q::QFactor, s::QSpace) = ⊗(QSpace(q),s)
⊗(s::QSpace, q::QFactor) = ⊗(s,QSpace(q))


###
# Lift to Joint Space
###

"""
    lift(q::QSpace, o::QOp, i::Int)

Lift an operator of a single factor into a joint space,
assuming that the factor is at position `i` of the tensor
product.

### Returns:
  - QOp : result of tensor product with appropriate identities
"""
@inline function lift(q::QSpace, o::QOp, i::Int)
    is = map(f -> f('i'), factors(q))
    splice!(is, i, [o])
    ⊗(is...)
end

###
# Partial Trace
###

"""
    ptr(s::QObj, o::QOp, n::Int...)

Partial trace of o, over subsystems in positions n, inside quantum space s.

### Returns:
  - (snew, onew)
  - snew : reduced system space
  - onew : reduced operator

"""
@inline ptr(s::QFactor, o::QOp, n::Int...) = tr(o)
@inline function ptr(s::QSpace, o::QOp, n::Int)
    v = subview(s,o)
    # Check that valid subsystem was specified
    dim = size(v)
    nsys = Int(fld(length(dim),2))
    (n > 0 && n <= nsys) || error("Partial trace over nonexistent system.")
    # Find second index matching n
    n2 = n+nsys
    # Helper function to truncate indices
    function trunc(l)
        l2 = [l...]
        deleteat!(l2, (n, n2))
        tuple(l2...)
    end

    # Remove indicated subsystem from system list
    sysind = [trunc(collect(1:length(dim)))[1:(nsys-1)]...]
    # Create new QSpace without indicated subsystem
    newsys = foldl(⊗, s.factors[sysind])
    if typeof(newsys)==QFactor
        newsys = QSpace(newsys)
    end
    # Create empty sparse matrix to hold result of trace
    m = spzeros(QComp, size(newsys(repeat("i",length(sysind))))...)
    # Create new view with easy indexing
    mv = subview(newsys, m)

    # Extract nonempty indices from original view
    vinds = CartesianIndices(v)
    oinds = LinearIndices(o)
    indices = []
    rows = rowvals(o)
    for c in 1:last(size(o))
        for r in nzrange(o, c)
            push!(indices, Tuple(vinds[oinds[rows[r], c]]))
        end
    end
    # Group indices that match when truncated
    for i in Set(map(trunc, indices)) # Set ensures no duplicates
        # Isolate diagonal elements within subgrouping, keep duplicates
        diagi = filter(e -> e[n]==e[n2], filter(j->trunc(j)==i, indices))
        if !isempty(diagi)
            # Set new sparse matrix elements as sums of diagonals
            mv[i...] = sum(map(j->v[j...], diagi))
        end
    end
    # Return both the new space, and the partially traced matrix
    newsys, m
end
@inline function ptr(s::QSpace, o::QOp, n::Int...)
    isempty(n) && return tr(o)
    n = sort([n...])
    s2, o2 = ptr(s,o,last(n))
    ptr(s2, o2, n[1:end-1]...)
end


#####################################################################
# Convenience functions for generating operators and states
###

###
# Sparse matrix basis generators by index
###

# Using (i) index for a 1D (vector) representation
Base.@propagate_inbounds Base.getindex(A::QFactor, i::Int) =
    sparsevec([i], [QComp(1)], length(A))
# Using (i,j) indices for the 2D representation
Base.@propagate_inbounds Base.getindex(A::QFactor, i::Int, j::Int) =
    sparse([i], [j], [QComp(1)], length(A), length(A))

# Using (i1,j1,k1,...,i2,j2,k2,...) indices for subsystem representation
@inline Base.@propagate_inbounds function Base.getindex(A::QSpace, inds::Int...)
    l = length(size(A))
    fs = factors(A)
    if length(inds) == l  # 1D (vector)
        ⊗([f[i] for (f,i) in zip(fs,inds)]...)
    elseif length(inds) == 2*l  # 2D (matrix)
        subinds = reshape(collect(inds), l, 2)
        ⊗([f[subinds[i,:]...] for (i,f) in enumerate(fs)]...)
    else
        error("Incorrect indexing.")
    end
end

# Overload function calls to match indexing from 0
@inline (q::QFactor)(i::Int...) = q[(1 .+ [i...])...]
@inline (q::QSpace)(i::Int...) = q[(1 .+ [i...])...]

###
# Ground states
###

"""
    groundvec(o::QObj)

Ground state vector of quantum space `o`.
"""
groundvec(o::QObj) = o[1]

"""
    groundvec(o::QObj)

Ground state density matrix of quantum space `o`.
"""
ground(o::QObj) = o[1,1]

"""
    projector(ψ::QKet)

Convert ket vector ψ into a projection operator.
"""
function projector(ψ::QKet)
    first(size(ψ)) == 1 && error("Not a column vector.")
    ρ = ψ * ψ'
    sparse(ρ / tr(ρ))
end
transition(ψl,ψr) = sparse(ψl * ψr') / (norm(ψl)*norm(ψr))

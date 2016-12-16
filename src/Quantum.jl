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
immutable QFactor <: QObj
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
QFactor(dim::Int, name::QName) = QFactor(dim, name, QOps("i" => speye(dim)))

# Helper functions
size(s :: QFactor) = s.dim
length(s :: QFactor) = s.dim
name(s :: QFactor) = s.name

# Overload function call of string to get operator
(q::QFactor)(arg::AbstractString) = q.ops[arg]

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
immutable QSpace <: QObj
    factors :: Vector{QFactor}
    ops     :: QOps
end
# Convenient constructors 
QSpace(f::QFactor) = QSpace(vec([f]), f.ops)
QSpace(dim::Int, name::QName) = QSpace(QFactor(dim, name))

# Helper functions
factors(s :: QSpace) = s.factors
@inline size(s :: QSpace) = tuple(map(size, factors(s))...)
length(s :: QSpace) = prod(size(s))
name(s :: QSpace) = join(map(name, factors(s)), " ⊗ ")

# Overload function call 
(q::QSpace)(arg::AbstractString) = q.ops[arg]


###
# Pretty printing
###
function show(io::IO, q::QObj)
    println(io, "$(typeof(q)): ", name(q))
    println(io, "Dims  : ", size(q))
    println(io, "Ops   : \"", join(keys(q.ops), "\", \""), '"')
end

###
# QView
###
"""
    QView{T,N,A<:AbstractArray,P} <: AbstractArray{T,N}

View on quantum operator that makes subsystem indices transparent.
Use the command ```subview``` to create a view.

**Example:**  ``T_{i;j} = a_{i}⊗b_{j}``
```julia
t = subview(sysa ⊗ sysb, a ⊗ b)
t[i,j] == a[i] * b[j]  # True
```

**Example:**  ``T_{ik;jl} = a_{ij}⊗b_{kl}``
```julia
t = subview(sysa ⊗ sysb, a ⊗ b)
t[i,k,j,l] == a[i,j] * b[k,l]  # True
```

### Fields:
  - data : reshaped array with raw subsystem indices
  - perm : permutation array to fix order of subsystem indices
"""
immutable QView{T,N,A<:AbstractArray,P} <: AbstractArray{T,N}
    data :: A
    perm :: P
end
# Suppress show output, since stored array differs from apparent indexing
Base.show(io::IO, A::QView) = println(io, "QView")
Base.showarray(io::IO, A::QView) = println(io, "QView")
Base.showarray(io::IO, A::QView, b::Bool) = println(io, "QView")
# Report size and length transparently, with correct permutation
Base.size(A::QView) = size(A.data)[A.perm]
Base.length(A::QView) = prod(size(A))
# Overload indexing so that the correctly permuted indices are actually called
@inline Base.@propagate_inbounds Base.getindex(A::QView, i::Int...) = 
    A.data[i[A.perm]...]
@inline Base.@propagate_inbounds Base.getindex(A::QView, i::Int) = 
    A.data[ind2sub(A.data,i)...]
@inline Base.@propagate_inbounds Base.setindex!(A::QView, v, i::Int...) =
    setindex!(A.data, v, i[A.perm]...)
@inline Base.sub2ind(A::QView, i::Int...) = sub2ind(A.data, i[A.perm]...)
@inline Base.ind2sub(A::QView, i::Int) = ind2sub(A.data, i)[A.perm]

"""
    unview(op::QView)

Extract original 2D sparse matrix from view.
"""
unview(op::QView) = op.data.parent


"""
    subview(s::QObj, op::QOp)

Return a QView, which is a view on quantum operator that makes 
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
        perm = reverse(1:d)
        op2 = reshape(op, rdim...)
    elseif l==2
        # Note: for a dimension list (2, 3, 4)
        #       the reshape should be: reshape(o, 4, 3, 2, 4, 3, 2)
        #       with permutation: [3, 2, 1, 6, 5 ,4]
        perm = vcat(reverse(1:d), reverse((d+1):(2*d)))
        op2 = reshape(op, rdim..., rdim...)
    else
        error("QView assumes an initial 1D or 2D array")
    end
    QView{eltype(op2),Int,typeof(op2),typeof(perm)}(op2, perm)
end
subview(s::QFactor, op::AbstractArray) = subview(QSpace(s),op)

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
    if :parent in fieldnames(op)
        op.parent
    else
        d = Int(sqrt(length(op)))
        d^2 == length(op) || error("Not a superket")
        reshape(op, d, d)
    end
end

"""
    superopl(op::QOp)

Convert `op` into a superoperator that performs left multiplication
by `op`. `superopl(op)*superket(v) == superket(op * v)`
"""
superopl(op::QOp) = speye(last(size(op))) ⊗ op

"""
    superopr(op::QOp)

Convert `op` into a superoperator that performs right multiplication
by `op`. `superopr(op)*superket(v) == superket(v * op)`
"""
superopr(op::QOp) = (op.') ⊗ speye(last(size(op)))

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
    j = QSpace(factors, Dict())
    # Compute all combinations of operator products
    for qq in product(map(s->s.ops, q)...) # Uses Iterators
        key = prod(map(first, qq))
        j.ops[key] = kron(map(last, qq)...)
    end
    j
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
    is = map(f -> f("i"), factors(q))
    splice!(is, i, [o])
    ⊗(is...)
end

###
# Partial Trace
###

"""
    ptrace(s::QObj, o::QView, n::Int...)
    ptrace(s::QObj, o::QOp, n::Int...)

Partial trace of o, over subsystems in positions n, inside quantum space s.

### Returns:
  - (snew, onew)
  - snew : reduced system space
  - onew : reduced operator 

"""
@inline ptrace(s::QFactor, v::QView, n::Int...) = trace(unview(v))
@inline ptrace(s::QFactor, v::QOp, n::Int...) = trace(v)
@inline ptrace(s::QSpace, o::AbstractArray, n::Int...) = ptrace(s,subview(s,o),n...)
@inline function ptrace(s::QSpace, v::QView, n::Int)
    # Check that valid subsystem was specified
    dim = size(v)
    nsys = Int(length(dim)/2)
    (n > 0 && n <= nsys) || error("Partial trace over nonexistent system.")
    # Find second index matching n
    n2 = n+nsys
    # Helper function to truncate indices
    function trunc(l)
        l2 = [l...]
        deleteat!(l2, (n, n2))
        (l2...)
    end

    # Remove indicated subsystem from system list
    sysind = [trunc(collect(1:length(dim)))[1:(nsys-1)]...]
    # Create new QSpace without indicated subsystem
    newsys = foldl(⊗, s.factors[sysind])
    # Create empty sparse matrix to hold result of trace
    m = spzeros(QComp, size(newsys(repeat("i",length(sysind))))...)
    # Create new view with easy indexing
    mv = subview(newsys, m)
    
    # Extract nonempty indices from original view
    vo = unview(v)
    psize = size(vo)
    rows = rowvals(vo)
    indices = Set() # Use a set to eliminate duplicates
    for c in 1:last(psize)
        for r in nzrange(vo, c)
            push!(indices, ind2sub(v, sub2ind(psize, rows[r], c)))
        end
    end
    # Group indices that match when truncated
    for i in map(trunc, indices) # Set ensures no duplicates
        # Isolate diagonal elements within subgrouping, keep duplicates
        diagi = filter(e -> e[n]==e[n2], filter(j->trunc(j)==i, collect(indices)))
        if !isempty(diagi)
            # Set new sparse matrix elements as sums of diagonals
            mv[i...] = sum(map(j->v[j...], diagi))
        end
    end
    # Return both the new space, and the partially traced matrix
    newsys, m
end
@inline function ptrace(s::QSpace, v::QView, n::Int...)
    isempty(n) && return trace(unview(v))
    n = sort([n...]) 
    s2, v2 = ptrace(s,v,last(n))
    ptrace(s2, v2, n[1:end-1]...)
end


#####################################################################
# Convenience functions for generating operators and states
###

###
# Sparse matrix basis generators by index
###

# Using (i) index for a 1D (vector) representation
Base.@propagate_inbounds Base.getindex(A::QObj, i::Int) = 
    sparsevec([i],[QComp(1)],length(A))
# Using (i,j) indices for the 2D representation
Base.@propagate_inbounds Base.getindex(A::QObj, i::Int, j::Int) = 
    sparse([i],[j],[QComp(1)],length(A),length(A))

# Using (i1,j1,k1,...,i2,j2,k2,...) indices for subsystem representation
@inline Base.@propagate_inbounds function Base.getindex(A::QObj, inds::Int...) 
    l = length(size(A))
    Ai = subview(A,A(repeat("i",l)))
    if length(inds) == 2*l
        i, j = ind2sub(unview(Ai),sub2ind(Ai, inds...))
        sparse([i],[j],[QComp(1)],length(A),length(A))
    elseif length(inds) == l
        i, j = ind2sub(unview(Ai),sub2ind(Ai, inds...,inds...))
        sparsevec([i],[QComp(1)],length(A))
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
    sparse(ρ / trace(ρ))
end
transition(ψl,ψr) = sparse(ψl * ψr') / (vecnorm(ψl)*vecnorm(ψr))


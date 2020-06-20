module QuantumBayesian

    using QuantumOptics
    using Interpolations
    using Distributed
    import Base: length, size, ndims, map

    #################################
    # Abstract types and type aliases
    abstract type QObj <: Any end

    # Note implementation is light-weight, using mostly type aliases
    # for the features in Base that do all the real work
    const Time = Float64
    const Efficiency = Float64
    const QOp = AbstractOperator
    const QKet = Ket
    ###

    ###
    # Trajectory
    #
    #   Enhanced Vector type to store temporal information for
    #   convenience, and auto-interpolate solution temporally
    ##

    struct Trajectory{T} <: AbstractVector{T}
        v :: AbstractVector{T}
        t :: AbstractRange{Time}
        interpolation
        # Input auto-interpolation
        function Trajectory{T}(t::AbstractRange{Time}, v::AbstractVector{T}) where T
            i = interpolate( (collect(t),), v, Gridded(Constant()) )
            new(v, t, i)
        end
    end
    @inline Base.@propagate_inbounds Base.getindex(T::Trajectory, i::Int...) = T.v[i...]
    (T::Trajectory)(t::Time) = T.interpolation[t]
    size(T::Trajectory) = size(T.v)
    length(T::Trajectory) = length(T.v)

    ###
    # Ensemble
    #
    #   Enhanced Array type to store temporal information for
    #   convenience, and easily give statistics over the collection
    ##

    struct Ensemble{T} <: AbstractArray{T, 2}
        a :: AbstractArray{T, 2}
        t :: AbstractRange{Time}
        n :: Integer
        function Ensemble{T}(t::AbstractRange{Time}, n::Integer, ::Type{T}) where {T}
            l = length(t)
            a = zeros(T, l, n)
            new(a, t, n)
        end
        function Ensemble{T}(t::AbstractRange{Time}, a::AbstractArray{T, 2}) where {T}
            length(t) == first(size(a)) || error("Incompatible array size")
            n = last(size(a))
            new(a, t, n)
        end
        function Ensemble{T}(tr::Trajectory{T}) where {T}
    	let e = Ensemble(tr.t, 1, T)
    	    e[1] = tr.v
                e
            end
        end
    end
    Ensemble(t::R, n::Integer, T::Type) where R<:AbstractRange{Time} = Ensemble{T}(t, n, T)
    Ensemble(t::R, a::AbstractVector{T}) where {T,R<:AbstractRange{Time}} = Ensemble{T}(t,reshape(a,length(a),1))
    Ensemble(t::R, a::AbstractArray{T,2}) where {T,R<:AbstractRange{Time}} = Ensemble{T}(t,a)
    Ensemble(tr::Trajectory{T}) where {T} = Ensemble{T}(tr)
    @inline Base.@propagate_inbounds Base.getindex(e::Ensemble, i::Int...) = e.a[i...]
    @inline Base.@propagate_inbounds Base.setindex!(e::Ensemble, v::AbstractVector, n::Int) = Base.setindex!(e.a,v,:,n)
    (e::Ensemble)(n::Integer) = Trajectory{eltype(e)}(e.t, e.a[:,n])
    size(e::Ensemble) = size(e.a)
    length(e::Ensemble) = length(e.a)
    ndims(e::Ensemble) = ndims(e.a)
    #Base.indices(e::Ensemble, d::Integer) = indices(e.a, d)

    function map(f::Function, e::Ensemble{T}) where {T}
        max = length(e.t)
        tr = zeros(T, max)
        for i in 1:max
            @inbounds tr[i] = f(e.a[i,:])
        end
        Trajectory{T}(e.t, tr)
    end
    mean(e::Ensemble) = map(mean, e)
    median(e::Ensemble) = map(median, e)
    std(e::Ensemble) = map(std, e)

    ###########################################################
    # Utilities
    ##

    """
        comm(a, b)

    Commutator product of `a` and `b`

    ### Returns:
      - Anti-Hermitian Operator: a * b - b * a
    """
    comm(a,b) = a * b' - b * a'

    """
        scomm(a)
    Superoperator for commutator with operator `a`.
    Assumes Hermitian superket.
    ### Returns:
      - Superoperator: scomm(a) * superket(b) == superket(a * b - b * a')
    """
    scomm(a) = spre(a) - spost(a')

    """
        acomm(a, b) = a ⊕ b
    Anti-commutator of operators `a` and `b`.
    ### Returns:
      - Hermitian operator: a * b' + b * a'
    """
    acomm(a, b) = a * b' + b * a'

    """
        sacomm(a)
    Superoperator for anticommutator with operator `a`.
    Assumes Hermitian superket.
    ### Returns:
      - Superoperator: scomm(a) * superket(b) == superket(a * b + b * a')
    """
    sacomm(a) = spre(a) + spost(a')

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
    ssand(a) = spre(a) * spost(a')

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

    ###########################################################
    # Simple propagators
    ##

    # Hamiltonian propagation
    """
        ham(dt::Time, H::QOp; ket=false)

    Return increment function for Hamiltonian evolution generated
    by `H` over a time step `dt`.

    Uses an exact (dense) matrix exponential, assuming no time-dependence.

    ### Returns:
      - ket=true  : (t::Time, ψ::QKet) -> u * ψ
      - ket=false : (t::Time, ρ::QOp)  -> u * ρ * u'

    """
    function ham(dt::Time, H::QOp; ket=false)
        u::QOp = SparseOperator(exp( -im * dt * DenseOperator(H)))
        ut = u'
        if ket
            (t::Time, ψ::QKet) -> u * ψ
        else
            (t::Time, ρ::QOp) -> u * ρ * ut
        end
    end
    function ham(dt::Time, H::Function; ket=false)
        (t::Time, state) -> ham(dt, H(t), ket=ket)(t, state)
    end

    # Superoperator Hamiltonian evolution
    """
        sham(dt::Time, H::QOp)
    Return increment function using a superoperator for Hamiltonian
    evolution generated by `H` over a time step `dt`.
    Uses an exact (dense) matrix exponential, assuming no time-dependence.
    ### Returns:
      - (t::Time, ρvec) -> u*ρvec : Evolution superoperator
    """
    function sham(dt::Time, H::QOp)
        u::QOp = SparseOperator(exp( -im * dt * DenseOperator(H)))
        l = spre(u) * spost(u')
        (t::Time, ρvec) -> l * ρvec
    end
    function sham(dt::Time, H::Function)
        (t::Time, ρvec) -> sham(dt, H(t))(t, ρvec)
    end

    # Runge-Kutta Hamiltonian evolution
    """
        rinc(t::Time, ρ::QOp, inc::Function, dt::Time)
    4th-order Runge-Kutta algorithm given increment `inc`.
    ### Returns:
      - dρ :: QOp
    """
    function rinc(t::Time, ρ::QOp, inc::Function, dt::Time)::QOp
        dρ1::QOp = inc(t, ρ)
        dρ2::QOp = inc(t + dt/2, ρ + dρ1*dt/ 2)
        dρ3::QOp = inc(t + dt/2, ρ + dρ2*dt/ 2)
        dρ4::QOp = inc(t + dt, ρ + dρ3*dt)
        dt*(dρ1 + 2*dρ2 + 2*dρ3 + dρ4)/6
    end

    """
        ham_rk4(dt::Time, H::QOp; ket=false)
    Return increment function for Hamiltonian evolution generated
    by Hamiltonian `H` over a time step `dt`.
    Uses a 4th-order Runge-Kutta integration method to construct the state
    increment from the first-order differential (master) equation.
    ### Returns:
      - ket=true  : (t::Time, ψ::QKet) -> ψnew
      - ket=false : (t::Time, ρ::QOp)  -> ρnew
    """
    function ham_rk4(dt::Time, H::Function; ket=false)
        inc(t::Time, ρ::QOp)::QOp = - im * comm(H(t),ρ)
        if ket
            inc(t::Time, ψ::QKet)::QKet = - im * H(t) * ψ
        end
        (t::Time, ρ) -> ρ + rinc(t, ρ, inc, dt)
    end
    function ham_rk4(dt::Time, H::QOp; ket=false)
        h(t) = H
        ham_rk4(dt, h, ket=ket)
    end

    # Jump-nojump Lindblad propagator
    """
        lind(dt::Time[, H]; clist=QOp[], flist=Function[])

    Return increment function over a time step `dt` for Lindblad dissipative
    evolution generated by (an optional) Hamiltonian `H` (operator or function),
    a list `clist` of constant dissipative operators, and a list `flist` of
    time-dependent functions returning dissipative operators.

    Uses the "jump no-jump" method to efficiently approximate the exact
    Lindblad propagator as a composition of Hamiltonian evolution, jumps,
    and no-jump informational backaction. Assumes small dt.
    [Physical Review A **92**, 052306 (2015)]

    ### Returns:
      - (t::Time, ρ(t)::QOp) -> ρ(t+dt)
    """
    function lind(dt::Time; clist=QOp[], flist=Function[])
        ns = Function[]
        ds = Function[]
        # Construct operations for constant operators
        if isempty(clist)
            push!(ns, (t, ρ) -> ρ)
        else
            Id = identityoperator(first(clist).basis_l)
            op = DenseOperator(Id - dt * mapreduce(a -> a' * a, +, clist))
            n::QOp = SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
            push!(ns, (t, ρ) -> n * ρ * n)
            push!(ds, (t, ρ) -> mapreduce(a -> a * ρ * a', +, clist) * dt)
        end
        # Construct operations for time-dependent operators
        if isempty(flist)
            push!(ns, (t, ρ) -> ρ)
        else
            function nf(t)
                Id = identityoperator(first(flist)(t).basis_l)
                op = DenseOperator(Id - dt * mapreduce(a -> a(t)' * a(t), +, flist))
                return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
            end
            push!(ns, (t, ρ) -> nf(t) * ρ * nf(t))
            push!(ds, (t, ρ) -> mapreduce(a -> a(t) * ρ * a(t)', +, flist) * dt)
        end
        push!(ds, (t, ρ) -> last(ns)(t, first(ns)(t, ρ)))
        (t::Time, ρ) -> mapreduce(f -> f(t, ρ), +, ds)
    end
    function lind(dt::Time, H; clist=QOp[], flist=Function[])
        # Rely on Hamiltonian to specify type of H
        h = ham(dt, H)
        # Apply Hamiltonian first, then the Lindblad increment
        (t::Time, ρ) -> lind(dt, clist=clist, flist=flist)(t, h(t, ρ))
    end

    # Runge-Kutta Lindblad propagator
    """
        lind_rk4(dt::Time[, H]; clist=QOp[], flist=Function[])
    Return increment function over a time step `dt` for Lindblad dissipative
    evolution generated by (an optional) Hamiltonian `H` (operator or function),
    a list `clist` of constant dissipative operators, and a list `flist` of
    time-dependent functions returning dissipative operators.
    Uses a 4th-order Runge-Kutta integration method to construct the state
    increment from the first-order Lindblad differential (master) equation.
    ### Returns:
      - (t::Time, ρ(t)::QOp) -> ρ(t) + dρ
    """
    lind_rk4(dt::Time, H) = ham_rk4(dt, H)
    function lind_rk4(dt::Time; clist=QOp[], flist=Function[])
        inc(t::Time, ρ::QOp)::QOp =
            mapreduce(a -> diss(a)(ρ), +, clist; init=zero(ρ)) +
            mapreduce(a -> diss(a(t))(ρ), +, flist; init=zero(ρ))
        (t::Time, ρ) -> ρ + rinc(t, ρ, inc, dt)
    end
    function lind_rk4(dt::Time, H::Function; clist=QOp[], flist=Function[])
        inc(t::Time, ρ::QOp)::QOp = - im * comm(H(t),ρ) +
            mapreduce(a -> diss(a)(ρ), +, clist; init=zero(ρ)) +
            mapreduce(a -> diss(a(t))(ρ), +, flist; init=zero(ρ))
        (t::Time, ρ) -> ρ + rinc(t, ρ, inc, dt)
    end
    function lind_rk4(dt::Time, H::QOp; clist=QOp[], flist=Function[])
        h(t) = H
        lind_rk4(dt, h, clist=clist, flist=flist)
    end

    # Superoperator Lindblad propagator
    """
        slind(dt::Time, H::QOp; clist=QOp[], flist=Function[])
    Return increment function over a time step `dt` for Lindblad dissipative
    evolution generated by (an optional) Hamiltonian `H` (operator or function),
    a list `clist` of constant dissipative operators, and a list `flist` of
    time-dependent functions returning dissipative operators.
    Uses direct matrix exponentiation of the total superoperator for
    the increment.
    ### Returns:
      - (t::Time, ρvec) -> u * ρvec : Superoperator for total evolution over dt
    """
    slind(dt::Time, H) = sham(dt, H)
    function slind(dt::Time, H::QOp; clist=QOp[], flist=Function[])
        h = -im*scomm(H)
        init = DenseSuperOperator(h.basis_l, h.basis_r)
        if isempty(flist)
            let l = DenseSuperOperator(h.basis_l, h.basis_r, (h + mapreduce(sdiss, +, clist; init=init)).data),
                u = exp(dt*l);
                (t::Time, ρvec) -> u * ρvec
            end
        else
            let l(t::Time) = DenseSuperOperator(h.basis_l, h.basis_r, (h + mapreduce(sdiss, +, clist; init=init).data) +
                                 mapreduce(f -> sdiss(f(t)), +, flist; init=init)),
                u(t::Time) = exp(dt*l(t))
                (t::Time, ρvec) -> u(t) * ρvec
            end
        end
    end
    function slind(dt::Time, H::Function; clist=QOp[], flist=Function[])
        (t::Time, ρvec) -> slind(dt, H(t), clist=clist, flist=flist)(t, ρvec)
    end


    ###
    # Diffusive stochastic evolution
    ###

    """
        meas(dt::Time, H; mclist=Tuple{QOp,Time,Float64}[],
                          mflist=Tuple{Function,Time,Float64}[],
                          clist=QOp[], flist=Function[])

    Return increment function over a time step `dt` for diffusive monitored
    evolution generated by Hamiltonian `H`, a list `mclist` of tuples specifying
    constant measurement operators, and a list `mflist` of tuples specifying
    time-dependent measurement operators. The tuples have the form `(m, τ, η)`,
    where `m` is a (generally non-Hermitian) operator (or time-dependent function)
    that specifies the measurement backaction (see below), `τ` is the timescale
    of the measurement collapse, and `η` is the measurement efficiency (bounded
    between 0 and 1).

    These quantities are related as follows:
      - Γ = 1/(2τη) : ensemble measurement dephasing rate
      - Γm = 1/(2τ) : dephasing rate produced by averaging the collected signal
      - γ = Γ - Γm = (1-η)/(2τη) : residual dephasing rate from signal loss
      - ``m = m_o - i m_p``
      - ``m_o = (m + m')/2``  : measured observable (Hermitian part of `m`)
      - ``m_p = i(m - m')/2`` : phase-backaction generator (anti-Hermitian part of `m`)

    The backaction takes the form of a three-step process:
      1. Sample `r` from a Gaussian distribution with mean ``⟨m_o⟩`` and variance τ/dt
      2. Conditionally update the state with a purity-preserving Kraus operator:
      - ``M_r = (dt/2πτ)^{1/4} e^{ -i m_p r dt/(2τ) - dt(r - m_o)^2/(4τ) }``
      3. Apply residual Lindblad dephasing evolution, including
      - `m` with rate ``γ = (1-η)/(2τη)``
      - the dephasing operators in `clist` and `flist`
      - natural Hamiltonian evolution `H`

    Uses the "jump no-jump" method to efficiently approximate the residual
    Lindblad dephasing as a composition of Hamiltonian evolution, jumps,
    and no-jump informational backaction. Assumes no time-dependence in `alist`,
    and small dt.  [Physical Review A **92**, 052306 (2015)]

    ### Returns:
      - (t::Time, ρ(t)::QOp) -> (ρ(t+dt)::QOp, rlist::Float64...)

    """
    function meas(dt, H; mclist=Tuple{QOp,Time,Efficiency}[], mflist=Tuple{Any,Time,Efficiency}[], clist=QOp[], flist=Function[])
        # Assemble readout generating functions and Kraus operators
        ros = Function[]
        gks = Function[]
        as = QOp[]
        afs = Function[]
        for (m, τ, η) in mclist
            println(readout)
            push!(ros, readout(dt, m, τ))
            push!(gks, gausskraus(dt, m, τ))
            # For inefficient measurements append to Lindblad dephasing
            if abs(1.0 - η) > 1e-4
                push!(as, sqrt((1.0-η)/(4*τ*η))*m)
            end
        end
        for (m, τ, η) in mflist
            push!(ros, readout(dt, m, τ))
            push!(gks, gausskraus(dt, m, τ))
            # For inefficient measurements append to Lindblad dephasing
            if abs(1.0 - η) > 1e-4
                push!(afs, t -> sqrt((1.0-η)/(4*τ*η))*m(t))
            end
        end
        # Assemble total Lindblad dephasing using jump/no-jump method
        l = lind(dt, H, clist=vcat(clist, as), flist=vcat(flist, afs))
        # Increment that samples each readout, applies all Kraus operators
        # then applies Lindblad dephasing (including Hamiltonian evolution)
        (t::Time, ρ) -> let rs = map(ro -> ro(t, ρ), ros),
                            gs = map(z -> z[2](t, z[1]), zip(rs,gks)),
                            ρ1 = foldr(sand, gs; init=ρ);
                        (l(t, ρ1/tr(ρ1)), rs...) end
    end

    """
        readout(dt::Time, m::QOp, τ::Time)
        readout(dt::Time, m::Function, τ::Time)

    Return a function that accepts a (time t, state ρ) tuple and returns a random
    number that is Gaussian-distributed about the mean ⟨(m + m')/2⟩ with standard
    deviation σ = sqrt(τ/dt). ```m``` may be a function of time t that returns a QOp.
    """
    function readout(dt::Time, m::QOp, τ::Time)
        mo = (m .+ m')/2
        σ = sqrt(τ/dt)
        (t::Time, ρ) -> σ*randn() + real(expect(ρ, mo))
    end
    function readout(dt::Time, m::Function, τ::Time)
        σ = sqrt(τ/dt)
        (t::Time, ρ) -> let mo = (m(t) .+ m(t)')/2;
                            σ*randn() + real(expect(ρ, mo)) end
    end

    """
        gausskraus(dt::Time, m::QOp, τ::Time)

    Return a function that accepts a Gaussian random sample `r` and returns a
    Kraus operator for the associated state back-action. Note that `m` is an
    operator of the form `m = m_o - i m_p`, where `m_o` is the Hermitian
    observable being monitored, and `m_p` is the generator for phase-backaction.

    The generated Kraus operator has the following Gaussian form:

      - ``M_r = (dt/2πτ)^(1/4) e^{ -i m_p r dt/(2τ) - dt(r - m_o)^2/(4τ) }``

    However, the constant factors have been eliminated in anticipation of
    subsequent state renormalization, to simply yield:

      - ``M_r = e^{ m r dt/(2τ) - m_o^2 dt/(4τ) }``

    """
    function gausskraus(dt::Time, m::QOp, τ::Time)
        mo = (m .+ m') / 2
        mo2 = mo^2 / 2
        v = dt/(2*τ)
        (t::Time, r) -> SparseOperator(exp((r*v)*m - v*mo2))
    end
    function gausskraus(dt::Time, m::Function, τ::Time)
        v = dt/(2*τ)
        (t::Time, r) -> let mo = (m(t) .+ m(t)') / 2
                            mo2 = mo^2 / 2
                            SparseOperator(exp((r*v) * m(t) - v*mo2)) end
    end

    ###
    # Crude trajectory integrator
    ###

    # Return trajectory array [f1(now), f2(now), ...]
    """
        trajectory(inc::Function, init,
                        tspan::Tuple{Time,Time}, fs::Function...;
                        dt::Time=1/10^4, points::Int=1000, verbose=true, readout=true)

    Compute time-stepped trajectory, starting from state `init`, incrementing with `inc`
    by time step `dt` for `t` in `tspan`, keeping `points` intermediate values
    of `f(ρ(t))` for each `f` in `fs`. If the trajectory is stochastic, optionally store
    the simulated readout as well.

    ### Returns:
      - Array of Trajectories for each function (and readout) value:
        [Trajectory(f(ρ(t)))..., Trajectory(r(t))...]

    """
    @inline function trajectory(inc::Function, init,
                        tspan::Tuple{Time,Time}, fs::Function...;
                        dt::Time=1/10^4, points::Int=1000, verbose=true, readout=true)
        # Simulate and collect data
        s = simulate(inc, init, tspan, fs...,
                     dt=dt, points=points, verbose=verbose, readout=readout)
        # Return interpolated trajectory objects
        ts = range(first(tspan), length=length(first(s)), stop=last(tspan))
        map(v -> Trajectory{typeof(first(v))}(ts, v), s)
    end

    """
        ensemble(n::Integer, inc::Function, init,
                 tspan::Tuple{Time,Time}, fs::Function...;
                 dt::Time=1/10^4, points::Int=1000, verbose=true, readout=true)

    Compute time-stepped ensemble of `n` trajectories, starting from state `init`,
    incrementing with `inc` by time step `dt` for `t` in `tspan`, keeping `points`
    intermediate values of `f(ρ(t))` for each `f` in `fs`. If the trajectory is
    stochastic, optionally store the simulated readouts as well.

    ### Returns:
      - Array of Ensembles for each function (and readout) value:
        [Ensemble(f(ρ(t)))..., Ensemble(r(t))...]

    """
    @inline function ensemble(n::Integer, inc::Function, init,
                      tspan::Tuple{Time,Time}, fs::Function...;
                      dt::Time=1/10^4, points::Int=1000, verbose=true, readout=true)
        # Information if desired
        if verbose
            pinit = inc(first(tspan), init)
            N = Int(fld(abs(last(tspan)-first(tspan)), dt))
            @info("Trajectories: ",n,", steps each: ",N,
                 ", points each: ",min(N,points),", values each = ", length(fs))
            readout && typeof(pinit) <: Tuple &&
                @info("Readouts: values each = ", length(pinit)-1)
        end

        inittime=time()
        # Simulate and collect data
        # Uses pmap to optionally support parallelization
        fn = _ -> simulate(inc, init, tspan, fs..., dt=dt, points=points, verbose=false, readout=readout)
        s = hcat(pmap(fn, 1:n)...)
        elapsed = time() - inittime

        # Return ensemble of each stored value
        ts = range(first(tspan), length=length(first(s)), stop=last(tspan))
        out = collect(Ensemble(ts, hcat(s[i,:]...)) for i in 1:first(size(s)))

        # Info if desired
        if verbose
            steps = length(first(out))
            @info("Time elapsed: ",elapsed," s, Steps: ",steps,
                 ", Steps per second: ",steps/elapsed)
        end
        out
    end

    function simulate(inc::Function, init,
                        tspan::Tuple{Time,Time}, fs::Function...;
                        dt::Time=1/10^4, points::Int=1000, verbose=true, readout=true)
        # Constants for integration
        t0 = first(tspan)              # Initial time
        tmax = last(tspan)             # Final time
        ts = t0:dt:tmax                      # simulated time range
        N = length(ts)                 # total num of points
        if points > N
            Ns = N                     # reset points if needed
        else
            Ns = points                # stored points
        end
        Nf = length(fs)                # stored f values per point
        Nl = Int(cld(N-1, points))     # steps per stored point
        Nldt = Nl*dt                   # time-step per stored point

        # Preallocate trajectory arrays for speed
        traj = [let fi=f(init); Array{typeof(fi)}(undef,Ns) end for f in fs]
        # Functions to update values
        function update!(i::Int, ρ)
            for k in 1:Nf
                @inbounds traj[k][i] = fs[k](ρ)
            end
        end
        next(t::Time, ρ) = inc(t, ρ)

        # Modify preallocation and updates for stochastic readout
        # Trial point to test for readout
        pinit = inc(first(tspan), init)
        if readout && typeof(pinit) <: Tuple
            # Stored readout variables per point
            Nr = length(pinit) - 1
            # Preallocate readout arrays for speed
    	rtraj = [Array{typeof(r)}(undef,Ns) for r in pinit[2:end]]
            append!(traj, rtraj)
    	# Specialize updates to handle readout tuples
            function update!(i::Int, t::Tuple)
                update!(i, first(t))
                for k in 1:Nr
                    @inbounds traj[Nf+k][i] = t[k+1]
                end
            end
            next(t::Time, tup::Tuple) = inc(t, first(tup))
        else
            Nr = 0
        end

        # Seed loop
        verbose && @info("Trajectory: steps = ",N-1,", points = ",Ns,", values = ",Nf)
        verbose && readout && Nr > 0 && @info("Readout: values = ",Nr)
        inittime = time()
        now = init
        tnow = t0
        update!(1, now)
        # loop
        for i in 2:Ns
            # inner loop without storage
            for _ in 1:Nl
                tnow += dt
                now = next(tnow, now)
            end
            # store point
            update!(i, now)
        end
        elapsed = time() - inittime
        # Performance summary
        verbose && @info("Time elapsed: ",elapsed," s, Steps per second: ",(N-1)/elapsed)
        traj
    end

    export simulate, ensemble, trajectory, gausskraus, readout, meas, lind, lind_rk4, ham, sham, rinc, ham_rk4, slind, map, mean, median, std

end

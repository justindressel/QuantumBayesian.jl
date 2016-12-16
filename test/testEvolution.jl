# Test Hamiltonian increment
h = ham(1e-2, 0*q("z"))
@test_approx_eq(h(0.0, q("x")), q("x"))
h = ham(1e-2, t->0*q("z"))
@test_approx_eq(h(0.0, q("x")), q("x"))
h = ham(π/2, q("x"))
@test_approx_eq(h(0.0, ground(q)), q(1,1))
h = ham(1e-2, 0*q("z"), ket=true)
@test_approx_eq(h(0.0, xplus), xplus)
h = ham(π/2, q("x"), ket=true)
@test_approx_eq(h(0.0, groundvec(q)), -im*q(1))

# Test superoperator Hamiltonian increment
h = sham(1e-2, 0*q("z"))
@test_approx_eq(h(0.0, superket(q("x"))), superket(q("x")))
h = sham(1e-2, t->0*q("z"))
@test_approx_eq(h(0.0, superket(q("x"))), superket(q("x")))
h = sham(π/2, q("x"))
@test_approx_eq(h(0.0, superket(ground(q))), superket(q(1,1)))

# Test RK4 Hamiltonian increment
h = ham_rk4(1e-2, 0*q("z"))
@test_approx_eq(h(0.0, q("x")), q("x"))
h = ham_rk4(1e-2, t->0*q("z"))
@test_approx_eq(h(0.0, q("x")), q("x"))
h = ham_rk4(1e-2, 0*q("z"), ket=true)
@test_approx_eq(h(0.0, xplus), xplus)

# Test jump/no-jump Lindblad increment
l = lind(1e-6, 0*q("z"), q("d"))
@test_approx_eq(l(0.0, q(1,1))[2,2], exp(-1e-6))

# Test RK4 Lindblad increment
l = lind_rk4(1e-6, 0*q("z"), q("d"))
@test_approx_eq(l(0.0, q(1,1))[2,2], exp(-1e-6))

# Test Superoperator Lindblad increment
l = slind(1e-6, 0*q("z"), q("d"))
@test_approx_eq(unsuperket(l(0.0, superket(q(1,1))))[2,2], exp(-1e-6))
l = slind(1e-6, t->0*q("z"), q("d"))
@test_approx_eq(unsuperket(l(0.0, superket(q(1,1))))[2,2], exp(-1e-6))

# Test trajectory code
t = trajectory(ham(0.0, q("z")), ground(q), (0.0, 1.0), ρ->ρ[1,1], dt=1e-3, points=100, verbose=false)
@test t[1].t == linspace(0.0,1.0,100)
@test_approx_eq(t[1][end], QComp(1))
@test_approx_eq(t[1](1.0), QComp(1))
@test size(t[1]) == size(t[1].v)
t = trajectory(ham(1e-4, (2π/4).*q("y")), ground(q), (0.0, 1.0), ρ->ρ[2,2], dt=1e-4, points=1000, verbose=false)
@test t[1].t == linspace(0.0,1.0,1000)
@test_approx_eq_eps(t[1][end], QComp(1), 1e-4)
@test_approx_eq_eps(t[1](1.0), QComp(1), 1e-4)

# Test trajectory point truncation
t = trajectory(ham(0.0, q("z")), ground(q), (0.0, 1.0), ρ->ρ[1,1], dt=1e-2, points=1000, verbose=false)
@test t[1].t == linspace(0.0,1.0,99)
@test_approx_eq(t[1][end], QComp(1))

# Test stochastic trajectory code
t = trajectory(meas(1e-3, q("z"), [(q("z"), 2.0, 1.0)], q("d")), ground(q), (0.0, 1.0), ρ->ρ[1,1], dt=1e-3, points=100, verbose=false)
@test t[1].t == linspace(0.0,1.0,100)
@test length(t) == 2

# Test stochastic ensemble code
e = Ensemble(t[1])
@test e(1).t == linspace(0.0,1.0,100)
@test size(e(1).v) == size(t[1].v)
@test size(e) == size(e.a)
@test length(e) == length(e.a)
@test ndims(e) == ndims(e.a)
e = Ensemble(t[1].t, 1, eltype(t[1].v))
e[1] = t[1].v
@test e[1][1] == t[1].v[1]
e = Ensemble(t[1].t, t[1].v) 
@test e[1][1] == t[1].v[1]
e = ensemble(2, meas(1e-3, q("z"), [(q("z"), 2.0, 1.0)], q("d")), ground(q), (0.0, 1.0), ρ->real(ρ[1,1]), dt=1e-3, points=100, verbose=false)
@test e[1].n == 2
@test length(mean(e[1])) == length(e[1](1))
@test length(std(e[1])) == length(e[1](1))
@test length(median(e[1])) == length(e[1](1))

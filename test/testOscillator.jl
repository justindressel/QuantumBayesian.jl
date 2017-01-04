# Test coherent state
o2 = osc(15)
v = coherentvec(o2, 1)
@test_approx_eq(v ⋅ v, 1)
@test_approx_eq_eps(expectvec(v, o2("n")), 1, 1e-11)
@test_approx_eq(expectvec(v, o2("n")), weakvaluevec(v, v, o2("n")))
@test_approx_eq_eps(expectvec(v, o2("d")), weakvaluevec(v, groundvec(o2), o2("d")), 1e-11)
@test_approx_eq(weakvaluevec(v, groundvec(o2), o2("n")), 0)

vp = projector(v)
@test_approx_eq(vp, coherent(o2, 1))

@test_approx_eq(trace(vp * vp), 1)
@test_approx_eq_eps(expect(vp, o2("n")), 1, 1e-11)
@test_approx_eq(expect(vp, o2("n")), weakvalue(vp, vp, o2("n")))
@test_approx_eq_eps(expect(vp, o2("d")), weakvalue(vp, ground(o2), o2("d")), 1e-11)
@test_approx_eq(weakvalue(vp, ground(o2), o2("n")), 0)


# Test coherent state
o2 = osc(15)
v = coherentvec(o2, 1)
@test v ⋅ v ≈ 1
@test expectvec(v, o2("n")) ≈ 1 atol=1e-11
@test expectvec(v, o2("n")) ≈ weakvaluevec(v, v, o2("n"))
@test expectvec(v, o2("d")) ≈ weakvaluevec(v, groundvec(o2), o2("d")) atol=1e-11
@test weakvaluevec(v, groundvec(o2), o2("n")) ≈ 0

vp = projector(v)
@test vp ≈ coherent(o2, 1)

@test trace(vp * vp) ≈ 1
@test expect(vp, o2("n")) ≈ 1 atol=1e-11
@test expect(vp, o2("n")) ≈ weakvalue(vp, vp, o2("n"))
@test expect(vp, o2("d")) ≈ weakvalue(vp, ground(o2), o2("d")) atol=1e-11
@test weakvalue(vp, ground(o2), o2("n")) ≈ 0

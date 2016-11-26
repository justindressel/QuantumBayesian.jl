using QuantumBayesian
using Base.Test

# Fixture definitions for tests
q = qubit();
@test size(q) == 2
@test length(q) == 2
@test name(q) == "Qubit"
@test q("z") == q.ops["z"]

f = QFactor(3, "Foo");
@test size(f) == 3
@test length(f) == 3
@test name(f) == "Foo"

o = osc(5);
@test size(o) == 5
@test length(o) == 5
@test name(o) == "Osc(5)"

s = ⊗(q,q,f)
@test size(s) == (2,2,3)
@test length(s) == 12
@test name(s) == "Qubit ⊗ Qubit ⊗ Foo"
@test factors(s) == Vector([q,q,f])

t = q ⊗ o ⊗ f ⊗ q
@test size(t) == (2,5,3,2)
@test length(t) == 60
@test name(t) == "Qubit ⊗ Osc(5) ⊗ Foo ⊗ Qubit"
@test factors(t) == Vector([q,o,f,q])

# Test tensor product subview indexing
a = ["a11" "a12" ; "a21" "a22"];
b = ["b11" "b12" ; "b21" "b22"];
c = ["c11" "c12" "c13" ; "c21" "c22" "c23"; "c31" "c32" "c33"];
abc = ⊗(a,b,c);
abcv = subview(s, abc);
@test abc[15] ==  abcv[15]
@test abcv[2,1,3,1,2,1] == a[2,1]*b[1,2]*c[3,1]
@test abcv.perm == [3,2,1,6,5,4]

tv = subview(t, t("yxii"))
@test tv.data.parent == t("yxii")
@test unview(tv) == t("yxii")
@test tv[1,1,3,1,1,2,1,2] == t("yxii")[sub2ind(tv, 1,1,3,1,1,2,1,2)]
@test (1,1,3,1,1,2,1,2) == ind2sub(tv, sub2ind(tv, 1,1,3,1,1,2,1,2))
t("yxii")[sub2ind(tv, 1,1,3,1,1,2,1,2)] = 9
@test t("yxii")[sub2ind(tv, 1,1,3,1,1,2,1,2)] == 9
@test tv[1,1,3,1,1,2,1,2] == 9

# Test sparse basis generation
e = s[1,2,3,2,2,1]
e2 = s(0,1,2,1,1,0)
@test e == e2
ev = subview(s, e)
@test ev[1,2,3,2,2,1] == QComp(1)
@test e[sub2ind(ev, 1,2,3,2,2,1)] == QComp(1)
e = s[1,2,3]
e2 = s(0,1,2)
@test e == e2
ev = subview(s, e)
@test ev[1,2,3] == QComp(1)
@test e[sub2ind(ev, 1,2,3)] == QComp(1)

# Test partial traces
@test ptrace(s, s("xii")) == trace(s("xii"))
@test ptrace(s, subview(s, s("xii"))) == trace(s("xii"))
pt, pm = ptrace(s, s("xii")./3, 3)
pt2, pm2 = ptrace(pt, pm./2, 2)
@test pm2 == q("x")
@test ptrace(s, s("xii")./6, 2, 3) == (q, q("x"))
@test ptrace(s, s("ixi")./6, 1, 3) == (q, q("x"))
@test ptrace(s, s("ixi")./6, 3, 1) == (q, q("x"))

pt, pm = ptrace(s, s("iii")./2, 2)
pt2, pm2 = ptrace(pt, pm./2, 1)
@test pm2 == f("i")
@test ptrace(s, s("iii")./4, 1, 2) == (f, f("i"))

# Test qubit relations
@test comm(q("x"),q("y")) == 2*im*q("z")
@test comm(q("y"),q("z")) == 2*im*q("x")
@test comm(q("z"),q("x")) == 2*im*q("y")
@test acomm(q("x"),q("y")) == spzeros(q.dim,q.dim)
@test acomm(q("y"),q("z")) == spzeros(q.dim,q.dim)
@test acomm(q("z"),q("x")) == spzeros(q.dim,q.dim)
@test acomm(q("x"),q("x")) == 2*q("i")
@test acomm(q("y"),q("y")) == 2*q("i")
@test acomm(q("z"),q("z")) == 2*q("i")
@test acomm(q("x"),q("y")) == spzeros(q.dim,q.dim)
@test acomm(q("y"),q("z")) == spzeros(q.dim,q.dim)
@test acomm(q("z"),q("x")) == spzeros(q.dim,q.dim)
xplus = (q[1] + q[2])/sqrt(2)
xplus2 = (q(0) + q(1))/sqrt(2)
@test xplus == xplus2
yplus = (q[1] + im*q[2])/sqrt(2)
@test weakvaluevec(xplus, yplus, q("z")) == QComp(-im)
xplusp = projector(xplus)
yplusp = projector(yplus)
@test weakvalue(xplusp, yplusp, q("z")) == QComp(-im)

# Test dissipator definition
@test diss(q("x"))(q("z")) == q("x") * q("z") * q("x")' - (1/2)*(q("x")'*q("x")*q("z") + q("z")*q("x")'*q("x"))
@test diss(o("x"))(o("n")) == o("x") * o("n") * o("x")' - (1/2)*(o("x")'*o("x")*o("n") + o("n")*o("x")'*o("x"))
@test diss(o("y"))(o("n")) == o("y") * o("n") * o("y")' - (1/2)*(o("y")'*o("y")*o("n") + o("n")*o("y")'*o("y"))

# Test innovation definition
@test inn(q("x"))(q("z")) == q("x") * q("z") +  q("z") * q("x")' - trace((q("x") + q("x")')*q("z"))*q("z")
@test inn(o("x"))(o("n")) == o("x") * o("n") +  o("n") * o("x")' - trace((o("x") + o("x")')*o("n"))*o("n")
@test inn(o("y"))(o("n")) == o("y") * o("n") +  o("n") * o("y")' - trace((o("y") + o("y")')*o("n"))*o("n")

# Test coherent state
o2 = osc(15)
v = coherentvec(o2, 1)
@test abs((v ⋅ v) - 1) < 1e-6
@test abs(expectvec(v, o2("n")) - 1) < 1e-6
@test abs(expectvec(v, o2("n")) - weakvaluevec(v, v, o2("n"))) < 1e-6
@test abs(expectvec(v, o2("d")) - weakvaluevec(v, groundvec(o2), o2("d"))) < 1e-6
@test abs(weakvaluevec(v, groundvec(o2), o2("n"))) < 1e-6
vp = projector(v)
@test abs((vp ⋅ vp) - 1) < 1e-6
@test abs(expect(vp, o2("n")) - 1) < 1e-6
@test abs(expect(vp, o2("n")) - weakvalue(vp, vp, o2("n"))) < 1e-6
@test abs(expect(vp, o2("d")) - weakvalue(vp, ground(o2), o2("d"))) < 1e-6
@test abs(weakvalue(vp, ground(o2), o2("n"))) < 1e-6


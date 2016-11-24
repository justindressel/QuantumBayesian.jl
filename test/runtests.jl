using QuantumBayesian
using Base.Test

# Fixture definitions for tests
q = qubit();
f = QFactor(5, "Foo");
o = osc(3);
s = ⊗(q,q,o)
t = q ⊗ f ⊗ o ⊗ q

# Test tensor product subview indexing
a = ["a11" "a12" ; "a21" "a22"];
b = ["b11" "b12" ; "b21" "b22"];
c = ["c11" "c12" "c13" ; "c21" "c22" "c23"; "c31" "c32" "c33"];
abc = ⊗(a,b,c);

abcv = subview(s, abc);
@test abc[15] ==  abcv[15]
@test abcv[2,1,3,1,2,1] == a[2,1]*b[1,2]*c[3,1]

tv = subview(t, t("yixi"))
@test tv[1,1,3,1,1,2,1,2] == t("yixi")[sub2ind(tv, 1,1,3,1,1,2,1,2)]
@test (1,1,3,1,1,2,1,2) == ind2sub(tv, sub2ind(tv, 1,1,3,1,1,2,1,2))
t("yixi")[sub2ind(tv, 1,1,3,1,1,2,1,2)] = 9
@test tv[1,1,3,1,1,2,1,2] == 9

# Test sparse basis generation
e = s[1,2,3,2,2,1]
ev = subview(s, e)
@test ev[1,2,3,2,2,1] == QComp(1)
@test e[sub2ind(ev, 1,2,3,2,2,1)] == QComp(1)
e = s[1,2,3]
ev = subview(s, e)
@test ev[1,2,3] == QComp(1)
@test e[sub2ind(ev, 1,2,3)] == QComp(1)

# Test partial traces
pt, pm = ptrace(s, 3, s("xii")./3)
pt2, pm2 = ptrace(pt, 2, pm./2)
@test pm2 == q("x")

pt, pm = ptrace(s, 2, s("iin")./2)
pt2, pm2 = ptrace(pt, 1, pm./2)
@test pm2 == o("n")

# Test coherent state
o2 = osc(15)
v = coherentvec(o2, 1)
@test abs((v ⋅ v) - 1) < 1e-6
@test abs(expectvec(v, o2("n")) - 1) < 1e-6
vp = projector(v)
@test abs((vp ⋅ vp) - 1) < 1e-6
@test abs(expect(vp, o2("n")) - 1) < 1e-6


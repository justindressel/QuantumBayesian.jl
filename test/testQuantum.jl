###
# Test quantum space definitions
#
# Recall:
# q = qubit();
# f = QFactor(3, "Foo");
# o = osc(5);
# s = ⊗(q,q,f)
# t = (q ⊗ (o ⊗ f)) ⊗ q
###

@test size(q) == 2
@test length(q) == 2
@test name(q) == "Qubit"
@test q('z') == q.ops['z']

@test size(f) == 3
@test length(f) == 3
@test name(f) == "Foo"

@test size(o) == 5
@test length(o) == 5
@test name(o) == "Osc(5)"

@test size(s) == (2,2,3)
@test length(s) == 12
@test name(s) == "Qubit ⊗ Qubit ⊗ Foo"
@test factors(s) == Vector([q,q,f])

@test size(t) == (2,5,3,2)
@test length(t) == 60
@test name(t) == "Qubit ⊗ Osc(5) ⊗ Foo ⊗ Qubit"
@test factors(t) == Vector([q,o,f,q])

f2 = QSpace(f)
f3 = QSpace(3, "Foo")
@test f2.factors == factors(f2)
@test f3.factors == factors(f3)
@test f2.ops["i"] == f2("i")
@test f3.ops["i"] == f3("i")

#show(f2)
#rdout = readavailable(rd)
#@test String(rdout) == "QuantumBayesian.QSpace: Foo\nDims  : (3,)\nOps   : \"i\"\n"

###
# Test tensor product subviews
###

a = ["a11" "a12" ; "a21" "a22"];
b = ["b11" "b12" ; "b21" "b22"];
c = ["c11" "c12" "c13" ; "c21" "c22" "c23"; "c31" "c32" "c33"];
abc = ⊗(a,b,c);
abcv = subview(s, abc);
# Index tests
@test abcv[2,1,3,1,2,1] == a[2,1]*b[1,2]*c[3,1]
@test length(abcv) == length(abc)
@test_throws ErrorException subview(s,abcv)

tv = subview(t, t("yxii"))
@test tv[1,4,3,2,2,5,3,2] == 2.0im

###
# Test sparse basis generation
###

e = s[1,2,3,2,2,1]
e2 = s(0,1,2,1,1,0)
@test e == e2
ev = subview(s, e)
@test ev[1,2,3,2,2,1] == QComp(1)
e = s[1,2,3]
e2 = s(0,1,2)
@test e == e2
ev = subview(s, e)
@test ev[1,2,3] == QComp(1)
@test_throws ErrorException s[1,2,3,4]

###
# Test partial trace functionality
###

@test ptr(q, q('i')) == tr(q('i'))
@test ptr(q, subview(q, q('i'))) == tr(q('i'))
@test ptr(s, s("xii")) == tr(s("xii"))
pt, pm = ptr(s, s("xii") ./ 3, 3)
pt2, pm2 = ptr(pt, pm ./ 2, 2)
@test pm2 == q('x')
@test ptr(s, s("xii")./6, 2, 3) == (q, q('x'))
@test ptr(s, s("ixi")./6, 1, 3) == (q, q('x'))
@test ptr(s, s("ixi")./6, 3, 1) == (q, q('x'))

pt, pm = ptr(s, s("iii")./2, 2)
pt2, pm2 = ptr(pt, pm./2, 1)
@test pm2 == f('i')
@test ptr(s, s("iii")./4, 1, 2) == (f, f('i'))

###
# Test lift operation
###

@test lift(s, q('x'), 1) == s("xii")
@test lift(s, q('x'), 2) == s("ixi")
@test lift(s, f('i'), 3) == s("iii")

###
# Test superoperator
###

@test unsuperket(superopl(o('d')) * superket(o('n'))) ≈ o('d') * o('n')
@test unsuperket(superopr(o('d')) * superket(o('n'))) ≈ o('n') * o('d')
@test unsuperket(ssand(o('d')) * superket(o('n'))) ≈ o('d') * o('n') * o('d')'
@test unsuperket(scomm(o('d')) * superket(o('n'))) ≈ o('d') * o('n') - o('n') * o('d')'

# Test dissipator definition
@test_approx_eq(diss(q("x"))(q("z")), q("x") * q("z") * q("x")' - (1/2)*(q("x")'*q("x")*q("z") + q("z")*q("x")'*q("x")))
@test_approx_eq(diss(o("x"))(o("n")), o("x") * o("n") * o("x")' - (1/2)*(o("x")'*o("x")*o("n") + o("n")*o("x")'*o("x")))
@test_approx_eq(diss(o("y"))(o("n")), o("y") * o("n") * o("y")' - (1/2)*(o("y")'*o("y")*o("n") + o("n")*o("y")'*o("y")))
@test_approx_eq(unsuperket(sdiss(o("y")) * superket(o("n"))), o("y") * o("n") * o("y")' - (1/2)*(o("y")'*o("y")*o("n") + o("n")*o("y")'*o("y")))

# Test innovation definition
@test inn(q("x"))(q("z")) == q("x") * q("z") +  q("z") * q("x")' - trace((q("x") + q("x")')*q("z"))*q("z")
@test inn(o("x"))(o("n")) == o("x") * o("n") +  o("n") * o("x")' - trace((o("x") + o("x")')*o("n"))*o("n")
@test inn(o("y"))(o("n")) == o("y") * o("n") +  o("n") * o("y")' - trace((o("y") + o("y")')*o("n"))*o("n")

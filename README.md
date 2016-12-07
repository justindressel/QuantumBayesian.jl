# QuantumBayesian

[![Build Status](https://travis-ci.org/justindressel/QuantumBayesian.svg?branch=master)](https://travis-ci.org/justindressel/QuantumBayesian)
[![codecov](https://codecov.io/gh/justindressel/QuantumBayesian/branch/master/graph/badge.svg)](https://codecov.io/gh/justindressel/QuantumBayesian)


Julia package for simulating quantum evolution efficiently, including the quantum Bayesian approach to stochastic measurement update and continuous readout.

# Example
```julia
using QuantumBayesian
using PyPlot

q = qubit()
```

    QuantumBayesian.QFactor: Qubit
    Dims  : 2
    Ops   : "d", "y", "x", "u", "z", "i"

```julia
Ω  = 2*π;        # Rabi frequency
τ = 3.0;         # Measurement collapse timescale
η = 1.0;         # Measurement efficiency
Γ = 1/(2*τ*η);   # Measurement dephasing rate
T = (0.0, 6*τ);  # Time duration of simulation;

dt = 1e-2;       # Simulation timestep (coarse to show method precision)

# Initial condition
init = ground(q)

# Time-dependent Hamiltonian
f(t) = 2*exp(-(t-3*τ)^2/2)/sqrt(2π)
H(t) = f(t)*(Ω/2)*q("y");

# Measurement dephasing
DM = sqrt(Γ/2)*q("z");

# Bloch coordinate expectation values 
fs = collect(ρ -> real(expect(ρ, q(l))) for l in ["x","y","z"])

# Integrate trajectory
@time out = trajectory(lind(dt, H, DM), init, T, fs..., dt=dt)

# PyPlot Whatnot
p = plot(out[1],out[2],label=L"$x$")
hold(true)
plot(out[1],out[3],label=L"$y$")
plot(out[1],out[4],label=L"$z$")
ax = gca()
ax[:set_ylim]([-1.1,1.1])
xlabel(L"$t (2\pi/\Omega)$")
ylabel("Bloch coordinates")
title("Lindblad Evolution: Jump/No-Jump")
legend()
hold(false)
show()
```

    INFO: Trajectory: steps = 1799, points = 1000, values = 3
    INFO: Time elapsed: 0.086969712 s, Steps per second: 20685.362278766657

![Lindblad Pulse Output](img/example_lindblad_pulse.png)

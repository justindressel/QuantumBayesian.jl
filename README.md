# QuantumBayesian

Julia package for simulating quantum evolution efficiently, including the quantum Bayesian approach to stochastic measurement update and continuous readout.

# Examples
```julia
using QuantumBayesian
using PyPlot

# Define global constants for evolution
const ωq = 2*π*5;
const Ω  = 2*π;
const dt = 1e-3;
const T1 = 100.0;
const T2 = 40.0;

# Create a qubit space
q = qubit()
```

    QuantumBayesian.QFactor: Qubit
    Dims  : 2
    Ops   : "d", "y", "x", "u", "z", "i"

```julia
function run(T; p=1000)
    # Hamiltonian
    const H = (ωq/2)*q("z")+(Ω/2)*(cos(π/6)*q("x")+sin(π/6)*q("y"));
    # Energy relaxation
    const D1 = sqrt(1/T1)*q("d");
    # Environmental dephasing
    const D2 = sqrt(1/T2)*q("z");
    # Lindblad increment
    inc = lind(dt, H, D1, D2);
    # Bloch coordinate expectations
    x(ρ) = real(expect(ρ, q("x")));
    y(ρ) = real(expect(ρ, q("y")));
    z(ρ) = real(expect(ρ, q("z")));
    # Initial ground state
    init = ground(q);
    # Compute trajectory
    trajectory(init, inc, dt, T, x, y, z, points=p)
end;

out = run(T2, p=5000)
```

    INFO: Trajectory: steps = 39999.0, points = 5000, values = 3
    INFO: Time elapsed: 1.454858859 s, Steps per second: 27493.388621555627


    (linspace(0.0,40.0,5000),
    [-0.025972; -0.05897; … ; -0.187372; -0.18605],
    
    [0.0353558; 0.0641149; … ; -0.0932612; -0.0997137],
    
    [-0.999037; -0.996197; … ; -0.908548; -0.908122])

```julia
figure()
plot(out[1],out[2])
plot(out[1],out[3])
plot(out[1],out[4])
```

![Example Output](tmp/example_lindblad.png)

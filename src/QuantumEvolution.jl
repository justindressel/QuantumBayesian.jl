
###########################################################
# Simple propagators
##

# Hamiltonian propagation
function ham(dt,H)
  const u = sparse(expm( -im * dt * full(H))); const ut = u'
  ρ -> u * ρ * ut
end

# Jump-nojump Lindblad propagator
function lind(dt,H,alist...)
  h = ham(dt,H)
  const n = sparse(sqrtm(eye(H) - dt * full(mapreduce(a -> a' * a, +, alist))))
  no(ρ) = n * ρ * n
  dec(ρ) = mapreduce(a -> a * ρ * a', +, alist) * dt
  ρ -> let ρu = h(ρ); no(ρu) + dec(ρu) end
end

# Runge-Kutta Lindblad propagator
function lind_runge(dt,H,alist...)
  inc(ρ) = - im * comm(H,ρ) * dt + sum(map(a -> diss(a)(ρ) * dt, alist))
  function rinc(ρ)
    local dρ1 = inc(ρ)
    local dρ2 = inc(ρ + dρ1 / 2)
    local dρ3 = inc(ρ + dρ2 / 2)
    local dρ4 = inc(ρ + dρ3)
    return (dρ1 + 2*dρ2 + 2*dρ3 + dρ4)/6
  end
  ρ -> ρ + rinc(ρ)
end

###
# Crude trajectory integrator
###

# Return trajectory array [f1(now), f2(now), ...] for each increment dt up to tmax
function trajectory(init, inc::Function, dt, tmax, fs::Function...; points=1000)
  # N steps, Nf stored values per step
  # p100 points per percentage, Nl steps per point
  const N = Int(floor(tmax / dt))
  const Nf = length(fs)
  const p100 = Int(floor(points/100))
  const Nl = Int(floor(N/points))
  # Preallocate trajectory arrays for speed
  traj = zeros((points+1,Nf)); ts = zeros((points+1,1))
  info("Trajectory: steps = $(N), points = $(points), values = $(Nf)")
  tic()
  now = init; tnow = 0
  for i in 1:100
    #@printf("\b\b\b\b\b\b\b\b\b\b\b\bpercent: %2d", i-1)
    offset = (i-1)*p100
    for j in 1:p100
      traj[j + offset, :] = [f(now) for f in fs]
      ts[j + offset] = tnow 
      for k in 1:Nl
        now = inc(now)
      end
      tnow += Nl*dt
    end
  end
  traj[points+1,:] = [f(now) for f in fs]
  ts[points+1] = tnow; 
  elapsed = toq()
  # Performance summary
  #print("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
  info("Time elapsed: $(elapsed) s, Steps per second: $(N / elapsed)")
  ts, traj
end

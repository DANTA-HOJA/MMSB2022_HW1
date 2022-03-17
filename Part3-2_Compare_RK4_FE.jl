# The ODE model. Exponential decay in this example
# The input/output format is compatible to Julia DiffEq ecosystem AND YOU SHOULD KEEP IT.
model(u, p, t) = cos(t)

# Forward Euler stepper 
step_euler(model, u, p, t, dt) = u .+ dt .* model(u, p, t)

# RK4 stepper 
function step_rk4(model, u, p, t, dt) # operator with prefix "." = element-wise
    k1 = dt .* model(u, p, t)
    k2 = dt .* model(u.+0.5*k1, p, t.+0.5*dt)
    k3 = dt .* model(u.+0.5*k2, p, t.+0.5*dt)
    k4 = dt .* model(u.+k3, p, t.+dt)
    return u .+ (k1 + 2*k2 + 2*k3 + k4)/6
end

# In house ODE solver
function mysolve(model, u0, tspan, p; dt=0.1, stepper)
    # Time points
    ts = tspan[1]:dt:tspan[end]
    # State variable at those time points
    us = zeros(length(ts), length(u0))
    # Initial conditions
    us[1, :] .= u0
    # Iterations
    for i in 1:length(ts)-1
        us[i+1, :] .= stepper(model, us[i, :], p, ts[i], dt)
    end
    # Results
    return (t = ts, u = us)
end

tspan = (0.0, 4.0π)
p = nothing
u0 = 0.0
dt_FE = 0.01
dt_RK4 = 0.1 # <===== 達到相同精度時，"RK4" 與相比 "FE", time_step大約可以好10倍

# Visualization
using Plots
Plots.gr(lw=2)

p1 = plot(title = "Part3-2", size=(600,400)) # 開一個空的圖，assign to p1

# Numerical solution: FE
sol_FE = mysolve(model, u0, tspan, p; dt=dt_FE, stepper=step_euler)
plot!(p1, sol_FE.t, sol_FE.u, label="FE method, dt=$dt_FE") # plot! = 疊在前一個圖上面

# Numerical solution: RK4
sol_RK4 = mysolve(model, u0, tspan, p; dt=dt_RK4, stepper=step_rk4)
plot!(p1, sol_RK4.t, sol_RK4.u, label="RK4 method, dt=$dt_RK4") # plot! = 疊在前一個圖上面

# tspan... = tspan[1], tspan[2], ....., tspan[n]
plot!(p1, sin, tspan..., label = "Analytical solution", linestyle=:dash, linecolor=:black) 
plot(p1) # show p1
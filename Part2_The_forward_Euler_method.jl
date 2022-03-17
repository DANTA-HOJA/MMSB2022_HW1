# The ODE model. Exponential decay in this example
# The input/output format is compatible to Julia DiffEq ecosystem AND YOU SHOULD KEEP IT.
model(u, p, t) = cos(t)

# Forward Euler stepper 
step_euler(model, u, p, t, dt) = u .+ dt .* model(u, p, t)

# In house ODE solver
function mysolve(model, u0, tspan, p; dt=0.1, stepper=step_euler)
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

# Visualization
using Plots
Plots.gr(lw=2)

dt_list = [0.01, 0.1, 0.2, 0.5, 0.7, 1]


p1 = plot(title = "Part2", size=(600,400)) # 開一個空的圖，assign to p1
for dt in dt_list
    sol = mysolve(model, u0, tspan, p; dt=dt, stepper=step_euler) # Numerical solution
    plot!(p1, sol.t, sol.u, label="FE method, dt=$dt") # plot! = 疊在前一個圖上面
end
# tspan... = tspan[1], tspan[2], ....., tspan[n]
plot!(p1, sin, tspan..., label = "Analytical solution", linestyle=:dash, linecolor=:black) 
plot(p1) # show p1
using Plots, DifferentialEquations
Plots.gr(lw=2) # lw = "line width" 

tspan = (0.0, 4.0π) # time range
p = nothing # parameter variable
u0 = 0.0 # state variable

model(u, p, t) = cos(t)

# Define a problem
prob = ODEProblem(model, u0, tspan, p)

# Solve the problem
sol = solve(prob)

# Visualize the solution
p1 = plot(title = "Part1", size=(600,400)) # 開一個空的圖，assign to p1
plot!(p1, sol, label="FE method") # Numerical solution
plot!(p1, sin, tspan..., label = "Analytical solution", linestyle=:dash) # plot! = 疊在前一個圖上面
plot(p1) # show p1
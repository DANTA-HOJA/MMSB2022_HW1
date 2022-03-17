# MMSB2022_HW1
110.2 BEBI5009_生物系統模擬

HW1 : https://ntumitolab.github.io/mmsb-bebi-5009/hw/hw-01.html


# Part 1: Julia’s ODE solver

> Please solve the ODE using DifferentialEquations.jl for  and plot the time series. Compare it to the analytical solution in one plot.

!["Forward Euler" compare to "Analytical solution"](/png/Part1_Julia_ODE_solver.png)

# Part 2: The forward Euler method

> Please try a range of dts to solve the ODE using the (home-made) forward Euler method for t = [0 , 4pi], plot the time series, and compare them to the analytical solution in one plot.

!["Forward Euler" in different dt"](/png/Part2_The_forward_Euler_method.png)

# Part 3: The RK4 method
1. Please try a range of dts to solve the ODE using the (home-made) fourth order Runge-Kutta (RK4) method for t = [0 , 4pi], plot the time series, and compare them to the analytical solution in one plot.

!["Forward Euler" in different dt"](/png/Part3-1_The_RK4_method.png)

2. Compared to the forward Eular method, which one is more efficient in terms of time step needed for the same accuracy? You could make a visual comparison by plotting the analytical and numerical solutions together.

![Compare efficient between "Forward Euler" and "Runge-Kutta (RK4)"](/png/Part3-2_Compare_RK4_FE.png)

> 達到相同精度時，"RK4" 與相比 "FE", "RK4" 的 time_step 大約可以好10倍


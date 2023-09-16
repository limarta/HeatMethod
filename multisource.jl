using DifferentialEquations
using Plots

function initialize_grid(dx,dy)
    grid = zeros(length(dx),length(dy))
    grid[1, 1] = 300.0
    grid[end, end] = 300.0
    grid
end

function laplacian(du,u,p,t)
    alpha, dx = p
    N,_ = size(du)
    limit(a) = (a == 0) ? 1 : (a == N+1 ? N : a)
    for I in CartesianIndices(du)
        x = I[1]
        y = I[2]
        dx1, dx2, dy1, dy2 = limit(x-1), limit(x+1), limit(y-1), limit(y+1)
        du[I] = alpha * (u[dx1, y] + u[dx2,y] + u[x,dy1] + u[x,dy2] - 4*u[I]) / dx^2
    end
end

N = 21
dx = range(-1,1,N)
dy = range(-1,1,N)
u0 = initialize_grid(dx,dy)
problem = ODEProblem(laplacian, u0,(0, 10.0), (10, step(dx)))
soln = solve(problem)
# heatmap(x,y,sol(t_viz), aspect_ratio=:equal, framestyle=:box, ticks=false, xlims=(-1,1), ylims=(-1,1), background_color_subplot=false, clims=cfunc)
begin
    heatmap(dx, dy, soln(0.1), aspect_ratio=:equal, xlims=(-1,1), ylims=(-1,1),legend=false)
    xx
    fig1 = scatter!([0.9, -0.2, -0.9],[0.9, -0.2, 0.9], markersize=5, color=:green, markerstrokewidth=0)
    fig2 = plot(bar(["a","b","c"], [1,2,3],  labels=["a","b","c"]), ylabel="Heat")
    plot(fig1, fig2, layout = (1, 2), legend = false)
end
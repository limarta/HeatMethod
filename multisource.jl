using DifferentialEquations
using Plots

function initialize_grid(dx,dy)
    grid = zeros(length(dx),length(dy))
    N = ceil(Int,length(dx)//2)
    grid[N, N] = 1.0
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

function varadhan_theorem(t, heat)
    sqrt.(-t * log.(heat)/4)
end

N = 21
dx = range(-1,1,N)
dy = range(-1,1,N)
u0 = initialize_grid(dx,dy)
problem = ODEProblem(laplacian, u0,(0, 1.0), (10, step(dx)))
soln = solve(problem
t = 0.01
dists = varadhan_theorem(t,soln(t))
minimum(dists)
begin
    fig1 = heatmap(dx, dy, soln(0.01), aspect_ratio=:equal, xlims=(-1,1), ylims=(-1,1),legend=false)
    fig2 = heatmap(dx, dy, dists, aspect_ratio=:equal, xlims=(-1,1), ylims=(-1,1),legend=false)
    plot(fig1, fig2, layout = (1, 2), legend = false)
end
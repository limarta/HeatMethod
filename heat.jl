using DifferentialEquations
using Plots

N = 20
dxy = range(-1,1, N)
x = dxy
y = dxy
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

function init_heat_map(dxy)
    N = length(dxy)
    u = zeros(N,N)
    mid = ceil(Int, N/2)
    u[mid,mid] =1.0
    u
end

u0 = init_heat_map(dxy)
p = (0.01, step(dxy))
prob = ODEProblem(laplacian, u0, (0.0, 20.0), p)
sol = solve(prob)
t_range = range(0, 20, 1000)
cfunc(x) = (0.0, maximum(x)+sqrt(100/(200-maximum(x)/2)))
anim = @animate for t in t_range
    heatmap(x,y,sol(t), aspect_ratio=:equal, framestyle=:box, ticks=false, xlims=(-1,1), ylims=(-1,1), background_color_subplot=false, clims=cfunc)
end every 5

gif(anim, "heat.gif", fps=20)


# animate(sol, every=4)
# p1 = heatmap(x,y,sol[2], aspect_ratio=:equal, framestyle=:box, ticks=false, xlims=(-1,1), ylims=(-1,1), ackground_color_subplot=false)
# plot(p1)
# savefig(p1, "posts/p2/heat.svg")

using Plots


using Eikonal
tsize = 1000
solver = FastMarching(tsize, tsize)
solver.v .= 1;

npoints = 4
for _ in 1:npoints
    (i, j) = rand(1:tsize, 2)
    init!(solver, (i, j))
end

march!(solver, verbose=true)

contour(1:(tsize+1), 1:(tsize+1), solver.t, levels=30,
        aspect_ratio=1, c=:coolwarm, size=(800, 600), ticks=false, framestyle=:box,
        title = "W", xlims=(1,(tsize+1)), ylims=(1,(tsize+1)))
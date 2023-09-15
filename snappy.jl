using WGLMakie, JSServe, Colors

using JSServe.DOM

begin
	markersize = JSServe.Slider(LinRange(1, 20, 100))
	App() do session
		return DOM.div(markersize, markersize.value)
	end
end

begin
	App() do session
        hue_slider = JSServe.Slider(LinRange(1, 360, 100))
        color = map(hue_slider) do hue
            HSV(hue, 0.5, 0.5)
        end
        positions = rand(Point3f, 10^6)
        fig = scatter(positions, markersize=markersize, color=color)
		return DOM.div(hue_slider, fig)
	end
end

App() do session
    slider = JSServe.Slider(5:10)
    color = map(hue_slider) do hue
		HSV(hue, 0.5, 0.5)
	end
    # n = slider.value.val
    # fig = Figure()
    
    # ax = Axis(fig[1,1], aspect=1.0)
    # x = range(0,1,n+1)
    # y = range(0,1,n+1)
    # xlims!(ax, (0,1))
    # ylims!(ax, (0,1))
    # hlines!(ax, y, color=:black, linewidth=7)
    # vlines!(ax, x, color=:black, linewidth=7)
    
    # scatter!(ax, 0, 0, marker=:circle, markersize=50, color=:red)
    # scatter!(ax, 1, 1, marker=:circle, markersize=50, color=:red)
    # xx = repeat(x, inner=2)[2:end]
    # yy = repeat(y, inner=2)[1:end-1]
    # lines!(ax, xx,yy, linewidth=10, color=:red)
    # hidedecorations!(ax)
    positions = rand(Point3f, 10^3)
    fig = scatter(positions, color=color)
    return DOM.div(slider, fig)
end

function fib(n::Int64)::Int64
    n <= 2 && return 1
    n == 3 && return 2
    fib(n - 1) + fib(n - 2)
end

@time fib(50);
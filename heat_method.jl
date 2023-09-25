### A Pluto.jl notebook ###
# v0.19.27

#> [frontmatter]
#> title = "Heat Method for Distance Computation"
#> date = "2023-09-01"
#> tags = ["blog", " "]
#> description = "Demonstration of the Heat Method"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ dab8582c-e8e2-443f-b52c-ac716b2ca12e
using PlutoUI

# ╔═╡ 358c1832-06c0-4918-b766-8c1de98c21d3
begin
	using DifferentialEquations , Plots, Eikonal, LinearAlgebra, SparseArrays
end

# ╔═╡ d79f5864-5193-4673-a593-1057ec15e927
begin
	using PlyIO
	using Arpack
	using ColorSchemes
	using WGLMakie
	using JSServe
	Page()
	WGLMakie.JSServe.SERVER_CONFIGURATION.listen_port[] = 8089
end


# ╔═╡ 7dc7f1eb-3e1f-4b75-8e4c-e3018b0259d6
# Note to self: What is Gridap.jl?

# ╔═╡ edce4ccf-7be7-4b0d-98aa-c572eac3d0ad
md"""
# Introduction

This notebook provides an overview of the Heat Method by Crane et al. used for geodesic computation [^1]. First we describe limitations with standard discretizations and shortest path algorithms and show the intuition behind the Heat Method. Then we provide some discretization details relevant to meshes and finally provide a Julia implementation. Originally this notebook was going to use pure Julia, but I think there is some value in highlighting some of the flagship packages in the ecosystem.
"""

# ╔═╡ efa0e0ba-8b30-4f69-9bc8-bdd89ca5f61a
PlutoUI.TableOfContents(depth=5)

# ╔═╡ 94ca9000-947d-44d2-8a6d-f9177c476345
md"""### What are Geodesics?

What is the shortest path between points $A$ and $B$? It's likely you have heard this question posed in the context of algorithms, and invariably is applicable to so many domains such as planning. For those of you in CS, you were most likely peppered with questions on finding shortest paths for graphs. The graph model turns out to be a perfect fit to solve shortest paths on because of the discrete nature of a graph. It is easy to store representations of vertices and edges as well as the paths along the graph by simply tabulating a collection of discrete objects. But can we ask the same question in other domains besides graphs, ones which are in fact continuous? 

Take a look below at several surfaces where you can imagine a small ant walking on. In all of them, intuitively the shortest path from $A$ and $B$ is the one denoted by the red curve. Sometimes it is simple to compute this shortest path. The flat sheet on the left is Euclidean sspace, and it is well known that the shortest path between any two points is the line crossing the two. What about the sphere on the right? Well the shortest paths there lie on the [great circle] (https://en.wikipedia.org/wiki/Great_circle) of the sphere (think equator). Like on graphs, there may be several shortest paths connecting the points. 
"""

# ╔═╡ ad43bc4b-3cc4-4960-bb46-eb6978620e61
let
# Sphere and bunny?
end

# ╔═╡ c235d925-3899-4cf7-9d9d-847ab81a7523
md"""
Now we see the problem is more general than graphs. When dealing with surfaces, the term [geodesic](https://en.wikipedia.org/wiki/Geodesic) has been coined to refer to shortest paths and is a central object in geometry. The question now is, given a surface and two points, what is the geodesic(s) connecting the two?
"""

# ╔═╡ ba6a8f2e-3e01-48c9-a07a-72fcde78e6f1
md"""
### Geodesics on Discretized Surfaces

There has been a concerted effort in the last two centuries to study geodesics on continuous, well-differentiable surfaces known as manifolds. While an interesting subject in its own right, there is not much focus on surfaces that arise naturally from computationally heavy domains such as graphics. By definition, our questions will be discrete in nature.

Often our surfaces can be represented as meshes or point clouds with a finite number of vertices and possibly edges connecting the vertices.
"""

# ╔═╡ 008abc67-4a20-4ba7-a7e0-a1922fd67387
# Add images of polygonal meshes + point clouds

# ╔═╡ 8c2aaf03-a48f-4241-9e78-3b7312ee71e2
md"""
### Issues with Discretization

Shortest path algorithms such as Dijkstra's seem to be the perfect fit for finding geodesics, right? At first this might seem true, but we have to remember that these algorithms answer a fundamentally different question - what is the shortest path along the edges of a graph? This means that all paths produced will necessarily follow the edges specified by the mesh. But why is this bad? 

Let us consider again the square surface above where each side is unit length. Clearly the shortest path from $A$ to $B$ is the diagonal with length $\sqrt{2}$. Now how would a shortest path algorithm do here? 

We need a graph to pass to our solver. So what if we simply use the four vertices and four sides of the square? Unfortunately we immediately run into problems - no matter how good the solver is, the optimal path from $A$ to $B$ is purportedly along the sides.

"""

# ╔═╡ 5a53c282-c0be-47b3-964f-2fedbfa25451
let
	x = [0, 1, 1, 0, 0]
	y = [0, 0, 1, 1, 0]
	highlighted_side_x = [0, 1]
	highlighted_side_y = [1, 1]
	# Create the plot
	Plots.plot(x, y, lw = 2, linecolor = :black, legend = false, xlims=(0,1), ylims=(0,1), size = (200,200), aspect_ratio=:equal, framestyle=:box, ticks=false, background_color=:white, margin=5*Plots.mm)
	Plots.plot!(highlighted_side_x, highlighted_side_y, seriestype = :shape, lw = 10, linecolor = :red)
	annotate!(-0.05,-0.05, "A", fontsize=5)
	annotate!(1.05,1.05, "B", fontsize=5)
	Plots.plot!([0,0], [0,1], lw=10, linecolor=:red)
end

# ╔═╡ a73d4efc-b9a7-4171-9bdc-c98e907dd2b7
md"""
Clearly our choice of the graph does not adequately describe the possible paths of a *continuous* surface, so what now? In many fields such as PDEs, numerical simulations commonly discretize continuous spaces into discrete ones, so we can try to use the same methods here. 

What should the discretization be though? It turns out there is not obvious choice of discretizing the space. The figure below shows one possible strategy that cuts the square into smaller triangles. 
"""

# ╔═╡ 6f2cdb0d-c00d-48ca-a2d0-ea462532895d
md"""Refine N: $(@bind refinement_n  PlutoUI.Slider(2:2:16,8,true))"""

# ╔═╡ 97460bf6-ce4f-4210-9664-8c6c43a9a382
let
	dx = range(0,1,refinement_n+1)
	Plots.plot()
	Plots.xlims!((0,1))
	Plots.ylims!((0,1))
	Plots.plot!(ticks=false,  framestyle=:box, size=(300, 300), legend=false, margin=5*Plots.mm)
	vline!(dx, line = :black, lw = 1)
	hline!(dx, line = :black, lw=1)

	xx = repeat(dx, inner=2)[2:end]
	yy = repeat(dx, inner=2)[1:end-1]
	Plots.plot!(xx,yy, color=:red, lw=3)
	for x in dx
		Plots.plot!([0,x], [x,0], color=:black)
	end
	for x in reverse(dx)
		Plots.plot!([x,1], [1,x], color=:black)
	end
	Plots.plot!([0,1], [1,0], color=:blue, lw=3)
	annotate!(-0.05, 1.05, "A", color = :black, fontsize = 5)
	annotate!(-0.05, -.05, "B", color = :black, fontsize = 5)
	annotate!(1.05, 1.05, "C", color = :black, fontsize = 5)
	annotate!(1.05, -0.05, "D", color = :black, fontsize = 5)
end

# ╔═╡ 0526ff5b-952c-4cba-842b-1059c67c12e1
md"""
Like before, the partitioned square can be represented as a graph and running a shortest path algorithm will produce the two paths shown in red and blue. This strategy successfully recovers the geodesic from $A$ to $D$ because it just so happens that the graph has the edges to describe this path. Despite our initial successes, this choice leads to the path in red. A careful look shows that this path *still* has length $2$ despite refining the graph! Further discretization will yield **no** improvement for the path from $B$ to $C$.

Adding an edge from $B$ to $C$ of length $\sqrt{2}$ can fix the problem for this particular path, but this fails to address the more general case where the start and ending points lie in arbitrary places on the square. Of course, there is always the option to add many irregular edges to the graph to express more paths, but it will never be enough to precisely answer every query. The downsides are that more edges means longer runtimes and possibly larger memory consumption to represent all of the visited points. 

Finally even if we were willing to take the hit in runtime and memory performance, this all assumes that we can compute the weights of the edges. In a square it is straightforward but asking the same question on more complicated surfaces becomes just as difficult as finding the geodesics themselves.  
"""

# ╔═╡ f0eb3973-f9c4-41fc-8f38-3bcb71e76c7d
let
	θ = range(0,2,100)
	c1 = Plots.Shape([(0.9i,0.9j) for (i,j) in zip(cospi.(θ), sinpi.(θ))])
	c2 = Plots.Shape([(0.5*i,0.5*j) for (i,j) in zip(cospi.(θ), sinpi.(θ))])
	T = 1.0
	Plots.plot(ticks=false,  framestyle=:none, size=(200, 200), legend=false, xlims=(-T,T), ylims=(-T,T))
	Plots.plot!(c1, color=:lightgray)
	Plots.plot!(c2, color=:white)
	Plots.scatter!([0.6], [0.1], markersize=5, markerstrokewidth=0,colormap=:red)
	Plots.scatter!([-0.7], [-0.4], markersize=5, markerstrokewidth=0, colormap=:red)
end

# ╔═╡ c8fc76f7-f900-460b-a6e3-a33a3386a8e0
Markdown.MD(Markdown.Admonition("warning", "Question", [md"""
Is there good discretization for this surface? Even if one could be found, would it work well with shortest path algorithms
"""]))

# ╔═╡ 19696716-d7ae-4ffd-8b73-399e5f02831b
md"""
# Intuition by Leveraging Heat

How do we use heat to produce a notion of distance? For now we will simplify the picture by working on square grids and later Let's continue with the running example of a metal square surface that is initially $0^{\circ}$F everywhere. Now let us take a red-hot needle and touch the center of the plate. At $t=0$ only the center of the plate is hot - say $100^\circ$ F. After a while, the heat flows away from the concentrated center to the boundaries, and this process is governed according to a diffusion equation known as the [heat equation](https://en.wikipedia.org/wiki/Heat_equation) [^2]. Time elapsed means that more of the plate becomes warmer while the center slowly cools down.
"""

# ╔═╡ 0350b1d6-4245-4b86-8b45-be9c00a16c77
md"""t =  $(@bind t_kernel PlutoUI.Slider(range(0,5,101), show_value=true, default=0.7))"""

# ╔═╡ ac43bbab-3014-4ece-b738-157e6367f734
md"""
Notably, the amount of heat at a given point on the plate is tied to how long it took for heat to diffuse from the center. Consider a point on the edge of the plate. We see that the temperature is consistently colder than the center, indicating that the point is quite far. Great. We have a notion of distance, but where do the geodesics come into play? By following where the heat flows away (e.g. decreasing temperatures), we can approximate the geodesics by finding paths that monotonically decrease in temperature. This is exactly what is happening when moving away from the center to the edge of the metal plate. 

The next example shows a more complicated setup with *two* source points in a insulated system. In an insulated system, the heat remains on the plate. In the multisource problem, the distance of the geodesics is the shortest path from any of the source points. The plot on the right shows heat measurements for three different locations. Can you map each dot to the correct bar?
"""

# ╔═╡ 3a23522e-d7d6-4461-bd33-878e1c823ea6
md"""
t = $(@bind t_multi PlutoUI.Slider(range(0,1,101), show_value=true, default=0.02))
"""

# ╔═╡ 7564b349-5d51-44a0-b78a-654cd1bbfb61
md"""
There are two interesting insights from the bar plot. The first is confirmation of our previous observation - the furthest point from the sources is consistently colder than points closer to them. The second is that the quality of the measurements is time sensitive. For large $t$, the plate has reached pretty close to equilibrium at around ~$0.04$, and so it is tougher to say which points are truly farther or whether it is a fluke with how the heat disperses. However, for $t$ very close to $0$, the difference in measurements is striking and better reflects the range of geodesic lengths. Notice that $t=0$ is uninformative, forcing us to use $t>0$.

It seems that to get high quality measurements, we need to operate in a middle ground. On one extreme, we must avoid running the simulation for too long or else the surface will reach an equilibrium. On the other extreme, we need to run the simulation for some small non-zero time to begin diffusion. The next section formalizes this intuition and presents our first algorithm to compute geodesics."""

# ╔═╡ 96fdedc9-3adf-4e46-9a50-d0b38bd38caf
md"""
### Varadhan's Approximation

The relationship between heat flow and distance is captured by an elegant expression known as Varadhan's formula.

$$\lim\limits_{t\rightarrow 0^+} t \log h(t,x,y) = -4d(x,y)^2$$

Here $d(x,y)$ is the geodesic distance between points $x$ and $y$. The function $h(t,x,y)$ is known as the *heat kernel*. The heat kernel describes how heat evolves on the surface if we initially placed one unit of heat at $x$. Then $h(t,x,y)$ refers to the amount of heat at $y$ at time $t$ when one unit of heat was placed at $x$. For example, Fig ? shows a heat kernel where $x=(0,0)$.

This formula gives us a direct way to estimate $d(x,y)$. First we place a unit of heat at at the source $x$ and simulate for a short period of time. Then we measure the temperature at $y$ to find the value of $h(t,x,y)$. For small $t$, we can do a substitution to obtain $d(x,y)$ to get a decent approximation.

First let's write the formula.

"""

# ╔═╡ c89ce545-5229-418e-a174-e2e4eddc1115
varadhan_formula(t, heat) = sqrt.(-4t * log.(heat))

# ╔═╡ 2dad5bf2-ac9c-4767-ac37-3abd252f338a
md"""
Below we have defined $\texttt{heat\_solution}$ to be the heat kernel with initial unit heat at pixel $(7,9)$. Evaluating $\texttt{heat\_solution}$ at a time $t$ gives a matrix representing $h(t,(7,9),y)$. The variable $\texttt{t\_varadhan}$ corresponds is the value of the slider.
"""

# ╔═╡ c814f0d2-cbb0-4db9-9499-993d51f42356
@bind t_varadhan PlutoUI.Slider(range(0,0.1, 101), show_value=true, default=0.005)

# ╔═╡ e1c19aea-d671-4c01-8bcf-119e7abb295f
md"""
Not great but not bad. We do see that the estimates get moderately better for smaller $t$ though the relative error is still significant near the source. Why do these artifacts appear? First, we cannot push $t$ arbitrary close to $0$. For extremely small values of $t$, the solver's numerical accuracy degrades (move the slider to the left to see numerical errors). This is not just a limitation of the differential equation solver. To simulate heat diffusion, the square was already discretized into an $N\times N$ grid with a single unit of heat placed at the cell $(7,9)$. ...
"""

# ╔═╡ 54c232ba-1175-40a3-b5a9-729450905e9f
md"""
### Eikonal Equation
Well Varadhan's formula was good first attempt, but unfortunately the fundamental inaccuracities due to numerical and discretization error puts a hard limit on this approach. Instead of refining results from Varadhan's formula, the paper instead takes an alternative approach. From Section ??? we discussed that geodesics correspond to find paths that cool down the further down you travel. Crane et al. leverage this observation by relating it to the Eikonal equation: 

$\vert \nabla u(x)\vert=1$ for $x\in M$ with boundary condition $u(s)=0$ for where $s$ is the source point. Here $\mathcal{M}$ represents the surface (e.g. square plate).


The Eikonal equation is a PDE used to describe shortest paths from a source point, and the solution $u$ represents the distance function from $s$. The boundary condition states that $u(x)$ must start at $0$ for the source point. More generally, the boundary condition can be expanded to include multiple source points or even a section of space. Below are some examples of solutions to the Eikonal equation with various sources [^3].
"""

# ╔═╡ dd03796a-0520-418a-88f3-f11547b05a19
let
	tsize = 1000
	solver = FastMarching(tsize, tsize)
	solver.v .= 1;
	
	npoints = 4
	for _ in 1:npoints
	    (i, j) = rand(1:tsize, 2)
	    init!(solver, (i, j))
	end
	
	march!(solver)
	
	Plots.contour(1:(tsize+1), 1:(tsize+1), solver.t, levels=30,
	        aspect_ratio=1, c=:coolwarm, size=(300, 300), ticks=false, framestyle=:box,
	        title = "W", xlims=(1,(tsize+1)), ylims=(1,(tsize+1)))
end

# ╔═╡ 217ac117-5415-4938-a543-8ebe5cca7898
md"""
Why does $u(x)$ describe the distance function? A helpful analogy is to think how spilt water moves across a surface. Tiny water particles flow away from the source, forming a puddle. The particle at the boundary, or wavefront, consist of particles that have traveled the furthest from the source and had to have traveled at the same rate to get to where they are at now. Thus the boundary represents all the particles that are equidistant from the source, in other words, the boundary is an *level set* mapping $u$ to the same distance. We can imagine letting this puddle continue growing and the wavefront expanding away. This of course represents $u$ increasing. 
"""

# ╔═╡ 17ca6c6a-d282-457d-961d-40275a01927a
md"""
So how does the analogy connect with the constraint $\vert \nabla u(x)\vert = 1$? Recall that $\nabla u(x)$ is the gradient and represents the direction of the largest change in distance at $x$. But since $u(x)$ is a distance function, it changes the most by increasing the distance it self, which means expanding the wavefront away from the origin as indicated by the arrows above. The speed of the wavefront is governed by the magnitude of $\nabla u(x)$. Notice that there is something funny going on with the units of $\nabla u(x)$. It basically represents the rate of change of the distance function with respect to the change in distance. Aha! That's where the $\nabla u(x)=1$ comes from. All it is trying to say is that "distance changes $1$ meter per meter"[^4]. Almost tautological...

So instead of evaluating Varadhan's formula, we evaluate the Eikonal equation and problem solved? Sorry to dissappoint again, but unfortunately the Eikonal equation is a bit too difficult to solve accurately, even by simulation. This isn't to say that people don't try to solve it, they do and have succeeded in getting good estimates. Fast marching is a whole class of approximation algorithms that takes the wave front analogy to the extreme. The problem is that many of these algorithms need to approximate the wave front which can be some what fickly. There are exact solvers, but even these incur large runtime costs, almost on the order of $O(n^3)$. Once again, the paper takes another turn to save the day.
"""

# ╔═╡ c912af64-147f-4ca5-a567-45f5c5e50303

md"""
# The Heat Method

Instead of solving the Eikonal equation directly, the Heat Method repurposes heat diffusion to solve it in three stages: simulating heat, computing gradients, and solving a Poisson equation. It makes use of the following observation on the solution to the heat equation.
"""


# ╔═╡ 3e920b70-d5b6-44f4-b257-e7568a458173
Markdown.MD(Markdown.Admonition("info", "Observation #3", [md"""
Let $h(x)$ be the heat kernel for source $s$ and small time $t$. Then $h(x)$ decreases approximately in the same directions as the distance function $u(x)$.
"""]))

# ╔═╡ d7ccc528-a811-4c31-8d64-fa2ce1e813cb
let
	# Picture of two vector fields. One for heat, one for eikonal.
end

# ╔═╡ e56e6e62-5f38-467d-83bd-daaf4a968044
md"""
Observation #3 says that $\nabla h(y)$ is roughly parallel to $\nabla u(x)$. But since $\vert\nabla u(x)\vert =1$, normalizing the vector field $\nabla h(x)$ gives a good estimate for $\nabla u(x)$. Although we have $\nabla u(x)$ instead of $u(x)$, this is a good first step so let's write out.
"""

# ╔═╡ 47e4da81-2c79-44bc-8d87-a9a0f2e45d4e
md"""
### Simulating Heat

On a flat surface, the heat flows according to the PDE [^5]

$$\Delta x = x_t$$

The operator $\Delta$ is the *Laplacian* and it usually takes the form $\Delta=\frac{\partial}{\partial x^2}+\frac{\partial}{\partial y^2}$ for $2$D surfaces. In the continuous setting, it is straightforward to compute the Laplacian of a function. In our setting, however, functions defined on the discretized structures (e.g. the grid or meshes) are not inherently continuous and instead are defined on a finite number of points. For example, the square grid defines functions on each of the $N\times N$ lattice points and can be represented as a $\mathbb{R}^{N\times N}$ vector. More generally discrete functions are $\mathbb{R}^{\vert V\vert}$ vectors where $\vert V \vert$ is the number of evaluable points. To adapt $\Delta$ for discrete functions, observe that $\Delta$ is a linear operator that maps scalar functions to new scalar functions. This means that finite version of $\Delta$ is simply a matrix mapping $\mathbb{R}^{\vert V\vert}$ to $\mathbb{R}^{\vert V\vert}$. Call this proposed matrix $L$.

Now defining the entries of $L$ is a whole ordeal on its own because it greatly depends on discretization and the domain (e.g. meshes, graphs, point clouds etc). Naturally, discretizing $\Delta$ means that $L$ loses some of the deseriable properties of $\Delta$. Lots of effort has been put in finding good and robust Laplacians for different domains, but luckily for us the square grid has a fairly simple $L$ that is intuitive and works great.
"""

# ╔═╡ 47b00e38-83d1-4888-baee-662bd716827c
md"""
##### Laplacian for a Grid

Recall that the square grid defines discrete functions for $N\times N$ points. Let $f\in \mathbb{R}^{N\times N}$ be one such function. The goal is to approxmimate 
$$\Delta=\frac{\partial }{\partial x^2}+\frac{\partial}{\partial y^2}$$, so we first discretize $$\frac{\partial}{\partial x^2}$$. Intuitively, the second derivative at a lattice point will depend both neighbors on the $x$ axis. Expanding out the Taylor series confirms this:

$\begin{align}
f(a+dx,b)&=f(a,b) + \frac{\partial}{\partial x}f(x,y)dx+\frac{1}{2}\frac{\partial}{\partial x^2}(dx)^2+o((dx)^2)\\
f(a-dx, b) &= f(a,b) - \frac{\partial}{\partial x} f(x,y) dx + \frac{1}{2}\frac{\partial}{\partial x^2} (dx)^2 + o((dx)^2)\\
\end{align}$
Summing both expressions yields

$\begin{align}
f(a+dx,b)+f(a-dx,b) = 2f(a,b) + \frac{\partial}{\partial x^2} f(x,y) (dx)^2 + o((dx)^4)
\end{align}$

or 

$$\frac{\partial}{\partial x^2} f(x,y) \approx\frac{[f(a+dx,b)-f(a,b)] + [f(a-dx,b) - f(a,b)]}{(dx)^2}$$

The same can be said for the derivative with respect to $y$:

$$\frac{\partial}{\partial y^2} f(x,y) \approx\frac{[f(x,b+dy)-f(a,b)] + [f(a,b-dy) - f(a,y)]}{(dy)^2}$$

Combining both derivatives yields our discrete Laplacian. Instead of using real-valued position $(a,b)$ we opt with indices representing the entries of a matrix. $L$ is then defined as

$\begin{align}
L_{ij} = &\frac{[f(i+1,j)-f(i,j)] + [f(i-1,j) - f(i,j)]}{(dx)^2}\\
&+\frac{[f(i,j+1)-f(i,j)]+[f(i,j-1)-f(i,j)]}{(dy)^2}\\
&=\frac{\sum_{k\in N(i,j)}[f(k)-f(i,j)]}{(dx)^2}
\end{align}$

where the last equality occurs when $dx=dy$. This is usually called the "five-point stencil" Laplacian since it requires a cross shape to evaluate.

To make sense of the boundary cells, the boundary conditions are used. There are two common ones, Dirichlet and Neumann. For Dirichlet, values at the boundary are fixed to $0$ (e.g. $f(s) = 0$) whereas for Neumann boundary conditions, the partial derivatives are fixed (e.g. $f_x(s) = c$).

First we define the Laplacian with Dirichlet boundary conditions.[^6]
"""

# ╔═╡ 387f2124-1802-405f-8d6e-3cfdcefe2f46
function laplacian_dirichlet(Δu,u,p,t)
	α, Δx,Δy = p
    N1, N2 = size(u)
	Δx² = Δx^2
	Δy² = Δy^2
    for j in 2:(N2-1)
        for i in 2:(N1-1)
            Δu[i, j] = 
				(u[i+1, j] + u[i-1, j] -2u[i,j]) / Δx² + 
				(u[i, j+1] + u[i, j-1] -2u[i,j])/ Δy²
        end
    end

	# Dirichlet condition enforced by dropping neighbors.
	
    # left/right edges
    for i in 2:(N1 - 1)
        Δu[i, 1] = (u[i+1, 1] + u[i-1, 1] -2u[i,1])/Δx² + (u[i, 2] - 2u[i, 1])/Δy²
        Δu[i, N2] = (u[i+1, N2] + u[i-1, N2] -2u[i,N2])/Δx² + (u[i, N2-1] - 2u[i, N2])/Δy²
    end

    # top/bottom edges
    for j in 2:(N2-1)
        Δu[1, j] = (u[1, j+1] + u[1, j-1] - 2u[1,j])/Δx² + (u[2, j] - 2u[1, j])/Δy²
        Δu[N1, j] = (u[N1, j+1] + u[N1, j-1] -2u[N1,j])/Δx² + (u[N1-1, j] - 2u[N1, j])/Δy²
    end

    # corners
    Δu[1, 1] = (u[2, 1]-2u[1,1])/Δx² + (u[1, 2]-2u[1, 1])/Δy²
    Δu[N1, 1] = (u[N1-1, 1]-2u[N1,1])/Δx² + (u[N1, 2]-2u[N1, 1])/Δy²
    Δu[1, N2] = (u[2, N2] -2u[1,N2])/Δx²+ (u[1, N2 - 1]-2u[1,N2])/Δy²
    Δu[N1, N2] = (u[N1 - 1, N2]-u[N1,N2])/Δx²+ (u[N1, N2-1]-u[N1, N2])/Δy²
	Δu .*= α
	Δu
end

# ╔═╡ 3ce304a9-3dfd-42f1-b3cf-308b6ce3cbac
md"""
The Laplacian with Neumann conditions is similar except now...
"""

# ╔═╡ ff09b1b2-3990-4bc7-8f6b-962e2ecfee3d
function laplacian_neumann(Δu, u, p, t)
	α, Δx, Δy = p
	Δx² = Δx^2
	Δy² = Δy^2
    n1, n2 = size(u)

    # internal nodes
    for j in 2:(n2 - 1)
        for i in 2:(n1 - 1)
            Δu[i, j] = (u[i+1, j]+u[i-1, j]-2u[i,j])/Δx² + (u[i,j+1] + u[i,j-1] - 2u[i, j])/Δy²
        end
    end

    # left/right edges
    for i in 2:(n1 - 1)
        Δu[i, 1] = (u[i+1, 1] + u[i-1, 1]-2u[i,1])/Δx² + (u[i, 2] - u[i, 1])/Δy²
        Δu[i, n2] = (u[i+1, n2]+u[i-1, n2]-2u[i,n2])/Δx² + (u[i, n2 - 1]-u[i, n2])/Δy²
    end

    # top/bottom edges
    for j in 2:(n2 - 1)
    	Δu[1, j] = (u[1, j+1]+u[1, j-1]-2u[1,j])/Δy² + (u[2, j]-u[1, j])/Δx²
        Δu[n1, j] = (u[n1, j+1]+u[n1, j-1]-2u[n1,j])/Δy² + (u[n1-1, j]-u[n1, j])/Δx²
    end

    # corners
    Δu[1, 1] = (u[2, 1]-u[1,1])/Δx² + (u[1, 2]-u[1, 1])/Δy²
    Δu[n1, 1] = (u[n1-1, 1]-u[n1,1])/Δx² + (u[n1, 2]-u[n1, 1])/Δy²
    Δu[1, n2] = (u[2, n2]-u[1,n2])/Δx² + (u[1, n2-1]-u[1, n2])/Δy²
    Δu[n1, n2]= (u[n1-1, n2]-u[n1,n2])/Δx² + (u[n1, n2 - 1]-u[n1, n2])/Δy²
	Δu *= α
end

# ╔═╡ 0bdff4e8-27ad-46ef-b9b3-a146d774cc6d
let
	N = 31
	u0 = zeros(N,N)
	u0[1,1] = 1.0
	u0[end,end] = 1.0
	L = 1.0
	dx = range(0,1,N)
	dy = range(0,1,N)
	p = (1, step(dx), step(dy))
	problem = ODEProblem(laplacian_neumann, u0,(0, 5.0), p)
	soln = solve(problem)

	Plots.heatmap(dx, dy, soln(t_multi), aspect_ratio=:equal, xlims=(0,L), ylims=(0,L))
	ii = [20, 9, 2]
	jj = [20, 9, 20]
	fig1 = Plots.scatter!(dx[ii],dy[jj], markersize=5, color=:green, markerstrokewidth=0, legend=false)
	measurements = [soln(t_multi)[i,j] for (i,j) in zip(ii,jj)]
	fig2 = Plots.plot(Plots.bar(["A","B","C"], measurements), ylabel="Heat", margin=5*Plots.mm)
	layout = @layout([a b{0.2w}])
    Plots.plot(fig1, fig2, layout = layout, legend = false)
	
end

# ╔═╡ ca26a33f-17c3-4a91-b95c-75b83409705e
md"""
To illustrate how DiffEq.jl works, let's write a function computing the heat kernel. Let's give the option as well to use either Dirichlet or Neumann boundary conditions.
"""

# ╔═╡ cffa4dc5-c9a8-4277-b87b-5dd3a5eff858
function heat_kernel(N, x, L=1.0, t=1.0, bc=:dirichlet)
	u0 = zeros(N,N)
	u0[x...] =1.0 # Initialize
	dxy = range(0, L, N)
	x = range(0,L,N)
	y = x
	p = (1.0, step(x), step(y))
	if bc == :dirichlet
		prob = ODEProblem(laplacian_dirichlet, u0, (0.0, t), p)
	elseif bc == :neumann
		prob = ODEProblem(laplacian_neumann, u0, (0.0, t), p)
	else
		error("BC not implemented")
	end
	return solve(prob)
end

# ╔═╡ e7080c15-ac7e-4106-8df4-65a668e39b83
let
	N = 31
	mid = ceil(Int, N/2)
	x = y = range(0,1,N)
	L = 1.0
	heat_sol = heat_kernel(N, (mid, mid), L, 5.0)
	cfunc(x) = (0.0, max(maximum(x), 0.001))
	Plots.heatmap(x,y,heat_sol(t_kernel), aspect_ratio=:equal, framestyle=:box, ticks=false, xlims=(0,L), ylims=(0,L), background_color_subplot=false, clims=cfunc, size=(400,400))
end

# ╔═╡ 2e4e7cf6-0034-4992-9ea0-f212b4111fc1
let
	N = 31
	L = 0.5
	mid = (7,9)
	heat_solution = heat_kernel(N, (7,9), L, 5.0)
	heat = heat_solution(t_varadhan)
	dist_estimates = varadhan_formula(t_varadhan, heat)

	# Compute true distances
	d(x,y) = sqrt((x-mid[1])^2 + (y-mid[2])^2)/N
	dist_true = d.((1:N)', 1:N)
	
	relative_error = abs.(dist_true - dist_estimates) ./ dist_true
	
	fig1 = Plots.heatmap(1:N,1:N, heat, aspect_ratio=:equal, framestyle=:box, ticks=false, xlims=(1,N), ylims=(1,N), background_color_subplot=false, title="Heat")
	fig2 = Plots.heatmap(1:N,1:N,dist_estimates, aspect_ratio=:equal, framestyle=:box, ticks=false, xlims=(1,N), ylims=(1,N), background_color_subplot=false, colormap=:viridis, title="Estimated Distance")
	fig3 = Plots.heatmap(1:N,1:N,dist_true, aspect_ratio=:equal, framestyle=:box, ticks=false, xlims=(1,N), ylims=(1,N), background_color_subplot=false, colormap=:viridis, title="True Distance")
	fig4 = Plots.heatmap(1:N,1:N,relative_error, aspect_ratio=:equal, framestyle=:box, ticks=false, xlims=(1,N), ylims=(1,N), background_color_subplot=false, colormap=:viridis, title="Relative Error")
	Plots.plot(fig1, fig2, fig3, fig4, label=(1,2,3,4))
end

# ╔═╡ 2a77d78c-1b8a-4690-9024-46a6794d8efd
md"""
### Computing Gradients
The next step is to take the gradient of the $h$. The discrete gradient is defined similarly.

$\begin{align}
\frac{\partial}{\partial x}h(a+dx,b) = h(a,b)+\frac{\partial}{\partial x}h(x,y)dx+o((dx)^2)\\
\frac{\partial}{\partial x}h(a-dx,b) = h(a,b)-\frac{\partial}{\partial x}h(x,y)dx+o((dx)^2)
\end{align}$

Combining terms yields

$$\frac{\partial}{\partial x}h(x,y) \approx \frac{h(a+dx,b)-h(a-dx,b)}{2dx}$$

Generally, the gradient field is a vector field and there are several equivalent ways of representing it for discrete structures. Here, the vector field is defined as a $\mathbb{R}^{\vert V\vert\times 2}$ matrix.
"""

# ╔═╡ 7a7a67a7-7958-43bb-bf54-36b80ecdf1de
function ∇(h::Matrix)
	N1, N2 = size(h)
	grad = zeros(N1+1, N2+1,2)
	for j in 1:N2
		for i in 1:N1
			if j == 1
				grad[i,j,1] = h[i,j]
			else
				grad[i,j,1] = h[i,j] - h[i,j-1]
			end

			if i == 1
				grad[i,j,2] = h[i,j]
			else
				grad[i,j,2] = h[i,j] - h[i-1,j]
			end
		end
	end
	for i in 1:N1
		grad[i,N2+1,1] = -h[i,N2]
	end
	for j in 1:N2
		grad[N1+1,j,2] = -h[N1,j]
	end
	return grad
end

# ╔═╡ b89464c0-baaa-4571-a744-91d5480c6ae1
function normalize_field(grad)
	Z = sqrt.(sum(grad.^2, dims=3))
	Z[Z .== 0] .= 1
	grad ./ Z
end

# ╔═╡ 160902db-e0b9-4d48-8fcb-41dbeaf429e8
begin
	N=17
	h = heat_kernel(N, (16,16), 1.0, 2.0, :dirichlet)
	nothing
end

# ╔═╡ 2415e55d-bed0-4050-893e-b6a7af00ef45
@bind t_grad PlutoUI.Slider(range(1e-6,2,101), show_value=true, default=0.5)

# ╔═╡ c0c94bdd-a7d0-46c3-a61e-6c0d40d8a3c9
begin
	heat_grid = h(t_grad)
	Plots.heatmap(1:N, 1:N, heat_grid, xlims=(1,N), ylims=(1,N),aspect_ratio=:equal, ticks=false)
	
	meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
	heat_grad_grid = ∇(heat_grid)
	heat_grad_normalized = normalize_field(heat_grad_grid)
	dd = (1:(N+1))
	Z =  2
	x, y = meshgrid((1:(N+1)).-0.5, (1:(N+1)))
	tx = vec(heat_grad_normalized[:,:,1]')/Z
	Plots.quiver!(x,y, quiver=(tx, zeros(size(tx))))
	ty = vec(heat_grad_normalized[:,:,2]')/Z
	x, y = meshgrid((1:(N+1)), (1:(N+1)).-0.5)
	Plots.quiver!(x,y, quiver=(zeros(size(ty)), ty))
end

# ╔═╡ 346c06a8-c91f-4c56-9840-67f2f02bd8c9
function div_grid(grad, dx, dy)
	N1, N2, _ = size(grad) .- (1,1,0)
	div = zeros(N1, N2)
	for j=1:N2
		for i=1:N1
			div[i,j] = (grad[i,j+1,1]-grad[i,j,1]) /(dx^2) +(grad[i+1,j,2]-grad[i,j,2])/(dy^2)
		end
	end
	div
end

# ╔═╡ 79d22285-1c69-495c-b05d-83d8c758ee46
function l(N,M, Δx, Δy)
	idx(i,j) = (j-1)*M + i
	Δx² = Δx^2
	Δy² = Δy^2
	Δ = spzeros(N*M, N*M)
	for j=2:(M-1)
		for i=2:(N-1)
			ii = idx(i,j)
			Δ[ii,ii] = -2/Δx² -2/Δy²
			Δ[ii,idx(i-1,j)] = 1/Δx²
			Δ[ii,idx(i+1,j)] = 1/Δx²
			Δ[ii,idx(i,j-1)] = 1/Δy²
			Δ[ii,idx(i,j+1)] = 1/Δy²
		end
	end
	for i in 2:(N - 1)
		ii = idx(i,1)
		jj = idx(i,M)
		Δ[jj, jj] = Δ[ii, ii] = -2/Δx²-2/Δy²

		Δ[ii,idx(i-1,1)] = Δ[ii,idx(i+1,1)] = Δ[jj, idx(i-1,M)] = Δ[jj, idx(i+1,M)] = 1/Δy²
		Δ[ii,idx(i,2)] = 1/Δx²
		Δ[jj,idx(i,M-1)] = 1/Δx²
	end
	for j in 2:(M-1)
		ii = idx(1,j)
		jj = idx(N,j)
		Δ[jj, jj] = Δ[ii, ii] = -2/Δx²-2/Δy²

		Δ[ii,idx(1,j-1)] = Δ[ii,idx(1,j+1)] = Δ[jj, idx(N,j-1)] = Δ[jj, idx(N,j+1)] = 1/Δx²
		Δ[ii,idx(2,j)] = 1/Δy²
		Δ[jj,idx(N-1,j)] = 1/Δy²
	end

    # corners
	ii = idx(1,1)
	Δ[ii,ii] = -2/Δx²-2/Δy²
	Δ[ii,idx(1,2)] = 1/Δx²
	Δ[ii,idx(2,1)] = 1/Δy²

	ii = idx(1,M)
	Δ[ii,ii] = -2/Δx²-2/Δy²
	Δ[ii,idx(1,M-1)] = 1/Δx²
	Δ[ii,idx(2,M)] = 1/Δy²

	ii = idx(N,1)
	Δ[ii,ii] = -2/Δx²-2/Δy²
	Δ[ii,idx(N,2)] = 1/Δx²
	Δ[ii,idx(N-1,1)] = 1/Δy²

	ii = idx(N,M)
	Δ[ii,ii] = -2/Δx²-2/Δy²
	Δ[ii,idx(N,M-1)] = 1/Δx²
	Δ[ii,idx(N-1,M)] = 1/Δy²
	Δ
end

# ╔═╡ a701e06d-553e-43b6-b36a-c68667bfd4b1
begin
	s = step(range(0,1,N))
	L_test = l(N,N, s, s)
	φ = div_grid(heat_grad_normalized, s,s)
	dist_grid = L_test \ vec(φ)

	Plots.heatmap(1:N, 1:N, reshape(heat_grid, N,N), xlims=(1,N), ylims=(1,N),aspect_ratio=:equal, ticks=false)
	# findall(x->x>0, φ)
end

# ╔═╡ 7c9b744e-72cd-449e-8228-a25b5c845233
md"""
### Solving the Poisson Equation

The final step is to transform the estimate of $\nabla u(x)$, $\vec{X}$, to an estimate of $u(x)$. To do so, we appeal to a general definition of the Laplacian using grad ($\nabla$) and div $(\nabla \cdot)$. The Laplacian is defined as

$$\Delta f= \nabla\cdot (\nabla f)$$

or "div grad of f". Verify that this definition is consistent when the Euclidean $\nabla$ and $\nabla \cdot$ are plugged in.

The right hand side of $\vec{X}=\nabla u(x)$ looks somewhat like the right hand side of the Laplacian. So applying $\nabla \cdot$ yields

$$\nabla\cdot \vec{X} = \nabla\cdot(\nabla u(x)) = \Delta u(x)$$

So solving this equation will yield $u(x)$, but is it any easier? In general no, but going back to the discrete setting simplifies this issue. Using our discrete operators so far, this amounts to solving

$$\nabla \cdot X = Lu$$

Notice that $\nabla\cdot X\in \mathbb{R}^{\vert V\vert}$, so as soon as we find a way to compute the "discrete div", the only thing left to do is to solve for $u$ as a system of linear equations!
"""

# ╔═╡ c3d64ff4-9d3e-44da-93f1-827b93042fc9
md"""
# Heat Method on Meshes
"""

# ╔═╡ 1bc67b2f-a1e9-4121-92cd-083b4ea9567b
Markdown.MD(Markdown.Admonition("danger", "", [md"""
**Warning**: I have tried my best but was unable to display the plots in this section for static viewing (what you are probably using right now). It seems that Pluto+WGLMakie is currently a bit brittle for static html files but hopefully this will be fixed soon. To my dismay, you will need to download and run Julia locally (or on Binder though this can take forever). For consolation, I left screenshots for those of you viewing online and a little suprise for the astute reader.
"""]))

# ╔═╡ 8a65dd84-2098-4765-8912-4ed6d32a9e0a
md"""
Up to this point, the examples were only of square plates, but the most interesting problems arise in 3D with meshes and point clouds. Luckily the Heat Method carries over to this regime. One thing we need to look out for is how how to do analogous calculations as those done one the square plate. For this, we need to propose discretizations for the operators.
"""

# ╔═╡ e48d51a3-debc-4339-84fe-20ee9613e808
let
	# Feature some meshes
end

# ╔═╡ 9a0aa371-9fbf-493f-ba4e-cb0801c2d5ef
md"""
### Discretizing ∇, ∇⋅, and Δ
"""

# ╔═╡ 861049ba-49d9-4f8c-a186-f1f95b282904
md"""

##### Laplace-Beltrami
The first order of business is to discretize the Laplacian $\Delta$ for meshes. You might recall that the Laplacian in $\mathbb{R}^3$ is defined as 

$$\Delta=\frac{\partial}{\partial x^2}+\frac{\partial}{\partial y^2}+\frac{\partial}{\partial z^2}$$

The problem is that this assumes that heat flows in a volume rather than a surface. To illustrate the difference, we plotted a sphere and a ball. The sphere is hollow and so the heat can only flow on its surface where as the ball conducts heat internally. 
"""

# ╔═╡ 83f7c9c7-cf37-44f0-8e51-ea6596d83605
let
	#simulation
end

# ╔═╡ 688712c4-b57f-49de-a1e9-3a3299eef60e
md"""
If the volumetric definition is used, then the heat values measured would correspond to geodesics that run *internal* to the mesh. Thus we need to use a different formulation of the Laplacian that handles surfaces of meshes well. For this, we introduce the [Laplace-Beltrami](https://en.wikipedia.org/wiki/Laplace–Beltrami_operator) operator. 

The Laplace-Beltrami (referred also as the Laplacian for short) operator is the 2D analog for surfaces embedded $\mathbb{R}^3$. Why is that? When considering meshes, they are often assumed to be manifolds, meaning that at any given point, the surface looks locally like a flat plane. So at a local scale, heat diffusion looks approximately the same as if we were simulating a flat surface. The generalization is defined as

How would the discrete analog look for a mesh? A mesh consists of vertices and faces, and a function can be defined over its vertices or faces. Here, a heat function $u(v)$ is defined over the vertices $v\in \mathcal{M}$ which can be represented compactly as a $\mathbb{R}^{\vert V\vert}$ vector. So the Laplacian must take in a scalar field $f(v)\in \mathbb{R}^{\vert V\vert}$ and produce a new scalar field over the vertices. This implies that the shape of $\Delta$ is represented as a matrix, more precisely a $|V|\times |V|$ matrix. Like [convolution filters](https://en.wikipedia.org/wiki/Kernel_(image_processing)#Details), there are many proposed Laplacian matrices that seek to approximate the continuous Laplace-Beltrami operator. It turns out that no discrete Laplacian is able to satisfy all of the properties of the continuous one, so each matrix has its own pitfalls. Here we use a common one called the **cotangent Laplacian**.

"""

# ╔═╡ 92b945aa-b19f-4732-963b-a3e8b42a6b02
Markdown.MD(Markdown.Admonition("info", "Note", [md"""
By far the most popular discretization for the Laplace-Beltrami operator is the cotangent Laplacian. You will see this choice in virtually most geometry processing papers as their first choice. It's main weakness is it is not robust even for rigid transformations and is sensitive to mesh changes. Moreover, there is not guarantee that the entries are all positive.
"""]))

# ╔═╡ d5b9c28b-e9fc-4e04-bf05-5b0b93da804e
md"""
The cotangent Laplacian, $L$, is defined as

$L_{ij}=\begin{cases}
w_{ij} & \text{if }i\neq j\\
-\sum\limits_{k\neq i} w_{ik} & \text{if }i=j
\end{cases}$
where $w_{ij}=\cot \alpha_i + \cot \alpha_j$. The figure below defines $\alpha_i$ and $\alpha_j$ for a pair of opposing vertices.
"""

# ╔═╡ 5bcf708c-8dc2-4a2d-a284-c17db2ea8b9a
md"""
This is actusally referred to as the *unweighted* version. The weighted version incorporates the vertex-based areas to weight the entries of $L$ [^7]. Vertices that amass a larger area need to be weighted down whereas vertices that cover little area are upweighted. Concretely, we define a mass matrix $A=\text{diag}(a)$ where $a\in\mathbb{R}^{\vert V\vert}$ are the vertex-based areas. The weighted cotangent Laplacian is then $A^{−1}L$.
"""

# ╔═╡ 3ed68b8f-db59-40fe-87ef-8df03f81f9df
md"""
##### Div and Grad
With the Laplacian defined, we now turn to finding manifold analogs of $\nabla$ and $\nabla\cdot$ while keeping faithful to Equation ?. It turns out that if we want to use the cotangent Laplacian for our choice of $\Delta$, then we need to be careful when defining $\nabla$ and $\nabla\cdot$. These operators need to satisfy the relationship $\nabla\cdot (\nabla u) = \Delta u$ (show yourself why it is true for Euclidean case).

Here we define a face-based gradient operator which takes a scalar field $\mathbb{R}^|V|$ and produces vectors for each face of the mesh. The vectors lie in the plane of the faces by design. For every triangle, the gradient is
$$(\nabla u)_f = \frac{1}{2A_f}\sum\limits_{j=1}^3 u_j (\textbf{N} \times \textbf{e}_j)$$
I'm slightly abusing the notation - the sum is over the three vertices of the triangle. $u_j$ are the scalar values of the vertices and $\textbf{e}_j$ is the edge opposite to the vertex. See the Figure below.

Explain what the definition of gradient is.

With $\nabla$ defined, the correspondencing choice of $\nabla\cdot$ should satisfy the relationship. The paper defines it as
$$(\nabla\cdot X)_v = \frac{1}{2}\sum\limits_{j} \cot\theta_1(e_2\cdot X_j) + \cot \theta_2(e_1\cdot X_j)$$
where $X$ is the vector field and $j$ is a sum over all the adjacent triangles to a vertex $v$. The angles $\theta_1$ and $\theta_2$ are the angles opposite of vertex $v$ on the triangle. The edges $e_1$ and $e_2$ eminate from $v$ and end at the opposite vertices. Take note of the ordering and see Figure .
"""

# ╔═╡ 3300a2ee-bc41-439a-8ec4-a981aab32a93
multicross(x,y) = reduce(hcat, cross.(eachcol(x), eachcol(y)))

# ╔═╡ 61f3e733-6527-43b9-97bd-08459e0878fc
vdot(x,y; dims=1) = sum(x .* y, dims=dims)

# ╔═╡ 9dd097b6-5f82-4fbe-b0d1-86756b7747d2
function cotlaplacian(V,F)
    nv = size(V,2)
    nf = size(F,2)
    #For all 3 shifts of the roles of triangle vertices
    #to compute different cotangent weights
    cots = zeros(nf, 3)
    for perm in [(1,2,3), (2,3,1), (3,1,2)]
        i, j, k = perm
        u = V[:,F[i,:]] - V[:, F[k,:]]
        v = V[:, F[j,:]] - V[:, F[k,:]]
        cotAlpha = vec(vdot(u,v; dims=1)) ./ norm.(eachcol(multicross(u,v)))
        cots[:,i] = cotAlpha
    end

    I = F[1,:]; J = F[2,:]; K = F[3,:];

    L = sparse([I;J;K], [J;K;I], [cots[:,1];cots[:,2];cots[:,3]], nv, nv)
    L = L + L'
    rowsums = vec(sum(L,dims=2))
    L = spdiagm(0 => rowsums) - L
    return 0.5 * L
end

# ╔═╡ f7914d06-c58e-4033-b895-9c069ec6eb4e
function face_area_normals(V,F)
    T = V[:,F]
    u = T[:,2,:] - T[:,1,:]
    v = T[:,3,:] - T[:,1,:]
    return 0.5 * multicross(u,v)
end

# ╔═╡ dc7eece7-942a-491a-9eaa-033c19112d32
function face_normals(V,F)
    A = face_area_normals(V,F)
    normalize!.(eachcol(A))
    A
end

# ╔═╡ 745e6b01-dfdc-4e41-a43b-8002cd5e8357
face_centers(V,F) = dropdims(sum(V[:,F], dims=2) ./ 3; dims=2)

# ╔═╡ 321527d2-068c-4e85-8947-f6d1d6fe4fd3
face_area(V,F) = norm.(eachcol(face_area_normals(V,F)))

# ╔═╡ 05d16cc7-3f0a-426c-954a-63d840708777
function vertex_area(V,F)
    B = zeros(size(V)[2])
    for f in eachcol(F)
        T = V[:,f]
        x = T[:,1] - T[:,2]
        y = T[:,3] - T[:,2]
        A = 0.5*(sum(cross(x,y).^2)^0.5)
        B[f] .+= A
    end
    B ./= 3
    return B
end

# ╔═╡ f324cbad-1f68-4359-9ea4-8eb8f13b27e1
function face_grad(V,F)
    A = face_area(V,F)
    N = face_normals(V,F)

    u = repeat(F[1,:], inner=3)
    v = repeat(F[2,:], inner=3)
    w = repeat(F[3,:], inner=3)
    uv = V[:,F[2,:]] - V[:,F[1,:]]
    vw = V[:,F[3,:]] - V[:,F[2,:]]
    wu = V[:,F[1,:]] - V[:,F[3,:]]
    J = 1:3*size(F,2)
    G2 = cross.(eachcol(N), eachcol(wu)) ./ A
    G1 = cross.(eachcol(N), eachcol(vw)) ./ A
    G3 = cross.(eachcol(N), eachcol(uv)) ./ A
    G1 = collect(Iterators.flatten(G1))
    G2 = collect(Iterators.flatten(G2))
    G3 = collect(Iterators.flatten(G3))
    g = sparse([J;J;J], [u;v;w], [G1; G2; G3], 3*size(F,2), size(V,2))
    g ./= 2
    return g
end

# ╔═╡ b7b393a9-e6f6-48ea-a4d5-d3e327f6bc18
function div(V,F)
    # ∇⋅ is |V|×3|F|
    uv = V[:,F[2,:]] - V[:,F[1,:]]
    vw = V[:,F[3,:]] - V[:,F[2,:]]
    wu = V[:,F[1,:]] - V[:,F[3,:]]
    cotan = -[sum(uv.*wu; dims=1); sum(vw.*uv; dims=1); sum(wu.*vw; dims=1)]
    s = [norm.(cross.(eachcol(uv), eachcol(wu)));;
        norm.(cross.(eachcol(vw), eachcol(uv)));;
        norm.(cross.(eachcol(wu), eachcol(vw)))]'
    cotan ./= s

    u = repeat(F[1,:], inner=3)
    v = repeat(F[2,:], inner=3)
    w = repeat(F[3,:], inner=3)
    J = 1:3*size(F,2)
    A = vec(-cotan[2,:]' .* wu + cotan[3,:]' .* uv)
    B = vec(cotan[1,:]' .* vw - cotan[3,:]' .* uv)
    C = vec(-cotan[1,:]' .* vw + cotan[2,:]' .* wu)
    ∇ = sparse([u;v;w], [J;J;J], [A;B;C], size(V,2), 3*size(F,2))
    ∇ ./= 2
end

# ╔═╡ 1e1b6bed-9224-4968-afb6-6bbc9d635191
md"""
### Examples
Now for some examples. We use the PlyIO package to read in meshes.
"""

# ╔═╡ 7d9647bb-b0a1-4a21-a496-b43f2c61e7fe
html"""
<head>
    <title>Toggle Online Image Button</title>
    <style>
        /* Initially hide the image */
        #toggleImage {
            display: none;
        }
    </style>
</head>
<body>
    <button id="toggleButton">?</button>
    <img id="toggleImage" src="https://media.tenor.com/6xwjsmMIAIoAAAAd/happy-happy-dog.gif" alt="Image">
    
    <script>
        // Get references to the button and the image
        const toggleButton = document.getElementById('toggleButton');
        const toggleImage = document.getElementById('toggleImage');

        // Add a click event listener to the button
        toggleButton.addEventListener('click', function() {
            // Toggle the image's visibility
            if (toggleImage.style.display === 'none') {
                toggleImage.style.display = 'block';
            } else {
                toggleImage.style.display = 'none';
            }
        });
    </script>
</body>
"""

# ╔═╡ b18c2afe-855c-4dd4-858d-f50a8cdd92fd
begin
# Download Stanford bunny and do some data maniuplation...
		bunny_file = download("https://raw.githubusercontent.com/naucoin/VTKData/master/Data/bunny.ply")
		bunny_ply = load_ply(bunny_file)
		V = stack(Array(bunny_ply["vertex"][i]) for i in ["x", "y", "z"])'
		F = stack(Array(bunny_ply["face"]["vertex_indices"])) .+ 1
		not_in = Set{Int}()
	n = size(V,2)
	for i=1:n
	    !(i in F) && push!(not_in, i)
	end
	V = V[:,setdiff(1:n, not_in)]
	n = size(V, 2)
	F = map(F) do i
	    !(i in not_in) && return i - sum(i .> not_in)
	    return i
	end
end

# ╔═╡ 860530d1-3f6c-4774-91be-01b7aec16f91
begin
	fig = Figure(resolution=(2500,900))
	ax1 = Axis3(fig[1,1], viewmode=:fit, aspect = (1, 1, 1), azimuth=0, elevation=0)
	hidespines!(ax1)
	hidedecorations!(ax1)
	ax2 = Axis3(fig[1,2], viewmode=:fit, aspect= (1,1,1), azimuth=0, elevation =0)
	hidespines!(ax2)
	hidedecorations!(ax2)
	L = cotlaplacian(V,F)
	A = vertex_area(V,F)
	u0 = zeros(size(L,1))
	u0[1] = 1.0
	heat = Observable(u0)
	∇_bunny = face_grad(V,F);
	Δ_bunny = div(V,F)
	anchor_points = face_centers(V,F)[:,1:20:end]
	show_gradients = Observable(true)
	heat_gradx = Observable(zeros(ceil(Int, size(F,2)/20)))
	heat_grady = Observable(zeros(ceil(Int, size(F,2)/20)))
	heat_gradz = Observable(zeros(ceil(Int, size(F,2)/20)))
	geodesics = Observable(zeros(size(V,2)))
end

# ╔═╡ ffdbb945-0e85-409a-9bd7-a5224f2724f9
md"""
Drag the sliders and bunnies.
"""

# ╔═╡ cbc43250-5168-4211-a92f-4f99209b07a0
md"""t $(@bind t_bunny PlutoUI.Slider(1:200, default=10, show_value=true)) v1 $(@bind v1 PlutoUI.Slider(1:1000, default=488, show_value=true))   v2 $(@bind v2 PlutoUI.Slider(1:1000, default=488, show_value=true)) v3 $(@bind v3 PlutoUI.Slider(1:1000, default=488, show_value=true))"""

# ╔═╡ 450ad839-16ce-406e-9269-665dba06937e
md"""Show Gradient Field$(@bind has_grad PlutoUI.CheckBox())"""

# ╔═╡ 273a2353-32c6-4509-aa39-d94e45907000
let
	mesh!(ax1,V[[3,1,2],:]',F[[2,1,3],:]', color = heat)
	arrows!(ax1, anchor_points[3,:],anchor_points[1,:],anchor_points[2,:], heat_gradz, heat_gradx, heat_grady; linewidth=0.001, arrowsize=0.001, lengthscale=0.001, arrowcolor=:red, linecolor=:red)
	mesh!(ax2,V[[3,1,2],:]',F[[2,1,3],:]', color = geodesics, colormap=:prism)
	fig
end

# ╔═╡ 92915d09-067f-4ea6-a6a9-0f519c9ea84d
md"""
# Appendix
"""

# ╔═╡ 7d9c0d8d-e52f-44b9-ae77-bbda953c498c
md"""
### Handwriting DiffEq.jl

There are two common ways to solve the heat equation if given the operator $L$. The first one is an iterative approach that uses the entries of $L$ to update the solution like how computing with a stencil works. The second method solves the equation in a manner similar to what a person would do to solve a PDE like the heat equation by hand - by finding the eigenfunctions.
"""

# ╔═╡ 070107d5-40e4-4f63-bdb9-9f315ebf18ba
md"""
#### Iterative Mode
Back in calculus, we used Euler's method to solve the solutions to a differential equation. The finer the resolution of the update - in our case time $t$, the more accurate our solutions were. Let $h_t(x)$ be the heat on the mesh at time $t$. We linearize the time dimension which yields roughly

$$-A^{-1}Lh_t \approx \frac{h_{t'}-h_t}{dt}$$
or

$$h_{t'} = A^{-1}(A-dtL)h_t$$

Equation (...) is referred to as foward-mode integration since we take the previously known value and use the derivative to update it. This works and is generally fast, but I want to point out a slight variant that happens to be more stable. It turns out that if we replace $Lh_t$ and $Lh_{t'}$ as our choice of the derivative, then this is referred to as semi-implicit integration. The step is instead

$$A h_t=(A+Ldt)h_{t'}$$
"""

# ╔═╡ 4513d343-80c2-48c9-b628-5df1fde04b76
function heat_implicit(L, A, init; dt=0.001, steps=1)
    M = spdiagm(A)
    D = cholesky(M+dt*L)
    heat = init
    for t=1:steps
        heat = D \ (M*heat)
    end
    return heat
end

# ╔═╡ 8c6eeb0b-54a6-44a7-9579-55a4de69e31d
begin
	_u0 = zeros(size(L,1))
	_u0[v1] = 1.0
	_u0[v2] = 1.0
	_u0[v3] = 1.0
	heat[] = heat_implicit(L, A, _u0; dt=0.0001, steps=t_bunny)
	show_gradients[] = has_grad
	heat_grad = reshape(∇_bunny * heat[], 3, :)
	heat_grad ./= sqrt.(sum(heat_grad .^ 2, dims=1))
	_heat_grad = heat_grad[:,1:20:end]
	heat_gradx[] = _heat_grad[1,:]
	heat_grady[] = _heat_grad[2,:]
	heat_gradz[] = _heat_grad[3,:]
	X = Δ_bunny * vec(heat_grad)
	dist = L \ X
	dist .-= minimum(dist)
	geodesics[] = dist
end

# ╔═╡ 948cee7b-2533-4693-a333-88231054ff83
md"""
#### Spectral Mode
So far computing heat flow has relied on simulations. Alternatively, we can appeal to theory to solve the heat equation. Normally, we need to find the *eigenfunctions* of the Laplacian subject to the boundary conditions. Specifically, $f$ is an eigenfunction if it satisfies the boundary conditions and $\Delta f = λf$. For heat diffusion with no source, the solutions to the heat equation are then $h(x,t) = e^{\lambda t}f(x)$. We can verify this by plugging it back into the PDE:

$$\lambda e^{\lambda t}f(x) = e^{\lambda t}f''(x) = \lambda e^{\lambda t}f(x)$$

The same logic applies to the discrete setting. Let $L$ be the Laplacian matrix. Then the eigenvalues and eigenvectors form a basis to the solutions to the PDE.
"""

# ╔═╡ 3ac61216-5029-47b5-85e4-3fc27f879e52
function heat_spectral(λ::Vector, ϕ::Matrix, init, t)
    """
    init - |V| or |V|×|C|
    note that init may be a single vector or length(t) vectors. If it is a single vector, then heat is diffused for each time t. If it is
    multiple vectors, then vector init[i] is diffused for time t[i]
    c = ϕ'*(A .* init) .* exp.(-λ * t')
    
    If A is defined, then inner products are treated as A-inner products
    """
    c = ϕ'*(init) .* exp.(-λ * t')
    heat = abs.(ϕ * c)
end

# ╔═╡ f798d3c7-278a-4c9b-aa89-1f8c2b94a938
function heat_spectral(λ, ϕ, A, init, t)
    c = ϕ' * (A*init) .* exp.(-λ * t')
    ϕ * c
end

# ╔═╡ 54f1f6fe-acc1-4308-b5f9-1694de5dab7f
md"""
### GPU Support
???
"""

# ╔═╡ 2a5ae2f3-f9a0-4ade-9307-f617a41e36ce
Markdown.MD(Markdown.Admonition("", "Extra: Shortest Paths on Graphs", [md"""
We can come full circle and apply ask the question can "we apply heat on a graph to compute the shortest path?". Really, so long as you have a way of taking Laplacians, gradients, and divergences then in principle you can repeat all steps and get a heat method for your domain of interest. Graphs are another type of setting with Laplacians. Laplacians defined for graphs are a popular tool in machine learning for normalization so their properties are well established. When I tried implementing this however, I encountered problems in defining a divergence discretization for graphs. After a quick google search, it was no surprise that I found a note by Keenan himself who explains how to do it. The set up is a bit janky, but works quite well. On your next software interview, consider using this method rather than BFS.
"""]))


# ╔═╡ a9dfc324-03ec-4884-a19e-4371ac069e1b
md"""
[^1]: Paper: https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/paperCACM.pdf
[^2]: We will address how to simulate heat in a later section.
[^3]: These plots were generated using the Fast Marching Method, another algorithm to approximate geodesics. The down side of the Fast Marching Method - and one of the original motiviations for the Heat Method - is that is is slow and not parallelizable. Code credit for the plots goes to [ffevotte](https://github.com/ffevotte). Check out [Eikonal.jl](https://github.com/triscale-innov/Eikonal.jl) for other nifty applications. 
[^4]: A friend of mine considers this to be the most based PDE of all PDEs.
[^5]: Thermal diffusivity is assumed to be $1$. For an introduction on the heat equation, check these [notes](https://web.stanford.edu/class/math220b/handouts/heateqn.pdf).
[^6]: Adapted from [this SciML tutorial](https://docs.sciml.ai/DiffEqDocs/stable/examples/beeler_reuter/).
[^7]: Wait vertex area? How does a vertex have an area? Here you can consider the barycentric area - the sum of the areas of the faces adjacent to the vertex divided by $3$.
"""

# ╔═╡ c1f1cb61-e802-4d72-8cc8-837450e9f35c
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Arpack = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Eikonal = "a6aab1ba-8f88-4217-b671-4d0788596809"
JSServe = "824d6782-a2ef-11e9-3a09-e5662e0c26f9"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlyIO = "42171d58-473b-503a-8d5f-782019eb09ec"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

[compat]
Arpack = "~0.5.4"
ColorSchemes = "~3.24.0"
DifferentialEquations = "~7.9.1"
Eikonal = "~0.1.1"
JSServe = "~2.2.10"
Plots = "~1.39.0"
PlutoUI = "~0.7.52"
PlyIO = "~1.1.2"
WGLMakie = "~0.8.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0-beta2"
manifest_format = "2.0"
project_hash = "92b21b001355236ecd687cecd838727252a04c65"

[[deps.ADTypes]]
git-tree-sha1 = "5d2e21d7b0d8c22f67483ef95ebdc39c0e6b6003"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.4"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractLattices]]
git-tree-sha1 = "f35684b7349da49fcc8a9e520e30e45dbb077166"
uuid = "398f06c4-4d28-53ec-89ca-5b2656b7603d"
version = "0.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "0d61921af2799487b80453a44abb57db7a0c1381"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.4.1"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "ef9997b3d5547c48b41c7bd8899e812a917b409d"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.4"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "0b816941273b5b162be122a6c94d706e3b3125ca"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.38"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["ArrayInterface", "BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NonlinearSolve", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "TruncatedStacktraces", "UnPack"]
git-tree-sha1 = "f7392ce20e6dafa8fee406142b1764de7d7cd911"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "4.0.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "601f7e7b3d36f18790e2caf83a882d88e9b71ff1"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.4"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "32abd86e3c2025db5172aa182b982debed519834"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "b86ac2c5543660d238957dbde5ac04520ae977a7"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "f9d7112bfff8a19a3a4ea4e03a8e6a91fe8456bf"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["DataStructures", "EnumX", "ExactPredicates", "Random", "SimpleGraphs"]
git-tree-sha1 = "bea7984f7e09aeb28a3b071c420a0186cb4fabad"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "0.8.8"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "89f3fbfe78f9d116d1ed0721d65b0b2cf9b36169"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.42.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cd6756e833c377e0ce9cd63fb97689a255f12323"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "1.33.4+0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DataStructures", "DocStringExtensions", "EnumX", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "6ece6f2956dea6380d8e5e6eadaab1bb9af20cce"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.129.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"
    DiffEqBaseZygoteExt = "Zygote"

    [deps.DiffEqBase.weakdeps]
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "Functors", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "42424e81924d4f463c6f8db8ce2978d51ba0aeaf"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.30.0"
weakdeps = ["OrdinaryDiffEq", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6b02e9c9d0d4cacf2b20f36c33710b8b415c5194"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.18.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "c3d11164d1b08c379bc3c6abae45fcd7250e8e35"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.9.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "b6def76ffad15143924a2199f72a5cd883a2e8a9"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.9"
weakdeps = ["SparseArrays"]

    [deps.Distances.extensions]
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "938fe2981db009f531b6332e31c58e9584a2f9bd"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.100"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Eikonal]]
deps = ["DataStructures", "Images", "LinearAlgebra", "PrecompileTools", "Printf"]
git-tree-sha1 = "ac89a6cf8c89a741448deb8692aaacba745ecee0"
uuid = "a6aab1ba-8f88-4217-b671-4d0788596809"
version = "0.1.1"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ErrorfreeArithmetic]]
git-tree-sha1 = "d6863c556f1142a061532e79f611aa46be201686"
uuid = "90fa49ef-747e-5e6f-a989-263ba693cf1a"
version = "0.5.2"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArraysCore", "Test"]
git-tree-sha1 = "276e83bc8b21589b79303b9985c321024ffdf59c"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.5"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "602e4585bcbd5a25bc06f514724593d13ff9e862"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.25.0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "aa9925a229d45fe3018715238956766fa21804d1"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.6"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "b12f05108e405dadcc2aff0008db7f831374e051"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.0"

[[deps.FastRounding]]
deps = ["ErrorfreeArithmetic", "LinearAlgebra"]
git-tree-sha1 = "6344aa18f654196be82e62816935225b3b9abe44"
uuid = "fa42c844-2597-5d31-933b-ebd51ab2693f"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "a20eaa3ad64254c61eeb5f230d9306e937405434"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.6.1"
weakdeps = ["SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "50351f83f95282cf903e968d7c6e8d44a5f83d0b"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9a68d75d466ccc1218d0552a8e1631151c569545"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "8e2d86e06ceb4580110d9e716be26658effc5bfd"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "da121cbdc95b065da07fbb93638367737969693f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.8+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "bb198ff907228523f3dee1070ceee63b9359b6ab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "424a5a6ce7c5d97cca7bcc4eac551b97294c54af"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.9"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "f57a64794b336d4990d90f80b147474b869b1bc4"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "3447781d4c80dbe6d71d239f7cfb1f8049d4c84f"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.6"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "b0b765ff0b4c3ee20ce6740d843be8dfce48487c"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.3.0"

[[deps.ImageMagick_jll]]
deps = ["JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1c0a2295cca535fabaf2029062912591e9b61987"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.10-12+3"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.ImageMorphology]]
deps = ["ImageCore", "LinearAlgebra", "Requires", "TiledIteration"]
git-tree-sha1 = "e7c68ab3df4a75511ba33fc5d8d9098007b579a8"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.3.2"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "44664eea5408828c03e5addb84fa4f916132fc26"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.8.1"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "8717482f4a2108c9358e5c3ca903d3a6113badc9"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.9.5"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "5fa9f92e1e2918d9d1243b1131abe623cdf98be7"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.25.3"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "be8e690c3973443bec584db3346ddc904d4884eb"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad37c091f7d7daf900963171600d7c1c5c3ede32"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "FastRounding", "LinearAlgebra", "Markdown", "Random", "RecipesBase", "RoundingEmulator", "SetRounding", "StaticArrays"]
git-tree-sha1 = "5ab7744289be503d76a944784bac3f2df7b809af"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.20.9"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "8e59ea773deee525c99a8018409f64f19fb719e6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.7"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "773125c999b4ebfe31e679593c8af7f43f401f1c"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.34"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSServe]]
deps = ["Base64", "CodecZlib", "Colors", "Dates", "Deno_jll", "HTTP", "Hyperscript", "LinearAlgebra", "Markdown", "MsgPack", "Observables", "RelocatableFolders", "SHA", "Sockets", "Tables", "Test", "ThreadPools", "URIs", "UUIDs", "WidgetsBase"]
git-tree-sha1 = "38be9964165e8693b63f2d5ba2b6154dfd69c3b1"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f9"
version = "2.2.10"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "327713faef2a3e5c80f96bf38d1fa26f7a6ae29e"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "61068b4df1e434c26ff8b876fbaf2be3e3e44d27"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.7.3"
weakdeps = ["FastBroadcast"]

    [deps.JumpProcesses.extensions]
    JumpProcessFastBroadcastExt = "FastBroadcast"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "884c2968c2e8e7e6bf5956af88cb46aa745c854b"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.1"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "90442c50e202a5cdf21a7899c66b240fdef14035"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.7"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "17e462054b42dcdda73e9a9ba0c67754170c88ae"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.4"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.0.1+1"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "Permutations", "Primes", "SimplePolynomials"]
git-tree-sha1 = "558a338f1eeabe933f9c2d4052aa7c2c707c3d52"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.1.12"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "Libdl", "LinearAlgebra", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "69cbd612e6e67ba2f8121bc8725bc9d04d803599"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.5.1"

    [deps.LinearSolve.extensions]
    LinearSolveCUDAExt = "CUDA"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMKLExt = "MKL_jll"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"

    [deps.LinearSolve.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    MKL_jll = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "0d097476b6c381ab7906460ef1ef1638fbce1d91"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.2"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "c88a4afe1703d731b1c4fdf4e3c7e77e3b176ea2"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.165"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "eb006abbd7041c28e0d16260e50a24f8f9104913"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Setfield", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "cf10f4b9d09da50f124ab7bcb530e57f700328f0"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.10"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "17d51182db2667962bc7e1d18b74881d0d0adbe6"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.7"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "8f52dbaa1351ce4cb847d95568cb29e62a307d93"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.6"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "1130dbe1d5276cb656f6e1094ce97466ed700e5a"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mods]]
git-tree-sha1 = "61be59e4daffff43a8cec04b5e0dc773cbb5db3a"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "1.3.3"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "fc8c15ca848b902015bd4a745d350f02cf791c2a"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.0"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.Multisets]]
git-tree-sha1 = "8d852646862c96e226367ad10c8af56099b4047e"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "EnumX", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "PrecompileTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "ee53089df81a6bdf3c06c17cf674e90931b10a73"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "1.10.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "963b004d15216f8129f6c0f7d187efa136570be0"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.7"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SciMLOperators", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "ba3ed480f991b846cf9a8118d3370d9752e7166d"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.55.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "9b02b27ac477cad98114584ff964e3052f656a0f"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.0"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "f9b1e033c2b1205cf30fd119f4e50881316c1923"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.1"
weakdeps = ["Requires", "TOML"]

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Permutations]]
deps = ["Combinatorics", "LinearAlgebra", "Random"]
git-tree-sha1 = "25e2bb0973689836bf164ecb960762f1bb8794dd"
uuid = "2ae35dd2-176d-5d53-8349-f30d82d94d4f"
version = "0.4.17"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PlyIO]]
git-tree-sha1 = "74619231a7aa262a76f82ae05c7385622d8a5945"
uuid = "42171d58-473b-503a-8d5f-782019eb09ec"
version = "1.1.2"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "d4c9ebdc6528a4aaf7cfcf43b482e927267b400d"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.6"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "ea78a2764f31715093de7ab495e12c0187f231d1"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.4"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "4c9f306e5d6603ae203c2000dd460d81a5251489"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "eeab25344bf9901146c0200a7ca64ea479f8bf5c"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.0"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "da095158bdc8eaccb7890f9884048555ab771019"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "552f30e847641591ba3f39fd1bed559b9deb0ef3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.1"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d7087c013e8a496ff396bae843b1e16d9a30ede8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.10"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "2b6d4a40339aa02655b1743f4cd7c03109f520c1"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.20"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.RingLists]]
deps = ["Random"]
git-tree-sha1 = "f39da63aa6d2d88e0c1bd20ed6a3ff9ea7171ada"
uuid = "286e9d63-9694-5540-9e3c-4e6708fa07b2"
version = "0.2.8"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "54ccb4dbab4b1f69beb255a2c0ca5f65a9c82f08"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.5.1"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "4b8586aece42bee682399c4c4aee95446aa5cd19"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.39"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "ChainRulesCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "6de099dba3c80e23bde1d19011161ea91a23ed6b"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.98.0"

    [deps.SciMLBase.extensions]
    ZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLNLSolve]]
deps = ["DiffEqBase", "LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "9dfc8e9e3d58c0c74f1a821c762b5349da13eccf"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.8"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "65c2e6ced6f62ea796af251eb292a0e131a3613b"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.6"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SetRounding]]
git-tree-sha1 = "d7a25e439d07a17b7cdf97eecee504c50fedf5f6"
uuid = "3cc68bcd-71a2-5612-b932-767ffbe40ab0"
version = "0.2.1"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "db0219befe4507878b1a90e07820fed3e62c289d"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleGraphs]]
deps = ["AbstractLattices", "Combinatorics", "DataStructures", "IterTools", "LightXML", "LinearAlgebra", "LinearAlgebraX", "Optim", "Primes", "Random", "RingLists", "SimplePartitions", "SimplePolynomials", "SimpleRandom", "SparseArrays", "Statistics"]
git-tree-sha1 = "b608903049d11cc557c45e03b3a53e9260579c19"
uuid = "55797a34-41de-5266-9ec1-32ac4eb504d3"
version = "0.8.4"

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "PackageExtensionCompat", "PrecompileTools", "Reexport", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "20aa9831d654bab67ed561e78917047143ecb9bf"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.19"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveNNlibExt = "NNlib"

    [deps.SimpleNonlinearSolve.weakdeps]
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"

[[deps.SimplePartitions]]
deps = ["AbstractLattices", "DataStructures", "Permutations"]
git-tree-sha1 = "dcc02923a53f316ab97da8ef3136e80b4543dbf1"
uuid = "ec83eff0-a5b5-5643-ae32-5cbf6eedec9d"
version = "0.3.0"

[[deps.SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "d537c31cf9995236166e3e9afc424a5a1c59ff9d"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.14"

[[deps.SimpleRandom]]
deps = ["Distributions", "LinearAlgebra", "Random"]
git-tree-sha1 = "3a6fb395e37afab81aeea85bae48a4db5cd7244a"
uuid = "a6525b86-64cd-54fa-8f65-62fc48bdc0e8"
version = "0.3.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "4b33e0e081a825dbfaf314decf58fa47e53d6acb"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.4.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "42d131931906bf4f0af97a7113c8456d0a8aff9d"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.6.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["Compat", "SHA", "Tables", "TupleTools"]
git-tree-sha1 = "19df33ca14f24a3ad2df9e89124bd5f5cc8467a2"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "1.0.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f295e0a1da4ca425659c57441bcb59abb035a4bc"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.8"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "03fec6800a986d191f64f5c0996b59ed526eda25"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.4.1"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "51621cca8651d9e334a659443a74ce50a3b6dfab"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.3"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "6e801d0da4c81d9cd6a05d97340404f9892fba85"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.16.0"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "b341540a647b39728b6d64eaeda82178e848f76e"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.62.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f02eb61eb5c97b48c153861c72fbbfdddc607e06"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.17"

[[deps.StructArrays]]
deps = ["Adapt", "ConstructionBase", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "0a3db38e4cce3c54fe7a71f831cd7b6194a54213"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.16"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.0+1"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "deea053391e5b352594030ac95bbc3be52855a69"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.19.5"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "ba4d38faeb62de7ef47155ed321dce40a549c305"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.2+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a1f34829d5ac0ef499f6d84428bd6b4c71f02ead"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "3c4535892eff963d14acee719df445287c2d8f98"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.5"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.TupleTools]]
git-tree-sha1 = "155515ed4c4236db30049ac1495e2969cc06be9d"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.4.3"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.WGLMakie]]
deps = ["Colors", "FileIO", "FreeTypeAbstraction", "GeometryBasics", "Hyperscript", "JSServe", "LinearAlgebra", "Makie", "Observables", "PNGFiles", "PrecompileTools", "RelocatableFolders", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "8a430666e4df430efe88a44f51eaca25f2dd293f"
uuid = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
version = "0.8.14"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "04a51d15436a572301b5abbb9d099713327e9fc4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.4+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─7dc7f1eb-3e1f-4b75-8e4c-e3018b0259d6
# ╟─edce4ccf-7be7-4b0d-98aa-c572eac3d0ad
# ╠═dab8582c-e8e2-443f-b52c-ac716b2ca12e
# ╠═efa0e0ba-8b30-4f69-9bc8-bdd89ca5f61a
# ╠═358c1832-06c0-4918-b766-8c1de98c21d3
# ╟─94ca9000-947d-44d2-8a6d-f9177c476345
# ╠═ad43bc4b-3cc4-4960-bb46-eb6978620e61
# ╟─c235d925-3899-4cf7-9d9d-847ab81a7523
# ╟─ba6a8f2e-3e01-48c9-a07a-72fcde78e6f1
# ╠═008abc67-4a20-4ba7-a7e0-a1922fd67387
# ╟─8c2aaf03-a48f-4241-9e78-3b7312ee71e2
# ╟─5a53c282-c0be-47b3-964f-2fedbfa25451
# ╟─a73d4efc-b9a7-4171-9bdc-c98e907dd2b7
# ╟─6f2cdb0d-c00d-48ca-a2d0-ea462532895d
# ╟─97460bf6-ce4f-4210-9664-8c6c43a9a382
# ╟─0526ff5b-952c-4cba-842b-1059c67c12e1
# ╟─f0eb3973-f9c4-41fc-8f38-3bcb71e76c7d
# ╟─c8fc76f7-f900-460b-a6e3-a33a3386a8e0
# ╟─19696716-d7ae-4ffd-8b73-399e5f02831b
# ╟─0350b1d6-4245-4b86-8b45-be9c00a16c77
# ╟─e7080c15-ac7e-4106-8df4-65a668e39b83
# ╟─ac43bbab-3014-4ece-b738-157e6367f734
# ╟─3a23522e-d7d6-4461-bd33-878e1c823ea6
# ╠═0bdff4e8-27ad-46ef-b9b3-a146d774cc6d
# ╟─7564b349-5d51-44a0-b78a-654cd1bbfb61
# ╟─96fdedc9-3adf-4e46-9a50-d0b38bd38caf
# ╠═c89ce545-5229-418e-a174-e2e4eddc1115
# ╟─2dad5bf2-ac9c-4767-ac37-3abd252f338a
# ╠═c814f0d2-cbb0-4db9-9499-993d51f42356
# ╠═2e4e7cf6-0034-4992-9ea0-f212b4111fc1
# ╟─e1c19aea-d671-4c01-8bcf-119e7abb295f
# ╟─54c232ba-1175-40a3-b5a9-729450905e9f
# ╠═dd03796a-0520-418a-88f3-f11547b05a19
# ╟─217ac117-5415-4938-a543-8ebe5cca7898
# ╟─17ca6c6a-d282-457d-961d-40275a01927a
# ╟─c912af64-147f-4ca5-a567-45f5c5e50303
# ╟─3e920b70-d5b6-44f4-b257-e7568a458173
# ╠═d7ccc528-a811-4c31-8d64-fa2ce1e813cb
# ╟─e56e6e62-5f38-467d-83bd-daaf4a968044
# ╟─47e4da81-2c79-44bc-8d87-a9a0f2e45d4e
# ╟─47b00e38-83d1-4888-baee-662bd716827c
# ╠═387f2124-1802-405f-8d6e-3cfdcefe2f46
# ╟─3ce304a9-3dfd-42f1-b3cf-308b6ce3cbac
# ╠═ff09b1b2-3990-4bc7-8f6b-962e2ecfee3d
# ╟─ca26a33f-17c3-4a91-b95c-75b83409705e
# ╠═cffa4dc5-c9a8-4277-b87b-5dd3a5eff858
# ╟─2a77d78c-1b8a-4690-9024-46a6794d8efd
# ╠═7a7a67a7-7958-43bb-bf54-36b80ecdf1de
# ╠═b89464c0-baaa-4571-a744-91d5480c6ae1
# ╠═160902db-e0b9-4d48-8fcb-41dbeaf429e8
# ╟─2415e55d-bed0-4050-893e-b6a7af00ef45
# ╟─c0c94bdd-a7d0-46c3-a61e-6c0d40d8a3c9
# ╠═346c06a8-c91f-4c56-9840-67f2f02bd8c9
# ╠═79d22285-1c69-495c-b05d-83d8c758ee46
# ╠═a701e06d-553e-43b6-b36a-c68667bfd4b1
# ╟─7c9b744e-72cd-449e-8228-a25b5c845233
# ╟─c3d64ff4-9d3e-44da-93f1-827b93042fc9
# ╟─1bc67b2f-a1e9-4121-92cd-083b4ea9567b
# ╟─8a65dd84-2098-4765-8912-4ed6d32a9e0a
# ╠═e48d51a3-debc-4339-84fe-20ee9613e808
# ╟─9a0aa371-9fbf-493f-ba4e-cb0801c2d5ef
# ╟─861049ba-49d9-4f8c-a186-f1f95b282904
# ╠═83f7c9c7-cf37-44f0-8e51-ea6596d83605
# ╟─688712c4-b57f-49de-a1e9-3a3299eef60e
# ╟─92b945aa-b19f-4732-963b-a3e8b42a6b02
# ╟─d5b9c28b-e9fc-4e04-bf05-5b0b93da804e
# ╠═9dd097b6-5f82-4fbe-b0d1-86756b7747d2
# ╟─5bcf708c-8dc2-4a2d-a284-c17db2ea8b9a
# ╟─3ed68b8f-db59-40fe-87ef-8df03f81f9df
# ╠═3300a2ee-bc41-439a-8ec4-a981aab32a93
# ╠═61f3e733-6527-43b9-97bd-08459e0878fc
# ╠═f7914d06-c58e-4033-b895-9c069ec6eb4e
# ╠═dc7eece7-942a-491a-9eaa-033c19112d32
# ╠═745e6b01-dfdc-4e41-a43b-8002cd5e8357
# ╠═321527d2-068c-4e85-8947-f6d1d6fe4fd3
# ╠═05d16cc7-3f0a-426c-954a-63d840708777
# ╠═f324cbad-1f68-4359-9ea4-8eb8f13b27e1
# ╠═b7b393a9-e6f6-48ea-a4d5-d3e327f6bc18
# ╟─1e1b6bed-9224-4968-afb6-6bbc9d635191
# ╠═d79f5864-5193-4673-a593-1057ec15e927
# ╟─7d9647bb-b0a1-4a21-a496-b43f2c61e7fe
# ╠═b18c2afe-855c-4dd4-858d-f50a8cdd92fd
# ╠═860530d1-3f6c-4774-91be-01b7aec16f91
# ╠═8c6eeb0b-54a6-44a7-9579-55a4de69e31d
# ╟─ffdbb945-0e85-409a-9bd7-a5224f2724f9
# ╟─cbc43250-5168-4211-a92f-4f99209b07a0
# ╟─450ad839-16ce-406e-9269-665dba06937e
# ╟─273a2353-32c6-4509-aa39-d94e45907000
# ╟─92915d09-067f-4ea6-a6a9-0f519c9ea84d
# ╟─7d9c0d8d-e52f-44b9-ae77-bbda953c498c
# ╟─070107d5-40e4-4f63-bdb9-9f315ebf18ba
# ╠═4513d343-80c2-48c9-b628-5df1fde04b76
# ╟─948cee7b-2533-4693-a333-88231054ff83
# ╠═3ac61216-5029-47b5-85e4-3fc27f879e52
# ╠═f798d3c7-278a-4c9b-aa89-1f8c2b94a938
# ╟─54f1f6fe-acc1-4308-b5f9-1694de5dab7f
# ╟─2a5ae2f3-f9a0-4ade-9307-f617a41e36ce
# ╟─a9dfc324-03ec-4884-a19e-4371ac069e1b
# ╟─c1f1cb61-e802-4d72-8cc8-837450e9f35c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

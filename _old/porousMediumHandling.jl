using ExtendableGrids
using VoronoiFVM
# using DifferentialEquations
using Statistics
using PyPlot

# Computation configuration
N_x    = 70
N_y    = 35
T_heat = 0.5
steadyState = true
diffeqFlag  = false
dt!    = 1
tend   = 1000

# Constants definition
const k       = 100;
const g       = [0 -1];
const c       = 0.001;
const a       = 0.01;
const λ       = 0.01;
const rho_ref = 1;
const T_ref   = 0;
const P_ref   = 0;

function porousMedium(Nx, Ny, Tₕₑₐₜ, steadyStateFlag, diffeqFlag = false, method = ImplicitEuler(), tend = 10.0, dt! = 1.0)

    # Grid definition
    X    = collect(range(0, 300, length = Nx))
    Y    = collect(range(0, 150, length = Ny))
    grid = simplexgrid(X,Y)	

    
    function flux!(f, u, edge)
        # f[1] = - rho_ref*k* (u[1,2]-u[1,1] - 
        #                      (rho_ref - a*(0.5*(u[2,1]+u[2,2]) - T_ref))*project(edge,g))
        f[1] = - rho_ref*k* (u[1,2]-u[1,1] - 
                             ( - a*(0.5*(u[2,1]+u[2,2])))*project(edge,g))
        f[2] = λ*(u[2,1]-u[2,2])
    end

    function storage!(f, u, node)
        f[1] = rho_ref
        f[2] = c*u[2]
    end

    # Boundary conditions definition
    function bcondition!(f, u, bnode)
        
        v = ramp(bnode.time,du=(0,Tₕₑₐₜ),dt=(0,1.0)) 
        
        if diffeqFlag
            boundary_dirichlet!(f,u,bnode,species=2,region=1,value= Tₕₑₐₜ)
        else
            boundary_dirichlet!(f,u,bnode,species=2,region=1,value= v)
        end
        
        boundary_dirichlet!(f,u,bnode,species=2,region=3,value= 0)
        boundary_dirichlet!(f,u,bnode,species=1,region=3,value= 0)
        
        boundary_neumann!(f,u,bnode,species=2,region=2,value= 0)
        boundary_neumann!(f,u,bnode,species=2,region=4,value= 0)
        boundary_neumann!(f,u,bnode,species=1,region=2,value= 0)
        boundary_neumann!(f,u,bnode,species=1,region=3,value= 0)
        boundary_neumann!(f,u,bnode,species=1,region=4,value= 0)
    end

    # Voronoi system definiion
    system = VoronoiFVM.System(grid; 
                               flux=flux!,  
                               bcondition = bcondition!,
                               storage = storage!,
                               species=[1,2])

    # Main solver
    if steadyStateFlag
        # Steady state solution
        sol = VoronoiFVM.solve(system)
        nf  = nodeflux(system,sol)
        
        return grid,sol,nf
    else
        # Transient solution
        if diffeqFlag
            inival  = unknowns(system,inival=0)
            problem = ODEProblem(system,inival,(0.0,tend))
            odesol  = DifferentialEquations.solve(problem, method, 
                      dt = dt!, adaptive = false)
            sol     = reshape(odesol,system)
            
        else
            sol = VoronoiFVM.solve(system,inival = 0,
                                   times = (0,tend), 
                                   Δu_opt = T_heat, 
                                   Δt_min=dt!*0.1, Δt = dt!)	
        end

        nf = []
        for t in 0:dt!:tend
            push!(nf,nodeflux(system,sol(t)))
        end
        
        return grid,sol,nf
    end
end

grid, sol, nf = porousMedium(N_x, N_y, T_heat, 
	                         steadyState, 
	                         diffeqFlag, 0, tend, dt!)


# plot

x    = collect(range(0, 300, length = N_x))
y    = collect(range(0, 150, length = N_y))

X    = x' .* ones(N_y)
Y    = ones(N_x)' .* y

PyPlot.clf()

ax = PyPlot.axes()
ax.set_aspect("equal")

PyPlot.contourf(X, Y, reshape(sol[1,:], (N_x, N_y))', levels = 20)
# PyPlot.axes().plot_surface(X, Y, reshape(sol[1,:], (N_x, N_y))', levels = 20)

PyPlot.set_cmap("RdBu_r")

PyPlot.suptitle("Steady state pressure distribution for T(heat)= " *string(T_heat) )

PyPlot.tick_params(direction= "in", which= "minor", length= 2, bottom= true, top= true, right= true, left= true)
PyPlot.tick_params(direction= "in", which= "major", length= 4, bottom= true, top= true, right= true, left=true)

PyPlot.xlim([0, 300])
PyPlot.ylim([0, 150])

PyPlot.xlabel("x")
PyPlot.ylabel("y")

PyPlot.colorbar()
PyPlot.gcf()

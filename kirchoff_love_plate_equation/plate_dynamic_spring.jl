####### CANTILEVER PLATE WITH SPRING BOUNDARY 

using StaticArrays 
using SparseArrays
using LinearAlgebra 
using BenchmarkTools
using Plots 
using DifferentialEquations


#Mesh generation
function gen_mesh(N::NTuple{2,Int}, Lx::Float64, Ly::Float64)
    mesh = (range(0, Lx, length=N[1]+1), range(0, Ly, length=N[2]+1))
    return mesh
end

#Stiffness matrix 1D
function gen_stiffmat1D(r::NTuple{1,AbstractRange})
    xmesh = r[1]
    Np1 = length(xmesh)
    N = Np1-1; h = abs(xmesh[end] - xmesh[1])/N; h2 = h*h
    e = ones(Np1) 
    A = spdiagm(-1 => -e[2:end], 0 => 2*e, 1 => -e[2:end])
    return A     
end

#Stiffness matrix 2D
function gen_stiffmat(r::NTuple{2,AbstractRange})
    xmesh = r[1]; ymesh = r[2]  
    Nxp1 = length(xmesh); Nyp1 = length(ymesh);   
    #..contribution in x-direction 
    A1Dxx = gen_stiffmat1D(tuple(r[1]))
    A2Dxx = kron(A1Dxx,I(Nyp1))
    #..contribution in y-direction
    A1Dyy = gen_stiffmat1D(tuple(r[2]))
    A2Dyy = kron(I(Nxp1),A1Dyy)
    #..sum of x-direction and y-direction contribution 
    A = A2Dxx + A2Dyy
    return A 
end

#Boundary Conditions
function applyBC!(A::SparseMatrixCSC, mesh, nx, ny, load, h_spacing, D)

    Nx = nx+1
    Ny = ny+1
    N = Nx*Ny
    k=10000            #k=0 means no spring so cantilever plate
    c=-k*2*h_spacing^3/D

    # South boundary (w=0)
    for i in 1:(Nx)
        A[i, :] .= 0.0
        A[i, i] = 1.0
        load[i] = 0.0
    end

    #############################

    i=Nx+1
    A[i,i] = 16.0; A[i,i+1] = -20.0; A[i,i+2] = 10.0
    A[i,i-Nx] = -6.0; A[i,i-Nx+1] = 4.0; A[i,i-Nx+2] = -2.0
    A[i,i+Nx] = -6.0; A[i,i+Nx+1] = 4.0; A[i,i+Nx+2] = -2.0
    A[i,i+2*Nx] = 2.0

    i=Nx+2
    A[i, i-1] = -7.0; A[i,i] = 21.0; A[i,i+1] = -9.0; A[i,i+2] = 1.0
    A[i,i+2*Nx] = 2.0

    for i in (Nx+3):(2*Nx-2)
        A[i,i+2*Nx] = 2.0
    end

    i=2*Nx-1
    A[i,i-2] = 1.0; A[i,i-1] = -9.0; A[i,i] = 21.0; A[i, i+1] = -7.0
    A[i,i+2*Nx] = 2.0

    i=2*Nx
    A[i,i-2] = 10.0; A[i,i-1] = -20.0; A[i,i] = 16.0
    A[i,i-Nx-2] = -2.0; A[i,i-Nx-1] = 4.0; A[i,i-Nx] = -6.0
    A[i,i+Nx-2] = -2.0; A[i,i+Nx-1] = 4.0; A[i,i+Nx] = -6.0
    A[i,i+2*Nx] = 2.0

    ##################

    for i in (2*Nx+1):Nx:(N-3*Nx+1)
        A[i,i] = 16.0; A[i,i+1] = -20.0; A[i,i+2] = 10.0
        A[i,i-Nx] = -6.0; A[i,i-Nx+1] = 4.0; A[i,i-Nx+2] = -2.0
        A[i,i+Nx] = -6.0; A[i,i+Nx+1] = 4.0; A[i,i+Nx+2] = -2.0
    end

    for i in (2*Nx+2):Nx:(N-3*Nx+2)
        A[i, i-1] = -7.0; A[i,i] = 21.0; A[i,i+1] = -9.0; A[i,i+2] = 1.0
    end
    
    for i in (3*Nx-1):Nx:(N-2*Nx-1)
        A[i,i-2] = 1.0; A[i,i-1] = -9.0; A[i,i] = 21.0; A[i, i+1] = -7.0
    end
    
    for i in (3*Nx):Nx:(N-2*Nx)
        A[i,i-2] = 10.0; A[i,i-1] = -20.0; A[i,i] = 16.0
        A[i,i-Nx-2] = -2.0; A[i,i-Nx-1] = 4.0; A[i,i-Nx] = -6.0
        A[i,i+Nx-2] = -2.0; A[i,i+Nx-1] = 4.0; A[i,i+Nx] = -6.0
    end

    #################

    i=N-2*Nx+1
    A[i,i] = 17.0; A[i,i+1] = -20.0; A[i,i+2] = 10.0
    A[i,i-Nx] = -7.0; A[i,i-Nx+1] = 4.0; A[i,i-Nx+2] = -2.0
    A[i,i+Nx] = c/2 - 5.0; A[i,i+Nx+1] = 4.0; A[i,i+Nx+2] = -2.0

    i=N-2*Nx+2
    A[i,i-1] = -7.0; A[i,i] = 22.0; A[i,i+1] = -9.0; A[i,i+2] = 1.0
    A[i,i-Nx] = -9.0;
    A[i,i+Nx] = c/2-7.0; 

    for i in (N-2*Nx+3):(N-Nx-2)
        A[i,i] = 21.0
        A[i,i-Nx] = -9.0
        A[i,i+Nx] = c/2-7.0
    end

    i=N-Nx-1
    A[i,i-2] = 1.0; A[i,i-1] = -9.0; A[i,i] = 22.0; A[i,i+1] = -7.0;
    A[i,i-Nx] = -9.0;
    A[i,i+Nx] = c/2-7.0; 

    i=N-Nx
    A[i,i-2] = 10.0; A[i,i-1] = -20.0; A[i,i] = 17.0
    A[i,i-Nx-2] = -2.0; A[i,i-Nx-1] = 4.0; A[i,i-Nx] = -7.0
    A[i,i+Nx-2] = -2.0; A[i,i+Nx-1] = 4.0; A[i,i+Nx] = c/2 - 5.0

    #################

    i=N-Nx+1
    A[i,i] = -4.0*c+14.0; A[i,i+1] = 2.0*c-16.0; A[i,i+2] = -c+8.0
    A[i,i-Nx] = -16.0; A[i,i-Nx+1] = 8.0; A[i,i-Nx+2] = -4.0
    A[i,i-2*Nx] = 8.0; A[i,i-2*Nx+1] = -4.0; A[i,i-2*Nx+2] = 2.0

    i=N-Nx+2
    A[i,i-1] = c-5.0; A[i,i] = -5.0*c+17.0; A[i,i+1] = c-7.0; A[i,i+2] = 1.0
    A[i,i-Nx-1] = 4.0; A[i,i-Nx] = -20.0; A[i,i-Nx+1] = 4.0
    A[i,i-2*Nx-1] = -2.0; A[i,i-2*Nx] = 10.0; A[i,i-2*Nx+1] = -2.0

    for i in (N-Nx+3):(N-2)
        A[i,i-1] = c-6.0; A[i,i] = -5.0*c+16.0; A[i,i+1] = c-6.0
        A[i,i-Nx-1] = 4.0; A[i,i-Nx] = -20.0; A[i,i-Nx+1] = 4.0
        A[i,i-2*Nx-1] = -2.0; A[i,i-2*Nx] = 10.0; A[i,i-2*Nx+1] = -2.0
    end

    i=N-1
    A[i,i-2] = 1.0; A[i,i-1] = c-7.0; A[i,i] = -5.0*c+17.0; A[i,i+1] = c-5.0
    A[i,i-Nx-1] = 4.0; A[i,i-Nx] = -20.0; A[i,i-Nx+1] = 4.0
    A[i,i-2*Nx-1] = -2.0; A[i,i-2*Nx] = 10.0; A[i,i-2*Nx+1] = -2.0

    i=N
    A[i,i-2] = -c+8.0; A[i,i-1] = 2.0*c-16.0; A[i,i] = -4.0*c+14.0
    A[i,i-Nx-2] = -4.0; A[i,i-Nx-1] = 8.0;A[i,i-Nx] = -16.0
    A[i,i-2*Nx-2] = 2.0; A[i,i-2*Nx-1] = -4.0; A[i,i-2*Nx] = 8.0

    return A, load
end

# Uniform load function for 2D plate
function uniform_load(n, q, D)
    nx, ny = n
    load = fill(q, (nx + 1) * (ny + 1))
    return -load
end

# Point load function for 2D plate
function point_load(n, F, D)
    nx, ny = n
    load = fill(0, (nx + 1) * (ny + 1))
    load[7*div((nx+1)*(ny+1),12)] = F
    load[div((nx+1)*(ny+1),12)] = -F
    return -load
end

# Partially Distributed Load function for 2D plate
function distributed_load(n, F, D, node_indices)
    nx, ny = n
    load = fill(0.0, (nx + 1) * (ny + 1))
    for idx in node_indices
        load[idx] = F/length(node_indices)
    end

    return -load  
end

# Parameters
E = 70e9       # Young's modulus (Pa)
h = 0.01        # Plate thickness (m)
ν = 0.3         # Poisson's ratio
q = 1000.0      # Uniform load (N/m^2)
nx, ny = 20, 20   # Number of intervals (keep them equal)
Lx, Ly = 2.0, 2.0  # Plate dimensions (m)
rho = 2700        # Density of plate (kg/m^3)
F = 1000.0       #Force function (N)

D = ( E * h^3) / (12 * (1 - ν^2))

# Generate mesh and matrices
N = (ny, nx)
mesh = gen_mesh(N, Lx, Ly)
h_spacing = Lx / nx  # Assuming uniform spacing in both directions
A = gen_stiffmat(mesh)
A = (A * A) 

# Load vector
f = point_load(N, F, D) #choose load

A, f = applyBC!(A, mesh, nx, ny, f, h_spacing, D)

A = A / h_spacing^4

# Calculate critical damping coefficient
λ_max = maximum(real(eigvals(Matrix(A))))
K_eff = D * λ_max
M = rho * h * Lx * Ly
c_crit = 2 * sqrt(K_eff * M)

function plate_dynamics!(du, u, p, t)
    D, A, f, rho, h = p
    N = length(f)
    
    du[1:N] .= u[N+1:end]
    
    du[N+1:end] .= (f - D * A * u[1:N]-0.0001*c_crit*u[N+1:end]) / (rho * h) #choose damping factor
end

# Initial conditions (zero displacement and velocity)
N = length(f)
u0 = zeros(2 * N)
u0[1:N] .= 0.0 
u0[N+1:end] .= 0.0  

tspan = (0.0, 2.0)

p = (D, A, f, rho, h)
prob = ODEProblem(plate_dynamics!, u0, tspan, p)
sol = solve(prob, Tsit5(), abstol = 10e-9, reltol = 10e-9)

println("Solution Ready")

#Plotting
nxp, nyp = nx + 1, ny + 1

displacement = reshape(hcat(sol.u...)[1:nxp * nyp, :], nxp, nyp, :)
displacement = displacement*10^3  #convert solution to mm
zlims = (minimum(displacement), maximum(displacement))

x = LinRange(0, Lx, nxp)
y = LinRange(0, Ly, nyp)
times = sol.t
frame_count = 500  # Desired number of frames
step_size = max(1, div(length(times), frame_count))
# Create an animated 3D surface plot
anim = @animate for i in 1:step_size:length(times)
    # surface(x, y, displacement[:, :, i]', title="Plate Deflection at t=$(round(times[i], digits=2)) s",
    #         xlabel="X-axis (m)", ylabel="Y-axis (m)", zlabel="Displacement (mm)",  
    #         zlims=zlims, aspect_ratio = Ly/Lx)
    wireframe(x, y, displacement[:, :, i]', title="Plate Deflection at t=$(round(times[i], digits=2)) s",
            xlabel="X-axis (m)", ylabel="Y-axis (m)", zlabel="Displacement (mm)",  
            zlims=zlims, aspect_ratio = Ly/Lx)
end

# Save animation
gif(anim, "dynamic_plate_spring.gif", fps=30)
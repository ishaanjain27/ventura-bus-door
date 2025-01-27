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
function applyBC!(A::SparseMatrixCSC, mesh, nx, ny, load)

    # South boundary (w=0)
    for i in 1:(nx + 1)
        A[i, :] .= 0.0
        A[i, i] = 1.0
        load[i] = 0.0
    end

    # North boundary (w=0)
    for i in 0:(nx)
        A[end-i, :] .= 0.0
        A[end-i, end-i] = 1.0
        load[end-i] = 0
    end

    # East and West boundary (w=0)
    for i in 2:(ny+1)
        idx = (i - 1) * (nx + 1) 
        A[idx, :] .= 0.0
        A[idx, idx] = 1.0
        load[idx] = 0.0
        if i == (ny+1)
            idx = i * (nx + 1)
            A[idx, :] .= 0.0
            A[idx, idx] = 1.0
            load[idx] = 0.0
        end
        idx = (i - 1) * (nx + 1) + 1  
        A[idx, :] .= 0.0
        A[idx, idx] = 1.0
        load[idx] = 0.0
    end

    # South West, South East, North West, North East corners (right next to the boundaries) => w'' = 0
    A[nx+3, nx+3] = A[2*(nx+1)-1, 2*(nx+1)-1] = A[end - nx - 2,end - nx - 2] = A[end - (2*(nx+1) - 2),end - (2*(nx+1) - 2)] = 18.0

    # One above south boundary => w'' = 0
    for i in (nx + 4) : (2*(nx+1)-2)
        A[i,i] = 19.0
    end

    # One below north boundary => w'' = 0
    for i in (nx+3) : (2*(nx+1)-3)
        A[end-i, end-i] = 19.0
    end 

    # One next to west boundary and one next to east boundary => w'' = 0
    for i in 1:(ny-3)
        idx = (i+1)*(nx+1)+2
        A[idx,idx] = 19.0
        idx += (nx-2)
        A[idx,idx] = 19.0
    end

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
N = (nx, ny)
mesh = gen_mesh(N, Lx, Ly)
h_spacing = Lx / nx  # Assuming uniform spacing in both directions
A = gen_stiffmat(mesh)
A = (A * A) 

# Load vector
f = point_load(N, F, D)

A, f = applyBC!(A, mesh, nx, ny, f)

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
    
    du[N+1:end] .= (f - D * A * u[1:N]-0.0001*c_crit*u[N+1:end]) / (rho * h)
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
    surface(x, y, displacement[:, :, i]', title="Plate Deflection at t=$(round(times[i], digits=2)) s",
            xlabel="X-axis (m)", ylabel="Y-axis (m)", zlabel="Displacement (mm)",  
            zlims=zlims)
end

# Save animation
gif(anim, "dynamic_plate_point_underdamped.gif", fps=30)
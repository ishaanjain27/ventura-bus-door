######## CANTILEVER PLATE WITH SPRING BOUNDARY 

using StaticArrays 
using SparseArrays
using LinearAlgebra 
using BenchmarkTools
using Plots 

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
    k=10000        #k=0 means no spring so cantilever plate
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

# Solve for displacement
function solve_plate(A, f)
    return A \ f
end

# Uniform load function for 2D plate
function uniform_load(n, q, D)
    nx, ny = n
    rhs = fill(q, (nx + 1) * (ny + 1))
    return -rhs / D
end

# Parameters 
E = 70e9       # Young's modulus (Pa)
h = 0.01        # Plate thickness (m)
ν = 0.3         # Poisson's ratio
q = 1000.0      # Uniform load (N/m^2)
nx, ny = 100, 100   # Number of intervals (keep spacing equal)
Lx, Ly = 2.0, 2.  # Plate dimensions (m)
h_spacing = Lx/nx

D = (E * h^3) / (12 * (1 - ν^2))
# Generate mesh and matrices
N = (ny, nx)
mesh = gen_mesh(N, Lx, Ly)
A = gen_stiffmat(mesh)
A = (A * A) 
# function plot_sparsity(A::SparseMatrixCSC; markersize=3, color=:black)
#     rows, cols, _ = findnz(A)  # Get the indices of non-zero entries
#     scatter(
#         cols, rows,
#         color=color,
#         markersize=markersize,
#         xlabel="Columns",
#         ylabel="Rows",
#         title="Sparsity Pattern",
#         legend=false,
#         aspect_ratio=:equal,  # Ensures square scaling of the axes
#         yflip=true  # Flips the rows for matrix-style orientation
#     )
# end

# Load vector
f = uniform_load(N, q, D)
A, f = applyBC!(A, mesh, nx, ny, f, h_spacing, D)
A = A / h_spacing^4

# plot_sparsity_fixed(A, markersize=3, color=:black)

# Plot results
function plot_displacement(w, nx, ny, Lx, Ly)
    x = range(0, Lx, length=nx + 1)
    y = range(0, Ly, length=ny + 1)
    z = reshape(w, nx + 1, ny + 1)'
    surface(
        x, y, z,
        title="Plate Deflection",
        xlabel="X (m)", ylabel="Y (m)", zlabel="Displacement (mm)",
        aspect_ratio=Ly/Lx,  
        camera=(60, 30)
    )
    ####### if you want wireframe plot comment the above and uncomment the below
    # wireframe(x, y, z, title="Plate Deflection",
    # xlabel="X (m)", ylabel="Y (m)", zlabel="Displacement (mm)", 
    # linewidth=0.5, linecolor=:black, aspect_ratio=Ly/Lx) 
    
end

w = solve_plate(A, f)
w = w*10^3
println("The maximum displacement is: $(-minimum(w))")
plot_displacement(w, nx, ny, Lx, Ly)
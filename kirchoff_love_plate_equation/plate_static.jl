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
nx, ny = 50, 25   # Number of intervals (keep spacing equal)
Lx, Ly = 2.0, 1.  # Plate dimensions (m)
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
A, f = applyBC!(A, mesh, nx, ny, f)
A = A / h_spacing^4

# plot_sparsity_fixed(A, markersize=3, color=:black)

# Plot results
function plot_displacement(w, nx, ny, Lx, Ly)
    x = range(0, Lx, length=nx + 1)
    y = range(0, Ly, length=ny + 1)
    z = reshape(w, nx + 1, ny + 1)'
    # surface(
    #     x, y, z,
    #     title="Plate Deflection",
    #     xlabel="X (m)", ylabel="Y (m)", zlabel="Displacement (mm)",
    #     aspect_ratio=Ly/Lx,  
    #     camera=(60, 30)
    # )
    ####### if you want wireframe plot comment the above and uncomment the below
    wireframe(x, y, z, title="Plate Deflection",
    xlabel="X (m)", ylabel="Y (m)", zlabel="Displacement (mm)", 
    linewidth=0.5, linecolor=:black, aspect_ratio=Ly/Lx) 
    
end

w = solve_plate(A, f)
w = w*10^3
println("The maximum displacement is: $(-minimum(w))")
plot_displacement(w, nx, ny, Lx, Ly)
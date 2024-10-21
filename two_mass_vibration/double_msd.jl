using Plots

#Input params.
m1 = 1.0
c1 = 0.2
k1 = 5.0

m2 = 1.0
c2 = 0.2
k2 = 5.0

deltat = 0.01


#IC's
x1_0 = 2.0
x2_0 = 1.0
v1_0 = 0.0
v2_0 = 0.0
T = 20
t0 = 0.0


#Creating Vectors
u = [x1_0; v1_0; x2_0; v2_0]
A = [0.0 1.0 0.0 0.0;
     -(k1+k2)/m1 -(c1 +c2)/m1 k2/m1 c2/m1;
     0.0 0.0 0.0 1.0;
    k2/m2 c2/m2 -k2/m2 -c2/m2]

#Initial lists
disp1 = [x1_0]
vel1 = [v1_0]
disp2 = [x2_0]
vel2 = [v2_0]

time = [t0]

#Iteration
t = t0
while t <= T
	global u,t
	uprime = A * u
	u_2 = u + deltat * uprime
	t += deltat
	
	push!(disp1,u_2[1])
	push!(disp2,u_2[3])
	push!(vel1,u_2[2])
	push!(vel2,u_2[4])
	push!(time, t)
	u = u_2
end	

plot1 = plot(time,disp1,ylimits=(-3,3),xlimits=(t0,T),title = "disp1")
plot2 = plot(time,disp2,ylimits=(-3,3),xlimits=(t0,T), title = "disp2")

plot(plot1, plot2, layout = (1, 2), legend = true)

savefig("sol_double_msd")

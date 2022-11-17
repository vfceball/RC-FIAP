def getSaT(EQ, dt, T, xi, npts):
	from math import pi, sqrt
	import numpy as np
	g = 9.81


	# -----------------------------------
	# Inputs:
	# -----------------------------------
	# EQ:   Filename which is a single column file in units g (e.g "Eq.txt")
	# dt:   Time step in seconds (e.g 0.01)
	# xi:   Elastic damping (e.g 0.05)
	# T:    Period in seconds (e.g 1.0)

	# -----------------------------------
	# Outputs:
	#------------------------------------
	# Sa:	Sa(T,%) - Pseudo-Spectral Acceleration in g
	# Sv: 	Sv(T,%) - Pseudo-Spectral Velocity in m/s
	# Sd: 	Sd(T,%) - Spectral Displacement in m
	# pga: 	Peak ground acceelration in g

	# Import the ground motion
	# with open(EQ, "r") as myfile:
	# 	accg = np.array(myfile.readlines().split()).astype(float)

	with open(EQ) as fp:
		accg = np.array([p for l in fp for p in l.split()]).astype(float)

	# with open(EQ) as f:
	# 	accg = f.readlines()
	# accg = np.array(accg)
	# accg = accg.astype(np.float)


	# with open(EQ, 'r', newline='\r') as f:
	# 	accg = np.genfromtxt(f, unpack=True)

	# accg = np.loadtxt(EQ, usecols=range(5))

	if T == 0.0:
		Sa, Sv, Sd, pga, amax = 0, 0, 0, 0, 0
		pga = max(max(accg), -min(accg))
	else:
		gamma = 0.5		 # Newmark terms (Set for Average Acceleration method)
		beta = 0.25		# Newmark terms (Set for Average Acceleration method)
		ms = 1.0 		# Set the mass to 1kg
		acc = accg*g  # Change the units of the record to m/s2
		p = -ms*acc  # Create the force in

		# Calculate the initial values
		k = ms*(2*pi/T)**2		# Stiffness in N/m (which will give T following assumption of mass)
		w = sqrt(k/ms)			# Circular frequency
		c = 2*xi*ms*w			# Damping coefficient
		a, v, u = np.zeros(npts), np.zeros(npts), np.zeros(npts) 	# Initialise some vectors
		a[0] = p[0]/ms			# Initial acceleration in m/s2
		a1 = gamma*c/beta/dt+ms/beta/dt/dt
		a2 = ms/beta/dt+(gamma/beta-1)*c
		a3 = (1/2/beta-1)*ms+dt*(gamma/beta-1)*c
		k_bar = k+a1  # Stiffness term (see Chopra book)
		# print('ms =', ms)
		# print('c =', c)
		# print('k =', k)
		# print('a1 =', a1)
		# print('kbar =', k_bar)
		# print('p =',p)
		# Loop the terms (This loop is to length minus 1)
		for ind in range(npts-1):
			u[ind+1] = (p[ind+1] + a1*u[ind] + a2*v[ind] + a3*a[ind])/k_bar
			v[ind+1] = gamma/beta/dt*(u[ind+1] - u[ind]) + (1-gamma/beta)*v[ind] + dt*(1-gamma/2/beta)*a[ind]
			a[ind+1] = 1/beta/dt**2*(u[ind+1] - u[ind]) - 1/beta/dt*v[ind] - (1/2/beta-1)*a[ind]

		# Calculate Spectral Values
		amax = max(max(a), -min(a))/g
		Sd = max(max(u), -min(u))
		Sv = Sd*w
		Sa = Sd*w**2/g
		pga = max(max(accg), -min(accg))

	return Sa, Sv, Sd, pga, amax, accg

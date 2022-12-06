from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np

def heatbdn(xl,xr,yb,yt,M,N,f):
	D = 1.									# diffusion coefficient
	h = float(xr-xl)/M; k=float(yt-yb)/N; m=M+1; n=N
	sigma = D*k/(h*h)
	a  = np.diag(1+2*sigma*np.ones(m)) + np.diag(-sigma*np.ones(m-1),1) 
	a += np.diag(-sigma*np.ones(m-1),-1)			# define matrix a
	a[  0,:] = np.hstack([-3, 4, -1, np.zeros(m-3)]) 	# Neumann conditions
	a[m-1,:] = np.hstack([np.zeros(m-3), -1, 4, -3]) 
	xvals = np.linspace(xl,xr,M+1)
	tvals = np.linspace(yb,yt,N+1)	
	w = np.zeros((m,N+1))   					# 2nd index is time index
	w[:,0] = f(xvals)						# initial conditions
	for j in range(n):
		b = w[:,j].copy(); b[0]=0; b[m-1]=0
		w[:,j+1] = np.linalg.solve(a,b) 
	print(w)
	[x,t] = np.meshgrid(xvals,tvals)			# 3-D plot of solution w
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_wireframe(x, t, w.T, rstride=1, cstride=1)
	ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('w')
	plt.show()
	return w



# test_heatbdn.py

def f(x):
	return np.sin(2*np.pi*x)**2

np.set_printoptions(precision=14,linewidth=250)

w = heatbdn(0,1,0,0.2,4,5,f)
print(w)
## problem 1
import sympy as sp
 
Ax,Ay,theta_ddot,phi_ddot = sp.symbols('Ax Ay theta_ddot phi_ddot')
m,L,g = sp.symbols('m L g')
theta = 30*sp.pi/180

I_bar_com = 1/12*m*L**2
I_disk_com = 1/2*(2*m)*(L/2)**2
I_disk_edge = I_disk_com+(2*m)*(L/2)**2

#sum of moments about disk edge
eq1 = Ax*L/2 == I_disk_edge*phi_ddot
#sum of moments about com of bar
eq2 = -Ay*sp.sin(theta)*L/2+Ax*sp.cos(theta)*L/2==I_bar_com*theta_ddot
#sum of forces in x dir
eq3 = Ax == L/2*phi_ddot+L/2*theta_ddot*sp.cos(theta)
#sum of forces in y dir
eq4 = Ay-m*g == -L/2*theta_ddot*sp.sin(theta)


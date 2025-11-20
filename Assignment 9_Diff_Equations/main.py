import my_lib as mb
import matplotlib.pyplot as plt
import math
from tabulate import tabulate

#INUPT---------------
filename = 'input.txt'
with open(filename) as data:
    arr = mb.read_data(data)
mb.floater(arr)
# ------------------

#Output-----------
res_file = open("output.txt",'w')
Res = []
# -----------------


#QUESTION -1 : DAMPED,FORCED SHO

#2nd order ODE -> Two coupled first order ODE
w, g = 1.2, 0.2

y_0, v_0, t_0 = 2, -1, 0
def dydt(y,v,t):
    return v
def dvdt(y,v,t):
    return 0.5*(2*math.cos(w*t) - 2*y - g*v)

sol = mb.RK4_coup(dydt, dvdt, y_0, v_0, t_0, dt = 0.05, tf = 100)

table = zip(sol[2],sol[0])
# print(tabulate(table, headers=["t", "y"], floatfmt=".4f"))

# plt.plot(sol[2],sol[0],".",markersize = 2,label="Damped and Forced SHO")
# plt.ylabel("y(t)")
# plt.xlabel("t")
# plt.legend()
# plt.savefig("plot1.png")
# plt.show()

#QUESTION -2:

#Boundary value problem with SHOOTING Method

T_0, T_L, dx, L = 40, 200, 0.05, 10

def dydx(y,v,x):
    return v
def dvdx(y,v,x):
    return 0.01*(y - 20)

#Undershooting
un_sht = 10
Y,V,X = mb.RK4_coup(dydx, dvdx, T_0,un_sht,0,dx,L)
y_under = Y[int(L/dx)]
# print(y_under)
# plt.plot(X,Y,'.',markersize=0.5,label="undershoot")

#Overshooting
ov_sht = 15
Y,V,X = mb.RK4_coup(dydx, dvdx,T_0,ov_sht,0,dx,L)
y_over = Y[int(L/dx)]
# print(y_over)
# plt.plot(X,Y,'.',markersize=0.5,label="overshoot")

v_0 = mb.interpol(200,un_sht,y_under,ov_sht,y_over)
# print(interpol)

y,v,x = mb.RK4_coup(dydx, dvdx, T_0,v_0,0,dx,L)
# print(y[int(L/dx)])


table = zip(x,y)
print(tabulate(table, headers=["x", "y"], floatfmt=".4f"))
#
# plt.plot(x,y,'.',markersize=2,label="Numerical Solution")
# plt.ylabel("Temperature, T(x)")
# plt.xlabel("Position, x")
# plt.legend()
# plt.savefig('Q2_plot.png')
# plt.show()


#QUESTION-3: PDE - Temperature profile of a conducting Rod.

lx, nx = 2, 100 #divsion of x domain
hx = lx/nx
# setting initital temperature profile and domain of the problem
X,profile,T = [],[],[]
for i in range(nx+1):
    if i == nx/2:
        profile.append(300) #profile = temperature of rod
    else:
        profile.append(0)
    X.append(i*hx)  #domain

profile = mb.PDE(X,profile, lx = 2, lt = 5, nx = 100, nt = 50000)

table = zip(X,profile)
# print(tabulate(table, headers=["X", "Temperature"], floatfmt=".4f"))

#QUESTION-4 Eigen-Value

# arr[0] is the matrix taken from input file...
EV = mb.EV_dom(arr[0],e = 0.001)

print('Eigen Value (Dominant, upto 3 decimals) = ', EV[0])
print('Iterations = ', EV[1])
print('Eigen vector (Normalized) = ', EV[2])





















for i in range(len(Res)):
    for j in range(len(Res[i])):
        res_file.writelines((Res[i][j]))
res_file.close()

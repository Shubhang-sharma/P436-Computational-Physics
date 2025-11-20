import my_lib as mb
import math
import matplotlib.pyplot as plt

# SHO
mu = 0.25
w2 = 2
def dvdt(x,v,t):
    return -mu*v - w2*x
def dxdt(x,v,t):
    return v

X,Y,Z = mb.RK4_coup(dxdt,dvdt,1,0,0,0.1,50)
print(X[10])
plt.plot(Z,X,'.',markersize=2)

# A,B,C = mb.RK4_coup(dxdt,dvdt,1,10,0,0.2,120)
# print(A[10])

# plt.plot(C,A,'.',markersize=2)
plt.show()

# ax = plt.axes(projection='3d')
#
# # Data for a three-dimensional line
#
# ax.plot3D(Z,X,Y,'gray')
# plt.show()


# # Lorentz Attractor
# sg = 10
# rho = 28
# b = 8/3
# def dxdt(x,y,z,t):
#     return sg*(y-x)
# def dydt(x,y,z,t):
#     return x*(rho-z)-y
# def dzdt(x,y,z,t):
#     return x*y - b*z
#
# X, Y, Z,T = mb.RK4_coup_dir(dxdt,dydt,dzdt,1,0,0,0,0.02)
#
# # plt.plot(Z,X,'.',markersize=4)
# # plt.show()
#
#
#
# # ax = plt.axes(projection='3d')
# cm = plt.get_cmap("Blues")
# col = [cm(float(i)/(len(T)-1)) for i in range(len(T))]
# # ax.plot3D(X,Y,Z)
# # plt.show()\
#
# plt.subplot(111,projection = '3d')
# plt.scatter(X,Y,Z,c=col)
# plt.show()

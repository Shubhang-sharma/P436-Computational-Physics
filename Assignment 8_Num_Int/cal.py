import my_lib as mb
import numpy as np
import matplotlib.pyplot as plt

filename = 'input.txt'
with open(filename) as data:
    arr = mb.read_data(data)
mb.floater(arr) #floats the data which was read as strings before

# END OF INPUT----------------------------


#OUTPUT---------------------------#
res_file = open('lineq_output.txt','w')
Res = []

def f(x):
    return float(1/x)
def g(x):
    return x*np.cos(x)
# for i in range(4,21,4):
#     print('midpoint ',i,'parts',mb.quadrature_midpoint(g,0,1.57,i),'\n')
#     print('trap',i,'parts',mb.quadrature_trapezoidal(g,0,1.57,i),'\n')
#     print('simpson',i,'parts',mb.quadrature_simpson(g,0,1.57,i),'\n')

x = [x*0.1 for x in range(0,17)]
y,z = [],[]

def b(ele):
    return (ele-3)*np.cos(ele) + np.sin(ele)
def c(ele):
    return -ele*np.cos(ele)-2*np.sin(ele)
for ele in x:
    y.append(c(ele))
    z.append(b(ele))
    if c(ele) > c(ele+0.01) and c(ele) > c(ele-0.01):
        print(c(ele))
    if b(ele) > b(ele+0.01) and b(ele) > b(ele-0.01):
        print(b(ele))
print(min(y))

# plt.plot(x,y)
# plt.plot(x,z)
# plt.show()




print(0.69314718,mb.quadrature_simpson(f, 1,2,6))
print(np.pi/2 - 1 ,mb.quadrature_simpson(g,0,np.pi/2,4))

import my_lib as mb
import math
import matplotlib.pyplot as plt



def p(a,b):
    return 1/(b-a)
def f(x):
    return 4/(1+(x*x))

# mb.quadrature_monte(f,p,0,1)


def g(x,y):
    return y + math.exp(x)
def h(x):
    return -math.cos(x) + math.exp(x)

# X,Y = [],[]
# mb.forward_euler(g,0,0,X,Y,0.05)
# plt.plot(X,Y,'+')
# mb.func_plot(h,0,3)
# plt.show()
#
# def g(x):
#     return math.sin(x) + math.exp(x)
# def Dg(x):
#     return math.cos(x) + math.exp(x)
# def h(x):
#     return -math.cos(x) + math.exp(x)
#
# X,Y = [],[]
# mb.backward_euler(g,Dg,0,0,X,Y,0.05)
# plt.plot(X,Y,'+')
# mb.func_plot(h,0,3)
# plt.show()

def dydx(y,x):
    return y - x*x + 1
def func(x):
    return x*x + 2*x + 1 - 0.5*math.exp(x)
    # return x*x

X,Y=[],[]
mb.RK4(dydx,0,0.5,X,Y,0.4)
plt.plot(X,Y,'.')
mb.func_plot(func,0,12)
plt.show()

import math
e = 0.0001
d = 0.000001
#
def f(x):
    return float(x - 2*math.cos(x))
def Df(x):
    return float(1 + 2*math.sin(x))
def g(x):
    return float(3*x + math.sin(x) - math.exp(x))
def Dg(x):
    return float(3+math.cos(x)-math.exp(x))
def h(x):
    return float(math.cos(x)-x**3)
def Dh(x):
    return float(-math.sin(x)-3*x**2)
def k(x):
    return float(x*math.exp(x) -2)
def Dk(x):
    return float(x*math.exp(x) + math.exp(x))
def H(x):
    return float(math.sin(x)/x)

def p(x,n,A,i=0):
    return A[i]

def deflate(A,r):
    for i in range(1,len(A)):
        A[i] = r*A[i-1] + A[i]
    return A

def D(A):
    D = []
    d = len(A)
    for i in range(0,d-1):
        D.append(A[i]*(d-i-1))
    return D
def D2(A):
    return D(D(A))
#
poly = [1,-1,-7,1,6]
sols = [1,-2,3,-1]

def p(poly,x):
    n = len(poly)
    sum = 0
    for i in range(0,n):
        sum+= poly[i]*(x**(n-i-1))
    return sum
print(p(poly,0))

e = 0.00001
d = 0.00001
def lag(poly,roots,e,b):
    n = len(poly)-1
    if len(roots) == n:
        return roots

    if abs(p(poly,b)) < e: #initial guess bang on!
        roots.append(b)
        return lag(deflate(poly,b),roots,e,b)

    g = p(D(poly),b)/p(poly,b)
    h = g**2 - (p(D2(poly),b)/p(poly,b))

    if g>0:
        a = n/(g+(((n-1)*(n*h-g*g))**0.5))
    elif g<0:
        a = n/(g-(((n-1)*(n*h-g*g))**0.5))
    if abs(a) < e and abs(p(poly,a)) < 0.00001:
        roots.append(b)
        return lag(deflate(poly,b),roots,e,b)
    else:
        return lag(poly,roots,e,b-a)

roots=[]
print(lag(poly,roots,0.0001,-3))

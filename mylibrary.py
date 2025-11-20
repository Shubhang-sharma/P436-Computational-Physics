#Library File, Computational Physics, Even Semester 2022, NISER.
# Palla Shubhang Sharma, 2011106.

#Eigen Value problem

def EV_dom(A,x0,e = 0.001, MAX = 30):
    ev = matrix_mult(A,x0)
    den = matrix_mult(ev,x0)
    prev,iter = 0,0
    while True or k<MAX:
        iter+=1
        ev = matrix_mult(A,ev)
        num = matrix_mult(ev,x0)
        ev_dom = num/den
        s = 0
        for j in range(len(ev)):
            s+= (ev[j][0]**2)
            evf = []
        for j in range(len(ev)):
            evf.append(round(ev[j][0]/(s**0.5),3))
        if abs(prev - ev_dom) < e:
            return round(ev_dom,3),iter,evf
        den = num
        prev = ev_dom

#PDE
def PDE(X,profile,lx,lt,nx,nt):
    import matplotlib.pyplot as plt
    hx, ht = lx/nx,lt/nt
    a = ht/(hx*hx)

    for j in range(0,nt+1):
        if j <= 5000:
            T = [j*ht for i in range(nx+1)]
            if j == 0 or j == 10 or j == 20 or j == 50 or j == 100 or j == 200 or j == 500 or j == 1000:
                # plt.subplot(111,projection = '3d')
                # plt.plot(X,T,profile)
                plt.plot(X,profile,label="t = "+str(j*ht))
        for k in range(0,nx+1): # matrix multiplication A*Vj
            if k==1:
                profile[k] = profile[k+1]*a + profile[k]*(1-2*a)
            if k==nx:
                profile[k] = profile[k-1]*a + profile[k]*(1-2*a)
            elif k!=nx and k!=0:
                profile[k] = profile[k-1]*a + profile[k]*(1-2*a) + profile[k+1]*a
    plt.ylabel("Temperature, T(x)")
    plt.xlabel("Position, x")
    plt.legend()
    plt.savefig("Q3_plot.png")
    plt.show()
    return profile

#ODE But BVP
def interpol(y,l,yl,h,yh):
    return l + ((h-l)/(yh-yl))*(y-yl)

#Ordinary DifferentiSal Equations----------- Second OrD
def RK4_coup_dir(dxdt,dydt,dzdt,x,y,z,t,dt,):
    X,Y,Z,T = [],[],[],[]
    for i in range(2000):

        X.append(x)
        Y.append(y)
        Z.append(z)
        T.append(t)

        k1x = dt*dxdt(x,y,z,t);
        k1y = dt*dydt(x,y,z,t);
        k1z = dt*dzdt(x,y,z,t);

        k2x = dt*dxdt(x+k1x/2,y+k1y/2,z+k1z/2,t+dt/2);
        k2y = dt*dydt(x+k1x/2,y+k1y/2,z+k1z/2,t+dt/2);
        k2z = dt*dzdt(x+k1x/2,y+k1y/2,z+k1z/2,t+dt/2);

        k3x = dt*dxdt(x+k2x/2,y+k2y/2,z+k2z/2,t+dt/2);
        k3y = dt*dydt(x+k2x/2,y+k2y/2,z+k2z/2,t+dt/2);
        k3z = dt*dzdt(x+k2x/2,y+k2y/2,z+k2z/2,t+dt/2);

        k4x = dt*dxdt(x+k3x,y+k3y,z+k3z,t+dt);
        k4y = dt*dydt(x+k3x,y+k3y,z+k3z,t+dt);
        k4z = dt*dzdt(x+k3x,y+k3y,z+k3z,t+dt);

        x += (k1x + 2*k2x + 2*k3x + k4x)/6
        y += (k1y + 2*k2y + 2*k3y + k4y)/6
        z += (k1z + 2*k2z + 2*k3z + k4z)/6
        t += dt

    return X, Y, Z, T

def RK4_coup(dxdt,dvdt,x,v,t,dt,tf):
    xol,vol,tol = [],[],[]

    for i in range(int(tf/dt)+1):

        tol.append(t)
        xol.append(x)
        vol.append(v)

        k1x = dt*dxdt(x,v,t)
        k1v = dt*dvdt(x,v,t)

        k2x = dt*dxdt(x+k1x/2,v+k1v/2,t+dt/2)
        k2v = dt*dvdt(x+k1x/2,v+k1v/2,t+dt/2)

        k3x = dt*dxdt(x+k2x/2,v+k2v/2,t+dt/2)
        k3v = dt*dvdt(x+k2x/2,v+k2v/2,t+dt/2)

        k4x = dt*dxdt(x+k3x/2,v+k3v/2,t+dt/2)
        k4v = dt*dvdt(x+k3x/2,v+k3v/2,t+dt/2)

        x += (k1x + 2*k2x + 2*k3x + k4x)/6
        v += (k1v + 2*k2v + 2*k3v + k4v)/6

        t += dt

    return [xol,vol,tol]


def func_plot(func,rangeA,rangeB,s = 'func'):
    import math
    import matplotlib.pyplot as plt
    funcX= [0.001*ele for ele in range(1000*rangeA,1000*rangeB)]
    funcY = [func(ele) for ele in funcX]
    plt.plot(funcX,funcY,label = s)
    plt.legend()


def forward_euler(g,initx,inity,x,sol,h,n=0):
    sol.append(inity)
    x.append(initx)
    inity=inity+h*g(initx,inity)
    initx+=h
    n+=1
    if n==3/h:
        return [x,sol]
    return forward_euler(g,initx,inity,x,sol,h,n)

def backward_euler(g,Dg,initx,inity,x,sol,h,n=0):
    sol.append(inity)
    x.append(initx)
    # inity=g(initx)+h*g(initx+h)
    def NR(yxh):
        return yxh - g(initx) + h*g(yxh)
    def DNR(yxh):
        return 1 + h*Dg(yxh)
    inity = root_newton_raphson(NR,DNR,0,0.0001)[0]
    initx+=h
    n+=1
    if n==3/h:
        return [x,sol]
    return backward_euler(g,Dg,initx,inity,x,sol,h,n)

def RK4(f,initx,inity,h,xf,x=[],sol=[]):
    N = int(xf/h)
    for i in range(N+1):
        x.append(initx)
        sol.append(inity)
        k1 = h*f(inity,initx)
        k2 = h*f (inity + k1/2, initx + h/2)
        k3 = h*f (inity + k2/2, initx + h/2)
        k4 = h*f(inity + k3,initx)
        inity+=(k1 + 2*k2 + 2*k3 + k4)/6
        initx+=h

    return x,sol




#Quadrature
def monte(f,a,b,n):
    arrF,arrN,arrSG=[],[],[]
    N = 10
    while N<=n:
        X = lcg_pRNG(10,N)
        FN,FN2 = 0,0
        for i in range(len(X)):
            X[i] = a + (b-a)*X[i]
            FN += f(X[i])
            FN2 += f(X[i])**2
        SG = ((FN2/N) - (FN/N)**2)**0.5
        FN = FN*(b-a)/N
        arrSG.append(SG)
        arrF.append(FN)
        arrN.append(N)
        N += 30
    return FN,SG,arrF,arrSG,arrN

def simpson(f,a,b,N):
    h = float(b-a)/N
    sum = f(a)+f(b)
    for i in range(1,N):
        sum+=(2**((i%2)+1)*f(a+i*h))
    return (h/3)*sum

def midpoint(f,a,b,N):
    h = float((b-a)/N)
    sum = 0
    for n in range(1,N+1):
        sum += f((2*a + (2*n-1)*h)/2)
    return sum*h

def trapezoidal(f,a,b,N):
    h = (b-a)/N
    sum = f(a)+f(b)
    for i in range(1,N):
        sum += 2*f(a+i*h)
    return (h/2)*sum


#Data Fitting

#Lagrange Interpolation

def chisqr_fit_linear(x,y,sg = 0):
    if sg == 0:
        sg = [1 for i in range(len(x))]
    s,sx,sy,sxx,syy,sxy = 0,0,0,0,0,0
    for i in range(len(x)):
        s   += 1/(sg[i]**2)
        sx  += x[i]/(sg[i]**2)
        sxy += x[i]*y[i]/(sg[i]**2)
        sy  += y[i]/(sg[i]**2)
        syy += (y[i]**2)/(sg[i]**2)
        sxx += (x[i]**2)/(sg[i]**2)
    delta = (s*sxx) - (sx**2)
    a1 = (sxx*sy-sx*sxy)/(delta)
    a2 = (sxy*s-sx*sy)/(delta)
    sga1 = (sxx/delta)**0.5
    sga2 = (s/delta)**0.5
    pearson_coeff = ((sxy)**2)/(sxx*syy)
    return [a1,a2,pearson_coeff,sga1,sga2]

def polfit(x,y,k):
    X = [[] for i in range(k+1)]
    Y = []
    sumx,sumy = [],[]

    for i in range(0,2*k+1):
        sx,sy = 0,0
        for j in range(len(x)):
            sx = sx + (x[j]**i)
            if i<k+1:
                sy = sy + ((x[j]**i)*y[j])
        sumx.append(sx)
        sumy.append(sy)

    for i in range(k+1):
        for j in range(k+1):
            if i<=j:
                X[i].append(sumx[i+j])
            if i>j:
                X[i][j] = X[j][i]
        Y.append(sumy[i])
    sol = cholesky(X,Y)

    tol = []
    for i in range(len(sol)):
        tol.append(sol[len(sol)-i-1])
    print(tol)
    return tol

def pol(a,b,c,x):
    return a*x**2+b*x+c




def poly_fit_init(x,y,k): #x-data, y-data, k-degree of polynomial
    def sumik(x,k):
        add=0
        for i in range(len(x)):
            add+=(x[i])**k
        return add

    def sumxyik(x,y,k):
        add=0
        for i in range(len(x)):
            add+=((x[i])**k)*y[i]
        return add
    X,Y=[],[]
    for i in range(k+1):  # Filling X,Y matrix with 0s accordingly
        rowx=[]
        rowy=[0]
        X.append(rowx)
        Y.append(rowy)
        for j in range(k+1):
            rowx.append(0)
    for i in range(len(X)):  # Filling X matrix with required elements
        for j in range(len(X)):
            if i==j:
                X[i][j]=sumik(x,2*i)

            if i!=j and i>j:
                X[i][j]=sumik(x,i+j)
            X[j][i]=X[i][j]
    for i in range(len(Y)): # Filling Y matrix with required elements
        Y[i][0]=sumxyik(x,y,i)
    return X,Y

def poly_fit(x,y,k):
    X,Yr = poly_fit_init(x,y,k)
    Y = []
    for i in range(len(Yr)):
        Y.append(Yr[i][0])
    sol = cholesky(X,Y)
    return sol
#roots------------------

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

def root_lag(poly,roots,e,b):
    n = len(poly)-1
    if len(roots) == n:
        return roots

    if abs(pol(poly,b)) < e: #initial guess bang on!
        roots.append(b)
        return root_lag(deflate(poly,b),roots,e,b)

    g = pol(D(poly),b)/pol(poly,b)
    h = g**2 - (pol(D2(poly),b)/pol(poly,b))

    if g>0:
        a = n/(g+(((n-1)*(n*h-g*g))**0.5))
    elif g<0:
        a = n/(g-(((n-1)*(n*h-g*g))**0.5))
    if abs(a) < e and abs(pol(poly,a)) < 0.00001:
        roots.append(b)
        return root_lag(deflate(poly,b),roots,e,b)
    else:
        return root_lag(poly,roots,e,b-a)

def poly(poly,x): #takes polynomial coefficient and x as input and produces the polynommial value at 'x'
    n = len(poly)
    sum = 0
    for i in range(n):
        sum+= poly[i]*(x**(i))
    return sum

def root_bracket(f,a,b):
    a = float(a)
    b = float(b)
    if f(a)*f(b) < 0:
        return [a,b]
    if f(a)*f(b) > 0:
        if abs(f(a)) < abs(f(b)):
            a = a - 0.2*(b-a)
        if abs(f(a)) > abs(f(b)):
            b = b + 0.2*(b-a)
        return root_bracket(f,a,b)


def root_bisection(f,a,b,e,out,i=0):
    if i==0:
        out.append(str(i) + '  ['+str(round(a,6)) +','+ str(round(b,6)) +']'+ '\n')
    i+=1
    if f(a)*f(b)>0: #bracket first
        a,b = root_bracket(f,a,b)[0],root_bracket(f,a,b)[1]
    out.append(str(i) + '  ['+str(round(a,6)) +','+ str(round(b,6)) +']'+ '\n')
    if abs(b-a) < e and f(a) < 0.0001: #if already is a solution just return
        return [a,i]
    c = float((a+b)/2)
    if f(c)*f(a) < 0: #bisect the interval to give a new interval
        b = c
    elif f(c)*f(b) < 0:
        a = c
    return root_bisection(f,a,b,e,out,i) #repeat the whole process for the new interval


def root_regula_falsi(f,a,b,e,out,i=0):
    if i==0:
        out.append(str(i) + '  ['+str(round(a,6)) +','+ str(round(b,6)) +']'+'\n')
    i+=1
    if f(a)*f(b)>0: #bracket first
        a,b = root_bracket(f,a,b)[0],root_bracket(f,a,b)[1]
    out.append(str(i) + '  ['+str(round(a,6)) +','+ str(round(b,6)) +']'+'\n')
    c = b - (b-a)*f(b)/(f(b)-f(a))
    if abs(b-a) < e and f(a) < 0.0001:
        return [a,i]
    elif abs(b-c) != 0 and abs(b-c) < e:
        return [c,i]
    elif abs(a-c) != 0 and abs(a-c) < e:
        return [c,i]
    if f(c)*f(a) < 0:
        b = c
    elif f(c)*f(b) < 0:
        a = c
    return root_regula_falsi(f,a,b,e,out,i)

def root_newton_raphson(f,Df,init,e,i=0,d=0.00001):
    init = float(init)
    temp = init
    init = init - (f(init)/Df(init))
    if abs(init-temp) < e and f(temp) < d:
        return temp,i
    i+=1
    return root_newton_raphson(f, Df, init, e, i)

#Lin-eq
#Gauss-seidal #Iterative Method--------------------
def seidel(A,B):
    iter = [0 for i in range(len(A))] #initial guess
    while True:
        c = 0  # converegence precision check artefact
        for i in range(len(A)): 
            sum = 0
            for j in range(len(A)):
                if j!=i:
                    sum += A[i][j]*iter[j]
            if  abs((1/A[i][i])*(B[i]-sum) - iter[i]) < 0.0000001: # convergence precision check with old iter
                c = 1
            else:
                c = 0
            iter[i] = round((1/A[i][i])*(B[i]-sum),7)  #result of a single iteration is stored in iter, i.e old elements replaced
        if c == 1: # convergence precision upto the given value and thus green signal to return teh solution
            return iter

#Jacobi #Iterative Method
def jacobi(A,B):
    iter = [0 for i in range(len(A))]
    while True:
        c = 0
        temp = iter #temp = old iter
        for i in range(len(A)):
            sum = 0
            for j in range(len(A)):
                if j!=i:
                    sum += A[i][j]*temp[j]
            iter[i] = (1/A[i][i])*(B[i]-sum)
            if  abs(temp[i] - iter[i]) < 0.0000001: #comparing temp i.e old iter with iter (new) #convergence precision check
                c = 1
            else:
                c = 0
        if c == 1: #green signal to return on checking the precision
            return iter

def rearrange_diag_dom(A,B):
    for i in range(len(A)):
        temp = 0
        for j in range(len(A[i])):
            if A[i][j] > A[i][temp]: #Searching the largest element, since it is the only prospective element to look at (satisfying the conditions for diagonal dominance)
                temp = j
        if temp != i: #swap if the row is in wrong position
            p = B[i]
            B[i] = B[temp]
            B[temp] = p
            for k in range(len(A[i])): #swapping with the row it should belong so that matrix becomes diagonal dominant
                t = A[i][k]
                A[i][k] = A[temp][k]
                A[temp][k] = t

    return A, B

def if_symmetric(A): #returns 1 if symmetric, i.e True
    c = 1
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i][j] != A[j][i]:
                c = 0
    return c

def cholesky(A,B):
    for i in range(len(A[0])):
        s = 0
        for k in range(i):
            s = s + (A[i][k]**2)
        A[i][i] = round((A[i][i] - s)**0.5,3)
        for j in range(i+1,len(A[0])):
            s = 0
            for k in range(i):
                s = s + (A[i][k]*A[k][j])
            A[j][i] = round((1/A[i][i])*(A[j][i] - s),3)
            A[i][j] = A[j][i]           # stores L and U in A which is essentially symmetric
    #forward sub result y
    y = []
    for i in range(0,len(A)):
        s = 0
        for j in range(i):
            s = s + A[i][j]*y[j]
        y.append((B[i] - s)/A[i][i])
    #to get the backward sub result x
    x = [round(y[len(A)-1]/A[len(A)-1][len(A)-1],3)]
    for i in reversed(range(len(A)-1)):
        s = 0
        for j in range(i+1,len(A)):
            s = s + A[i][j]*x[j-i-1]
        x.insert(0,round((1/A[i][i])*(y[i] - s),3))
    return x


def ludecomp(A,B):
    for j in range(0,len(A)):
        for i in range(1,j+1): #Decomposition-L and U storing in A
        # with diagonal elements as that of U. since elements of L are esentially known and = 1
            s = 0
            for k in range(0,i-1):
                s = s + A[i][k]*A[k][j]
            A[i][j] = A[i][j] - s
        for i in range(j+1,len(A)):
            t = 0
            for k in range(0,j-1):
                t = t + A[i][k]*A[k][j]
            A[i][j] = round((A[i][j] - t)/A[j][j],3)

    #forward-backward sub result
    y = []
    for i in range(0,len(A)):
        s = 0
        for j in range(i):
            s = s + A[i][j]*y[j]
        y.append((B[i] - s)/A[i][i])
    #backward sub
    x = [round(y[len(A)-1]/A[len(A)-1][len(A)-1],3)]
    for i in reversed(range(len(A)-1)):
        s = 0
        for j in range(i+1,len(A)):
            s = s + A[i][j]*x[j-i-1]
        x.insert(0,round((1/A[i][i])*(y[i] - s),3))
    return x

def gauss(M):
    for i in range(0,len(M)):
        if M[i][i] == 0:
            #-------------find the largest leading term......
            temp = i
            for j in range(i+1,len(M)):
                if abs(M[j][i]) > abs(M[temp][i]):
                    row_large = j
                    temp = row_large
                else:
                    row_large = temp
            #---------------- swap the rows.....
            for k in range(0,len(M[0])):
                t = M[row_large][k]
                M[row_large][k] = M[i][k]
                M[i][k] = t
            # ----------------------- time to make the leading term 1.....
        d = M[i][i]
        if d == 0:
            return
        for l in range (0,len(M[0])):
            M[i][l] = float(M[i][l])/float(d)
            # ----------------------------- substracting the rows to make'em 0......
        for m in range(0,len(M)):
            if m != i :
                s = M[m][i]
                for n in range(0,len(M[0])):
                    M[m][n] = round(M[m][n] - s*M[i][n],3)
    return M


#Random Number Generators
seed = 10

def rwalk(N):

    import matplotlib.pyplot as plt
    import numpy as np
    x,y = [],[]
    for i in range(N):
        x.append(myrand())
        y.append(myrand())

    X, Y = [0], [0]

    rm = 0
    for i in range(N):
        x[i] = 2*x[i] - 1 # numbers now in -1 to 1 range
        y[i] = 2*y[i] - 1

        X.append(X[i]+x[i]) #Simulating walk via incrementing the cordinates by random numbers(already generated).
        Y.append(Y[i]+y[i])

        rm = rm + ((X[i+1]-X[i])**2 + (Y[i+1]-Y[i])**2)
    rms = round((rm/N)**0.5,3) #RMS

    s = round(((X[N])**2 + (Y[N])**2)**0.5,3) #Displacement
    #Adding details to the plot itslef
    plt.title("N =" +  str(N) + "\n" + "Displacement =" + " " + str(s) + " " + "RMS =" + " " + str(rms))
    plt.plot(X,Y,lw = 0.7)
    plt.scatter(X[0],Y[0],c = 'magenta',label='Start')
    plt.scatter(X[-1],Y[-1],c = 'purple',label='Finish')
    plt.legend()
    plt.grid(True)


def myrand(a=1103515245,c = 12345,m = 32768): # returns random numbers in range(0,1)
    global seed
    seed = ((a*seed + c)%m)/m
    return seed

def matrix_mult(A,B,M=[]):
    R = []
    S = []
    if len(A[0]) == len(B): #The form (mxn)(nxk) is only allowed.
        for i in range(len(A)):
            s = 0
            for k in range(len(B[0])):
                for j in range(len(B)):
                    s = round(s + (A[i][j]*B[j][k]),2) #the definition of matrix multiplcation
                S.append(s)
                s = 0
            R.append(S)
            # R.append('\n')
            S = []
        M.append(R)
        M.append('\n')

        return(R)
    elif len(A) == len(B) and len(A[0]) ==1 and len(B[0]) == 1: #Dot product, i.e. coloumn matrix
        DP = 0
        for i in range(len(A)):
            DP = round(DP + (A[i][0]*B[i][0]),2) #Definition of dot product#
        M.append(str(DP) + '\n')
        M.append('\n')

        return(DP)
def dot_product(a,b):
    DP = 0
    if len(a) == len(b):
        for i in range(len(a)):
            DP += round(a[i]*b[i],3)
    return DP


# Input/Output-----------
def read_data(data):
    A,B,C,D,E,F,G,H,I,J,K,L = [],[],[],[],[],[],[],[],[],[],[],[]
    arr = [A,B,C,D,E,F,G,H,I,J,K,L] # arr is the list which stores the input data
    j = 0
    for line in data:
        line = line.strip().split() #stripping line by line and splitting based on whitespace
        # print(line)
        if len(line) != 0 :
            arr[j].append(line)
        elif len(line) == 0 :
            j += 1
    return arr

def print_array(M):
    for i in range(0,len(M)):
        print(M[i])
        print()

#floats all the strings read from input file
def floater(data_list):
    for i in range(len(data_list)):
        if len(data_list[i])!=0:
            for j in range(len(data_list[i])):
                for k in range(len(data_list[i][j])):
                    data_list[i][j][k] = round(float(data_list[i][j][k]),4)
#takes the table outputs two lists
def table_input(A):
    x,y=[],[]
    for i in range(len(A)):
            x.append(A[i][0])
            y.append(A[i][1])
    return x, y

def stringer(data_list): #STRINGS ALL THE ELEMENTS OF A 1D LIST
    for i in range(len(data_list)):
        data_list[i] = str(data_list[i]) + '  '
    return data_list

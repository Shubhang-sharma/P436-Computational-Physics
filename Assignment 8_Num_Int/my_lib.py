
# Quadrature ~ Numerical Integration Methods

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
# --------------------------------------------------------------------------
#DATA FITTINGS

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

# Polynomial fit
def polfit(x,y,k):
    X = [[] for i in range(k+1)]
    Y = []
    sumx,sumy = [],[]
    for i in range(0,2*k+1):
        sx,sy = 0,0
        for j in range(len(x)):
            sx = sx + (x[j]**i)
            sy = sy + ((x[j]**i)*y[j])
        sumx.append(sx)
        sumy.append(sy)
    for i in range(k+1):
        for j in range(k+1):
            X[i].append(sumx[i+j])
        Y.append(sumy[i])
    M = [X,Y]
    return M
# --------------------------------------------------------------------------------------------------------------------------

#roots------------------

def p(poly,x): #takes polynomial coefficient and x as input and produces the polynommial value at 'x'
    n = len(poly)
    sum = 0
    for i in range(0,n):
        sum+= poly[i]*(x**(n-i-1))
    return sum

def bracket(f,a,b):
    a = float(a)
    b = float(b)
    if f(a)*f(b) < 0:
        A.append(a)
        A.append(b)
        return A
    if f(a)*f(b) > 0:
        if abs(f(a)) < abs(f(b)):
            a = a - 1.5*(b-a)
        if abs(f(a)) > abs(f(b)):
            b = b + 1.5*(b-a)
        return bracket(f,a,b)

def bisect(f,a,b):
    if abs(b-a) < e and f(a) < d:
        return a
    c = float((a+b)/2)
    if f(c)*f(a) < 0:
        b = c
    elif f(c)*f(b) < 0:
        a = c
    return bisect(f,a,b)

def regula_falsi(f,a,b):
    a = float(a)
    b = float(b)
    c = b - (b-a)*f(b)/(f(b)-f(a))
    if abs(b-a) < e and f(a) < d:
        return a
    elif abs(b-c) != 0 and abs(b-c) < e:
        return c
    elif abs(a-c) != 0 and abs(a-c) < e:
        return c
    if f(c)*f(a) < 0:
        b = c
    elif f(c)*f(b) < 0:
        a = c
    return regula_falsi(f,a,b)

def newton_raphson(f,Df,init,e,d,i=0):
    init = float(init)
    temp = init
    init = init - (f(init)/Df(init))
    if abs(init-temp) < e and f(temp) < d:
        return temp,i
    i+=1
    return newton_raphson(f, Df, init,e,d,i)

# def lag(poly,roots,i=0,b=0,e,d):
#     i+=1
#     n = len(poly)-1
#     g = p(D(poly),b)/p(poly,b)
#     h = g**2 - (p(D2(poly),b)/p(poly,b))
#
#     if g>0:
#         a = n/(g+(((n-1)*(n*h-g*g))**0.5))
#     elif g<0:
#         a = n/(g-(((n-1)*(n*h-g*g))**0.5))
#
#     if abs(p(poly,b)) < e:
#         roots.append(b)
#         lag(deflate(poly,b),roots,i)
#     elif abs((b-a)) < e and abs(p(poly,a)) < d:
#         roots.append(a)
#         lag(deflate(poly,b),roots,i)
#     else:
#         lag(poly,roots,i,a)

#LINEAR EQUATIONS----------
def if_symmetric(A):
    c = 1
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i][j] != A[j][i]:
                c = 0
    return c
#Gauss Jordan Elimination
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
    # --------------------
# LU Decomposition
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

# Cholesky Decomposition
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

#Gauss-seidal #Iterative Method--------------------
def seidel(A,B,iter,p,ITERN=0):
    while True:
        c = 0  # converegence precision check artefact
        for i in range(len(A)):
            sum = 0
            for j in range(len(A)):
                if j!=i:
                    sum += A[i][j]*iter[j]
            if  abs((1/A[i][i])*(B[i]-sum) - iter[i]) < 10**(-p): # convergence precision check with old iter
                c+=1
            iter[i] = round((1/A[i][i])*(B[i]-sum),p)  #result of a single iteration is stored in iter, i.e old elements replaced
        ITERN+=1
        if c == len(A): # convergence precision upto the given value and thus green signal to return teh solution
            return iter,ITERN

#Jacobi #Iterative Method

def jacobi(A,B,init,ITERN=0):
    while True:
        c=0
        init_new = []
        for i in range(len(A)):
            sum = 0
            for j in range(len(A)):
                if j!=i:
                    sum += A[i][j]*init[j]
            if  abs((1/A[i][i])*(B[i]-sum) - init[i]) <= 0.000001: #comparing temp i.e old iter with iter (new) #convergence precision check
                c+=1
            init_new.append((1/A[i][i])*(B[i]-sum))
        init = init_new
        ITERN+=1
        if c == len(A): #green signal to return on checking the precision
            return init_new,ITERN

def rearrange_diag_dom(A,B):
    for i in range(len(A)):
        temp = 0
        for j in range(len(A[i])):
            if A[i][j] > A[i][temp]: #Searching the largest element, since it is the only prospective element to look at (satisfying the conditions for diagonal dominance)
                temp = j
        if temp != i: #swap if the row is in wrong position
            for k in range(len(A[i])): #swapping with the row it should belong so that matrix becomes diagonal dominant
                t,s = A[i][k],B[i]
                A[i][k],B[i] = A[temp][k],B[temp]
                A[temp][k],B[temp] = t,s
    return A
#---------------------------
#Random Number Generators-----------------------
def quicky_pRNG(seed,num,c=4.1):
    #c = 3.5, 4.01#
    nums = []
    x = seed
    for i in range(num):
        nums.append(x)
        x = c*x*(1-x)
    return nums

def lcg_pRNG(seed,num,a=1103515245,c = 12345,m = 32768): # returns random numbers in range(0,1)
    x = seed
    Y = []
    for i in range(num+1):
        Y.append(x/m)
        x = ((a*x + c)%m)
    return Y







#WARMUP FUNCTIONS----------------------------
def sumON(N,R):
    s = 0
    for i in range(int(N)):
        s = s + (2*i + 1) # iterative addition
    R.append(str(s) + '\n')
    R.append('\n')
    return(s)
# Factorial of a natural number "N" .
def fac(N,R):
    s = 1
    for i in range(1,int(N)+1):
        s = s*i # iterative product
    R.append(str(s) + '\n')
    R.append('\n')
    return(s)

# sum of terms in "AP"
def sumAP(init,n,d,R): #input = initial value, "n" for sum upto, common difference
    a = init
    s = init
    for i in range(0,int(n) -1):
        a = a + d
        s = s + a
    R.append(str(s) + '\n')
    R.append('\n')
    return(s)
# Sum of terms in "GP"
def sumGP(init,n,r,R): #input = initial value, "n" for sum upto, common ratio
    a = init
    s = init
    for i in range(0,int(n)-1):
        a = a*r
        s = s + a
    R.append(str(s) + '\n')
    R.append('\n')
    return(s)
# Sum of terms in "HP"
def sumHP(init,n,d,R): #input = initial value, "n" for sum upto, common difference
    a = init
    s = init
    for i in range(0,int(n)-1):
        a = 1/((1/a) + d)
        s = round(s + a,3)
    R.append(str(s) + '\n')
    R.append('\n')
    return(s)
# Series Sum
def sumSERIES(b,R): #takes input to specify b - accurate upto b decimal places.
    s = 0
    n = 1
    X = []
    Y = []
    while True:
        if (round((((-1)**(n+1))/(2**n) + s),int(b)) - round(s,int(b))) != 0:
            s = s + ((-1)**(n+1))/(2**n)
            Y.append(s)
            X.append(n)
            n += 1
        else:
            break
    import matplotlib.pyplot as mp
    mp.plot(X,Y)
    mp.show()
    R.append(str(round(s,4)) + '\n')
    R.append('\n')
    return(s)

def matrix_mult(A,B,M):
    R = []
    S = []
    if len(A[0]) == len(B): #The form (mxn)(nxk) is only allowed.
        for i in range(len(A)):
            s = 0
            for k in range(len(B[0])):
                for j in range(len(B)):
                    s = round(s + (A[i][j]*B[j][k]),2) #the definition of matrix multiplcation
                S.append(str(s) + "   ")
                s = 0
            R.append(S)
            R.append('\n')
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

#CLass
class myComplex:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def sum(A,B,M): #takes two complex numbers and an output list as input#
        R = str((A.x + B.x)) +' '+ '+' +' '+ str((A.y + B.y)) + 'i'
        M.append(str(R) + '\n')
        M.append('\n')

        return(R)
    def mult(A,B,M):
        R = str((A.x*B.x)-(A.y*B.y)) +' '+ '+' +' '+ str((A.x*B.y)+(A.y*B.x)) + 'i'
        M.append(str(R) + '\n')
        M.append('\n')

        return(R)
    def mod(A,M):
        R = round((A.x**2 + A.y**2)**(0.5),3)
        M.append(str(R) + '\n')
        M.append('\n')
        return(R)

#INPUT AND OUTPUT
def read_data(data):
    arr = [[] for i in range(11)] # arr is the list which stores the input data
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
                    data_list[i][j][k] = round(float(data_list[i][j][k]),2)
def stringer(data_list): #STRINGS ALL THE ELEMENTS OF A 1D LIST
    for i in range(len(data_list)):
        data_list[i] = str(data_list[i]) + '  '
    return data_list

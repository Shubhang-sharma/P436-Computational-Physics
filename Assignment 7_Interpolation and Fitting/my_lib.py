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
            iter[i] = round(    (1/A[i][i])*(B[i]-sum),7)  #result of a single iteration is stored in iter, i.e old elements replaced
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

def rearrange_diag_dom(A):
    for i in range(len(A)):
        temp = 0
        for j in range(len(A[i])):
            if A[i][j] > A[i][temp]: #Searching the largest element, since it is the only prospective element to look at (satisfying the conditions for diagonal dominance)
                temp = j
        if temp != i: #swap if the row is in wrong position
            for k in range(len(A[i])): #swapping with the row it should belong so that matrix becomes diagonal dominant
                t = A[i][k]
                A[i][k] = A[temp][k]
                A[temp][k] = t
    return A
#---------------------------
#FUNCTIONS - INPUT AND OUTPUT
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

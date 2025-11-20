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
# Sum of first N "Odd" Natural Numbers
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

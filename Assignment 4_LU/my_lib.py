
def if_symmetric(A):
    c = 1
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i][j] != A[j][i]:
                c = 0
    return c

# LU Decomposition
def ludecomp(A,B):
    for j in range(0,len(A)):
        for i in reversed(range(2,j+1)): #Decomposition-L and U storing in A
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
    y = [B[0]]
    for i in range(1,len(A)):
        s = 0
        for j in range(i):
            s = s + A[i][j]*y[j]
        y.append(B[i] - s)
    x = [round(y[len(A)-1]/A[len(A)-1][len(A)-1],3)]
    for i in reversed(range(len(A)-1)):
        s = 0
        for j in range(i+1,len(A)):
            s = s + A[i][j]*x[j-i-1]
        x.insert(0,round((1/A[i][i])*(y[i] - s),3))
    print('solution',x)

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
            A[i][j] = A[j][i]
    y = []
    for i in range(0,len(A)):
        s = 0
        for j in range(i):
            s = s + A[i][j]*y[j]
        y.append(B[i] - s)
    x = [round(y[len(A)-1]/A[len(A)-1][len(A)-1],3)]
    for i in reversed(range(len(A)-1)):
        s = 0
        for j in range(i+1,len(A)):
            s = s + A[i][j]*x[j-i-1]
        x.insert(0,round((1/A[i][i])*(y[i] - s),3))
    print('solution',x)
    print()

def print_array(M):
    for i in range(0,len(M)):
        print(M[i])
        print()

def floater(data_list):
    for i in range(len(data_list)):
        if len(data_list[i])!=0:
            for j in range(len(data_list[i])):
                for k in range(len(data_list[i][j])):
                    data_list[i][j][k] = round(float(data_list[i][j][k]),2)

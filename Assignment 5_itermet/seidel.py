import my_lib as mb



A = [[4,1,-1,1],[1,4,-1,-1],[-1,-1,5,1],[1,-1,1,3]]
B = [-2,-1,0,1]
def seidel(A,B):
    iter = [0,0,0,0]
    norm = 0
    while [ele < 0.001 for ele in NM]:
        for i in range(len(A)):
            sumI,sumP = 0,0
            for j in range(len(A)):
                if i!=j:
                    sumP += A[i][j]*iter[j]
            nV = iter[i]
            iter[i] = (1/A[i][i])*(B[i]-sumI-sumP)
            NM.append(nv-iter[i])
    return iter

print(seidel(A,B))

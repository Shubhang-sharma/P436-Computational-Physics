def parray(M):
    for i in range(0,4):
        print(M[i])
        print()


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
        for l in range (0,5):
            M[i][l] = float(M[i][l])/float(d)
            # ----------------------------- substracting the rows to make'em 0......
        for m in range(0,4):
            if m != i :
                s = M[m][i]
                for n in range(0,5):
                    M[m][n] = M[m][n] - s*M[i][n]
    # --------------------

import my_lib as mb
import matplotlib.pyplot as mp

#INUPT---------------
filename = 'input.txt'
with open(filename) as data:
    arr = mb.read_data(data)
mb.floater(arr) #floats the data which was read as strings before

# END OF INPUT----------------------------



#OUTPUT---------------------------#
res_file = open('warmup_output.txt','w')
Res = []
#----------------------------------

#QUESTION 1--------------------------#

mb.sumON(arr[0][0][0],Res)

mb.fac(arr[1][0][0],Res)


#QUESTION 2-------------------------#
data = arr[2][0]
mb.sumAP(data[0],data[1],data[2],Res)
data = arr[3][0]
mb.sumGP(data[0],data[1],data[2],Res)
data = arr[4][0]
mb.sumHP(data[0],data[1],data[2],Res)

#QUESTION 3-----------------------#
#Prints accurate upto required number of decimal places
mb.sumSERIES(arr[5][0][0],Res) #takes input to plot sum upto "n" terms.

#QUESTION 4------------------------#


mb.matrix_mult(arr[6], arr[7],Res) #AB

mb.matrix_mult(arr[8], arr[9],Res) #D.C

mb.matrix_mult(arr[7], arr[8],Res) #BC

#QUESTION 5--------------------------------------#

Z = mb.myComplex(arr[10][0][0],arr[10][0][1])
W = mb.myComplex(arr[11][0][0],arr[11][0][1])

mb.myComplex.sum(Z,W,Res)

mb.myComplex.mult(Z,W,Res)

mb.myComplex.mod(Z,Res)

mb.myComplex.mod(W,Res)

#OUTPUT-----------------------------------------#
for i in range(len(Res)):
    for j in range(len(Res[i])):
        res_file.writelines((Res[i][j]))
res_file.close()

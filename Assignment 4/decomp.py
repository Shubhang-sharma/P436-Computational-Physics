import my_lib as  mb
import math
# Cholesky Decomposition

file = 'input.txt'
with open(file) as data:
    A,B = [],[]
    arr = [A,B] # arr is the list which stores the input data
    j = 0
    for line in data:
        line = line.strip().split() #stripping line by line and splitting based on whitespace
        # print(line)
        if len(line) != 0 :
            arr[j].append(line)
        elif len(line) == 0 :
            j += 1
mb.floater(arr)


M = arr[0]
N = arr[1][0]
mb.print_array(M)

print()

if mb.if_symmetric(M):
    mb.cholesky(M,N)

mb.print_array(M)

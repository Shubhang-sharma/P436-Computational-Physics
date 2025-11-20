import my_lib as mb
import matplotlib.pyplot as plt

#INUPT---------------
filename = 'input.txt'
with open(filename) as data:
    arr = mb.read_data(data)
mb.floater(arr) #floats the data which was read as strings before

# END OF INPUT----------------------------


#OUTPUT---------------------------#
res_file = open('lineq_output.txt','w')
Res = []
#----------------------------------
x = arr[0][0]
y = arr[1][0]

X = mb.polfit(x,y,2)

for i in range(len(X[1])): #Construction of Augmented Matrix
    X[0][i].append(X[1][i])
S = mb.gauss(X[0])
mb.print_array(S)
sol = []
for i in range(len(S)): #Extraction of solution
    sol.append(S[i][len(S[i])-1])
print('The least square polynomial fit coefficients are ',sol)
def f(x):
    a = sol[0]
    b = sol[1]
    c = sol[2]
    return (c*x*x) + (b*x) + a
func = [f(x) for x in range(0,101,1)]

fx = [x for x in range(0,101,1)]
plt.plot(fx,func)
plt.plot(x,y,'.')

for i in range(0,101,1):
    if f(i) > f(i + 0.1) and f(i) > f(i - 0.1):
        hx = i
        hy = f(i)
plt.plot(hx,hy,'bo',markersize=4)
print('max at r =',hx,'with h = ',hy)
plt.show()

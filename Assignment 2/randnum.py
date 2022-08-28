import my_lib as mb
import matplotlib.pyplot as plt
filename = "input.txt"
with open(filename) as data:
    arr = mb.read_data(data)
mb.floater(arr)

res_file = open('randnum_output.txt','w')
Res = []
# Question 1

Res.append(mb.stringer(mb.my_lcg(arr[0][0][0],int(arr[0][0][1]),0)))
Res.append('\n \n')

# Question 2
Res.append(str(mb.vol_sphere()) + '\n')

# Question 3
mb.rwalk(1,5,300)
plt.show()
mb.rwalk(2,8,600)
plt.show()
mb.rwalk(7,3,900)
plt.show()

for i in range(len(Res)):
    for j in range(len(Res[i])):
        res_file.writelines((Res[i][j]))
res_file.close()

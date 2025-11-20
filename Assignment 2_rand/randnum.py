import my_lib as mb
import matplotlib.pyplot as plt
filename = "input.txt"
with open(filename) as data:
    arr = mb.read_data(data)
mb.floater(arr)

res_file = open('randnum_output.txt','w')
Res = []
# Question 1
d = arr[0][0]
Res.append(mb.stringer(mb.my_lcg(d[0],int(d[1]))))
Res.append('\n \n')

# Question 2
Res.append(str(mb.vol_sphere()) + '\n')

# Question 3
d = arr[1][0]
mb.rwalk_lcg(d[0],d[1],int(d[2])) #takes two seeds  and Steps as input for the random walk
plt.savefig('rwalk_300.png')
d = arr[2][0]
mb.rwalk_lcg(d[0],d[1],int(d[2]))
plt.savefig('rwalk_600.png')
d = arr[3][0]
mb.rwalk_lcg(d[0],d[1],int(d[2]))
plt.savefig('rwalk_900.png')


for i in range(len(Res)):
    for j in range(len(Res[i])):
        res_file.writelines((Res[i][j]))
res_file.close()

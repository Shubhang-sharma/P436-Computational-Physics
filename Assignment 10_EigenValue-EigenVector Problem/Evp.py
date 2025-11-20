import my_lib as mb

#INUPT---------------
filename = 'input.txt'
with open(filename) as data:
    arr = mb.read_data(data)
mb.floater(arr) #floats the data which was read as strings before


A = arr[0]
x0 = arr[1]



def EV_dom(A,x0 = [1,2,3],MAX = 30):
    ev = mb.matrix_mult(A,x0)
    den = mb.matrix_mult(ev,x0)
    prev,iter = 0,0
    while True or k<MAX:
        ev = mb.matrix_mult(A,ev)
        num = mb.matrix_mult(ev,x0)
        ev_dom = num/den
        if abs(prev - ev_dom) < e:
            return ev_dom,iter,ev
        den = num
        prev = ev_dom
        iter+=1

print(EV_dom(A))

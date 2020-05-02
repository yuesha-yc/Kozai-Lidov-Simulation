def numsum(x):
    for i in range(0,x):
        if i == 0:
            print(0)
        else:
            print(str(i+i-1))

n = int(input("Input number"))
numsum(n)
import random

def run(t,n=10):
    for i in range(t):
        if random.random() < 0.5:
            n -= 1
        else:
            n += 1
        if n == 0:
            return 0
        if n == 11:
            return 11
    return n
        

import sys
import matplotlib.pyplot as plt

with open(sys.argv[1], 'r') as fin:
    line = fin.readlines()[int(sys.argv[2])]
    numbers = [float(s) for s in line.strip().split(' ')]

    plt.hist(numbers, 10000, normed=10, facecolor='green', alpha=0.75)
    plt.show()
    #plt.plot()

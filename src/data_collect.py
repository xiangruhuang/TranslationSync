import sys

def parse(filename):
    fin = open(filename, 'r')
    lines = fin.readlines()

    iters = []
    times = []
    loss = []

    for line in lines:
        tokens = line.split(', ')
        for token in tokens:
            lr = token.split('=')
            l = lr[0]
            r = lr[1]
            if l.startswith('iter'):
                iters.append(int(r))
            if l.startswith('elapsed_time'):
                times.append(float(r))
            if l.startswith('min_loss'):
                loss.append(float(r))

    fin.close()

    return iters, times, loss

iters, times, loss = parse(sys.argv[1])
print(loss)

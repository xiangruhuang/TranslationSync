import sys

with open(sys.argv[1], 'r') as fin:
    lines = fin.readlines()
    strings = {0:[], 1:[], 2:[], 3:[]}
    for count, line in enumerate(lines):
        print line
        a = [float(token.split('=')[-1].strip().lstrip('(').rstrip(')')) for token in line.split(',')] 
        p = 1.0-a[0]
        if count >= 8: 
            sigma = 0.04
        else:
            sigma = 0.01
        ID = (count % 8) / 2
        if ID == 0:
            name = '$G_{dr}$'
        if ID == 1:
            name = '$G_{di}$'
        if ID == 2:
            name = '$G_{sr}$'
        if ID == 3:
            name = '$G_{si}$'
        ans = ''
        ans += '%s & %.1f & %.2f ' % (name, p, sigma, )
        for i in range(4, 7):
            if a[i] < a[i-3]:
                ans += '& \\textbf{%.2fe-2}' % (a[i]*100)
            else:
                ans += '& %.2fe-2' % (a[i]*100)
        if a[8] < a[7]:
            ans += '& \\textbf{%.3fs}' % a[8]
        else:
            ans += '& %.3fs' % a[8]
        for i in range(1, 4):
            if a[i] < a[i+3]:
                ans += '& \\textbf{%.2fe-2}' % (a[i]*100)
            else:
                ans += '& %.2fe-2' % (a[i]*100)
        if a[7] < a[8]:
            ans += '& \\textbf{%.3fs}' % a[7]
        else:
            ans += '& %.3fs' % a[7]
        ans += '\\\\'
        strings[ID].append(ans)

    for j in range(4):
        string = strings[j]
        print(string[1])
        print(string[3])
        print(string[0])
        print(string[2])
        print('\\hline')

t = genfromtxt('/Users/tom/K2TargetManagement/Field0/K2Campaign0targets.csv',delimiter=',',usecols=[0,3],skip_header=1)

s = genfromtxt('cdpp.txt',dtype=None)

epicarr = np.zeros(len(t))
magarr = np.zeros(len(t))
cdpparr = np.zeros(len(t))
for i,line in enumerate(s):
    epic, mag = [x for x in t if str(int(x[0])) in line[0]][0]
    epicarr[i] = epic
    magarr[i] = mag
    cdpparr[i] = line[1]

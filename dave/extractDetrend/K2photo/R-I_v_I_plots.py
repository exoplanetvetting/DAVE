import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

tbl10 = Table.read('/Users/bryanmann/Desktop/10_table', format='ascii')
tbl50 = Table.read('/Users/bryanmann/Desktop/50_table', format='ascii')
tbl90 = Table.read('/Users/bryanmann/Desktop/90_table', format='ascii')

fontsize = 12

"""plt.scatter(tbl10['R-I'],tbl10['rmag_1'], c='grey')
plt.xlabel('R-I', fontsize=fontsize)
plt.ylabel('R', fontsize=fontsize)
plt.title('R v R-I', fontsize=18)"""
plt.ylim(17,13)
#plt.hold(True)

plt.scatter(tbl50['R-I'], tbl50['rmag_1'], c='yellow')
plt.xlabel('R-I', fontsize=fontsize)
plt.ylabel('R', fontsize=fontsize)
plt.title('R-I v I', fontsize=18)
plt.legend('50% Membership Probability', label=
plt.hold(True)

plt.scatter(tbl90['R-I'], tbl90['rmag_1'], c='b')
plt.xlabel('R-I', fontsize=fontsize)
plt.ylabel('R', fontsize=fontsize)
plt.title('R-I v I', fontsize=18)
plt.hold(True)
plt.show()
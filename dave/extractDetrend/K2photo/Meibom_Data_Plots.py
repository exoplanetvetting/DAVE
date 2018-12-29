import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

Meibom_Data = Table.read('/Users/bryanmann/Documents/NASA_Kepler_2.0/Compiled Lists/Meibom_with_Kepler_Rot', format='ascii')

fontsize = 12

plt.scatter(Meibom_Data['(B-V)0_1'], Meibom_Data['col2'], c='k')
plt.axis(fontsize=fontsize)
plt.xlabel('B-V', fontsize=fontsize)
plt.ylabel('Rotation Period [d]', fontsize=fontsize)
plt.title('Color vs Period', fontsize=18)
plt.hold(True)
plt.tight_layout()
plt.show()
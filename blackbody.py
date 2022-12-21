import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

# constants

k = 1.381e-23
h = 6.626e-34
c = 2.998e8
h_divided_by_k = 4.80e-11

def RJlaw(nu, T):
    return 8.0 * np.pi * (nu ** 2) * (k * T) / (c ** 3)
def Planck(nu, T):
    ret = (nu ** 3) / (np.exp(h_divided_by_k * nu / T) - 1)
    return 8.0 * np.pi * h / (c ** 3) * ret

nu_all = np.linspace(1, 2e15, 10000)
nu_RJ = np.linspace(1, 2.0e14, 200)

temp = np.linspace(4000, 6000, 10)
color = np.linspace(0, 0.7, 10)

sm = plt.cm.ScalarMappable(cmap=plt.cm.hot, norm=plt.Normalize(vmin=4000, vmax=6000))
plt.clf()

plt.plot(nu_RJ/1e14, RJlaw(nu_RJ, 4000)*1e15, 'k--', linewidth=1.0, label='LJ law, T=4000K')

for T, C in zip(temp, color):
    plt.plot(nu_all/1e14, Planck(nu_all, T)*1e15, c=plt.cm.hot(C))
    
plt.legend()
cbar = plt.colorbar(sm)
cbar.set_label('Temperature')
plt.xlabel('frequency / 1e14')
plt.ylabel('Intensity')
plt.title('Blackbody radiation')
plt.show()
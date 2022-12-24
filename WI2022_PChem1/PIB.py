import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import warnings

warnings.filterwarnings("ignore")

# Electron in a one-dimensional box
# of length 10 nm

m = 9.109e-31
a = 1.0e-8
h = 6.626e-34

def Energylevel(n):
    return (n**2)*(h**2)/(8.*m*(a**2))
def waveftn(x, n):
    return np.sqrt(2./a)*np.sin(n*np.pi*x/a)

fig = plt.figure(figsize=(5,8))
x = np.arange(0, 1.01e-8, 1.0e-10)
ones = np.ones(101)

for i in range(1, 6):
    
    E = ones * Energylevel(i) * 1.0e21
    plt.plot(x, E, 'k-')
    plt.plot(x, waveftn(x, i) / 2.0e4 + E, label=f'n={i}')
    
plt.xlabel('x')
plt.ylabel('E (scaled)')
plt.ylim(0, 20)
plt.title('Energy level and wavefunction')
plt.legend(loc='upper right')
plt.show()
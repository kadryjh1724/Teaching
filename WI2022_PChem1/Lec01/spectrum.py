import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import warnings

warnings.filterwarnings("ignore")

# constants

m_e = 9.10938e-31
e = 1.60218e-19
eo = 8.85419e-12
h = 6.626e-34
hc = 1.98645e-25

def E(n):
    return -m_e*(e**4)/(8.*(eo**2)*(h**2)*(n**2))
def deltaE(m, n):
    return abs(E(m) - E(n))
def lamda(m, n):
    return hc / deltaE(m, n)

def wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an 
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B,1)

clim=(350,780)
norm = plt.Normalize(*clim)
wl = np.arange(clim[0],clim[1]+1,2)
colorlist = list(zip(norm(wl),[wavelength_to_rgb(w) for w in wl]))
spectralmap = matplotlib.colors.LinearSegmentedColormap.from_list("spectrum", colorlist)

ones = np.ones(10)
y = np.linspace(0, 3, 10)

fig = plt.figure(figsize = (15, 3))

# start: m, end: n

for n in range(1,4):
    
    for m in range(n+1, n+8):
        
        wavelength = lamda(m, n) * 1e9
        norm_wavelength = (wavelength - 350) / 430
        plt.plot(wavelength*ones, y, c=spectralmap(norm_wavelength))
        
plt.xlabel('wavelength (nm)')
plt.title('Hydrogen atomic spectrum')
plt.show()
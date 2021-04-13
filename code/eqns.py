import numpy as np


def MM09e2(B_V, age):
    ''' 
    Eqn 2 
    http://adsabs.harvard.edu/abs/2009ApJ...695..679M
    '''
    a = 0.50
    b = 0.15
    P = np.sqrt(age) * (np.sqrt(B_V - a)) - b * (B_V - a)
    return P


def MM09e3(B_V, age):
    ''' 
    Eqn 3 
    http://adsabs.harvard.edu/abs/2009ApJ...695..679M
    '''
    c = 0.77
    d = 0.40
    f = 0.60
    P = age**0.52 * (c * (B_V - d)**f)
    return P


def Angus2015(B_V, age):
    '''
    Compute the rotation period expected for a star of a given color (temp) and age

    NOTE: - input Age is in MYr
          - output Period is in days

    Eqn 15 from Angus+2015
    http://adsabs.harvard.edu/abs/2015MNRAS.450.1787A

    '''
    P = (age ** 0.55) * 0.4 * ((B_V - 0.45) ** 0.31)

    return P


def Angus2015_age(B_V, P):
    '''
    invert the above eqn
    '''
    age = np.power(P / (0.4 * ((B_V - 0.45) ** 0.31)), 1. / 0.55)
    return age


def Noyes1984_eqn4(B_V):
    '''
    Eqn 4 from the classic Noyes et al. (1984)
    http://adsabs.harvard.edu/abs/1984ApJ...279..763N
    
    empirical definition of Tau (convective turnover timescale) vs B-V color
    '''
    x = 1 - B_V
    
    # actuall the eqn for log tau_c
    tau = np.piecewise(x, [(x < 0),(x >= 0)], 
                       [lambda y: (1.362 - 0.14 * y), 
                        lambda y: (1.362 - 0.166 * y + 0.025 * y**2 - 5.323 * y**3)])
    
    return 10**tau


def Wright2011_tau(BV):
    '''
    http://adsabs.harvard.edu/abs/2011ApJ...743...48W
    '''
    # values from their Table 2
    BV0 = [0.46,0.61,0.76,0.92,1.13,1.32,1.41,1.50,1.55,1.61]
    BV1 = [0.61,0.75,0.92,1.12,1.31,1.41,1.49,1.55,1.60,1.95]
    logtau = [1.01, 1.08, 1.18, 1.32, 1.41, 1.49, 1.71, 1.94, 1.97, 2.12]

    clr = (np.array(BV0) + np.array(BV1)) / 2.
    
    # fit parameters I created
    ff = [0.51138488, -0.24907552, 1.00734295]
    tau = 10**np.polyval(ff, BV)
    
    return tau


def _gaus(x, a, b, x0, sigma):
    """
    Simple Gaussian function
    Parameters
    ----------
    x : float or 1-d numpy array
        The data to evaluate the Gaussian over
    a : float
        the amplitude
    b : float
        the constant offset
    x0 : float
        the center of the Gaussian
    sigma : float
        the width of the Gaussian
    Returns
    -------
    Array or float of same type as input (x).
    """
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + b


def CranmerSaar2011_eqn36(Teff):
    """
    https://ui.adsabs.harvard.edu/abs/2011ApJ...741...54C
    
    Teff in Kelvin
    
    Returns Tau in days
    """
    tau = 314.24 * np.exp(-1*(Teff / 1952.5) - ((Teff/6250)**18)) + 0.002
    return tau


def teff2bv(teff, logg=4.3, feh=0):
    """
    Relation from Sekiguchi & Fukugita (2000)
    https://arxiv.org/abs/astro-ph/9904299
    """    
    # from their Tbl 4:
    t = [-813.3175, 684.4585, -189.923, 17.40875]
    f = [1.2136, 0.0209]
    d1, g1, e1 = -0.294, -1.166, 0.3125
    BV = t[0] + t[1]*np.log10(teff) + t[2]*(np.log10(teff))**2 + \
            t[3]*(np.log10(teff))**3 + f[0]*feh + f[1]*feh**2 \
            + d1*feh*np.log10(teff) + g1*logg + e1*logg*np.log10(teff)
    return BV


def bv2teff(BV, logg=4.3, feh=0):
    """
    Relation from Sekiguchi & Fukugita (2000)
    https://arxiv.org/abs/astro-ph/9904299
    """
    # Full Sample, Tbl 3
    c = np.array([3.939654, -0.395361, 0.2082113, -0.0604097])
    f1, f2, g1, h1 = 0.027153, 0.005036, 0.007367, -0.01069
    
    logTeff = c[0] + c[1]*BV + c[2]*(BV**2) + c[3]*(BV**3) + \
               f1*feh + f2*(feh**2) + \
               g1*logg + h1*BV*logg
    return 10**logTeff
    
 
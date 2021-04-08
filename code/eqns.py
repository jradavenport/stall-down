import numpy as np


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
    # P = (age ** 0.55) * 0.4 * ((B_V - 0.45) ** 0.31)
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
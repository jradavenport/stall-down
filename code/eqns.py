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
    Compute the rotation period expected for a star of a given color (temp) and age

    NOTE: - output Age is in MYr
          - input Period is in days

    Eqn 15 from Angus+2015
    http://adsabs.harvard.edu/abs/2015MNRAS.450.1787A

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

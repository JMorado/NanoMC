def fugacity_coeff(temperature, pressure):
    """
    Temperature in Kelvin and pressure in atm.
    c3 is wrong in the original article.
    """

    import numpy as np

    c1 = np.exp(-3.8402*temperature**(1.0/8.0) + 0.5410)
    c2 = np.exp(-0.1263*temperature**(1.0/2.0) - 15.980)
    c3 = 300.0*np.exp(-0.011901*temperature - 5.941)

    fugacity = c1*pressure - c2*pressure*pressure + c3*(np.exp(-pressure/300.0)-1.0)
    return np.exp(fugacity)


print(fugacity_coeff(300.0,1513.0)*1513.0)

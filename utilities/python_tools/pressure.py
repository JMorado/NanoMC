import math

pressure = float(input("Pressure (bar):"))
T = float(input("Temperature (K):"))


kb = 1.38064852e-23
planck = 6.62607004e-34
mass = 2.0 * 1.6737236e-27
beta = 1.0 / (kb * T)

thermal_wav = planck * planck * beta / (2 * math.pi * mass)
thermal_wav = thermal_wav ** (1/2.)

chemical_pot = math.log(pressure * thermal_wav ** 3 * beta * 1e5) * T

print("The chemical potential is: " + str(chemical_pot))

import math
import os


# Generate series of uvt files
# General parameters
seed = "nanotube_10_10"
T = 300.0

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


def chemical_pot(pressure_ideal, temperature):
    """
    Compute the chemical potential for a given pressure.
    """
    # Compute real pressure ( a.k.a. fugacity )
    bar_to_atm = 0.986923
    pressure = pressure_ideal * fugacity_coeff(T,pressure_ideal * bar_to_atm)

    # Determine the chemcial potential for a given pressure
    kb = 1.38064852e-23
    planck = 6.62607004e-34
    mass = 2.0 * 1.6737236e-27
    beta = 1.0 / (kb * T)
    thermal_wav = planck * planck * beta / (2 * math.pi * mass)
    thermal_wav = thermal_wav ** (1/2.)

    chemical_pot = math.log(pressure * thermal_wav ** 3 * beta * 1e5) * T

    return chemical_pot, pressure


def write_input(u, T,x,name,nanotube_file):
    """
    Subroutine that writes input file to
    """

    ensemble = "uvt"
    radius = 6.780

    f = open(name,'w')
    f.write("atomistic                ! swcnt model \n")
    f.write( str(ensemble) + "        ! ensemble \n")
    f.write(str(u) + " \n")
    f.write(nanotube_file + " \n")
    f.write("482                      ! Number of MC sweeps \n")
    f.write("2000000                  ! number of MC sweeps \n")
    f.write(str(T) + "                ! temperature \n")
    f.write("491.90242935             !  \n")
    f.write(str(radius) + "                     ! swcnt radius \n")
    f.write("12.0                     ! fluid_cut_off \n")
    f.write("12.0                     ! fluid_cnt_cut_off \n")
    f.write("2.96                     ! sigma LJ 12-6 fluid \n")
    f.write("34.2                     ! epsilon LJ 12-6 fluid \n")
    f.write("3.18                     ! sigma fluid-C \n")
    f.write("30.95	              ! epsilon fluid-CNT \n")
    f.write("uvt_tra                  ! outputseed \n")
    f.write("995293                   ! random generator seed \n")
    f.write("100                      ! frequency of properties writting \n")
    f.write("100                      ! frequency of trajectory writting \n")
    f.write("random                   ! mode to generate initial configuration \n")
    f.write(str(x) + " \n")

    return f.close()


# Create main folder and change to it

main_folder = str(seed) + "_temp_" + str(T)
os.mkdir(main_folder)
os.chdir(main_folder)

# Generate files for different pressures
P0 = 50.0
dP = 100.0
nanotube_file = "nanotube-10-10-200.xyz"

for i in range(10):
    P = P0 + i * dP
    u, preal = chemical_pot(P, T)
    # Create sub folder
    sub_folder = "p_" + '{0:03f}'.format(preal)
    os.mkdir( sub_folder )
    os.chdir( sub_folder )

     # TODO: copy new.xyz, nanotube-10-10.xyz, nanomc_nvt.exe, input_fil
    os.system("cp ../../tmp/nanomc_uvt.exe .")
    os.system("cp ../../tmp/" + str(nanotube_file) +  " .")
    os.system("cp ../../tmp/input_file .")

    print(u, preal, P)
    # Write input file
    name = "input_uvt_atomistic.nanomc"
    write_input(u, T,0.5,name)

    # TODO: copy nanotube-10-10.xyz, nanomc_uvt.exe, input_uvt

    # Go back to main folder
    os.chdir("..")



import math
import os


seed = "nanotube_10_10"
pressure = 10.0
T = 300.0


# Determine the chemcial potential for a given pressure
# TODO: include fugacity
kb = 1.38064852e-23
planck = 6.62607004e-34
mass = 2.0 * 1.6737236e-27
beta = 1.0 / (kb * T)

thermal_wav = planck * planck * beta / (2 * math.pi * mass)
thermal_wav = thermal_wav ** (1/2.)

chemical_pot = math.log(pressure * thermal_wav ** 3 * beta * 1e5) * T

print("The chemical potential is: " + str(chemical_pot))




def write_input(T,x,name):
    """
    Subroutine that writes input file to
    """

    ensemble = "nvt"
    nanotube_file = "nanotube-10-10-200.xyz"

    f = open(name,'w')
    f.write("atomistic                ! swcnt model \n")
    f.write( str(ensemble) + "                      ! ensemble \n")
    f.write(nanotube_file + " \n")
    f.write("482                      ! Number of MC sweeps \n")
    f.write("5000000                  ! number of MC sweeps \n")
    f.write(str(T) + "                    ! temperature \n")
    f.write("491.90242935             !  \n")
    f.write("6.78                     ! swcnt radius \n")
    f.write("12.0                     ! fluid_cut_off \n")
    f.write("12.0                     ! fluid_cnt_cut_off \n")
    f.write("2.96                     ! sigma LJ 12-6 fluid \n")
    f.write("34.2                     ! epsilon LJ 12-6 fluid \n")
    f.write("3.18                     ! sigma fluid-C \n")
    f.write("30.95	              ! epsilon fluid-CNT \n")
    f.write("mjoao                    ! outputseed \n")
    f.write("995293                   ! random generator seed \n")
    f.write("100000                   ! frequency of trajectory writting \n")
    f.write("100000                   ! frequency of properties writting \n")
    f.write("random                   ! mode to generate initial configuration \n")
    f.write(str(x) + "\n")

    return f.close()


# Create main folder and change to it

main_folder = str(seed) + "_pressure_" + str(pressure) + "_temp_" + str(T)
os.mkdir(main_folder)
os.chdir(main_folder)

# Generate files for different displacements
x0 = 0.05
dx = 0.075
for i in range(15):

    x = x0 + i * dx

    # Create sub folder
    sub_folder = "x_" + '{0:03f}'.format(x)
    os.mkdir( sub_folder )
    os.chdir( sub_folder )

    # TODO: copy new.xyz, nanotube-10-10.xyz, nanomc_nvt.exe, input_fil
    os.system("cp ../../tmp/nanotube*xyz .")
    os.system("cp ../../tmp/nanotube*xyz .")
    os.system("cp ../../tmp/input_file .")

    # Write input file
    name = "input_nvt_atomistic.nanomc"
    write_input(T,x,name)



    # Go back to main folder
    os.chdir("..")



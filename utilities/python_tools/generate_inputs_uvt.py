from fugacity import *
from chemical_potential import *
import math
import os



def write_input(file_name, chemical_pot, T,nanotube_file):
    """
    Subroutine that writes input file to
    """

    ensemble = "uvt"
    radius = 6.780
    length= 491.90242935

    f = open(name,'w')
    f.write("atomistic                ! swcnt model \n")
    f.write("uniform                  ! displacement \n")
    f.write(str(ensemble) + "         ! ensemble \n")
    f.write(str(chemical_pot) + "     ! chemical potential \n")
    f.write(nanotube_file + "         ! nanotube file \n")
    f.write("482                      ! particles \n")
    f.write("10000000                 ! number of MC sweeps \n")
    f.write(str(T) + "                ! temperature \n")
    f.write(str(length) + "           ! length \n")  # 8 dc
    f.write(str(radius) + "           ! swcnt radius \n")
    f.write(str(radius) + " " + str(radius) + " " + str(length) + " 90.0 90.0 90.0    ! cell \n ")
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

    return f.close()


# General parameters
seed = "nanotube_10_10"
T = 293.15

# Create main folder and change to it
main_folder = str(seed) + "_temp_" + str(T)
os.mkdir(main_folder)
os.chdir(main_folder)

bar_to_atm = 0.986923
bar_to_pa = 100000.0

# Generate files for different pressures (bar)
pressure_range = [5, 10, 20, 40, 80, 160, 320, 640, 800]
nanotube_file = "nanotube-10-10-200.xyz"


print("T (K), Pideal (Pa), Peff (Pa), chemical_pot (K)")
for P in pressure_range:
    fug_coeff = fugacity_coeff(T,P*bar_to_atm) # (adimensional)
    p_eff = fug_coeff * P * bar_to_pa # In Pascal (converting from bar to Pa)

    chemical_pot = chemical_potential(T, p_eff)

    # Create sub folder
    sub_folder = "p_" + '{0:03f}'.format(p_eff)
    os.mkdir( sub_folder )
    os.chdir( sub_folder )

     # TODO: copy new.xyz, nanotube-10-10.xyz, nanomc_nvt.exe, input_fil
    os.system("cp ../../tmp/nanomc_uvt.exe .")
    os.system("cp ../../tmp/" + str(nanotube_file) +  " .")
    os.system("cp ../../tmp/input_file .")

    print(T, P*bar_to_pa, p_eff, chemical_pot, fug_coeff)
    # Write input file
    name = "input_uvt_atomistic.nanomc"
    write_input(name, chemical_pot, T, nanotube_file)

    # Go back to main folder
    os.chdir("..")



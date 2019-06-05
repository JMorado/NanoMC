####*****INPUT FILE****####
#All the units are in ANGSTROMs
#If a line begins with a sharp (#), it is considered a comment and, thus, it is ignored
#If you don't want a parameter to be taken in consideration write None
#To run the simulation write >>> python run_simulation.py in the command line


#####*****Simulation Properties*****####
model = "molecular" #options:"molecular";"atomic" #the atomic option is not implemented yet
n_cycles = 100
n_particles = 500
temperature = 77 #Kelvin
initial_position = 'random' #options: "origin"; "random"; "file"
file_name = None #choose None if the initial position is not definided by coordinates in a file. Otherwise write 'filename.txt'
output_file_name = None #if None, the file is not created. If you want to create a file write the name of it with the '.txt' extension in the end.


#####*****Cell/Geometry*****######
cell_type = "swcnt" #options: "vacuum"; "swnct"

####*****Cell Properties*****##### Only used in the swcnt algorithm.
height = 100
radius = 10
xy_boundary_conditions = "reflexion" #options: "reflexion"; "hard_wall"
z_boundary_conditions = "periodicity" #options: "periodicity"; "hard_wall"

#####*****Particle Potentials*****##### They work for all systems.
inter_particle = "lennard_jones_12_6" #options: "lennard_jones_12_6"; None
inter_part_cut_off = 8

#####*****Cell Potentials*****##### Only used in the swcnt algorithm. If this was not the cell type choosen, this section has no meaning.
particle_cell = "lennard_jones_12_6" #options: "lennard_jones_12_6"; None
#There is not cut_off for this potential. However, it may be necessary to implement if the swcnt is too wide.


#####****Box Parameteres****#####
box = "no" #If parameter is 'yes', there will be no boundary conditions in the z axis for the swcnt
#####*****--!--BOX---PARAMETERS--!--*****##### Only used if the box parameter was marked as 'yes'.
box_height = 10
box_width = 10





#####NOT USED RIGHT NOW####
'''
#####*****Box Potentials*****###### Only used if the box parameter was marked as 'yes'.
box_inter_particle = "no"
box_inter_particle_cut_off = 0
box_particle = "no"
box_particle_cut_off = 0
'''

#!/usr/bin/env python

import funcoes_randomwalk as rw
import random
random.seed(3) #seed escolhida
import math
import numpy as np



def random_walk_simulation_molecular():
    #Creates the simulation object
    simulation = rw.Simulation()
    #Opens the file object
    file_io = rw.Files()
    file_io.screen_init(simulation.dict_parameters)
    

    #Defines the functions that will be used in the simulation
    energy_inter = simulation.energy_inter()
    energy_cell_part = simulation.energy_cell_part()
    xy_boundary_conditions = simulation.xy_boundary_conditions_f()
    z_boundary_conditions = simulation.z_boundary_conditions_f()
    cell_algorithm = simulation.cell_properties()
    
    
    #Create the matrix with the initial position of the particles
    matrix = simulation.create_matrix(file_io) #file_object in case the simulation is initialized by a file
    time = rw.Time()
    difusion = rw.Difusion()
    displacement = rw.Displacements()
    coefs = []
    
    for mov in range(simulation.n_cycles):
        if mov in range(0, simulation.n_cycles, 100):
            #print mov
            pass
        displacement.d_ciclo = 0
        
        for p in range(simulation.n_particles):
            ####### RANDOM-WALK STEP OBJECT #######
            RWstep = rw.RandomWalk(matrix, p, simulation)
            
        
            ####### CELL PROPERTIES #######
            deltaE = cell_algorithm(RWstep, xy_boundary_conditions, z_boundary_conditions, energy_inter, energy_cell_part)

            ###### PROBABILITY DETERMINATION #####
            if deltaE <= 0:
                p_accept = 1.0
            else:
                p_accept = math.e ** ( - (deltaE) / simulation.temperature)
        
        
            ###### MOVEMENT ACCEPTANCE ######
            if p_accept > random.uniform(0,1):
                new_position = RWstep.get_new_position()
                matrix[p][0] = new_position[0]
                matrix[p][1] = new_position[1]
                matrix[p][2] = new_position[2]

                displacement.update_displacements_ciclo(RWstep.displacement) #updates the displacement
                
                
             
        displacement.update_displacements()
        time.refresh_time(displacement.get_cycle_displacement()) #updates time
	
	difusion.self_difusion(displacement.get_total_displacement(), time.get_total_time())
        
        coefs.append(difusion.self_difusion(displacement.get_total_displacement(), time.get_total_time()))
        print(difusion.self_difusion(displacement.get_total_displacement(), time.get_total_time()))
        file_io.write_file('coord/cord_conc_const_np_{0}_{1}.txt'.format(simulation.n_particles, mov), matrix) #guarda coordenadas         
    file_io.write_coef('deff_conc_const_np_{0}.txt'.format(simulation.n_particles), coefs) #guarda Deff

    print "#### THE SIMULATION HAS FINISHED ####"
    return matrix
    
    




random_walk_simulation_molecular()





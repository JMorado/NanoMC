#!/usr/bin/env python
import math
import random
import numpy as np
import parameters #importa input file
pi = math.pi



class RandomWalk(object):
    def __init__(self, matrix, index, simulation_object):
        self.index = index
        self.simulation_object = simulation_object
        self.matrix = np.array(matrix)
        self.x = matrix[index][0]
        self.y = matrix[index][1]
        self.z = matrix[index][2]
        self.T = simulation_object.temperature
        self.h = float(simulation_object.h)
        self.r = float(simulation_object.r)
        self.r_corte_inter = int(simulation_object.inter_part_cut_off)
        self.old_position = np.array(matrix[index])
        self.displacement = self.deslocamento()
        self.new_position = np.array(self.matrix[index] + np.array(self.displacement))

    def deslocamento(self):
        dx = - math.log(random.uniform(0, 1))
        teta = random.uniform(0, 2*pi)
        psi = math.acos(2 * random.uniform(0,1) - 1)

        self.d = [dx*math.cos(teta)*math.sin(psi), dx*math.sin(teta)*math.sin(psi), dx*math.cos(psi)]

        return self.d


    def periodicity_z(self):
        if abs(self.new_position[2]) > self.h:
            excess = abs(self.new_position[2]) - self.h

            if self.new_position[2] < 0:
                self.new_position[2] = self.h - excess

            else:
                self.new_position[2] = - self.h + excess

            return self.new_position[2]

        else:
            return self.new_position[2]


    def swcnt_particle(self):
        z_old = self.r - math.hypot(self.old_position[0], self.old_position[1])
        z_new = self.r - math.hypot(self.new_position[0], self.new_position[1])

        X12 = (z_old ** 12 - z_new ** 12)
        Y12 = ((z_new ** 12) * (z_old ** 12))
        X6 = (z_old ** 6 - z_new ** 6)
        Y6 = ((z_new ** 6) * (z_old ** 6))

        deltaE = 4 * 300 * ( ((3 ** 12) * (X12 / Y12)) - ((3 ** 6) * (X6 / Y6)))

        return deltaE

    def swcnt_particle_old(self):
        z_old = self.r - math.hypot(self.old_position[0], self.old_position[1])

        energia_old = (4 * 30 * (( 3 / z_old) ** 24 - ( 3 / z_old) ** 12))

        return energia_old

    def swcnt_particle_new(self):
        z_new = self.r - math.hypot(self.new_position[0], self.new_position[1])

        energia_new = (4 * 30 * (( 3 / z_new) ** 24 - ( 3 / z_new) ** 12))

        return energia_new

    def box_algorithm(self, reflexao, periodicity_z, energy_inter, energy_cell_part):
        if abs(self.old_position[2]) <= self.h:
                reflexao(self)
                deltaE = energy_inter(self) + energy_cell_part(self)
                return deltaE
        elif abs(self.old_position[2]) > self.h:
            periodic_box = self.periodic_box()
            if abs(self.new_position[2]) <= self.h:
                reflexao(self)
                deltaE_inter = self.inter_molecular_new()
                deltaE_cell = self.swcnt_particle_new()
                deltaE = deltaE_inter + deltaE_cell
                return deltaE
            else:
                return 0


    def regular_algorithm(self, reflexao, periodicity_z, energy_inter, energy_cell_part):
        reflexao(self)
        periodicity_z(self)
        deltaE = energy_inter(self) + energy_cell_part(self)
        return deltaE

    def vacuum_algorithm(self, reflexao, periodicity_z, energy_inter, energy_cell_part):
        deltaE = energy_inter(self)
        return deltaE

    def periodic_box(self):
        width = float(self.simulation_object.box_width)
        height = float(self.simulation_object.box_height)

        if  abs(self.new_position[2]) > self.h and abs(self.new_position[2]) < self.h + 2*width:
            #we are inside the box
            if abs(self.new_position[1]) > width:
                excess = abs(self.new_position[1]) - width
                self.new_position[1] = math.copysign(width, - self.new_position[1]) + excess

            if abs(self.new_position[0]) > height:
                excess = abs(self.new_position[0]) - height
                self.new_position[0] = math.copysign(height, - self.new_position[0]) + excess

            if abs(self.new_position[0]) < height and abs(self.new_position[1]):
                #in case that we don't need to use the periodic conditions
                return self.new_position

        elif abs(self.new_position[2]) > self.h + 2*width:
            #we are outside the box
            self.new_position = self.old_position
            return self.new_position

        elif abs(self.new_position[2]) < self.h and abs(self.old_position[2]) > self.h:
            #we are moving from the box to the nanotube
            distance_to_origin = math.hypot(self.new_position[0], self.new_position[1])
            if distance_to_origin >= self.r:
                self.new_position = self.old_position
                return self.new_position
            else:
                return self.new_position
        else:
            #we are inside the cell
            return self.new_position


    def reflexao(self):
        desloc_x = self.displacement[0]
        desloc_y = self.displacement[1]
        while (self.r) < math.hypot(self.x + desloc_x, self.y + desloc_y):
            #Grava as posicoes iniciais das particulas que serao necessaras para o calculo de zeta_inicial.
            posicao_i_x = self.x
            posicao_i_y = self.y


            #Calculo do deslocamento necessario para a particula se mover para cima da circunferencia.
            vectord = math.hypot(desloc_x, desloc_y)
            vector_unitario_dx = desloc_x / vectord
            vector_unitario_dy = desloc_y / vectord

            a = 1.0
            b = 2.0 * ( (self.x*desloc_x + desloc_y*self.y) / vectord)
            c = (self.x * self.x + self.y * self.y) - (self.r*self.r)

            if (c  > 0):
                print(a,b,c)
                print(self.r, self.x, self.y)

            d = (- b + math.sqrt(b*b - 4*a*c)) / 2*a

            #print(a,b,c,d,self.x,self.y)
            #Mover a particula para a circunferencia
            self.x = self.x + (desloc_x / vectord) * d
            self.y = self.y + (desloc_y / vectord) * d


            #Calculo do deslocamento que ainda falta percorrer.
            d_restox = (desloc_x - vector_unitario_dx*d)
            d_restoy = (desloc_y - vector_unitario_dy*d)

            #Calculo do angulo de incidencia (que sera == ao angulo de reflexao). Determina-se tambem se a reflexao
            #e feita no sentido horario ou anti-horario.
            modulo_i = math.hypot(self.x, self.y)


            alpha = math.acos(round((desloc_x*self.x + desloc_y*self.y) / (modulo_i * vectord), 4))

            angulo_i = (pi / 2) - alpha


            zeta = math.atan(self.y / self.x)

            try:
                zeta_i = math.atan(posicao_i_x / posicao_i_y)
            except:
                print(posicao_i_x, posicao_i_y)
                zeta_i =  zeta

            if zeta_i > zeta and zeta_i < zeta + pi:
                beta = - 2 * angulo_i
            elif zeta_i < zeta and zeta_i > zeta + pi:
                beta = 2 * angulo_i
            else:
                beta = pi


            #Rotacao do vector deslocamento segundo beta radianos.
            desloc_x = (d_restox*math.cos(beta) - d_restoy*math.sin(beta))
            desloc_y = (d_restox*math.sin(beta) + d_restoy*math.cos(beta))


            print(self.x,desloc_x)
            print(self.y,desloc_y)
        self.new_position[0] = self.x + desloc_x
        self.new_position[1] = self.y + desloc_y

        print(self.new_position,"****")

        return self.x, self.y, desloc_x, desloc_y


    def z_hard_wall(self):
        if (self.h) < abs(self.z + self.displacement[2]):
            self.new_position = self.old_position

        return self.new_position

    def xy_hard_wall(self):
        if (self.r) < math.hypot(self.x + self.displacement[0], self.y + self.displacement[1]):
            self.new_position = self.old_position

        return self.new_position

    def inter_molecular_old(self):
        #determina distancia absoluta entre dois pontos sem ter em consideracao as condicoes periodicas
        dif_old = (self.matrix - self.old_position) ** 2

        #determina a distancia entre dois pontos tendo em consideracao as condicoes periodicas
        dif_old_pz = (2 * self.h - (abs(self.matrix[:,2:] - self.old_position[2])))

        dif_old[:,2:] = np.minimum(dif_old[:,2:], dif_old_pz)

        #obtem-se quadrado da distancia entre 2 pontos no espaco tridimensional
        dist_old = dif_old.sum(1)

        #cut off
        dist_old[dist_old > self.r_corte_inter*self.r_corte_inter] = 0


        energia_old = (4 * 30 * (( 3 / dist_old[dist_old != 0.]) ** 24 - ( 3 / dist_old[dist_old != 0.]) ** 12))
        energia_old_total = np.sum(energia_old)

        return energia_old_total


    def inter_molecular_new(self):
        #determina distancia absoluta entre dois pontos sem ter em consideracao as condicoes periodicas
        dif_new = (self.matrix - self.new_position) ** 2

        #determina a distancia entre dois pontos tendo em consideracao as condicoes periodicas
        dif_new_pz = (2 * self.h - (abs(self.matrix[:,2:] - self.new_position[2])))

        dif_new[:,2:] = np.minimum(dif_new[:,2:], dif_new_pz)

        #obtem-se quadrado da distancia entre 2 pontos no espaco tridimensional
        dist_new = dif_new.sum(1)

        #cut off
        dist_new[dist_new > self.r_corte_inter*self.r_corte_inter] = 0


        energia_new = (4 * 30 * (( 3 / dist_new[dist_new != 0.]) ** 24 - ( 3 / dist_new[dist_new != 0.]) ** 12))
        energia_new_total = np.sum(energia_new)

        return energia_new_total


    def inter_molecular(self):

        self.matrix = np.delete(self.matrix, self.index, 0) #elimina a linha com o index da matriz
        #determina distancia absoluta entre dois pontos sem ter em consideracao as condicoes periodicas
        dif_old = (self.matrix - self.old_position) ** 2
        dif_new = (self.matrix - self.new_position) ** 2

        #CONDICOES PERIODICAS
        if self.old_position[2] < 0:
            periodic_old = self.old_position[2] + self.h * 2
        else:
            periodic_old = self.old_position[2] - self.h * 2
        if self.new_position[2] < 0:
            periodic_new = self.new_position[2] + self.h * 2
        else:
            periodic_new = self.new_position[2] - self.h * 2

        dif_old_pz = (self.matrix[:,2] - periodic_old) ** 2
        dif_new_pz = (self.matrix[:,2] - periodic_new) ** 2

        dif_old[:,2] = np.minimum(dif_old[:,2], dif_old_pz) #escolhe o valor minimo
        dif_new[:,2] = np.minimum(dif_new[:,2], dif_new_pz)

        #obtem-se quadrado da distancia entre 2 pontos no espaco tridimensional
        dist_old = dif_old.sum(1)
        dist_new = dif_new.sum(1)

        #cut off
        dist_old[dist_old > self.r_corte_inter*self.r_corte_inter] = 0
        dist_new[dist_new > self.r_corte_inter*self.r_corte_inter] = 0


        energia_old = (4 * 30 * (( 3*3 / dist_old[dist_old != 0]) ** 6 - ( 3*3 / dist_old[dist_old != 0]) ** 3))
        energia_old_total = np.sum(energia_old)

        energia_new = (4 * 30 * (( 3*3 / dist_new[dist_new != 0]) ** 6 - ( 3*3 / dist_new[dist_new != 0]) ** 3))
        energia_new_total = np.sum(energia_new)


        deltaE = energia_new_total - energia_old_total


        return deltaE



    def inter_molecular_OLD(self):
        u_old_total = 0
        u_new_total = 0

        for item in self.matrix:
            if item[0] != self.old_position[0] and item[1] != self.old_position[1] and item[2] != self.old_position[2]:
                distancia_inter_particula_old = math.sqrt( (self.old_position[0] - item[0])*(self.old_position[0] - item[0])
                        + (self.old_position[1] - item[1])*(self.old_position[1] - item[1])
                        + (self.old_position[2] - item[2])*(self.old_position[2] - item[2]) )

                distancia_inter_particula_new = math.sqrt( ((self.new_position[0]) - item[0])*((self.new_position[0]) - item[0])
                        + ((self.new_position[1]) - item[1])*((self.new_position[1]) - item[1])
                        + ((self.new_position[2]) - item[2])*((self.new_position[2]) - item[2]) )

                if distancia_inter_particula_old < self.r_corte_inter:
                    u_old = (4 * 30 * (( 3 / distancia_inter_particula_old ) ** 12 - ( 3 / distancia_inter_particula_old) ** 6))
                    u_old_total += u_old
                if distancia_inter_particula_new < self.r_corte_inter:
                    u_new = (4 * 30 * (( 3 / distancia_inter_particula_new ) ** 12 - ( 3 / distancia_inter_particula_new) ** 6))
                    u_new_total += u_new


                if (self.old_position[2] > (self.h - self.r_corte_inter) and self.old_position[2] <= self.h) or (self.old_position[2] < (-self.h + self.r_corte_inter) and self.old_position[2] > - self.h):
                    if self.old_position[2] > 0:
                        distancia_inter_particula_old_2 =  math.sqrt( (self.old_position[0] - item[0])*(self.old_position[0] - item[0])
                        + (self.old_position[1] - item[1])*(self.old_position[1] - item[1])
                        + ((self.old_position[2] - 2*self.h) - item[2])*((self.old_position[2] - 2*self.h) - item[2]) )

                    else:
                        distancia_inter_particula_old_2 =  math.sqrt( (self.old_position[0] - item[0])*(self.old_position[0] - item[0])
                        + (self.old_position[1] - item[1])*(self.old_position[1] - item[1])
                        + ((self.old_position[2] + 2*self.h) - item[2])*((self.old_position[2] + 2*self.h) - item[2]) )

                    if distancia_inter_particula_old_2 < self.r_corte_inter:
                        u_old = (4 * 30 * (( 3 / distancia_inter_particula_old_2 ) ** 12 - ( 3 / distancia_inter_particula_old_2) ** 6))
                        u_old_total += u_old



                if (self.z) > (self.h - self.r_corte_inter) and (self.z) <= self.h or ((self.z) < (-self.h + self.r_corte_inter) and (self.z) > - self.h):
                    if (self.z) > 0:
                        distancia_inter_particula_new_2 =  math.sqrt( ((self.new_position[0]) - item[0])*((self.new_position[0]) - item[0])
                        + ((self.new_position[1]) - item[1])*((self.new_position[1]) - item[1])
                        + (((self.new_position[2]) - 2*self.h) - item[2])*(((self.new_position[2]) - 2*self.h) - item[2]) )

                    else:
                        distancia_inter_particula_new_2 =  math.sqrt( ((self.new_position[0]) - item[0])*((self.new_position[0]) - item[0])
                        + ((self.new_position[1]) - item[1])*((self.new_position[1]) - item[1])
                        + (((self.new_position[2]) + 2*self.h) - item[2])*(((self.new_position[2]) + 2*self.h) - item[2]) )


                    if distancia_inter_particula_new_2 < self.r_corte_inter:
                        u_new = (4 * 30 * (( 3 / distancia_inter_particula_new_2 ) ** 12 - ( 3 / distancia_inter_particula_new_2) ** 6))
                        u_new_total += u_new

        deltaE = u_new_total - u_old_total
        return deltaE


    def get_new_position(self):
        return self.new_position


class Simulation(RandomWalk):
    def __init__(self):
        self.dict_parameters = {'file_name': parameters.file_name, 'output_file_name': parameters.output_file_name,
                      'model': parameters.model, 'n_particles' : parameters.n_particles, 'n_cycles': parameters.n_cycles,
                      'temperature': parameters.temperature, 'initial_position': parameters.initial_position,
                      'cell_type': parameters.cell_type, 'inter_particle': parameters.inter_particle,
                      'inter_particle_cut_off': parameters.inter_part_cut_off, 'particle_cell': parameters.particle_cell,
                      'radius': parameters.radius, 'height': parameters.height,
                      'xy_boundary_conditions': parameters.xy_boundary_conditions, 'z_boundary_conditions': parameters.z_boundary_conditions, 'box': parameters.box,
                      'box_width': parameters.box_width, 'box_height': parameters.box_height}


        self.output_file_name = parameters.output_file_name
        self.file_name = parameters.file_name
        self.model = parameters.model
        self.n_particles = parameters.n_particles
        self.n_cycles = parameters.n_cycles
        self.temperature = parameters.temperature
        self.initial_position = parameters.initial_position
        self.cell_type = parameters.cell_type
        self.inter_particle = parameters.inter_particle
        self.inter_part_cut_off = parameters.inter_part_cut_off
        self.particle_cell = parameters.particle_cell
        self.r = parameters.radius
        self.h = parameters.height
        self.xy_boundary_conditions = parameters.xy_boundary_conditions
        self.z_boundary_conditions = parameters.z_boundary_conditions
        self.box = parameters.box
        self.box_height = parameters.box_height
        self.box_width = parameters.box_width



    def energy_inter(self):
        if self.inter_particle == None:
            return self.nonDefined
        else:
            return RandomWalk.inter_molecular



    def energy_cell_part(self):
        if self.particle_cell == None:
            return self.nonDefined
        elif self.particle_cell == "lennard_jones_12_6":
            return RandomWalk.swcnt_particle



    def nonDefined(self, *arg):
        '''
        Dummy function that returns zero
        '''
        return 0


    def xy_boundary_conditions_f(self):
        if self.xy_boundary_conditions == "reflexion":
            return RandomWalk.reflexao
        elif self.xy_boundary_conditions == "hard_wall":
            return RandomWalk.xy_hard_wall
        else:
            return self.nonDefined


    def difusion(self):
        if self.difusion == "self":
            return Difusion.self_difusion
        elif self.difusion == "transport":
            return Difusion.transport_difusion
        else:
            return self.nonDefined

    def z_boundary_conditions_f(self):
        if self.z_boundary_conditions == "periodicity":
            return RandomWalk.periodicity_z
        elif self.z_boundary_conditions == "hard_wall":
            return RandomWalk.z_hard_wall
        else:
            return self.nonDefined


    def cell_properties(self):
        if self.cell_type == 'vacuum':
            return RandomWalk.vacuum_algorithm

        elif self.cell_type == "swcnt":
            if self.box == 'yes':
                return RandomWalk.box_algorithm
            elif self.box == 'no':
                return RandomWalk.regular_algorithm




    def create_matrix(self, file_object):
        if self.initial_position == 'origin':
            return [[0.,0.,0.] for p in range(self.n_particles)] #[0.,0.,-110.]
            #return [[0.,0.,-99.], [0,0,99.]]
        elif self.initial_position == 'random':
            matrix = []
            while len(matrix) != parameters.n_particles:
                x = random.uniform(-self.r, self.r)
                y = random.uniform(-self.r, self.r)
                if (self.r*self.r - x*x) > y*y:
                    z = random.uniform(-self.h, self.h)
                    matrix.append([x,y,z])
            return matrix
        elif self.initial_position == 'file':
            return file_object.open_file(self.file_name)


    def write_file(self, file_object, matrix):
        if self.output_file_name != None:
            file_object.write_file(self.output_file_name, matrix)




class Files(object):
    def __init__(self):
        pass

    def screen_init(self, parameters):
        print "################## MONTE CARLO SIMULATION ###################"
        print "##### Your simulation has the following paramemeters: #######"
        print "\n"
        for key in parameters:
            print "##### {0}: {1}".format(key, parameters[key])
        print "\n"
        print "#############################################################"
        print "Number of steps done:"



    def open_file(self, name):
        proto_matriz = []
        matriz_posicao = []

        for line in open(name, 'r'):
            items = line.rstrip('\n').split()
            proto_matriz.append(items)

        for item in proto_matriz:
            matriz_posicao.append([float(item[0]), float(item[1]), float(item[2])])

        return matriz_posicao


    def write_file(self, name, matrix):
        file = open(name, 'w')
        for item in matrix:
            x = str(item[0])
            y = str(item[1])
            z = str(item[2])
            file.write(x + ' ' + y + ' ' + z + "\n")

        file.close()

    def write_coef(self, name, coefs):
        file = open(name, 'w')
        for c in coefs:
            file.write(str(c) + "\n")
        file.close()



class Difusion:
    def __init__(self):
        self.flux = 0

    def self_difusion(self, disp, time):
        if time != 0:
            return disp / (6 * time) / 6.02e23


    def transport_difusion(self):
        return - self.flux * (parameters.height / - 1000)

    def flux(self, old_pos, new_pos):
        if old_pos[2] < 0 and new_pos[2] > 0: #caso passe na seccao para a esq
            self.flux += 1

        elif old_pos[2] > 0 and new_pos[2] < 0: #caso passe na seccao para a dir
            self.flux -= 1




class Time: #class de tempo
    def __init__(self):
        self.v_quadrado = ((3 * 6.02e23 * 1.3806488e-23 * parameters.temperature) / (2.000)) * 1000
        self.tempo_total = 0


    def refresh_time(self, deslocamento_total):
        tempo = math.sqrt(((deslocamento_total / parameters.n_particles) / self.v_quadrado)) / 10e10
        self.tempo_total += tempo

    def get_total_time(self):
        return self.tempo_total


class Displacements: #class de deslocamentos
    def __init__(self):
        self.d_total = 0
        self.d_ciclo = 0


    def update_displacements_ciclo(self, displacement):
        self.diclo = 0
        for disp in displacement:
            self.d_ciclo += disp**2



    def update_displacements(self):
        self.d_total += self.d_ciclo


    def get_cycle_displacement(self):
        return self.d_ciclo

    def get_total_displacement(self):
        return self.d_total











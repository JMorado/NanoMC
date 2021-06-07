PROGRAM nanomc_uvt 
  USE io_module
  USE simulation_module
  USE energy_module
  USE pbc_displacement_module
  IMPLICIT NONE
  !----------------------------------------------------------------------!
  ! Ad-hoc variables for this program
  !----------------------------------------------------------------------!
  CHARACTER(LEN=100)                  :: input_file
  TYPE(Simulation)                    :: sim
  INTEGER(4)                          :: i, j


  ! Potential related variables
  INTEGER(4)                          :: npoints
  REAL(8)                             :: r(3), potential, r2

  npoints = 200 

  ! Input file name
  READ(5,*) input_file

  ! Restart the random number generator with a given seed
  CALL srand(sim%seed)

  ! Read the input_file
  CALL read_input(sim, input_file)

  ! Determine and use reduced units from this point
  ! TODO: Joao Morado 17.12.2018  
  ! TODO: not correctly implemented yet
  ! CALL sim%generate_reduced_units()

  ! Allocate the position array and generate initial random distribution
  ALLOCATE( sim%coord(3,1))
  sim%coord=0
  sim%coord(3,1) = sim%length/2.0

  ! Choose correct energy functions
  SELECT CASE (sim%model)
  CASE("atomistic")
     CALL read_cnt_xyz(sim)
  END SELECT


  ! Fluid-CNT contribution to the energy
  DO i=1,npoints
     sim%coord(1,1) =  sim%coord(1,1) + 0.001 * i
     potential = 0
     DO j=1,sim%ncnt
        ! Compute distance r
        r = sim%coord(:,1)-sim%coord_cnt(:,j)
        ! Apply PBC to r
        CALL pbc_z_distance(r, sim%length)
        ! Compute r squared of the periodic distance
        r2 = SUM(r*r)

        potential = potential + lj_fluid_cnt_atomistic(sim,r2)
     END DO
     WRITE(6,*) SQRT(SUM(sim%coord(1:2,1)**2)), potential, sim%coord(:,1)
  END DO


  STOP
END PROGRAM nanomc_uvt

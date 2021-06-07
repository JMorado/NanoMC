PROGRAM bd_nvt
  USE brownian_dynamics_module
  USE io_module
  USE simulation_module
  USE std_output_module
  IMPLICIT NONE

  !----------------------------------------------------------------------!
  ! Ad-hoc variables for this program
  !----------------------------------------------------------------------!
  ! Simulation and time objects
  TYPE(Simulation)                    :: sim

  ! Iteration over MC sweeps (i) and particles (o)
  INTEGER(4)                          :: i, o, na

  ! File names
  CHARACTER(LEN=100)                  :: input_file
  CHARACTER(LEN=100)                  :: traj_file_name, prop_file_name
  CHARACTER(LEN=100)                  :: fff

  ! Time
  REAL(8)                             :: t, dt

  t = 0.0d0
  dt = 1e-12

  CALL std_output_initialize()


  ! Input file name
  READ(5,*) input_file

  WRITE(6,'(A)') " * Reading input file"

  ! Read the input_file
  CALL read_input(sim, input_file)

  ! Trajectory and property file name
  traj_file_name = trim(sim%output_seed) // ".xyz"
  prop_file_name = trim(sim%output_seed) // ".dat"


  ! Restart the random number generator with a given seed
  CALL srand(sim%seed)

  ! Write simulations details to standard output
  CALL STD_OUTPUT_SIMULATION_DETAILS(sim)
  ! Determine and use reduced units from this point
  ! TODO: Joao Morado 17.12.2018  
  ! TODO: not correctly implemented yet
  ! CALL sim%generate_reduced_units()

  ! Allocate the position array and generate initial random distribution
  WRITE(6,'(A)') " * Allocating hydrogen coordinate matrix."
  ALLOCATE( sim%coord(3,sim%npart))
  sim%coord=0

  WRITE(6,'(A)') " * Preparing initial system configuration"
  ! Generate random distribution of initial particles
  CALL generate_random_distribution(sim)


  ! TODO:
  fff =  "new.xyz"
  CALL read_initial_xyz(sim, fff)

  ! Allocate the initial coordinate array and array with real z coordinates
  ALLOCATE( sim%coord_init(3,sim%npart), sim%coord_eff(3,sim%npart))
  sim%coord_eff  = sim%coord
  sim%coord_init = sim%coord

  CALL read_cnt_xyz(sim)


  ! Initialize time module
  sim%fluid_mass = 2.0d0*1.6737236d-27
  CALL STD_OUTPUT_STARTING_SIM()

  !--------------------------------------------------------------------------------
  !                              Start BD simulation
  !--------------------------------------------------------------------------------

  t = 0
  DO i=1,sim%nsweeps
     CALL brownian_dynamics_integrator(sim, dt)
     t = t + dt
     !WRITE(6,*) t
     !square_displacement = tim%compute_SD(sim%coord_init, sim%coord_eff, sim%npart, 1)
     !--------------------------------------------------------------------------------
     ! Write properties to output
     IF (MOD(i,sim%ntwx) .eq. 0) THEN
       !  WRITE(6,'(A,I10.2,I10.4)') " MC sweep number: ", i, sim%npart
        ! TODO: implement write_properties function in the io_module
        !CALL write_properties()

     END IF

     ! Write trajectory to output
     IF (MOD(i,sim%ntpr) .eq. 0) THEN
        CALL write_trajectory(sim , traj_file_name, i)
     END IF
     !--------------------------------------------------------------------------------

  END DO


  !--------------------------------------------------------------------------------
  !                                End NVT simulation
  !--------------------------------------------------------------------------------

  CALL STD_OUTPUT_END_SIM()
  STOP
END PROGRAM bd_nvt



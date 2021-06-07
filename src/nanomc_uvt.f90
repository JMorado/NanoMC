PROGRAM nanomc_uvt 
  USE time_module
  USE io_module
  USE simulation_module
  USE energy_module
  USE pbc_module
  USE displacement_module
  USE monte_carlo_module
  USE std_output_module
  IMPLICIT NONE
  !----------------------------------------------------------------------!
  ! Ad-hoc variables for this program
  !----------------------------------------------------------------------!
  CHARACTER(LEN=100)                  :: input_file
  TYPE(Simulation)                    :: sim
  LOGICAL(1)                          :: nvt_accepted, uvt_accepted
  TYPE(Time)                          :: tim
  INTEGER(4)                          :: i, j
  REAL(8), PARAMETER                  :: prob_disp = 1.0 / 3.0
  REAL(8)                             :: n_acc, n_tot, frac_acc
  PROCEDURE(energy_function), POINTER :: energy_func_ptr => null ()
  CHARACTER(LEN=100)                  :: traj_file_name, prop_file_name
  INTEGER(4) :: o


  CALL std_output_initialize()

  ! Number of accepted steps, number of total steps
  n_acc    = 0
  n_tot    = 0

  ! Input file name
  READ(5,*) input_file

  WRITE(6,'(A)') " * Reading input file"
  ! Read the input_file
  CALL read_input(sim, input_file)


  traj_file_name = trim(sim%output_seed) // ".xyz"
  prop_file_name = trim(sim%output_seed) // ".dat"


  ! Restart the random number generator with a given seed
  CALL srand(sim%seed)

  ! Write simulations details to standard output
  CALL STD_OUTPUT_SIMULATION_DETAILS(sim)

  ! Allocate the position array and generate initial random distribution
  WRITE(6,'(A)') " * Allocating hydrogen coordinate matrix."
  ALLOCATE( sim%coord(3,sim%npart+1000))
  sim%coord=0

  WRITE(6,'(A)') " * Preparing initial system configuration"
  ! Generate random distribution of initial particles
  CALL generate_random_distribution(sim)

  ! Allocate the initial coordinate array and array with real z coordinates
  ALLOCATE( sim%coord_init(3,sim%npart+1000), sim%coord_eff(3,sim%npart+1000))
  sim%coord_eff = sim%coord

  ! Choose correct energy functions
  SELECT CASE (sim%model)
  CASE("atomistic")
     ! Read CNT coordinates
     WRITE(6,'(A)') " * Branched into CNT atomistic model."
     CALL read_cnt_xyz(sim)
     energy_func_ptr => particle_energy_atomistic
  CASE("continuum")
     WRITE(6,'(A)') " * Branched into CNT continuum model."
     energy_func_ptr => particle_energy_continuum
  END SELECT

  WRITE(6,'(A,F16.8)') " * Reservoir pressure: ", sim%pressure_reservoir
  CALL STD_OUTPUT_STARTING_SIM()
  WRITE(6,'(A20,A20)') adjustl(" MC sweep number"), adjustl("Number of particles")


  DO i=1,sim%nsweeps
     IF (RAND() .lt. prob_disp) THEN
        ! Select a particle at random
        o = INT(sim%npart*RAND())+1
        CALL mc_move(sim, energy_func_ptr, nvt_accepted,o)
     ELSE
        IF (sim%npart .gt. 0) THEN
           CALL mc_exchange(sim, energy_func_ptr)
        END IF
     END IF

     ! Write properties to output
     IF (MOD(i,sim%ntwx) .eq. 0) THEN
        WRITE(6,'(I20,I20)') i, sim%npart
        ! TODO: implement write_properties function in the io_module
        !CALL write_properties()
     END IF

     ! Write trajectory to output
     IF (MOD(i,sim%ntpr) .eq. 0) THEN
        CALL write_trajectory(sim , traj_file_name, i, .true.)
     END IF

  END DO

  CALL STD_OUTPUT_END_SIM()

  ! frac_acc = n_acc / n_tot * 100.0
  ! WRITE(6,'(A,F8.4)') "Percentage of accepted trial moves: ", frac_acc

  !CALL write_xyz(sim)
  STOP
END PROGRAM nanomc_uvt

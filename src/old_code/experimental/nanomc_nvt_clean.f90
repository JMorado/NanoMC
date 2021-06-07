PROGRAM nanomc_nvt 
  USE time_module
  USE io_module
  USE simulation_module
  USE energy_module
  USE pbc_displacement_module
  USE monte_carlo_module
  USE std_output_module
  IMPLICIT NONE

  !----------------------------------------------------------------------!
  ! Ad-hoc variables for this program
  !----------------------------------------------------------------------!
  ! Simulation and time objects
  TYPE(Simulation)                    :: sim
  TYPE(Time)                          :: tim

  ! Iteration over MC sweeps (i) and particles (o)
  INTEGER(4)                          :: i, o, na

  ! Pointer for the energy function
  PROCEDURE(energy_function), POINTER :: energy_func_ptr => null ()

  ! File names
  CHARACTER(LEN=100)                  :: input_file
  CHARACTER(LEN=100)                  :: traj_file_name, prop_file_name
  CHARACTER(LEN=100)                  :: fff

  ! MC acceptance statistics
  REAL(8)                             :: n_acc, n_tot, frac_acc
  LOGICAL(1)                          :: accepted

  ! Diffusion
  REAL(8)                             :: sd, dif

  CALL std_output_initialize()

  ! Number of accepted steps, number of total steps
  n_acc    = 0.0d0
  n_tot    = 0.0d0

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


  fff =  "new.xyz"
  CALL read_initial_xyz(sim, fff)

  ! Allocate the initial coordinate array and array with real z coordinates
  ALLOCATE( sim%coord_init(3,sim%npart), sim%coord_eff(3,sim%npart))
  sim%coord_eff  = sim%coord
  sim%coord_init = sim%coord

  ! Choose correct energy functions
  SELECT CASE (sim%model)
  CASE("atomistic")
     ! Read CNT coordinates
     WRITE(6,'(A)') " * Branched into CNT atomistic model."
     CALL read_cnt_xyz(sim)
     energy_func_ptr => particle_energy_atomistic
  CASE("continuum")
     WRITE(6,'(A)') " * Branched into CNT continuum model."
     energy_func_ptr => particle_energy
  END SELECT

  ! Initialize time module
  sim%fluid_mass = 2.0d0*1.6737236d-27
  CALL tim%initialization(sim%fluid_mass, sim%temperature)
  CALL STD_OUTPUT_STARTING_SIM()

  !--------------------------------------------------------------------------------
  !                              Start NVT simulation
  !--------------------------------------------------------------------------------
  DO i=1,sim%nsweeps
     DO o=1,sim%npart
        !TODO:put o back
        CALL mc_move(sim, energy_func_ptr, accepted, o)

        IF (accepted) THEN
           n_acc = n_acc + 1
        END IF
        n_tot = n_tot + 1
     END DO

     frac_acc = n_acc / n_tot
     CALL tim%update_time(1.0d0)
     sd = tim%compute_SD(sim%coord_init,sim%coord_eff,sim%npart, 1)
     dif = tim%self_diffusion(sim%npart, 1, sd,tim%t )
     WRITE(6,*) i, frac_acc, tim%t, sd, dif
     !--------------------------------------------------------------------------------
     ! Write properties to output
     IF (MOD(i,sim%ntwx) .eq. 0) THEN
         WRITE(6,'(A,I10.2,I10.4,F10.4)') " MC sweep number: ", i, sim%npart, frac_acc
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
  n_tot = sim%nsweeps * sim %npart
  frac_acc = n_acc / n_tot * 100.0
  WRITE(6,'(A,F8.4)') "Percentage of accepted trial moves: ", frac_acc
  STOP
END PROGRAM nanomc_nvt



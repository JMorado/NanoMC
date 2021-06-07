PROGRAM nanomc_nvt 
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

  ! Ad-hoc time variablesn 
  REAL(8) :: t1, t2, t3, t4, t5, t6
  REAL(8) :: t1sweep, t2sweep, t3sweep, t4sweep, t5sweep, t6sweep
  ! DIF COEFF
  REAL(8) :: d1, d2, d3, d4, d5, d6
  REAL(8) :: square_displacement
  REAL(8) :: pideal, preal, lambda_real, lambda_ideal, lambda_virial
  REAL(8) :: density, density0
  REAL(8) :: avg_virial, knudsen

  CALL std_output_initialize()

  ! Number of accepted steps, number of total steps
  n_acc    = 0
  n_tot    = 0

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
     energy_func_ptr => particle_energy_continuum
  END SELECT

  ! Initialize time module
  sim%fluid_mass = 2.0d0*1.6737236d-27
  CALL tim%initialization(sim%fluid_mass, sim%temperature)
  CALL STD_OUTPUT_STARTING_SIM()

  sim%disp = 3.0
  !--------------------------------------------------------------------------------
  !                              Start NVT simulation
  !--------------------------------------------------------------------------------


 ! density0 = (sim%npart / 6.022140d23) / (sim%volume * 1d-30)
  
 ! density = density0 + density0 * density0 * 1.4232327032795003d-05

 ! lambda_virial = 1.0d0 / (SQRT(2d0)* 3.14159265359 * (sim%sigma_fluid*1d-10) ** 2 * density * 6.022140d23 )
 ! lambda_virial = lambda_virial * 1e10 !m to A 
 !   preal = (sim%npart * sim%temperature * 1.381d-23 + virial(sim) * 1.381d-23) / (sim%volume*1d-30) 
 ! preal = (virial(sim) + sim%npart * sim%temperature) * 1.381d-23 / (sim%volume*1d-30)
 ! lambda_real = 1.38064852d-23 * sim%temperature / (SQRT(2d0)* 3.14159265359 * sim%sigma_fluid ** 2 * preal * 1d-30)

 ! pideal = sim%npart * sim%temperature * 1.381d-23 / (sim%volume*1d-30)



 !TODO  lambda_ideal = 1.38064852d-23 * sim%temperature / ( sqrt(2) * 3.14159265359 * sim%sigma_fluid ** 2 * pideal * 1d-30)
 



 ! write(6,*) lambda_ideal, lambda_real, lambda_virial
  
  !sim%disp= lambda_viria
  
 ! WRITE(6,*) pideal*1e-5, preal*1e-5, lambda_real, lambda_ideal


  t1 = 0
  t2 = 0
  t3 = 0
  t4 = 0
  t5 = 0
  t6 = 0

  avg_virial = 0
  DO i=1,sim%nsweeps
     t1sweep = 0
     t2sweep = 0
     t3sweep = 0
     t4sweep = 0
     t5sweep = 0
     t6sweep = 0


     !TODO:COMMENT IF NOT NECESSARY

    !avg_virial = avg_virial + virial(sim)
    !preal = (sim%npart * sim%temperature * 1.381d-23 + virial(sim)) / (sim%volume*1d-30)
    !preal = (virial(sim) / sim%npart + sim%npart * sim%temperature) * 1.381d-23 / (sim%volume*1d-30)
    !lambda_real = 1.38064852d-23 * sim%temperature / (SQRT(1d0)* 3.14159265359 * sim%sigma_fluid ** 2 * preal * 1d-30)

     !knudsen =  lambda_ideal / (2*sim%radius)
     !WRITE(6,*) "KNUDSEN", knudsen
    !sim%disp = ( 1.0d0 / (1 + exp(-1.0*(knudsen-0.5d0))) ) * (1d0 / (2*sim%radius)) &
    !     + ( 1.0d0 / (1 + exp(1.0*(knudsen-0.5d0)) )) * (1.0d0/lambda_ideal)
    !sim%disp = 1.0d0/sim%disp

    !sim%disp =  (1.0/lambda_ideal)
    !sim%disp = 1.0d0 / sim%disp
    !sim%disp =  lambda_real

    !sim%disp= lambda_ideal!2*sim%radius!lambda_ideal!((lambda_ideal* 2*sim%radius)) /  (lambda_ideal+ 2*sim%radius)
    na = 0
     DO o=1,sim%npart
        !TODO:put o back
        CALL mc_move(sim, energy_func_ptr, accepted, o)

        IF (accepted) THEN
           ! MC move was accepted
           t1sweep = t1sweep +  SQRT(SUM(sim%dr**2)) / tim%v_rrels
           t2sweep = t2sweep + SQRT(SUM(sim%dr**2)) / tim%v_mean
           t3sweep = t3sweep + SQRT(SUM(sim%dr**2)) / tim%v_rms

           t4sweep = t4sweep + SUM(sim%dr**2) 
           t5sweep = t5sweep + SUM(sim%dr**2)
           t6sweep = t6sweep + SUM(sim%dr**2)

           n_acc = n_acc + 1
           na = na + 1
        END IF
     END DO

     ! Determine tsweep
     if (na .ne. 0) then
     na = sim%npart
     t1sweep = t1sweep / na
     t2sweep = t2sweep / na
     t3sweep = t3sweep / na


     t4sweep = SQRT(t4sweep / na) / tim%v_rrels
     t5sweep = SQRT(t5sweep / na) / tim%v_mean
     t6sweep = SQRT(t6sweep / na) / tim%v_rms



  else
     t1sweep=0
     t2sweep=0
     t3sweep=0
     t4sweep=0
     t5sweep=0
     t6sweep=0

     end if


     !t4sweep = SQRT(t4sweep / sim%npart) / tim%v_rrels
     !t5sweep = SQRT(t5sweep / sim%npart) / tim%v_mean
     !t6sweep = SQRT(t6sweep / sim%npart) / tim%v_rms

     ! Update total time
     t1 = t1 + t1sweep
     t2 = t2 + t2sweep
     t3 = t3 + t3sweep
     t4 = t4 + t4sweep
     t5 = t5 + t5sweep
     t6 = t6 + t6sweep



     ! TODO: CALL tim%add_time(tcycle)
     square_displacement = tim%compute_SD(sim%coord_init, sim%coord_eff, sim%npart, 1)

     d1 = tim%self_diffusion(sim%npart, 1, square_displacement, t1)
     d2 = tim%self_diffusion(sim%npart, 1, square_displacement, t2)
     d3 = tim%self_diffusion(sim%npart, 1, square_displacement, t3)
     d4 = tim%self_diffusion(sim%npart, 1, square_displacement, t4)
     d5 = tim%self_diffusion(sim%npart, 1, square_displacement, t5)
     d6 = tim%self_diffusion(sim%npart, 1, square_displacement, t6)

     WRITE(6,*) "--- MC SWEEP", i, sim%disp, lambda_real, lambda_ideal, lambda_ideal / (2*sim%radius)
     WRITE(6,*) "1-Ds", t1, square_displacement, d1
     WRITE(6,*) "2-Ds", t2, square_displacement, d2
     WRITE(6,*) "3-Ds", t3, square_displacement, d3
     WRITE(6,*) "4-Ds", t4, square_displacement, d4
     WRITE(6,*) "5-Ds", t5, square_displacement, d5
     WRITE(6,*) "6-Ds", t6, square_displacement, d6


     !--------------------------------------------------------------------------------
     ! Write properties to output
     IF (MOD(i,sim%ntwx) .eq. 0) THEN
         WRITE(6,'(A,I10.2,I10.4)') " MC sweep number: ", i, sim%npart
        ! TODO: implement write_properties function in the io_module
        !CALL write_properties()
     END IF

     ! Write trajectory to output
     IF (MOD(i,sim%ntpr) .eq. 0) THEN
        CALL write_trajectory(sim , traj_file_name, i, .true.)
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



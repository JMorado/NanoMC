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
    REAL(8)                             :: n_acc_nvt, n_acc_uvt, n_steps_nvt, n_steps_uvt
    INTEGER(4) :: o


    CALL std_output_initialize()

    ! Number of accepted steps, number of total steps
    n_acc_nvt = 0.0
    n_acc_uvt = 0.0
    n_steps_nvt = 0.0
    n_steps_uvt = 0.0

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
    WRITE(6,'(A15,A15,A15,A15,A15,A15,A15)') " MC sweep", "# particles", "E (kJ/mol)",&
            "Acc. rate NVT", "Acc. rate uVT", "Frac. NVT", "Frac. uVT"


    DO i=1,sim%nsweeps
        IF (RAND() .lt. prob_disp) THEN
            ! Select a particle at random
            o = INT(sim%npart*RAND())+1
            CALL mc_move(sim, energy_func_ptr, nvt_accepted,o)

            IF (nvt_accepted) THEN
                n_acc_nvt = n_acc_nvt + 1.0
            END IF
            n_steps_nvt = n_steps_nvt + 1.0
        ELSE
            IF (sim%npart .gt. 0) THEN
                CALL mc_exchange(sim, energy_func_ptr, uvt_accepted)

                IF (uvt_accepted) THEN
                    n_acc_uvt = n_acc_uvt + 1.0
                END IF
                n_steps_uvt = n_steps_uvt + 1.0
            END IF
        END IF

        ! Write properties to output
        IF (MOD(i,sim%ntwx) .eq. 0) THEN
            WRITE(6,'(I15,I15,E15.6, F15.4, F15.4, F15.4, F15.4)') i, sim%npart, total_system_energy(sim),&
                    n_acc_nvt/n_steps_nvt, n_acc_uvt/n_steps_uvt, n_steps_nvt/i, n_steps_uvt/i
        END IF

        ! Write trajectory to output
        IF (MOD(i,sim%ntpr) .eq. 0) THEN
            CALL write_trajectory(sim , traj_file_name, i, .true.)
        END IF

    END DO

    CALL write_xyz(sim, "output_final.xyz")
    CALL STD_OUTPUT_END_SIM()

    STOP
END PROGRAM nanomc_uvt

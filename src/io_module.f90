MODULE io_module
    USE simulation_module
    USE constants_module
    USE cell_module
    IMPLICIT NONE

CONTAINS
    SUBROUTINE read_input(simulation_instance, input_file)
        !========================================================================!
        ! This functions read the input NanoMC file                              !
        !------------------------------------------------------------------------!
        ! simulation_instance (in) : Simulation instance                         !
        ! input_file          (in) : name of the NanoMC input file               !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation),  INTENT(INOUT) :: simulation_instance
        CHARACTER(LEN=50), INTENT(IN)    :: input_file
        INTEGER(4)                       :: input_seed

        input_seed = 100
        OPEN(input_seed,file=input_file)

        ! Simulation variables
        READ(input_seed,*) simulation_instance%model
        READ(input_seed,*) simulation_instance%displacement_type
        READ(input_seed,*) simulation_instance%ensemble

        ! Ensemble and Model specific output part
        SELECT CASE (simulation_instance%ensemble)
        CASE ("nvt")
            SELECT CASE (simulation_instance%model)
            CASE ("atomistic")
                READ(input_seed,*) simulation_instance%cnt_file
            END SELECT

        CASE ("uvt")
            READ(input_seed,*) simulation_instance%chemical_pot              ! Note that NanoMC read chemical potentials in units of Kelvin (they are divided by kB already)
            SELECT CASE (simulation_instance%model)
            CASE ("atomistic")
                READ(input_seed,*) simulation_instance%cnt_file
            END SELECT
        END SELECT

        READ(input_seed,*) simulation_instance%npart
        READ(input_seed,*) simulation_instance%nsweeps
        READ(input_seed,*) simulation_instance%temperature

        ! CNT parameters
        READ(input_seed,*) simulation_instance%length
        READ(input_seed,*) simulation_instance%radius

        ! Cell
        READ(input_seed,*) simulation_instance%cell_instance%a, simulation_instance%cell_instance%b, &
                simulation_instance%cell_instance%c, simulation_instance%cell_instance%alpha, &
                simulation_instance%cell_instance%beta, simulation_instance%cell_instance%gamma

        ! Lennard-Jones 12-6 parameters
        READ(input_seed,*) simulation_instance%fluid_cut_off
        READ(input_seed,*) simulation_instance%cnt_cut_off
        READ(input_seed,*) simulation_instance%sigma_fluid
        READ(input_seed,*) simulation_instance%eps_fluid
        READ(input_seed,*) simulation_instance%sigma_cnt
        READ(input_seed,*) simulation_instance%eps_cnt

        ! Seed variables
        READ(input_seed,*) simulation_instance%output_seed
        READ(input_seed,*) simulation_instance%seed

        ! Output variables
        READ(input_seed,*) simulation_instance%ntwx
        READ(input_seed,*) simulation_instance%ntpr

        ! Input configuration variables
        READ(input_seed,*) simulation_instance%initial_config_type

        ! Read input file name if initial_config_mode is file
        SELECT CASE (simulation_instance%initial_config_type)
        CASE("file")
            READ(input_seed,*) simulation_instance%initial_xyz_file
        END SELECT

        ! Close file
        CLOSE(input_seed)

        ! Compute input dependent variable
        simulation_instance%beta = 1.0d0 / (kb*simulation_instance%temperature)

        simulation_instance%volume = pi * simulation_instance%radius * simulation_instance%radius * simulation_instance%length
        simulation_instance%volume_eff = pi * (simulation_instance%radius-simulation_instance%sigma_cnt)**2 &
                * simulation_instance%length

        ! Density calculation
        simulation_instance%density_fuid = simulation_instance%npart / simulation_instance%volume
        simulation_instance%density_cnt = simulation_instance%ncnt / simulation_instance%volume


        SELECT CASE (simulation_instance%ensemble)
        CASE("uvt")
            simulation_instance%volume = pi * simulation_instance%radius * simulation_instance%radius &
                    * simulation_instance%length
            simulation_instance%thermal_wav = (planck*planck*simulation_instance%beta) / (2*pi*2*1.00794*atomicmass_to_kg)
            simulation_instance%thermal_wav = (simulation_instance%thermal_wav) ** (1/2.)

            simulation_instance%activity = exp(simulation_instance%chemical_pot / simulation_instance%temperature) &
                    / (simulation_instance%thermal_wav ** 3.0)

            simulation_instance%pressure_reservoir = exp(simulation_instance%chemical_pot/simulation_instance%temperature) &
                    / (simulation_instance%beta * simulation_instance%thermal_wav ** 3.0) 
        END SELECT
    END SUBROUTINE read_input

    SUBROUTINE read_initial_xyz(simulation_instance, input_file)
        !========================================================================!
        ! This subroutines reads the initial configuration from a .xyz file      !
        ! It accepts .xyz files which contain both CNT and H atoms               !
        ! with the requirement that the CNT comes first.                         !
        !------------------------------------------------------------------------!
        ! simulation_instance (in) : Simulation instance                         !
        ! input_file          (in) : name of the xyz file                        !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation), INTENT(INOUT)  :: simulation_instance
        CHARACTER(LEN=50), INTENT(IN)    :: input_file
        INTEGER(4)                       :: k, l
        INTEGER(4)                       :: input_seed, xyz_nat
        CHARACTER(LEN=2)                 :: dummy
        REAL(8)                          :: tmp(3)


        input_seed = 101

        OPEN(input_seed, file=input_file)

        READ(input_seed,*) xyz_nat
        READ(input_seed,*)

        IF (xyz_nat .eq. simulation_instance%npart) THEN
            ! File only contains hydrogen atoms
            DO k=1,simulation_instance%npart
                READ(input_seed,*) dummy, (simulation_instance%coord(l,k),l=1,3)
            END DO
        ELSE
            ! .xyz file contains the CNT, ignore it and just read the hydrogen molecules c.o.m.
            READ(input_seed,*) dummy, (tmp(k),k=1,3)
            IF (dummy .ne. "C") THEN
                WRITE(6,*) "Number of hydrogen molecules not equal to the number of atoms in &
                        the .xyz file and first atom type is not C."
                WRITE(6,*) "If the CNT is present it has to come before the hydrogen molecules."
            ELSE
                DO k=2,(xyz_nat-simulation_instance%npart)
                    READ(input_seed,*) dummy, (tmp(l),l=1,3)
                END DO

                DO k=1,simulation_instance%npart
                    READ(input_seed,*) dummy, (simulation_instance%coord(l,k),l=1,3)
                END DO
            END IF
        END IF
        CLOSE(input_seed)
    END SUBROUTINE read_initial_xyz


    SUBROUTINE read_cnt_xyz(simulation_instance)
        !========================================================================!
        ! This subroutines reads the CNT xyz file n from a .xyz file.            !
        !------------------------------------------------------------------------!
        ! simulation_instance (in) : Simulation instance                         !
        !========================================================================!
        IMPLICIT NONE

        TYPE(Simulation), INTENT(INOUT) :: simulation_instance
        CHARACTER(LEN=20) :: dummy
        INTEGER(4) :: k,l

        OPEN(101,file=simulation_instance%cnt_file)
        READ(101,*) simulation_instance%ncnt
        READ(101,*) dummy ! Dummy
        ALLOCATE(simulation_instance%coord_cnt(3,simulation_instance%ncnt))
        DO k=1,simulation_instance%ncnt
            READ(101,*) dummy, (simulation_instance%coord_cnt(l,k),l=1,3)
        END DO
        CLOSE(101)
    END SUBROUTINE read_cnt_xyz

    SUBROUTINE write_matrix(matrix, output_file)
        !========================================================================!
        ! This subroutine write a 2-dimensional matrix to output_file.           !
        !                                                                        !
        !   output_file options are:                                             !
        ! - std_out : write matrix to standard output                            !
        ! - if present and not std_out: write matrix to output_file              !
        ! - if not present: write matrix to "output.matrix"                      !
        !------------------------------------------------------------------------!
        ! matrix(:,:)                   (in) : matrix                            !
        ! output_file                   (in) : described above                   !
        !========================================================================!
        IMPLICIT NONE

        REAL(8), INTENT(IN)                            :: matrix(:,:)
        CHARACTER(LEN=100), INTENT(INOUT), OPTIONAL    :: output_file

        INTEGER(4)                                     :: k,l, nr, nc
        INTEGER(4)                                     :: output_seed

        nr = SIZE(matrix,1)
        nc = SIZE(matrix,2)

        IF (PRESENT(output_file)) THEN
            IF (output_file .eq. "std_out") THEN
                output_seed = 6
                DO k=1,nr
                    WRITE(output_seed,*) (matrix(l,k),l=1,nc)
                END DO
                return
            END IF
        ELSE
            output_file = "output.matrix"
        END IF

        output_seed = 101
        OPEN(output_seed,file="output.matrix")
        DO k=1,nr
            WRITE(output_seed,*) (matrix(l,k),l=1,nc)
        END DO
        CLOSE(output_seed)

        return
    END SUBROUTINE write_matrix


    SUBROUTINE write_xyz(simulation_instance, output_file_name)
        !========================================================================!
        ! This subroutines reads the a xyz file from output_file_name            !
        !------------------------------------------------------------------------!
        ! simulation_instance (in) : Simulation instance                         !
        ! output_file_name    (in) : name of the xyz file                        !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation), INTENT(IN)                   :: simulation_instance
        CHARACTER(LEN=100), INTENT(INOUT), OPTIONAL    :: output_file_name
        CHARACTER(LEN=100)                             :: output_file
        INTEGER(4)                                     :: output_seed
        INTEGER(4)                                     :: k,l

        output_seed = 101

        IF (.not. PRESENT(output_file_name)) THEN
            output_file = "cnt_hydrogen.xyz"
        ELSE
            output_file = output_file_name
        END IF

        OPEN(output_seed,file=output_file)

        ! Choose correct energy functions
        SELECT CASE (simulation_instance%model)
        CASE("atomistic")
            WRITE(output_seed,*) int(simulation_instance%ncnt+simulation_instance%npart)
            WRITE(output_seed,*) "CNT atomistic + Hydrogen"
            DO k=1,simulation_instance%ncnt
                WRITE(output_seed,*) "C", (simulation_instance%coord_cnt(l,k),l=1,3)
            END DO
        CASE("continuum")
            WRITE(output_seed,*) int(simulation_instance%npart)
            WRITE(output_seed,*) "CNT continuum + Hydrogen"
        END SELECT

        DO k=1,simulation_instance%npart
            WRITE(output_seed,*) "H", (simulation_instance%coord(l,k),l=1,3)
        END DO

        CLOSE(output_seed)
        return
    END SUBROUTINE write_xyz

    SUBROUTINE write_trajectory(simulation_instance, output_file, step, no_cnt)
        !========================================================================!
        ! This subroutines writes the simulation trajectory                      !
        !------------------------------------------------------------------------!
        ! simulation_instance (in) : Simulation instance                         !
        ! output_file         (in) : name of the xyz file                        !
        ! step                (in) : current simulation step                     !
        ! no_cnt              (in) : flag that signals if CNT is to be written   !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation),   INTENT(IN)                 :: simulation_instance
        CHARACTER(LEN=100), INTENT(IN)                 :: output_file
        INTEGER(4),         INTENT(IN)                 :: step
        LOGICAL(4),         INTENT(IN)                 :: no_cnt
        CHARACTER(LEN=100)                             :: comment
        INTEGER(4)                                     :: output_seed, k, l
        LOGICAL                                        :: exist

        output_seed = 101

        INQUIRE(file=output_file, exist=exist)
        IF (exist) THEN
            OPEN(output_seed, file=output_file, status="old", position="append", action="write")
        ELSE
            OPEN(output_seed, file=output_file, status="new", action="write")
        END IF

        IF (no_cnt) THEN
            WRITE (comment, "(I10.4)") step
            WRITE(output_seed,*) (simulation_instance%npart)
            WRITE(output_seed,*) comment

            DO k=1,simulation_instance%npart
                WRITE(output_seed, *) "H ", (simulation_instance%coord(l,k),l=1,3)
            END DO
        ELSE
            WRITE (comment, "(A,I10.4)") "Sweep number", step
            WRITE(output_seed,*) (simulation_instance%npart+simulation_instance%ncnt)
            WRITE(output_seed,*) comment

            DO k=1,simulation_instance%ncnt
                WRITE(output_seed, *) "C ", (simulation_instance%coord_cnt(l,k),l=1,3)
            END DO

            DO k=1,simulation_instance%npart
                WRITE(output_seed, *) "H ", (simulation_instance%coord(l,k),l=1,3)
            END DO
        END IF

        CLOSE(output_seed)
    END SUBROUTINE write_trajectory
END MODULE io_module

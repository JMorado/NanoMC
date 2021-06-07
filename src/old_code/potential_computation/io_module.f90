MODULE io_module
  USE simulation_module
  IMPLICIT NONE

  REAL(8), PARAMETER :: kb = 1.38064852d-23          ! m^2 kg s^-2 K-1
  REAL(8), PARAMETER :: planck = 6.62607004d-34      ! m^2 kg s^-1
  REAL(8), PARAMETER :: pi= 4.0d0 * ATAN(1.0d0)


  CONTAINS
    SUBROUTINE read_input(this, input_file)
      IMPLICIT NONE
      !------------------------------------------------------------------------
      !
      !
      !-------------------------------------------------------------------------
      CLASS(Simulation), INTENT(INOUT) :: this
      CHARACTER(LEN=50), INTENT(IN)    :: input_file
      INTEGER(4)                       :: input_seed


      input_seed = 100
      OPEN(input_seed,file=input_file)

      ! Simulation variables
      READ(input_seed,*) this%model
      READ(input_seed,*) this%ensemble

      ! Ensemble and Model specific output part
      SELECT CASE (this%ensemble)
         CASE ("nvt")
            SELECT CASE (this%model)
            CASE ("atomistic")
               READ(input_seed,*) this%cnt_file
            END SELECT

         CASE ("uvt")
            READ(input_seed,*) this%chemical_pot
            SELECT CASE (this%model)
            CASE ("atomistic")
               READ(input_seed,*) this%cnt_file
            END SELECT
      END SELECT

      READ(input_seed,*) this%npart
      READ(input_seed,*) this%nsweeps
      READ(input_seed,*) this%temperature

      ! TODO: Joao Morado 17.12.2018
      ! TODO: include xy and z boundary conditions
      ! TODO: include different initial distribution generators
      ! CNT parameters
      READ(input_seed,*) this%length
      READ(input_seed,*) this%radius

      ! Lennard-Jones 12-6 parameters
      READ(input_seed,*) this%fluid_cut_off
      READ(input_seed,*) this%cnt_cut_off
      READ(input_seed,*) this%sigma_fluid
      READ(input_seed,*) this%eps_fluid
      READ(input_seed,*) this%sigma_cnt
      READ(input_seed,*) this%eps_cnt

      ! Seed variables
      READ(input_seed,*) this%output_seed
      READ(input_seed,*) this%seed

      ! Output varialbes
      READ(input_seed,*) this%ntwx
      READ(input_seed,*) this%NTPR

      ! Input configuration variables
      READ(input_seed,*) initial_config_mode


      ! Read input file name if initial_config_mode is file
      SELECT CASE (this%initial_config_mode)
      CASE("file")
         READ(input_seed,*) initial_xyz_file
      END SELECT

      ! Close file
      CLOSE(input_seed)

      ! Compute input dependent variable
      this%beta = 1.0d0 / (kb*this%temperature)


      SELECT CASE (this%ensemble)
      CASE("uvt")
         this%volume = pi * (this%radius)**2 * this%length
         this%volume_eff = pi * (this%radius-this%sigma_cnt)**2 * this%length
         this%thermal_wav = (planck*planck*this%beta) / (2*pi*2*1.6737236d-27)
         this%thermal_wav = (this%thermal_wav) ** (1/2.)
         ! Activity z
         this%activity = exp(this%chemical_pot / this%temperature) / (this%thermal_wav ** 3.0)
         this%pressure = exp(this%chemical_pot/this%temperature) / (this%beta * this%thermal_wav ** 3.0)
      END SELECT


   
    END SUBROUTINE read_input

    SUBROUTINE read_initial_xyz(this, input_file)
      IMPLICIT NONE
      !------------------------------------------------------------------------
      ! This subroutines read the initial configuration from a .xyz file.
      ! It accepts .xyz files which contain both CNT and H atoms
      ! with the requirement that the CNT comes first.
      !------------------------------------------------------------------------
      CLASS(Simulation), INTENT(INOUT) :: this
      CHARACTER(LEN=50), INTENT(IN)    :: input_file
      INTEGER(4)                       :: k, l
      INTEGER(4)                       :: input_seed, xyz_nat
      CHARACTER(LEN=2)                 :: dummy
      REAL(8)                          :: tmp(3)

      
      input_seed = 101

      OPEN(input_seed, file=input_file)

      READ(input_seed,*) xyz_nat
      READ(input_seed,*)

      IF (xyz_nat .eq. this%nat) THEN
         ! File only contains hydrogen atoms
         DO k=1,this%nat
            READ(input_seed,*) dummy, (this%coord(l,k),l=1,3)
         END DO
      ELSE
         ! .xyz file contains the CNT, ignore it and just read the hydrogen molecules c.o.m.
         READ(input_seed,*) dummy, (tmp(k),k=1,3)
         IF (dummy .ne. "C") THEN
            WRITE(6,*) "Number of hydrogen molecules not equal to the number of atoms in the .xyz file and first atom type is not C."
            WRITE(6,*) "If the CNT is present it has to come before the hydrogen molecules."
         ELSE
            DO k=2,(xyz_nat-this%nat)
               READ(input_seed,*) dummy, (tmp(l),l=1,3)
            END DO

            DO k=1,this%nat
               READ(input_seed,*) dummy, (this%coord(l,k),l=1,3)
            END DO
         END IF
      END IF
      CLOSE(input_seed)

    END SUBROUTINE read_initial_xyz


    SUBROUTINE read_cnt_xyz(this)
      IMPLICIT NONE
      !------------------------------------------------------------------------
      ! This subroutines read the CNT .xyz file.
      ! The .xyz file name is stored in the this%cnt file variable.
      ! The CNT atom number is stored in the this%cnt variable.
      ! The CNT coordinates are stored in the this%coord_cnt(3,this%nct) array.
      !------------------------------------------------------------------------
      CLASS(Simulation), INTENT(INOUT) :: this
      CHARACTER(LEN=20) :: dummy
      INTEGER(4) :: k,l

      OPEN(101,file=this%cnt_file)
      READ(101,*) this%ncnt
      READ(101,*) dummy ! Dummy
      ALLOCATE(this%coord_cnt(3,this%ncnt))
      DO k=1,this%ncnt
         READ(101,*) dummy, (this%coord_cnt(l,k),l=1,3)
      END DO
      CLOSE(101)
    END SUBROUTINE read_cnt_xyz

    SUBROUTINE write_matrix(matrix, output_file)
      IMPLICIT NONE
      !------------------------------------------------------------------------
      ! This subroutines writes matrix to output_file.
      ! output_file options:
      ! - std_out : write matrix to standard output
      ! - if present and not std_out: write matrix to output_file
      ! - if not present: write matrix to "output.matrix"
      !------------------------------------------------------------------------
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
    END SUBROUTINE


    SUBROUTINE write_xyz(this, output_f)
      IMPLICIT NONE
      !------------------------------------------------------------------------
      ! This subroutine write an .xyz file with name "output.xyz"
      ! The writte .xyz file will first contain the H atoms
      ! and then the SWCNT C atoms.
      !------------------------------------------------------------------------
      CLASS(Simulation), INTENT(IN)                  :: this
      CHARACTER(LEN=100), INTENT(INOUT), OPTIONAL    :: output_f
      CHARACTER(LEN=100)                             :: output_file
      INTEGER(4)                                     :: output_seed
      INTEGER(4)                                     :: k,l

      output_seed = 101

      IF (.not. PRESENT(output_f)) THEN
         output_file = "cnt_hydrogen.xyz"
      ELSE
         output_file = output_f
      END IF

      OPEN(output_seed,file=output_file)

      ! Choose correct energy functions
      SELECT CASE (this%model)
      CASE("atomistic")
         WRITE(output_seed,*) int(this%ncnt+this%npart)
         WRITE(output_seed,*) "CNT atomistic + Hydrogen"
         DO k=1,this%ncnt
            WRITE(output_seed,*) "C", (this%coord_cnt(l,k),l=1,3)
         END DO
      CASE("continuum")
         WRITE(output_seed,*) int(this%npart)
         WRITE(output_seed,*) "CNT continuum + Hydrogen"
      END SELECT

      DO k=1,this%npart
         WRITE(output_seed,*) "H", (this%coord(l,k),l=1,3)
      END DO
    
      CLOSE(output_seed)
      return
    END SUBROUTINE write_xyz
END MODULE io_module

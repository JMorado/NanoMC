MODULE simulation_module
  IMPLICIT NONE

  TYPE, PUBLIC :: Simulation
     ! Simulation variables
     CHARACTER(LEN=100)    :: output_seed                      ! Output file seed
     CHARACTER(LEN=100)    :: input_file                       ! Output file seed
     CHARACTER(LEN=20)     :: model                            ! Simulation model
     CHARACTER(LEN=3)      :: ensemble                         ! Ensemble
     INTEGER(4)            :: nsweeps                          ! Number of MC sweeps
     INTEGER(4)            :: npart                            ! Number of particles
     INTEGER(4)            :: seed                             ! Random generator seed
     REAL(8)               :: temperature                      ! Temperature

     ! CNT parameters
     REAL(8)               :: length                           ! SWCNT length
     REAL(8)               :: radius                           ! SWCNT radius

     ! Lennard-Jones 12-6 parameters 
     REAL(8)               :: fluid_cut_off                    ! h2-h2 rcut_off 
     REAL(8)               :: cnt_cut_off                      ! h2-swcnt rcut_off
     REAL(8)               :: sigma_fluid                      ! sigma h2-h2 
     REAL(8)               :: eps_fluid                        ! epsilon h2-h2
     REAL(8)               :: sigma_cnt                        ! sigma h2-swcnt
     REAL(8)               :: eps_cnt                          ! epsiolon h2-swcnt

     ! Grand canonical variables
     REAL(8)               :: beta                             ! 1/(kb*T)
     REAL(8)               :: pressure                         ! pressure
     REAL(8)               :: pressure_reservoir
     REAL(8)               :: volume                           ! volume
     REAL(8)               :: volume_eff
     REAL(8)               :: thermal_wav                      ! De Broglie thermal wavelength
     REAL(8)               :: chemical_pot                     ! Chemical potential
     REAL(8)               :: activity                         ! z = exp(beta*chemical_pot) / thermal_wav ** 3


     !
     REAL(8)               :: density
     REAL(8)               :: disp
     REAL(8)               :: fluid_mass


     ! Atomist model variables
     CHARACTER(LEN=100)    :: cnt_file                         ! CNT xyz input file
     INTEGER(4)            :: ncnt                             ! SWCNT number of atoms


     ! Write to output frequency variables
     INTEGER(4)            :: ntwx                             ! ntwx: rate of output of trajectory to .xzy file
     INTEGER(4)            :: ntpr                             ! ntpr: rate of output of properties to output file

     ! Initial configuration variables
     CHARACTER(LEN=50)     :: initial_config_type              ! can be either random or file
     CHARACTER(LEN=100)    :: initial_xyz_file                 ! initial .xyz file





     ! Coordinate arrays
     REAL(8), DIMENSION(:,:), ALLOCATABLE :: coord
     REAL(8), DIMENSION(:,:), ALLOCATABLE :: coord_eff
     REAL(8), DIMENSION(:,:), ALLOCATABLE :: coord_init
     REAL(8), DIMENSION(:,:), ALLOCATABLE :: coord_cnt


     ! Variables of a given MC step
     REAL(8)               :: deltaE                           ! deltaE = E_new-E_old
     REAL(8)               :: dr(3)

     CONTAINS
       PROCEDURE :: generate_random_distribution =>  generate_random_distribution
       PROCEDURE :: generate_reduced_units => generate_reduced_units
  END TYPE Simulation


  CONTAINS
    SUBROUTINE generate_random_distribution(this)
      IMPLICIT NONE
      !------------------------------------------------------------------------
      ! This subroutines generates an initial random distribution of particles
      ! inside the CNT.
      ! It starts by generating the random x coordinates and then generates the
      ! y coordinate with the constraint the the x**2+y**2 <= r**2.
      !------------------------------------------------------------------------
      CLASS(Simulation), INTENT(INOUT) :: this
      REAL(8) :: r2, x2 
      INTEGER(4) :: k


      this%coord = 0
      r2 = this%radius * this%radius                              ! Radius squared

      DO k=1,this%npart
         ! Generate x coordinate
         this%coord(1,k) = RAND()*2*this%radius-this%radius
         x2 = this%coord(1,k) * this%coord(1,k)
         ! Generate y coordinate
         this%coord(2,k) = RAND()*2*SQRT(r2-x2) - SQRT(r2 - x2)
         this%coord(3,k) = RAND()*this%length
      END DO
    END SUBROUTINE generate_random_distribution

    SUBROUTINE generate_reduced_units(this)
      IMPLICIT NONE
      !------------------------------------------------------------------------
      ! This subroutines makes the program to work with reduced unit by 
      ! transforming variables like the temperature, SWCNT length, radius
      ! and LJ 12-6 parameters into reduced variables.
      ! TODO: confirm that everything is working properly. I think it is not
      ! TODO: JoaoMorado. It was not 17.12.2018
      !------------------------------------------------------------------------
      CLASS(Simulation), INTENT(INOUT) :: this

      WRITE(6,*) "Old variables:"
      WRITE(6,*) "T=", this%temperature
      WRITE(6,*) "L=", this%length
      WRITE(6,*) "r=", this%radius
      WRITE(6,*) "sigma_cnt=", this%sigma_cnt
      WRITE(6,*) "eps_cnt=", this%eps_cnt


      ! Technically we do not multiply by the Boltzmann
      ! because it cancels in the acceptance probability

      this%temperature = this%temperature  / this%eps_fluid        
      this%length      = this%length / this%sigma_fluid
      this%radius      = this%radius / this%sigma_fluid

      ! Scale CNT interaction parameters
      this%sigma_cnt   = this%sigma_cnt / this%sigma_fluid
      this%eps_cnt     = this%eps_cnt   / this%eps_fluid



      WRITE(6,*) "New reduced variables:"
      WRITE(6,*) "T*=", this%temperature
      WRITE(6,*) "L*=", this%length
      WRITE(6,*) "r*=", this%radius
      WRITE(6,*) "sigma_cnt*=", this%sigma_cnt
      WRITE(6,*) "eps_cnt*=", this%eps_cnt
    END SUBROUTINE generate_reduced_units
END MODULE simulation_module

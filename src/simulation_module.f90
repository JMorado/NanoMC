MODULE simulation_module
    USE cell_module
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
        CHARACTER(LEN=100)    :: displacement_type                ! gaussian or uniform


        ! Cell parameters
        TYPE(Cell)            :: cell_instance                    ! Cell instance

        ! CNT parameters
        REAL(8)               :: length                           ! SWCNT length
        REAL(8)               :: radius                           ! SWCNT radius

        ! Lennard-Jones 12-6 parameters 
        REAL(8)               :: fluid_cut_off                    ! h2-h2 rcut_off 
        REAL(8)               :: cnt_cut_off                      ! h2-swcnt rcut_off
        REAL(8)               :: sigma_fluid                      ! sigma h2-h2 
        REAL(8)               :: eps_fluid                        ! epsilon h2-h2
        REAL(8)               :: sigma_cnt                        ! sigma h2-swcnt
        REAL(8)               :: eps_cnt                          ! epsilon h2-swcnt

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
    END TYPE Simulation


CONTAINS
    SUBROUTINE generate_random_distribution(self)
        !========================================================================!
        ! This subroutines generates an initial random distribution of particles !
        ! inside the CNT                                                         !
        ! It starts by generating the random x coordinates and then generates the!
        ! y coordinate with the constraint the the x**2+y**2 <= r**2.            !
        !------------------------------------------------------------------------!
        ! self (in) : Simulation instance                                        !
        !========================================================================!
        IMPLICIT NONE
        CLASS(Simulation), INTENT(INOUT) :: self
        REAL(8) :: r2, x2
        INTEGER(4) :: k

        self%coord = 0
        r2 = self%radius * self%radius

        DO k=1,self%npart
            ! Generate x coordinate
            self%coord(1,k) = RAND()*2*self%radius-self%radius
            x2 = self%coord(1,k) * self%coord(1,k)
            ! Generate y coordinate
            self%coord(2,k) = RAND()*2*SQRT(r2-x2) - SQRT(r2 - x2)
            self%coord(3,k) = RAND()*self%length
        END DO
    END SUBROUTINE generate_random_distribution
END MODULE simulation_module

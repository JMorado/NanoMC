MODULE time_module
  IMPLICIT NONE

  TYPE, PUBLIC :: Time
     ! Time parameters
     REAL(8) :: t                    ! Time
 
     ! Velocity parameters
     REAL(8) :: v_rms                ! Root mean square speed
     REAL(8) :: v_mean               ! Mean speed
     REAL(8) :: v_rrels              ! Root mean square relative velocity 

     ! Displacements parameterr
     !REAL(8) :: msd                  ! Mean square deviation

   CONTAINS
     PROCEDURE     :: initialization => initialization
     PROCEDURE     :: update_time => update_time
     PROCEDURE     :: self_diffusion => self_diffusion
     PROCEDURE     :: compute_SD => compute_SD
  END TYPE Time


CONTAINS
  SUBROUTINE initialization(this, mass, temperature)
    IMPLICIT NONE
    !-----------------------------------------------------------------
    ! This subroutines initializes a time object by setting initial
    ! time equal to 0 and by computing the root mean square velocity
    !
    ! The average velocity of gas particles is found using the root
    ! mean square velocity formula:
    !
    ! v_rms = SQRT(3*kB*T/mass)
    !
    !
    ! The expected value of the Maxwell-Boltzman distribution is the
    ! mean speed:
    !
    ! v_mean = SQRT(8*kB*T/(pi*mass))
    !
    !
    ! The root mean square relative velocity corresponds to the
    ! average relative speed of the molecules:
    !
    ! v_rrels =  SQRT(6*kB*T/mass)
    !-----------------------------------------------------------------
    ! SOURCE:
    ! Maxwell-Boltzmann distribution
    ! https://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution
    !
    ! Kinetic Theory, David Tong
    ! http://www.damtp.cam.ac.uk/user/tong/kintheory/kt.pdf
    !-----------------------------------------------------------------
    Class(Time), INTENT(INOUT) :: this
    REAL(8), INTENT(IN)        :: mass, temperature

    ! Initiate time
    this%t = 0

    ! Compute different velocities
    this%v_rrels = SQRT((6*1.381d-23*temperature)/mass) * 1d10
    this%v_rms   = SQRT((3*1.381d-23*temperature)/mass) * 1d10
    this%v_mean  = SQRT((8*1.381d-23*temperature)/(3.14159265359*mass)) * 1d10

  END SUBROUTINE initialization

  SUBROUTINE update_time(this, dt)
    IMPLICIT NONE
    !-----------------------------------------------------------------
    ! This subroutines updates the simulation time given dt.
    !-----------------------------------------------------------------
    Class(Time), INTENT(INOUT)   :: this
    REAL(8),  INTENT(IN)         :: dt

    this%t = this%t + dt
  END SUBROUTINE update_time

  FUNCTION compute_SD(this, initial_coord, effective_coord, npart, dimen) RESULT(SD)
    IMPLICIT NONE
    !-----------------------------------------------------------------
    ! This subroutine computes the squared displacement.
    !-----------------------------------------------------------------
    Class(Time), INTENT(IN)    :: this
    INTEGER(4) , INTENT(IN)    :: npart, dimen
    REAL(8)    , INTENT(IN)    :: initial_coord(:,:), effective_coord(:,:)

    REAL(8)                    :: SD
    INTEGER(4)                 :: k


    ! Compute MSD
    SD = 0
    DO k=1,npart
       SD = SD + ( effective_coord(3,k)-initial_coord(3,k) )**2
       !SD = SD + SUM(( effective_coord(:,k)-initial_coord(:,k) )**2) ! 3D mean square
    END DO

  END FUNCTION compute_SD

  FUNCTION self_diffusion(this, npart, dimen, square_displacement, total_t) RESULT(Ds)
    IMPLICIT NONE
    !-----------------------------------------------------------------
    ! This function calculates the diffusion coefficients along the
    !  using the Einstein formula:
    !
    ! Deff = (Total_DeltaR)^2 / (N*2*dim*time)
    !-----------------------------------------------------------------
    ! SOURCE:
    !
    !
    !-----------------------------------------------------------------
    Class(Time), INTENT(INOUT) :: this
    INTEGER(4) , INTENT(IN)    :: npart, dimen
    REAL(8), INTENT(IN)        :: square_displacement
    REAL(8), OPTIONAL          :: total_t

    REAL(8)                    :: Ds

    IF ( .NOT. PRESENT(total_t)) THEN
       total_t =  this%t
    END IF


    ! Determine the self-diffusion coefficient
    Ds = square_displacement / (2.0d0 * dimen * total_t * npart)
    Ds = Ds * 1.0d-20

    !WRITE(6,*) t, square_displacement, Ds
  END FUNCTION self_diffusion

END MODULE time_module

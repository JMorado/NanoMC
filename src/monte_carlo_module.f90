MODULE monte_carlo_module
  USE simulation_module
  USE energy_module
  USE pbc_displacement_module
  

  IMPLICIT NONE
  ABSTRACT INTERFACE
     FUNCTION energy_function(this, index, coord)
       IMPORT                            :: Simulation
       REAL(8)                           :: energy_function
       Class(Simulation), INTENT(IN)     :: this
       REAL(8),           INTENT(IN)     :: coord(:)
       INTEGER(4),        INTENT(IN)     :: index
     END FUNCTION energy_function
  END INTERFACE
 
 
CONTAINS


  SUBROUTINE MC_MOVE(this, energy_func, accepted, o)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! This subroutine attempts to perform a stsandard MC move
    ! by displacing a particle.
    !
    ! TODO: 
    !   ! Select a particle to be displaced
    !o=INT(this%npart*RAND())+1
    !
    !
    ! SOURCE:
    ! "Understanding Molecular Simulation:
    ! From Algorithms To Applications"
    ! Daan Frenkel, Berend Smit
    ! pp. 32-33
    !----------------------------------------------------------------------! 
    CLASS(Simulation), INTENT(INOUT)    :: this
    LOGICAL(1), INTENT(OUT)             :: accepted
    INTEGER(4), INTENT(IN)              :: o
    REAL(8)                             :: deltaE, e_old, e_new, prob_acceptance
    REAL(8)                             :: coord_tmp(3)
    REAL(8)                             :: dr(3)
    PROCEDURE(energy_function), POINTER :: energy_func



    accepted = .FALSE.
    ! Energy old configuration
    e_old = energy_func(this, o, this%coord(:,o))

 
    dr =  displacement( this%disp )
    coord_tmp = this%coord(:,o) + dr(:)
    DO WHILE (SQRT(SUM(coord_tmp(1:2)**2)) .ge. this%radius)
       dr =  displacement( this%disp )
       coord_tmp = this%coord(:,o) + dr
    END DO
    CALL pbc_z(coord_tmp, this%length)

    ! Energy new configuration
    e_new = energy_func(this, o, coord_tmp(:))

    ! Energy difference
    deltaE = e_new-e_old
    ! Metropolis criteria
    IF (deltaE .le. 0) THEN
       this%dr = dr
       this%deltaE = deltaE

       this%coord(:,o) = coord_tmp
       this%coord_eff(:,o) =  this%coord_eff(:,o) + dr

       accepted = .TRUE.
    ELSE IF ((deltaE / this%temperature) .le. 13) THEN

       prob_acceptance = min(exp(- deltaE / this%temperature),1.0)

       IF (prob_acceptance .ge. RAND()) THEN
          this%dr = dr
          this%deltaE = deltaE

          this%coord(:,o) = coord_tmp
          this%coord_eff(:,o) =  this%coord_eff(:,o) + dr

          accepted = .TRUE.
       END IF
    END IF
  END SUBROUTINE MC_MOVE

  SUBROUTINE MC_EXCHANGE(this, energy_func)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! This subroutine attempts to perform a grand canonical MC move
    ! by exchanging a particle with a reservoir.
    ! 
    ! SOURCE:
    ! "Understanding Molecular Simulation:
    ! From Algorithms To Applications"
    ! Daan Frenkel, Berend Smit
    ! pp. 32-33
    !----------------------------------------------------------------------!
    CLASS(Simulation), INTENT(INOUT) :: this
    REAL(8)                          :: energy, prob_acceptance, arg
    REAL(8)                          :: r2, x2
    REAL(8)                          :: coord_tmp(3)
    INTEGER(4)                       :: o
    PROCEDURE(energy_function), POINTER :: energy_func


    IF (RAND() .lt. 0.5) THEN
       IF (this%npart .eq. 0) THEN
          ! There are no particles
          return
       END IF

       ! TODO: Joao Morado 17.12.2018  
       ! TODO: Test for overlap of particle upon insertion a displacement
       o = INT(this%npart*RAND())+1

       ! Energy of particle to be removed
       energy = energy_func(this, o, this%coord(:,o))

       ! OLD ACCEPTANCE CRITERIA (Frenkel and Smith)
       !arg =  exp(energy/this%temperature) * (this%npart * 1d0) / (this%pressure * this%beta * this%volume * 1d-30)

       ! New acceptance criteria (Tildesley)
       arg = energy/this%temperature
       arg = arg + LOG( REAL(this%npart) / (this%activity*this%volume*1d-30)) 

       IF (arg .gt. 0) THEN
          this%coord(:,o) = this%coord(:,this%npart)
          this%npart = this%npart - 1
       ELSE IF (arg .gt. -30.0) THEN
          prob_acceptance = min(1.0d0,EXP(arg))
          IF ( RAND() .lt. prob_acceptance ) THEN
          this%coord(:,o) = this%coord(:,this%npart)
          this%npart = this%npart - 1
          END IF
       END IF
    ELSE
       ! New particle at random poistion
       r2 = this%radius * this%radius
       coord_tmp(1) = RAND()*2*(this%radius-this%sigma_cnt)-(this%radius-this%sigma_cnt)
       x2 = coord_tmp(1) * coord_tmp(1)
       coord_tmp(2) = RAND()*2*SQRT(r2-x2) - SQRT(r2 - x2)
       coord_tmp(3) = RAND()*this%length

       ! Energy of particle to be inserted
       energy = energy_func(this, this%npart+1 , coord_tmp)

       ! New acceptance criteria (Tildesley)
       arg = - energy/this%temperature
       arg = arg + LOG(this%activity*this%volume*1d-30/REAL(this%npart+1))

       IF (arg .gt. 0) THEN
          this%coord(:,this%npart+1) = coord_tmp
          this%npart = this%npart + 1
       ELSE IF (arg .gt. -30.0) THEN
          prob_acceptance = min(1.0d0,EXP(arg))
          IF ( RAND() .lt. prob_acceptance ) THEN
          this%coord(:,this%npart+1) = coord_tmp
          this%npart = this%npart + 1
          END IF
       END IF
    END IF

  END SUBROUTINE mc_exchange
END MODULE monte_carlo_module

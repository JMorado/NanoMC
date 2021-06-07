MODULE brownian_dynamics_module
  USE simulation_module
  USE pbc_displacement_module
CONTAINS
  SUBROUTINE compute_forces(this, forces, Epot)
    IMPLICIT NONE
    !----------------------------------------------------------------------------
    ! Computes the forces 
    !----------------------------------------------------------------------------  
    CLASS(Simulation), INTENT(IN)  :: this
    REAL(8),           INTENT(OUT) :: forces(:,:)
    REAL(8),           INTENT(OUT) :: Epot
    INTEGER(4)                     :: k, l
    REAL(8)                        :: dx, dy, dz, r2, r(3)
    REAL(8)                        :: sigma2fluid, sigma2cnt
    REAL(8)                        :: fr2, fr6, fxi, fyi, fzi, fpr


    forces = 0.0d0
    Epot = 0.0d0
    sigma2fluid = this%sigma_fluid * this%sigma_fluid
    sigma2cnt = this%sigma_cnt * this%sigma_cnt

    DO k=1,this%npart
       DO l=k+1,this%npart
          dx = this%coord(1,k)-this%coord(1,l)
          dy = this%coord(2,k)-this%coord(2,l)
          dz = this%coord(3,k)-this%coord(3,l)

          ! TODO:make this better
          r = (/dx, dy, dz/)
          CALL pbc_z_distance(r, this%length)
          dz = r(3)


          r2 = dx*dx + dy*dy + dz*dz


          IF (r2 .lt. this%fluid_cut_off * this%fluid_cut_off) THEN
             fr2 = sigma2fluid / r2
             fr6 = fr2 * fr2 * fr2
             fpr = 48.d0 * this%eps_fluid * fr6 * ( fr6 - 0.5d0 ) / r2

             fxi = fpr * dx
             fyi = fpr * dy
             fzi = fpr * dz

             forces(1,k) = forces(1,k) + fxi
             forces(2,k) = forces(2,k) + fyi
             forces(3,k) = forces(3,k) + fzi

             forces(1,l) = forces(1,l) - fxi
             forces(2,l) = forces(2,l) - fyi
             forces(3,l) = forces(3,l) - fzi

             Epot = Epot + 4.d0 * this%eps_fluid * fr6 * (fr6 - 1.d0)
          END IF

       END DO
       DO l=1,this%ncnt
          dx = this%coord(k,1)-this%coord_cnt(l,1)
          dy = this%coord(k,2)-this%coord_cnt(l,2)
          dz = this%coord(k,3)-this%coord_cnt(l,3)

          ! TODO:make this better
          r = (/dx, dy, dz/)
          CALL pbc_z_distance(r, this%length)
          dz = r(3)
          r2 = dx*dx + dy*dy + dz*dz

          IF (r2 .lt. this%cnt_cut_off * this%cnt_cut_off) THEN
             fr2 = sigma2cnt / r2
             fr6 = fr2 * fr2 * fr2
             fpr = 48.d0 * this%eps_cnt * fr6 * ( fr6 - 0.5d0 ) / r2

             fxi = fpr * dx
             fyi = fpr * dy
             fzi = fpr * dz

             forces(1,k) = forces(1,k) + fxi
             forces(2,k) = forces(2,k) + fyi
             forces(3,k) = forces(3,k) + fzi

             Epot = Epot + 4.d0 * this%eps_cnt * fr6 * (fr6 - 1.d0)
          END IF
       END DO
    END DO
  END SUBROUTINE compute_forces



  SUBROUTINE brownian_dynamics_integrator(this, dt)
    IMPLICIT NONE
    !----------------------------------------------------------------------------
    ! Propagates equation of motion
    !
    !
    ! Source
    ! Allen, M. P.; Tildesley, D. J.
    ! Computer Simulation of Liquids, Second edition.
    ! Oxford University Press: Oxford, United Kingdom, 2017.
    ! pp. 383-387
    !----------------------------------------------------------------------------
    CLASS(Simulation), INTENT(INOUT)         :: this
    INTEGER(4)                               :: k,j
    REAL(8), INTENT(IN)                   :: dt
    REAL(8), ALLOCATABLE                  :: forces(:,:),  G(:,:)
    REAL(8)                               :: Epot, D0, dr(3)
    REAL(8)                               :: mean, std_dev, fric, coord_tmp(3)


   


    D0 = 4.759896405273463d-08     ! m^2 s^-1
    D0 = 4759896405273.463         ! A^2 s^-1 

    ! D0 =

    fric = 1.0d0

    mean = 0.0d0
    std_dev = SQRT(2*1.381d-23*this%temperature*dt/fric)


    ALLOCATE(forces(3,this%npart))
    ALLOCATE(G(3,this%npart))

    ! Compute forces

    CALL compute_forces(this,forces,Epot)
    DO k=1,this%npart
       WRITE(6,*) (forces(j,k),j=1,3)
     end do

    CALL random_displacement(G, this%npart, mean, std_dev)

    DO k=1,this%npart
       !this%coord(:,k) = this%coord(:,k) + forces(:,k) * D0 * dt / this%temperature +  G(:,k) * SQRT(2.0d0*D0*dt)

       dr = forces(:,k) * dt / fric +  G(:,k)

       CALL reflexion(this,dr,this%coord(:,k))
       this%coord(3,k) = this%coord(3,k) + dr(3)


              ! fric = kb*T / D


       CALL pbc_z(this%coord(:,k), this%length)
    END DO

    WRITE(6,*) Epot


    !WRITE(6,*) this%coord


    DEALLOCATE(forces, G)

  END SUBROUTINE BROWNIAN_DYNAMICS_INTEGRATOR


  SUBROUTINE random_displacement(G, nat, mean, std_dev)
    IMPLICIT NONE
    REAL(8), INTENT(OUT)      :: G(:,:)
    INTEGER(4), INTENT(IN)    :: nat
    REAL(8), INTENT(OUT)      ::  mean, std_dev
    INTEGER(4)                :: j, k
    
    G = 0
    DO j=1,nat
       DO k=1,3
          G(k,j) = rand_normal(mean,std_dev)
       END DO
    END DO

    return
  END SUBROUTINE random_displacement



  FUNCTION rand_normal(mean,stdev) RESULT(c)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: mean, stdev
    REAL(8)             :: c,temp(2),theta,r
    IF(stdev <= 0.0d0) THEN
       WRITE(*,*) "Standard Deviation must be +ve"
    ELSE
       CALL RANDOM_NUMBER(temp)
       r=(-2.0d0*log(temp(1)))**0.5
       theta = 2.0d0*3.1415*temp(2)
       c= mean+stdev*r*sin(theta)
    END IF
  END FUNCTION rand_normal


END MODULE brownian_dynamics_module

!----------------------------------------------------------------------! 
! Modules create at some point but not necessary anymore.
! Joao Morado
!----------------------------------------------------------------------! 
!!! This module contains the following loose modules:
!!! TODO: move type
!!! TODO: virial
!!! TODO: virial_local
!!! TODO: impulsive_pressure
!!! TODO: reflexion
!!! TODO: fpath_z / fpath_xy
!!! TODO: compute_vel_tst
!!! TODO: PES2
!!! TODO: rand_nroma
!!! TODO: max_boltz_vel
!!! TODO: PES
!!! TODO: mc_exchange_cavity_biased
!!! TODO: is_in_cavity

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

FUNCTION virial_local(simulation_instance, k) RESULT(virial)
    !-------------------------------------------------------------
    ! Computes the virial
    !-----------------------------------------------------------
    IMPLICIT NONE
    ! Input class
    TYPE(Simulation), INTENT(IN) :: simulation_instance
    ! Intermediate variables
    INTEGER(4), INTENT(IN) :: k
    INTEGER(4) :: l
    REAL(8) :: r(3), r2
    REAL(8) :: attract_term, rep_term
    REAL(8) :: sigmaf2, sigmacnt2
    ! Output variable
    REAL(8) :: virial


    sigmaf2 = simulation_instance%sigma_fluid*simulation_instance%sigma_fluid          ! Sigma fluid squared
    sigmacnt2 = simulation_instance%sigma_cnt*simulation_instance%sigma_cnt            ! Sigma CNT squared

    virial = 0

    DO l=k+1,simulation_instance%npart
        ! Compute distance r
        r = simulation_instance%coord(:,k)-simulation_instance%coord(:,l)
        ! Apply PBC to r
        CALL pbc_z_distance(r, simulation_instance%length)

        ! Compute r squared of the periodic distance
        r2 = SUM(r*r)


        IF (r2 .lt. simulation_instance%fluid_cut_off * simulation_instance%fluid_cut_off) THEN
            attract_term = (sigmaf2 / r2) * (sigmaf2 / r2) * (sigmaf2 / r2)
            rep_term = attract_term * attract_term
            virial = virial + 48 * simulation_instance%eps_fluid * (rep_term - 0.5 * attract_term) !/ SQRT(r2)
        END IF
    END DO


    DO l=1,simulation_instance%ncnt
        ! Compute distance r
        r = simulation_instance%coord(:,k)-simulation_instance%coord_cnt(:,l)
        ! Apply PBC to r
        CALL pbc_z_distance(r, simulation_instance%length)

        ! Compute r squared of the periodic distance
        r2 = SUM(r*r)

        ! If within cut off
        IF (r2 .lt. simulation_instance%cnt_cut_off * simulation_instance%cnt_cut_off) THEN
            attract_term = (sigmacnt2 / r2) * (sigmacnt2 / r2) * (sigmacnt2 / r2)
            rep_term = attract_term * attract_term
            virial = virial + 48 * simulation_instance%eps_cnt * (rep_term - 0.5 * attract_term) ! / SQRT(r2)!
        END IF
    END DO

    virial = (1.0d0/6.0d0) * virial

    RETURN

END FUNCTION virial_local


  SUBROUTINE move_type(this, o, dr)
    IMPLICIT NONE
    CLASS(Simulation), INTENT(INOUT)    :: this
    INTEGER(4), INTENT(IN)              :: o
    REAL(8), INTENT(INOUT)              :: dr(3)
    INTEGER(4)                          :: j
    REAL(8)                             :: nbulk, nsurf, area, vol, dx, s


    nbulk = 0
    nsurf = 0
    

    s = 1.2
    ! Area and volume
    area = 2*3.14159265359*(this%radius-this%sigma_cnt)*this%length
    vol = 3.14159265359 * (this%radius-this%sigma_cnt)**2 * this%length

    DO j=1,this%npart
       IF (SQRT(SUM(this%coord(1:2,j)**2)) .lt. this%radius-this%sigma_cnt-s) THEN
          nbulk = nbulk + 1d0
       ELSE
          nsurf = nsurf + 1d0
       END IF
    END DO


    

    IF (SQRT(SUM(this%coord(1:2,o)**2)) .lt. this%radius-this%sigma_cnt-s) THEN
       dx = sqrt(2.0d0) * 3.14159265359 * (nbulk / vol) * this%sigma_fluid**2
       dx = 1.0d0 / dx
       dr = displacement( dx )

    ELSE
       dx = sqrt(2.0d0) * 3.14159265359 * (nsurf / area) * this%sigma_fluid
       dx = 1.0d0 / dx
       dr = displacement( dx)
    END IF

    !WRITE(6,*) nbulk,nsurf,dx,SQRT(SUM(this%coord(1:2,o)**2))
    return
  END SUBROUTINE MOVE_TYPE


  FUNCTION virial(this)
    !-------------------------------------------------------------
    ! Computes the virial
    !----------------------------------------------------------- 
    IMPLICIT NONE
    ! Input class
    CLASS(Simulation), INTENT(IN) :: this
    ! Intermediate variables
    INTEGER(4) :: k, l
    REAL(8) :: r(3), r2
    REAL(8) :: attract_term, rep_term
    REAL(8) :: sigmaf2, sigmacnt2
    ! Output variable
    REAL(8) :: virial


    sigmaf2 = this%sigma_fluid*this%sigma_fluid          ! Sigma fluid squared
    sigmacnt2 = this%sigma_cnt*this%sigma_cnt            ! Sigma CNT squared


    virial = 0
    DO k=1,this%npart
       DO l=k+1,this%npart
          ! Compute distance r
          r = this%coord(:,k)-this%coord(:,l)
          ! Apply PBC to r
          CALL pbc_z_distance(r, this%length)

          ! Compute r squared of the periodic distance
          r2 = SUM(r*r)

          IF (r2 .lt. this%fluid_cut_off * this%fluid_cut_off) THEN
             attract_term = (sigmaf2 / r2) * (sigmaf2 / r2) * (sigmaf2 / r2)
             rep_term = attract_term * attract_term
             virial = virial + 48 * this%eps_fluid * (rep_term - 0.5 * attract_term) 
          END IF
       END DO

       DO l=1,this%ncnt
          ! Compute distance r
          r = this%coord(:,k)-this%coord_cnt(:,l)
          ! Apply PBC to r
          CALL pbc_z_distance(r, this%length)
      
          ! Compute r squared of the periodic distance
          r2 = SUM(r*r)
        
          ! If within cut off
          IF (r2 .lt. this%cnt_cut_off * this%cnt_cut_off) THEN
             attract_term = (sigmacnt2 / r2) * (sigmacnt2 / r2) * (sigmacnt2 / r2)
             rep_term = attract_term * attract_term
             virial = virial + 48 * this%eps_cnt * (rep_term - 0.5 * attract_term)
          END IF
      END DO

       ! Fluid-CNT contribution to the eergy
       ! r2 = this%radius - SQRT(SUM( this%coord(1:2,k)*this%coord(1:2,k) ))
       ! r2 = r2*r2
       ! attract_term = (sigmacnt2 / r2) * (sigmacnt2 / r2) * (sigmacnt2 / r2)
       ! rep_term = attract_term * attract_term
       ! virial = virial + 48 * this%eps_cnt * (rep_term - 0.5 * attract_term) / SQRT(r2)
    END DO

    virial = (1.0d0/3.0d0) * virial 
    RETURN
  END FUNCTION virial

    FUNCTION virial_local(this, k) RESULT(virial)
    !-------------------------------------------------------------
    ! Computes the virial
    !----------------------------------------------------------- 
    IMPLICIT NONE
    ! Input class
    CLASS(Simulation), INTENT(IN) :: this
    ! Intermediate variables
    INTEGER(4), INTENT(IN) :: k
    INTEGER(4) :: l
    REAL(8) :: r(3), r2
    REAL(8) :: attract_term, rep_term
    REAL(8) :: sigmaf2, sigmacnt2
    ! Output variable
    REAL(8) :: virial


    sigmaf2 = this%sigma_fluid*this%sigma_fluid          ! Sigma fluid squared
    sigmacnt2 = this%sigma_cnt*this%sigma_cnt            ! Sigma CNT squared

    virial = 0
    
       DO l=k+1,this%npart
          ! Compute distance r
          r = this%coord(:,k)-this%coord(:,l)
          ! Apply PBC to r
          CALL pbc_z_distance(r, this%length)

          ! Compute r squared of the periodic distance
          r2 = SUM(r*r)


          IF (r2 .lt. this%fluid_cut_off * this%fluid_cut_off) THEN
             attract_term = (sigmaf2 / r2) * (sigmaf2 / r2) * (sigmaf2 / r2)
             rep_term = attract_term * attract_term
             virial = virial + 48 * this%eps_fluid * (rep_term - 0.5 * attract_term) !/ SQRT(r2)
          END IF
       END DO


       DO l=1,this%ncnt
        ! Compute distance r
          r = this%coord(:,k)-this%coord_cnt(:,l)
          ! Apply PBC to r
          CALL pbc_z_distance(r, this%length)
      
          ! Compute r squared of the periodic distance
          r2 = SUM(r*r)
        
          ! If within cut off
          IF (r2 .lt. this%cnt_cut_off * this%cnt_cut_off) THEN
             attract_term = (sigmacnt2 / r2) * (sigmacnt2 / r2) * (sigmacnt2 / r2)
             rep_term = attract_term * attract_term
             virial = virial + 48 * this%eps_cnt * (rep_term - 0.5 * attract_term) ! / SQRT(r2)!
          END IF
    END DO
    
    virial = (1.0d0/6.0d0) * virial
    
    RETURN

  END FUNCTION virial_local




  FUNCTION impulsive_pressure(this,sigma,eps, cut_off)
    IMPLICIT NONE
    CLASS(Simulation), INTENT(IN) :: this
    REAL(8), INTENT(IN)           :: sigma, eps, cut_off
    REAL(8)                       :: impulsive_pressure


    impulsive_pressure = 8 * 3.14159265359 * this%density ** 2 * eps * sigma ** 3 / 3.0d0
    impulsive_pressure = impulsive_pressure * ( (sigma/cut_off)**9 - (sigma/cut_off)**3 )
  END FUNCTION IMPULSIVE_PRESSURE


  FUNCTION fpath_z(this,o) RESULT(lambda)
    IMPLICIT NONE
    CLASS(Simulation), INTENT(INOUT)    :: this
    INTEGER(4), INTENT(IN)              :: o
    REAL(8)                             :: lambda, r(3), volume, density
    INTEGER(4)                          :: k, ncount

    ncount = 0
    DO k=1,this%npart
       IF (k .eq. o) THEN
          CYCLE
       END IF
       r = this%coord(:,o)-this%coord(:,k)
       CALL pbc_z_distance(r, this%length)

       IF (SQRT(SUM(r(1:2)**2)) .le. this%sigma_fluid) THEN
          ncount = ncount + 1
       END IF
    END DO


    IF (ncount .eq. 0) then
       lambda = this%length * (-LOG(RAND()))
       return
    END IF
    volume = 3.1415 * (this%sigma_fluid/2.0d0) ** 2 * this%length
    density = ncount / volume

    lambda = 1.0d0
    lambda = lambda / (sqrt(2.0d0) * 3.1415 * this%sigma_fluid ** 2 * density)
    lambda = lambda * (-LOG(RAND()))


    return
  END FUNCTION FPATH_z


 SUBROUTINE reflexion(this,dr,coord)
    IMPLICIT NONE
    CLASS(Simulation), INTENT(INOUT)    :: this

    REAL(8), INTENT(IN) :: dr(3)
    REAL(8), INTENT(INOUT) :: coord(3)

    REAL(8) :: a, b, c, d, hypoth,hypothd, dx, dy, initialx, initialy
    real(8) :: drestx, dresty, moduloi, alpha,angulo_i, zeta,beta,zeta_i, deslocx,deslocy
    

    !RITE(6,*) "********************************", coord(1),coord(2),coord(3)

    deslocx = dr(1)
    deslocy = dr(2)
    hypoth = SQRT( (coord(1)+deslocx)**2 + (coord(2)+deslocy)**2)

    !RITE(6,*) "COORD",coord, hypoth
    !RITE(6,*) deslocx,deslocy
    DO WHILE (this%radius < hypoth)
       initialx = coord(1)
       initialy = coord(2)

       hypothd = SQRT( deslocx**2 + deslocy**2)
       dx = deslocx / hypothd
       dy = deslocy / hypothd

       a = 1.0d0
       b = 2.0d0 * ( (coord(1)*deslocx + coord(2)*deslocy ) / hypothd )
       c = (coord(1)*coord(1) + coord(2)*coord(2)) - (this%radius*this%radius)

       d = (-b + sqrt(b*b - 4*a*c)) / (2*a)

       
       ! Move particle to the circunferenc
       !IF (c .gt. 0) then
       !   WRITE(6,*) "CCCC",c
       !   WRITE(6,*) coord, this%radius,coord(1)*coord(1),coord(2)*coord(2),this%radius**2
       !   STOP
       !END IF
       coord(1) = coord(1) +  (deslocx/hypothd)*d
       coord(2) = coord(2) +  (deslocy/hypothd)*d


       ! Compute displacement left
       drestx = (deslocx-dx*d)
       dresty = (deslocy-dy*d)

       ! Compute incidence angle which will be equal to reflexion angle
       moduloi = sQRT(coord(1)**2 + coord(2)**2)


       alpha = ACOS( (deslocx*coord(1)+deslocy*coord(2)) / (moduloi * hypothd))
       angulo_i = 3.14159265359 / 2.0d0 - alpha


       
       zeta = ATAN(coord(2) / coord(1))

     
       zeta_i = ATAN(initialx / initialy)
       !WRITE(6,*)"ZETA", zeta,zeta_i, coord(2),coord(1)


       IF ((zeta_i .gt. zeta) .and. (zeta_i .lt. (zeta +3.14159265359 ))) THEN
          beta = - 2.0d0 * angulo_i
       ELSE IF ((zeta_i .lt. zeta) .and. (zeta_i .gt. (zeta +3.14159265359 ))) THEN
          beta = 2.0d0 * angulo_i
       ELSE
          beta = 3.14159265359
       END IF

       deslocx=(drestx*cos(beta)-dresty*sin(beta))
       deslocy=(drestx*sin(beta)+dresty*cos(beta))
       hypoth = SQRT( (coord(1)+deslocx)**2 + (coord(2)+deslocy)**2)

    END DO

    coord(1) = coord(1) + deslocx
    coord(2) = coord(2) + deslocy

    return

  END SUBROUTINE reflexion
  FUNCTION fpath_xy(this,o) RESULT(lambda)
    IMPLICIT NONE
    CLASS(Simulation), INTENT(INOUT)    :: this
    INTEGER(4), INTENT(IN)              :: o
    REAL(8)                             :: lambda, r(3), volume, density
    INTEGER(4)                          :: k, ncount

    ncount = 0
    DO k=1,this%npart
       IF (k .eq. o) THEN
          CYCLE
       END IF
       r = this%coord(:,o)-this%coord(:,k)
       CALL pbc_z_distance(r, this%length)

       IF (abs(r(3)) .le. this%sigma_fluid) THEN
          ncount = ncount + 1
       END IF
    END DO


    IF (ncount .eq. 0) then
       lambda = 2 * this%radius !2*(-LOG(RAND()))
       return
    END IF
    volume = 3.1415 * this%radius ** 2 * this%sigma_fluid
    density = ncount / volume

    lambda = 1.0d0
    lambda = lambda / (sqrt(2.0d0) * 3.1415 * this%sigma_fluid ** 2 * density)
    lambda = lambda * (-LOG(RAND()))


    return
  END FUNCTION FPATH_XY



  FUNCTION compute_vel_tst(this, deltaE, temperature)
    IMPLICIT NONE
    !-----------------------------------------------------------------
    ! This subroutines computes the velocity of a given process
    ! accordingly to the equation (in the spirit of TST):
    !
    ! v = <v^2> exp(-deltaE/T)
    !
    ! It it assumed that deltaE is given in Kelvin units.
    !-----------------------------------------------------------------
    ! SOURCE:
    !
    !-----------------------------------------------------------------
    Class(Time), INTENT(IN)    :: this
    REAL(8), INTENT(IN)        :: deltaE, temperature
    REAL(8)                    :: compute_vel_tst


   ! if (deltaE .lt. 0) THEN
   !    compute_vel_tst  = 0
   ! else
       compute_vel_tst = 1e-13 * exp(-deltaE/temperature)  !this%v_rms !* EXP(-deltaE / temperature)
   ! end if

    !WRITE(6,*) compute_vel_tst, deltaE
  END FUNCTION compute_vel_tst

  SUBROUTINE PES2(this,energy_func,o)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! This subroutine attempts to perform a standard MC move
    ! by displacing a particle.
    !
    ! SOURCE:
    ! "Understanding Molecular Simulation:
    ! From Algorithms To Applications"
    ! Daan Frenkel, Berend Smit
    ! pp. 32-33
    !----------------------------------------------------------------------! 
    CLASS(Simulation), INTENT(INOUT)    :: this
    INTEGER(4), INTENT(IN)           :: o 
    REAL(8)                          :: coord_tmp(3), dir(3)
    REAL(8)                          :: energy_points(0:20)
    INTEGER(4)                       :: j,i
    PROCEDURE(energy_function), POINTER :: energy_func
    REAL(8) :: minim, maxim, minimloc,maximloc 

  energy_points(0) = energy_func(this, o, coord_tmp)
  minimloc = 0
  maximloc = 0
  minim = energy_points(0)
  maxim = energy_points(0)
  j = 20
  DO i=1,j 
     coord_tmp = this%coord(:,o) + this%dr*(i*1.0d0/j)
     energy_points(i) = energy_func(this, o, coord_tmp)

     IF (energy_points(i) .lt. minim) THEN
        ! New minimum was found
        minim = energy_points(i)
        minimloc = i
     END IF

     IF (energy_points(i) .gt. maxim) THEN
        ! New maximum was found
        maxim = energy_points(i)
        maximloc = i
     END IF
     !     write(6,*) i, energy_points(i), this%deltaE, (i*1d0/j)
  END DO

  IF (maximloc .gt. minimloc) THEN
     ! deltaE > 0
     this%deltaE = maxim-minim
  ELSE
     ! deltaE < 0
     this%deltaE = minim-maxim
  END IF
END SUBROUTINE PES2

FUNCTION rand_normal(this,mean,stdev) RESULT(c)
  IMPLICIT NONE
  Class(Time), INTENT(IN)    :: this

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

FUNCTION max_boltz_vel(this, temperature) RESULT(vel)
  IMPLICIT NONE
  Class(Time), INTENT(IN)    :: this
  REAL(8), INTENT(IN) :: temperature
  INTEGER(4) :: k
  REAL(8)    :: vel(3), mass

  mass = 2*1.6737236d-27   ! Kg
  vel = 0
  DO k=1,3
     ! Variance 1 Mean 0
     vel(k) = SQRT(1.381d-23*temperature / mass) * rand_normal(this,0d0, 1d0)
  END DO
  vel = vel * 1e10 ! convert to Angstrom/s
END FUNCTION max_boltz_vel


  SUBROUTINE PES(this,energy_func,o)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! This subroutine attempts to perform a standard MC move
    ! by displacing a particle.
    !
    ! SOURCE:
    ! "Understanding Molecular Simulation:
    ! From Algorithms To Applications"
    ! Daan Frenkel, Berend Smit
    ! pp. 32-33
    !----------------------------------------------------------------------! 
    CLASS(Simulation), INTENT(INOUT)    :: this
    INTEGER(4), INTENT(IN)           :: o 
    REAL(8)                          :: coord_tmp(3), dir(3)
    REAL(8)                          :: energy_points(0:100), en,r2
    INTEGER(4)                       :: j,i
    PROCEDURE(energy_function), POINTER :: energy_func
    REAL(8) :: minim, maxim, minimloc,maximloc 




    dir = this%dr / SQRT(DOT_PRODUCT(this%dr,this%dr))
    coord_tmp = this%coord(:,o)
    energy_points(0) = energy_func(this, o, coord_tmp)
    energy_points(1) = energy_func(this, o, coord_tmp + 0.01 * dir)

    WRITE(6,*) "Starting...."
    WRITE(6,*)0,energy_points(0)
    write(6,*)1,energy_points(1)
    if (energy_points(1) .gt. energy_points(0)) then
       ! Ascending the hill
       DO i=2,100
          energy_points(i) = energy_func(this, o, coord_tmp + i * 0.01 * dir)
          write(6,*)i,energy_points(i)

          IF (energy_points(i) .lt. energy_points(i-1)) THEN
             ! maxima found
             this%deltaE = energy_points(i-1)-energy_points(0)
             WRITE(6,*) "Maxima found on ascending hill"
             WRITE(6,*) "deltaE", this%deltaE

             return
          END IF
       END DO
       WRITE(6,*) "No maxima found on initial ascending hill."
    else
       ! Descending the hill
       DO i=2,100
          energy_points(i) = energy_func(this, o, coord_tmp + i * 0.01 * dir)

          write(6,*)i,energy_points(i)
          IF (energy_points(i) .gt. energy_points(i-1)) THEN
             !
             minim = energy_points(i-1)
             DO j=i+1,100
                ! Now we are ascending
                energy_points(j) = energy_func(this, o, coord_tmp + j * 0.01 * dir)

                write(6,*)j, energy_points(j)

                IF (energy_points(j) .lt. energy_points(j-1)) THEN
                   ! maxima found
                   this%deltaE = energy_points(j-1)-minim
                   WRITE(6,*) "Maxima found on descending hill"
                   WRITE(6,*) "deltaE", this%deltaE

                   return
                END IF
             END DO

             this%deltaE = minim-energy_points(0)
             WRITE(6,*) "No maxima found on initial descending hill. But minima found."
             WRITE(6,*) "deltaE", this%deltaE
             return
          END IF
       END DO
       WRITE(6,*) "No maxima found on initial descending hill."
       WRITE(6,*) "deltaE", this%deltaE
    end if



    return
    DO i=2,100
        coord_tmp = coord_tmp + (i * 0.01) * dir 
        en = energy_func(this, o, coord_tmp)
        r2= SUM(coord_tmp(1:2)**2) 
        if (r2 > (this%radius-2.0)**2) THEN
            return       
        end if
        write(6,*) "iiii",i, r2, en
    END DO

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    return

    dir = this%dr / SQRT(DOT_PRODUCT(this%dr,this%dr))
    coord_tmp = this%coord(:,o)
    
    DO i=1,20
        coord_tmp = coord_tmp + (i/20.0) * dir 
        en = energy_func(this, o, coord_tmp)
        r2= SUM(coord_tmp(1:2)**2) 
        if (r2 > (this%radius-2.0)**2) THEN
            return       
        end if
        write(6,*) "iiii",i, r2, en
    END DO
    return



   
    energy_points(0) = energy_func(this, o, coord_tmp)
    minimloc = 0
    maximloc = 0
    minim = energy_points(0)
    maxim = energy_points(0)
    j = 20
    DO i=1,j 
        coord_tmp = this%coord(:,o) + this%dr*(i*1.0d0/j)
        energy_points(i) = energy_func(this, o, coord_tmp)

        IF (energy_points(i) .lt. minim) THEN
           ! New minimum was found
           minim = energy_points(i)
           minimloc = i
        END IF

        IF (energy_points(i) .gt. maxim) THEN
           ! New maximum was found
           maxim = energy_points(i)
           maximloc = i
        END IF
  !     write(6,*) i, energy_points(i), this%deltaE, (i*1d0/j)
     END DO

     IF (maximloc .gt. minimloc) THEN
        ! deltaE > 0
        this%deltaE = maxim-minim
     ELSE
        ! deltaE < 0
        this%deltaE = minim-maxim
     END IF
  END SUBROUTINE PES

  SUBROUTINE mc_exchange_cavity_biased(this, pcn)
    !TODO
    !TODO test more this subroutine
    !TODO
    !TODO
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! Attempt to exchange a particle with a reservoir
    ! Cavity biased approach
    !----------------------------------------------------------------------!
    CLASS(Simulation), INTENT(INOUT) :: this
    REAL(8)          , INTENT(INOUT) :: pcn(:,:)
    REAL(8)                          :: energy, prob_acceptance, arg
    REAL(8)                          :: r2, x2
    REAL(8)                          :: coord_tmp(3)
    INTEGER(4)                       :: o, nt, k, n_cavity
    REAL(8), ALLOCATABLE             :: test_points(:,:)
    LOGICAL(1), ALLOCATABLE          :: in_cavity(:)
    REAL(8)                          :: radius, theta, pcorr
    
    IF (RAND() .lt. 0.5) THEN
       !---------------------------------------------------------------------------!
       !-------------------------------- Insertion --------------------------------!
       !---------------------------------------------------------------------------!
       nt=20
       ALLOCATE(test_points(3,nt), in_cavity(nt))

       n_cavity=0
       DO k=1,nt
          theta = 2*3.1415*RAND()
          radius = (this%radius-this%sigma_cnt)*sqrt(RAND())
          test_points(1,k) = radius * cos(theta)
          test_points(2,k) = radius * sin(theta)
          test_points(3,k) = this%length*RAND()
          in_cavity(k) = is_in_cavity(this, test_points(:,k))
          if (in_cavity(k) .eqv. .true.) THEN
            n_cavity = n_cavity + 1
          END IF
       END DO
 
       IF (n_cavity .eq. 0) THEN
          ! normal insertion attempt
          ! New particle at random poistion
          r2 = (this%radius-this%sigma_cnt) * (this%radius-this%sigma_cnt)                           ! Radius squared
          coord_tmp(1) = RAND()*(2*this%radius-2*this%sigma_cnt)-(this%radius-this%sigma_cnt)
          x2 = coord_tmp(1) * coord_tmp(1)                         ! Coordinate x squared
          coord_tmp(2) = RAND()*2*SQRT(r2-x2) - SQRT(r2 - x2)
          coord_tmp(3) = RAND()*this%length

          ! Energy new particle
          energy = particle_energy(this, this%npart+1 , coord_tmp)


          ! New acceptance criteria (Tildesley)
          arg = - energy/this%temperature
          arg = arg + LOG(this%activity*this%volume_eff*1d-30/REAL(this%npart+1))

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
       ELSE
          ! Insert a particle at one of the points which was found to be in a cavity
          DO k=1,nt
             IF (in_cavity(k) .eqv. .true.) THEN
                EXIT
             END IF
          END DO

          ! Energy new particle
          energy = particle_energy(this, this%npart+1 , test_points(:,k))

          pcn(1,this%npart+1) = pcn(1,this%npart+1) + (real(n_cavity) / real(nt))
          pcn(2,this%npart+1) = pcn(2,this%npart+1) + 1
          pcorr = pcn(1,this%npart+1) /  pcn(2,this%npart+1)


          arg = - energy/this%temperature
          arg = arg + LOG(pcorr*this%activity*this%volume_eff*1d-30/REAL(this%npart+1))
          !write(6,*) pcorr, arg

          IF (arg .gt. 0) THEN
             this%coord(:,this%npart+1) = test_points(:,k)
             this%npart = this%npart + 1
          ELSE IF (arg .gt. -30.0) THEN
             prob_acceptance = min(1.0d0,EXP(arg))
             IF ( RAND() .lt. prob_acceptance ) THEN
                this%coord(:,this%npart+1) = test_points(:,k)
                this%npart = this%npart + 1
             END IF
          END IF  

       END IF


       DEALLOCATE(test_points, in_cavity)

   ELSE

       !---------------------------------------------------------------------------!
       !--------------------------------- Deletion --------------------------------!
       !---------------------------------------------------------------------------!
       IF (this%npart .eq. 0) THEN
          ! There are no particles
          return
       END IF


       !TODO: test for overlap of particle upon insertion a displacement

       ! Choose a particle randomly
       o = INT(this%npart*RAND())+1

     
       IF (pcn(1,this%npart-1) .eq. 0) THEN
          pcorr = pcn(1,this%npart) /  pcn(2,this%npart)
       ELSE
          pcorr = pcn(1,this%npart-1) /  pcn(2,this%npart-1)
       END IF

       IF (RAND() .lt. (1-pcorr)**20) THEN
           write(6,*) "deletion normal" 
          ! Energy of particle to be removed
          energy = particle_energy(this, o, this%coord(:,o))


          ! OLD ACCEPTANCE CRITERIA (Frenkel and Smith)
          !arg =  exp(energy/this%temperature) * (this%npart * 1d0) / (this%pressure * this%beta * this%volume * 1d-30)
          ! New acceptance criteria (Tildesley)
          arg = energy/this%temperature
          arg = arg + LOG( REAL(this%npart) / (this%activity*this%volume_eff*1d-30)) 

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
  write(6,*) "deletion biased" 
          energy = particle_energy(this, o, this%coord(:,o))
          arg = energy/this%temperature
          arg = arg + LOG( pcorr * REAL(this%npart) / (this%activity*this%volume_eff*1d-30)) 

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
       END IF

    END IF



  END SUBROUTINE mc_exchange_cavity_biased



  FUNCTION is_in_cavity(this, test_point)
    CLASS(Simulation), INTENT(IN) :: this
    REAL(8),           INTENT(IN) :: test_point(3)
    INTEGER(4)                    :: l
    REAL(8)                       :: r
    LOGICAL(1)                    :: is_in_cavity


    is_in_cavity = .true.
    DO l=1,this%npart
       r = SQRT(SUM( (test_point-this%coord(:,l))**2 ))

       IF (r .lt. 3.2) THEN
          is_in_cavity = .false.
          return
       END IF
    END DO


  END FUNCTION IS_IN_CAVITY

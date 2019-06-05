MODULE pbc_displacement_module
  USE simulation_module
  IMPLICIT NONE


CONTAINS
  SUBROUTINE pbc_x(r, length)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! Applies period boundary conditions to the inter-particle distance
    ! r along the x direction.
    !----------------------------------------------------------------------!
    REAL(8), INTENT(INOUT) :: r(3)
    REAL(8), INTENT(IN)    :: length


    DO WHILE (r(1) .gt. length)
       r(1) = r(1) - length
    END DO

    DO WHILE(r(1) .lt. 0)
       r(1) = r(1) + length
    END DO
  END SUBROUTINE pbc_x

  SUBROUTINE pbc_y(r, length)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! Applies period boundary conditions to the inter-particle distance
    ! r along the y direction.
    !----------------------------------------------------------------------!
    REAL(8), INTENT(INOUT) :: r(3)
    REAL(8), INTENT(IN)    :: length


    DO WHILE (r(2) .gt. length)
       r(2) = r(2) - length
    END DO

    DO WHILE(r(2) .lt. 0)
       r(2) = r(2) + length
    END DO
  END SUBROUTINE pbc_y

  SUBROUTINE pbc_z(r, length)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! Applies period boundary conditions to the inter-particle distance
    ! r along the z direction.
    !----------------------------------------------------------------------!
    REAL(8), INTENT(INOUT) :: r(3)
    REAL(8), INTENT(IN)    :: length


    DO WHILE (r(3) .gt. length)
       r(3) = r(3) - length
    END DO

    DO WHILE(r(3) .lt. 0)
       r(3) = r(3) + length
    END DO
  END SUBROUTINE PBC_Z


  SUBROUTINE pbc_z_distance(r, length)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! Applies period boundary conditions to the inter-particle distance
    ! r along the z direction.
    !----------------------------------------------------------------------!
    ! SOURCE OF PBC:
    ! "Understanding Molecular Simulation:
    ! From Algorithms To Applications"
    ! Daan Frenkel, Berend Smit
    ! p. 68
    !
    ! PBC expression: xr=xi-xj
    ! PBC expression: xr=xr-box*nint(xr/box)
    !----------------------------------------------------------------------!
    REAL(8), INTENT(INOUT) :: r(3)
    REAL(8), INTENT(IN)    :: length

    r(3) = r(3) - length*NINT(r(3)/length)
  END SUBROUTINE pbc_z_distance

  FUNCTION new_pos(old_pos, dx)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! Returns a new position instead of a displacement
    ! (see FUNCTION displacement)
    !----------------------------------------------------------------------!
    ! SOURCE:
    ! "Understanding Molecular Simulation:
    ! From Algorithms To Applications"
    ! Daan Frenkel, Berend Smit
    ! p. 43
    !----------------------------------------------------------------------!
    REAL(8), INTENT(IN)    :: old_pos(:), dx
    REAL(8)                :: new_pos(3)

    new_pos(1) =  old_pos(1) + ( RAND()*2*dx - dx )
    new_pos(2) =  old_pos(2) + ( RAND()*2*dx - dx )
    new_pos(3) =  old_pos(3) + ( RAND()*2*dx - dx )
  END FUNCTION new_pos


  FUNCTION displacement(dx)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    ! This function returns a random displacement so that the
    ! x, y and z coordinates are in the range [-dx,dx].
    !----------------------------------------------------------------------!
    ! SOURCE:
    ! "Understanding Molecular Simulation:
    ! From Algorithms To Applications"
    ! Daan Frenkel, Berend Smit
    ! p. 43
    !----------------------------------------------------------------------!
    REAL(8), INTENT(IN)       :: dx
    REAL(8)                   :: displacement(3), dr, theta, psi,dt,ds, factor


    !ds = 5.2d-7
    !dt = 1d-20


    !factor = 2 * ds * dt
    factor = 2.0d0 *  dx ** 2
    dr = rand_normal(0.0d0,sqrt(factor)*1d0) !* 1d10

    displacement = (/RAND()*2d0 - 1.0d0,RAND()*2d0 - 1d0,RAND()*2d0 - 1d0/)
    displacement = displacement / SQRT(DOT_PRODUCT(displacement,displacement))
    displacement = displacement * dr


    !displacement(1) = dr*cos(theta)*sin(psi)
    !displacement(2) = dr*sin(theta)*sin(psi)
    !displacement(3) = dr*cos(psi)

    return 
  END FUNCTION displacement


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
END MODULE pbc_displacement_module

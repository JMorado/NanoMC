MODULE pbc_module
    USE simulation_module
    IMPLICIT NONE

CONTAINS
    SUBROUTINE pbc_z(r, cell)
        !========================================================================!
        ! Applies Periodic Boundary Conditions (PBC)                             !
        !                                                                        !
        ! Reference:                                                             !
        ! Understanding Molecular Simulation                                     !
        ! Smit B and Frenkel D                                                   !
        ! Elsevier, 2002. https://doi.org/10.1016/B978-0-12-267351-1.X5000-7.    !
        ! p. 68                                                                  !
        !                                                                        !
        ! PBC expression: xr=xi-xj                                               !
        ! PBC expression: xr=xr-box*nint(xr/box)                                 !
        !------------------------------------------------------------------------!
        ! cell(3)                   (in) : cell vectors                          !
        ! r(N,3)                 (inout) : array containing the positions        !
        !========================================================================!
        REAL(8), INTENT(INOUT) :: r(3)
        REAL(8), INTENT(IN)    :: cell(3)
        INTEGER(4)             :: dim

        ! xy
        !DO dim=1,2
        !    IF (ABS(cell(dim)) .gt. 1d-8) THEN
        !        r(dim) = r(dim) - cell(dim)*NINT(r(dim)/cell(dim))
        !    END IF
        !END DO

        ! z

        DO WHILE (r(3) .gt. cell(3))
            r(3) = r(3) - cell(3)
        END DO

        DO WHILE(r(3) .lt. 0)
            r(3) = r(3) + cell(3)
        END DO

        !DO WHILE (SQRT(SUM(coord_tmp(1:2)**2)) .ge. self%radius)
        !   dr =  displacement( self%disp )
        !   coord_tmp = self%coord(:,o) + dr
        !END DO
    END SUBROUTINE pbc

    SUBROUTINE pbc_z_distance(r, length)
        IMPLICIT NONE
        !========================================================================!
        ! Applies Periodic Boundary Conditions (PBC) along z direction           !
        !                                                                        !
        ! Reference:                                                             !
        ! Understanding Molecular Simulation                                     !
        ! Smit B and Frenkel D                                                   !
        ! Elsevier, 2002. https://doi.org/10.1016/B978-0-12-267351-1.X5000-7.    !
        ! p. 68                                                                  !
        !                                                                        !
        ! PBC expression: xr=xi-xj                                               !
        ! PBC expression: xr=xr-box*nint(xr/box)                                 !
        !------------------------------------------------------------------------!
        ! cell(3)                   (in) : cell vectors                          !
        ! r(N,3)                 (inout) : array containing the positions        !
        !========================================================================!
        REAL(8), INTENT(INOUT) :: r(3)
        REAL(8), INTENT(IN)    :: length

        r(3) = r(3) - length*NINT(r(3)/length)
    END SUBROUTINE pbc_z_distance
    
    
    SUBROUTINE reflection(simulation_instance,coord,dr)
        IMPLICIT NONE
        TYPE(Simulation), INTENT(INOUT)    :: simulation_instance

        REAL(8), INTENT(IN) :: dr(3)
        REAL(8), INTENT(INOUT) :: coord(3)

        REAL(8) :: a, b, c, d, hypoth,hypothd, dx, dy, initialx, initialy
        real(8) :: drestx, dresty, moduloi, alpha,angulo_i, zeta,beta,zeta_i, deslocx,deslocy

        
        deslocx = dr(1)
        deslocy = dr(2)
        hypoth = SQRT( (coord(1)+deslocx)**2 + (coord(2)+deslocy)**2)

        DO WHILE (simulation_instance%radius < hypoth)
            initialx = coord(1)
            initialy = coord(2)

            hypothd = SQRT( deslocx**2 + deslocy**2)
            dx = deslocx / hypothd
            dy = deslocy / hypothd

            a = 1.0d0
            b = 2.0d0 * ( (coord(1)*deslocx + coord(2)*deslocy ) / hypothd )
            c = (coord(1)*coord(1) + coord(2)*coord(2)) - (simulation_instance%radius*simulation_instance%radius)
            d = (-b + sqrt(b*b - 4*a*c)) / (2*a)

            ! Move particle to the circunferenc
            coord(1) = coord(1) +  (deslocx/hypothd)*d
            coord(2) = coord(2) +  (deslocy/hypothd)*d


            ! Compute displacement left
            drestx = (deslocx-dx*d)
            dresty = (deslocy-dy*d)

            ! Compute incidence angle which will be equal to reflexion angle
            moduloi = SQRT(coord(1)**2 + coord(2)**2)

            alpha = ACOS( (deslocx*coord(1)+deslocy*coord(2)) / (moduloi * hypothd))
            angulo_i = 3.14159265359 / 2.0d0 - alpha

            zeta = ATAN2(coord(2), coord(1))
            zeta_i = ATAN2(initialy, initialx)

            IF ((zeta_i .gt. zeta) .and. (zeta_i .lt. (zeta+3.14159265359))) THEN
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
        RETURN

    END SUBROUTINE reflection
END MODULE pbc_module

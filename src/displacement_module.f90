MODULE displacement_module
    IMPLICIT NONE

CONTAINS
    FUNCTION get_displacement_uniform(dr) RESULT(displacement)
        !========================================================================!
        ! This function returns a unifrom random displacement so that the        !
        ! x, y and z coordinates are in the range [-dr,dr]                       !
        !                                                                        !
        ! Reference:                                                             !
        ! Understanding Molecular Simulation                                     !
        ! Smit B and Frenkel D                                                   !
        ! Elsevier, 2002. https://doi.org/10.1016/B978-0-12-267351-1.X5000-7.    !
        ! pp. 43                                                                 !
        !------------------------------------------------------------------------!
        ! dr                           (in) : displacement size                  !
        !========================================================================!
        IMPLICIT NONE
        REAL(8), INTENT(IN)            :: dr
        REAL(8)                        :: displacement(3)
        REAL(8)                        :: theta, psi

        displacement = (/RAND()*2d0 - 1.0d0,RAND()*2d0 - 1d0,RAND()*2d0 - 1d0/)
        displacement = displacement / SQRT(DOT_PRODUCT(displacement,displacement))
        displacement = displacement * dr

        ! Using polar coordinates
        !displacement(1) = dr*cos(theta)*sin(psi)
        !displacement(2) = dr*sin(theta)*sin(psi)
        !displacement(3) = dr*cos(psi)

        RETURN
    END FUNCTION get_displacement_uniform

    FUNCTION get_displacement_gaussian(dx) RESULT(displacement)
        IMPLICIT NONE
        !========================================================================!
        ! This function returns a gaussian random displacement                   !
        !------------------------------------------------------------------------!
        ! dx                           (in) : displacement size                  !
        !========================================================================!
        REAL(8), INTENT(IN)       :: dx
        REAL(8)                   :: displacement(3)
        REAL(8)                   :: dr, theta, psi,dt,ds, factor

        factor = 2.0d0 * dx ** 2
        dr = rand_normal(0.0d0,sqrt(factor)*1d0) !* 1d10

        displacement = (/RAND()*2d0 - 1.0d0,RAND()*2d0 - 1d0,RAND()*2d0 - 1d0/)
        displacement = displacement / SQRT(DOT_PRODUCT(displacement,displacement))
        displacement = displacement * dr

        RETURN
    END FUNCTION get_displacement_gaussian

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
END MODULE displacement_module
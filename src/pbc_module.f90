MODULE pbc_module
    IMPLICIT NONE

CONTAINS
    SUBROUTINE pbc(r, cell)
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

        DO dim=1,3
            IF (ABS(cell(dim)) .gt. 1d-8) THEN
                r(dim) = r(dim) - cell(dim)*NINT(r(dim)/cell(dim))
            END IF
        END DO

        !dr =  displacement( self%disp )
        !coord_tmp = self%coord(:,o) + dr(:)
        !DO WHILE (SQRT(SUM(coord_tmp(1:2)**2)) .ge. self%radius)
        !   dr =  displacement( self%disp )
        !   coord_tmp = self%coord(:,o) + dr
        !END DO
        !CALL pbc_z(coord_tmp, self%length)
    END SUBROUTINE pbc
END MODULE pbc_module

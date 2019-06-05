MODULE pbc_displacement_module
  USE simulation_module
  IMPLICIT NONE


CONTAINS
  SUBROUTINE pbc_z(r, length)
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
    ! 
    ! WARNING! Frenkel Smit expression can cause problem because NINT(-0.2)=0
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
    ! Give particle random displacement
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
    REAL(8), INTENT(IN)    :: dx
    REAL(8)                :: displacement(3)


    displacement = (/RAND()*2*dx - dx,RAND()*2*dx - dx,RAND()*2*dx - dx/)
  END FUNCTION displacement
END MODULE pbc_displacement_module

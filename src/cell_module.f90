MODULE cell_module
  IMPLICIT none

  TYPE, PUBLIC :: Cell
     ! Axial Lengths and Angles
     REAL(8)  :: a, b, c                           ! Cell axial lenghts
     REAL(8)  :: alpha, beta, gamma                ! Cell angles

   CONTAINS
     PROCEDURE :: cell_volume => cell_volume

  END TYPE CELL

CONTAINS
  FUNCTION cell_volume(self)
    !========================================================================!
    ! This function calculates the volume of the simulation cell             !
    !                                                                        !
    ! Reference:                                                             !
    ! http://webmineral.com/help/CellDimensions.shtml#.YLhfqCYo_Io           !
    !------------------------------------------------------------------------!
    CLASS(Cell), INTENT(IN)        :: self
    REAL(8)                        :: cell_volume

    cell_volume = self%a*self%b*self%c*(1- &
            COS(self%alpha)*COS(self%alpha)-COS(self%beta)*COS(self%beta)-COS(self%gamma)*COS(self%gamma))

    cell_volume = cell_volume + 2*(COS(self%alpha)*COS(self%beta)*COS(self%gamma)) ** (1/2.)

    RETURN
  END FUNCTION cell_volume

END MODULE cell_module

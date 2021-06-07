MODULE cell_module
  IMPLICIT none

  TYPE, PUBLIC :: Cell
     ! CNT parameters
     REAL(8)               :: length                           ! SWCNT length
     REAL(8)               :: radius                           ! SWCNT radius


   CONTAINS
     PROCEDURE :: Cell

  END TYPE CELL

CONTAINS
  SUBROUTINE cell(this)
    !-------------------------------------------------------
    ! Cell type object constructor
    !------------------------------------------------------
    TYPE(Cell)            :: this

    this%length = 10
    this%radius = 10

  END SUBROUTINE init

  SUBROUTINE str(this)
    TYPE(Cell)            :: this

    WRITE(6,*) "LOL"
  END SUBROUTINE str

END MODULE cell_module

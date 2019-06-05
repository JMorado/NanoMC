

MODULE test
    USE FRUIT
    IMPLICIT NONE


CONTAINS                          !fortran 95 limits subroutine name to 31 char.
  subroutine test_calculator_should_produce_4_when_2_and_2_are_inputs
    use calculator, only: add
    integer:: result

    call add (2,2,result)
    call assert_equals (4, result)
  end subroutine test_calculator_should_produce_4_when_2_and_2_are_inputs
END MODULE test


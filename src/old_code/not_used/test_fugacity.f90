PROGRAM test
  USE virial
  IMPLICIT NONE


  REAL(8) :: press, temp, fug


  temp = 200d0
  press = 50.0                  !bar
  fug = log_fugacity(press,temp)
  write(6,*)fug, exp(fug)



END PROGRAM test

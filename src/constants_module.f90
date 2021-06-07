MODULE constants_module
  IMPLICIT NONE

  REAL(8), PARAMETER :: pi = ATAN(1.0)*4.0d0
  REAL(8), PARAMETER :: kb = 1.38064852d-23                  ! J K^-1
  REAL(8), PARAMETER :: atomicmass_to_kg = 1.66054d-27       ! amu -> kg
  REAL(8), PARAMETER :: planck = 6.62607004d-34              ! m^2 kg s^-1

END MODULE constants_module

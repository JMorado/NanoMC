MODULE virial
  IMPLICIT NONE


CONTAINS
  FUNCTION LOG_FUGACITY(pressure, temperature) RESULT(fugacity)
    IMPLICIT NONE
    !-----------------------------------------------------------!
    ! Fugacity coefficients for hydrogen gas between
    ! 0 and 1000 centigrades for pressures up to 3000 atm
    !-----------------------------------------------------------!
    REAL(8), INTENT(IN) :: pressure, temperature
    REAL(8)             :: fugacity
    REAL(8)             :: c1, c2, c3

    c1 = EXP(-3.8402*temperature**(1d0/8d0) + 0.5410)
    c2 = EXP(-0.1263*temperature**(1d0/2d0) -15.980)
    c3 = 300d0*EXP(-0.11901*temperature - 5.941)


    fugacity = c1*pressure - c2*pressure*pressure + c3*(EXP(-pressure/300d0)-1d0)
  END FUNCTION LOG_FUGACITY


  FUNCTION factorial(n)
    IMPLICIT NONE
    !----------------------------------------------------
    ! Computes the factorial of n.
    !----------------------------------------------------
    INTEGER(4), INTENT(IN)    :: n
    INTEGER(4)                :: i
    REAL(8)                   :: factorial

    IF (n .eq. 0) THEN
       factorial = 1.0d0
       return
    ELSE
       factorial = 1.0d0
       DO i=1,n
          factorial = factorial * i*1d0
       END DO
    END IF
  END FUNCTION factorial

  FUNCTION b2_coeff(nterms,sigma,epsilon,temperature) RESULT(b2)
    IMPLICIT NONE
    !----------------------------------------------------
    ! Computes the second virial coefficient for the
    ! Lennard-Jones 12-6 potential.
    !
    ! Units : A^3
    !
    ! SOURCE:
    ! http://demonstrations.wolfram.com/SecondVirialCoefficientsUsingTheLennardJonesPotential/
    !----------------------------------------------------
    INTEGER(4), INTENT(IN)    :: nterms
    REAL(8),    INTENT(IN)    :: sigma, epsilon, temperature
    REAL(8)                   :: sum, sum_local
    REAL(8)                   :: b2
    INTEGER(4)                :: n


    b2  = -3.14159265359 * (2d0/3d0) * sigma ** 3

    sum = 0
    DO n=0,nterms
       sum_local = 2 ** (n-3d0/2d0) / factorial(n)
       sum_local = sum_local * GAMMA((2d0*n-1d0)/4d0)
       sum_local = sum_local * (epsilon / temperature) ** ((2d0*n+1d0)/4d0)


       sum = sum + sum_local
    END DO

    b2 = b2 * sum

  END FUNCTION b2_coeff



END MODULE VIRIAL

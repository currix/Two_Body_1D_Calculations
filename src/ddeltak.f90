FUNCTION ddeltak(n1, n2)
!     
!     Kronecker delta function
!
!     
!     by Currix TM
!
  USE nrtype
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: n1, n2
  REAL(KIND = DP) :: ddeltak
  !
  ddeltak = 0.0_DP
  !
  IF (n1 == n2)  ddeltak = 1.0_DP
  !
END FUNCTION ddeltak


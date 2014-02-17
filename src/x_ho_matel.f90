FUNCTION xhomatel(m, n)
  !     
  !     Harmonic Oscillator Matrix Element < m | x | n >
  !
  !     
  !     by Currix TM
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: m, n
  REAL(KIND = DP) :: xhomatel
  !
  SELECT CASE (m-n)
     CASE (1)
        xhomatel = SQRT((REAL(n,DP) + 1.0_DP)/2.0_DP)
     CASE (-1)
        xhomatel = SQRT(REAL(n, DP)/2.0_DP)
     CASE DEFAULT
        xhomatel = 0.0_DP
  END SELECT
  !
END FUNCTION xhomatel

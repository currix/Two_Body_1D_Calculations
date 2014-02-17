FUNCTION x2homatel(m, n)
  !     
  !     Harmonic Oscillator Matrix Element < m | x^2 | n >
  !
  !     
  !     by Currix TM
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: m, n
  REAL(KIND = DP) :: x2homatel
  !
  SELECT CASE (m-n)
     CASE (2)
        x2homatel = SQRT( (REAL(n,DP) + 1.0_DP) * (REAL(n,DP) + 2.0_DP) )/2.0_DP
     CASE (-2)
        x2homatel = SQRT( (REAL(n,DP) - 1.0_DP) * REAL(n,DP) )/2.0_DP
     CASE (0)
        x2homatel = (REAL(n,DP) + 0.5_DP)
     CASE DEFAULT
        x2homatel = 0.0_DP
  END SELECT
  !
END FUNCTION x2homatel

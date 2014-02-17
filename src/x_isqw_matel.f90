FUNCTION x_box_matel(m, n)
  !     
  USE nrtype
  !
  !     ISQW Matrix Element < m | x | n >
  !
  !     m, n = 1,2,3,...
  !     
  !     by Currix TM
  !
  IMPLICIT NONE
  INTEGER(KIND = I4B), INTENT(IN) :: m, n
  REAL(KIND = DP) :: x_box_matel
  !
  x_box_matel = 0.0_DP
  !
  IF (MOD(m,2) == MOD(n,2)) THEN
     !     SAME PARITY
     x_box_matel = 0.0_DP
  ELSE 
     !     DIFFERENT PARITY
     x_box_matel = (-1.0_DP)**((m-n-1)/2)/(m-n)**2 + &
          (-1.0_DP)**((m+n-1)/2)/(m+n)**2 
  ENDIF
  !
  RETURN
  !
END FUNCTION x_box_matel

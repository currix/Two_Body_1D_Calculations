FUNCTION x2_box_matel(m, n)
  !     
  !     ISQW Matrix Element < m | x^2 | n > except a factor of 8 x_b^2/Pi^2
  !
  !     m, n = 1,2,3,...
  !     
  !     by Currix TM
  !
  USE nrtype
  USE constants
  !
  IMPLICIT NONE
  INTEGER(KIND = I4B) :: m, n
  !
  REAL(KIND = DP) :: x2_box_matel
  !
  x2_box_matel = 0.0_DP
  !
  IF (MOD(m,2) /= MOD(n,2)) THEN
     !     DIFFERENT PARITY
     x2_box_matel = 0.0_DP
  ELSE 
     !     SAME PARITY
     !Same parity case (odd-odd)
     IF ( MOD(n,2) /= 0.0_dp ) THEN  
        IF ( m == n ) THEN
           x2_box_matel = 2.0_DP * &
                ( &
                1.0_DP/6.0_DP - 1.0_DP/((PI_D*REAL(n,DP))**2) &
                )    
        ELSE 
           x2_box_matel = ( 8.0_dp /(PI_D**2)) * &
                ( &
                ((-1.0_DP)**((m-n)/2)) / (REAL(m-n,DP)**2)  &
                + &
                ((-1.0_DP)**((n+m)/2)) / (REAL(n+m,DP)**2) &
                )
        ENDIF
        !Same parity case (even-even)
     ELSE 
        IF ( m == n ) THEN
           x2_box_matel = 2.0_DP * &
                ( &
                1.0_DP/6.0_DP - 1.0_DP/((PI_D*REAL(n,DP))**2) &
                )            
        ELSE
           x2_box_matel = ( 8.0_DP/(PI_D**2)) * &
                ( &
                ((-1.0_DP)**((m-n)/2)) / (REAL(m-n,DP)**2) &
                - &
                ((-1.0_DP)**((m+n)/2)) / (REAL(n+m,DP)**2) &
                )
        ENDIF
     ENDIF
  ENDIF
  !
  RETURN
  !
END FUNCTION x2_box_matel

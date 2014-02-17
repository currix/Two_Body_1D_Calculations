FUNCTION F(X,Y)
  ! for D01DAF
  !
  USE nrtype
  USE vardef_ri
  USE module_2body
  !
  IMPLICIT NONE
  !
  ! Arguments 
  REAL(KIND = DP), INTENT(IN) :: X, Y
  REAL(KIND = DP) :: F
  !
  INTEGER(KIND = I4B) :: X1,X2,I,J
  ! F = PSI_i*(x1+x2)*PSI_j COME FA A CONOSCERE INDEX_I E INDEX_J????????
  !
  ! a partire da X e Y devo calcolare il posto più vicino in x grid
  X1 = NINT((X-X_min)/delta_X)+1
  X2 = NINT((Y-X_min)/delta_X)+1
  I = F_index_i
  J = F_index_j
  !
  F = avec_2b_X(x1,x2,I)*(X_grid(x1)+X_grid(Y))*avec_2b_X(x1,x2,J)
  !
  print*, X1, X2, F
  !
  RETURN
END FUNCTION F
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION PHI1(Y)
  ! for D01DAF
  !
  USE nrtype
  USE vardef_ri
  !
  IMPLICIT NONE
  !
  ! Arguments 
  REAL(KIND = DP), INTENT(IN) :: Y
  !
  REAL(KIND = DP) :: PHI1
  !
  PHI1 = X_min
  !
END FUNCTION PHI1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION PHI2(Y)
  ! for D01DAF
  !
  USE nrtype
  USE vardef_ri
  !
  IMPLICIT NONE
  !
  ! Arguments 
  REAL(KIND = DP), INTENT(IN) :: Y
  !
  REAL(KIND = DP) :: PHI2
  !
  PHI2 = X_max
  !
END FUNCTION PHI2

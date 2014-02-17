MODULE vardef_ri
  !
  ! $Id: var_definition.f90,v 1.1 2013/07/01 15:54:04 curro Exp $
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  !
  ! XGRID DIMENSION
  INTEGER(KIND = I4B) :: dim_X
  ! Xmin Xmax and delta_X step
  REAL(KIND = DP) :: X_min, X_max, Delta_X
  ! XGRID
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: X_Grid      
  !
  ! for the numerical E1 and E2 calculation
  INTEGER(KIND = I4B) :: F_index_I, F_index_J
  !
END MODULE vardef_ri

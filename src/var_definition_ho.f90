MODULE vardef_ri
  !
  ! $Id: var_definition_ho.f90,v 1.1 2013/06/07 07:26:35 curro Exp $
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
  ! One body variables
  !
  ! DIMENSION OF THE 1body HO BASIS
  INTEGER(KIND = I4B) :: dim_1b_HO
  !
  ! HARMONIC BASIS, EIGENVALUES AND EIGENVECTORS
  REAL(KIND = DP),  DIMENSION(:), ALLOCATABLE :: Aval_1b_Har
  REAL(KIND = DP),  DIMENSION(:,:), ALLOCATABLE :: Har_1b_Bas, Avec_1b_Har, Avec_1b_Har_X
  ! EIGENVECTOR DERIVATIVES
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Avec_1b_Har_Der_X
  !
  !
  !
END MODULE vardef_ri

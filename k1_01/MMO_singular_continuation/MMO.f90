!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   lor :     The Lorenz Equations
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
	  
	  !set variable type for parameters and state variables
	  DOUBLE PRECISION a,b,x,y,Period,alpha,mu,kappa,zeta,epsilon_sq,epsilon_b,delta,k1,k2,k3,k4,k5,k6,k7,k8,kn7,artificial

	  a = U(1)
	  b = U(2)
	  x = U(3)
	  y = U(4)
  
	  !define parameters
    
	  Period= PAR(11)
  
	  mu= PAR(1)
	  alpha= PAR(2)
	  epsilon_b= PAR(3)
	  kappa= PAR(4)
	  epsilon_sq= PAR(5)
	  zeta= PAR(6)
	  delta= PAR(7)
	  artificial=PAR(8)
  
	  F(1)= (mu - alpha*a - a*b*y)/artificial 
	  F(2)= epsilon_b*(1 - b*x - a*b*y)
	  F(3)= (b*x-x**2 + 3*a*b*y - zeta*x + delta)/(epsilon_sq)
	  F(4)= kappa*(x**2 - y - a*b*y)/(epsilon_sq)
  

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

! Parameter values for the starting orbit in lor.dat (original parameter values in heteroclinic paper, check eigenvectors.m for param values with k_1=0.1) :
     	PAR(1:2) = (/ 0.96969696969696969696969696969697, 0.091226247860007630974266135126527 /) 
		PAR(3:4) = (/ 0.0037056277056277056277056277056277, 3.7962796283345614511972598847933 /) 
		PAR(5:6) = (/ 0.053971962616822429906542056074766, 0.98473192783466186187392528100194 /) 
		PAR(7) = 0.000012121212121212121212121212121212
		PAR(8) = 1.0

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS

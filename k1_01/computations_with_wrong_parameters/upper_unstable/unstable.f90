!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   one_slow :  slow manifold of system corresponding to third bullet
!	point in Kuhen and Szmolyan 2015
!   Use method from Farjami, Kirk, Osinga 2017 (?) to find stable and 
!	unstable manifolds of slow manifolds.
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION a, b, x, y, Period, alpha, mu, kappa, zeta, epsilon_sq, epsilon_b, delta, k1, k2, k3, k4, k5, k6,k7,k8,kn7
	  
	  !system variables
	   a= U(1)
	   b= U(2)
	   x= U(3)
	   y= U(4)
	  
	   !system parameters
	   
	   k1= PAR(1)
	   k2= PAR(2)
	   k3= PAR(3)
	   k4= PAR(4)
	   k5= PAR(5)
	   k6= PAR(6)
	   k7= PAR(7)
	   kn7= PAR(8)
	   k8= PAR(9)
	      
	   Period= PAR(11)
	   
	   mu= k7/k8
	   alpha= (k1*k5*kn7)/(k3*k8*sqrt(2*k2*k8))
	   epsilon_b= ((k1**2)*k5)/(2*k2*k3*k8)
	   kappa= (sqrt(2*k2*k8))/k5
	   epsilon_sq= (k3*k8)/(k1*k5)
	   zeta= k4/(sqrt(2*k2*k8))
	   delta= k6/k8

       F(1)= Period*(mu - alpha*a - a*b*y) 
       F(2)= Period*epsilon_b*(1 - b*x - a*b*y)
       F(3)= Period*(b*x-x**2 + 3*a*b*y - zeta*x + delta)/(epsilon_sq)
	   F(4)= Period*kappa*(x**2 - y - a*b*y)/(epsilon_sq)

      END SUBROUTINE FUNC
!---------------------------------------------------------------------- 
      SUBROUTINE STPNT(NDIM,U,PAR,T) 
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

	  PAR(1:24)= (/ 0.1, 250.0, 0.035, 20.0, 5.35, 1E-4, 0.8, 0.1, 0.825,& !PAR(1:9) table 1 parameters
	  0.0, 0.0,& !PAR(10:11) PAR(11) is where AUTO stores the value of the period so we want to avoid it
	  6.8333504443159869005123251424624, 0.28, 0.72616596871012748190465871734031, 0.18100096576321164240402462599352,& !PAR(12:15) equilibrium of fast subsystem at b=0.3 startpoint
	  6.8333504443159869005123251424624, 0.28, 0.72616596871012748190465871734031, 0.18100096576321164240402462599352,& !PAR(16:19) equilibrium of fast subsystem at b=0.3 endpoint
	  10.399327486524577796943974841419, 0.0, & !PAR(20:21) unstable stpt eigenvalue and radius
	  -0.077757481143788528554766091127121, 0.90024085357178045422514716243999, 0.42839255325755315283349487867513 /) !PAR (22:24) unstable stpt eigenvector
	  
	  U(1)= PAR(12)
	  U(2)= PAR(13)
	  U(3)= PAR(14)
	  U(4)= PAR(15)
	  
      END SUBROUTINE STPNT
!---------------------------------------------------------------------- 
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
	  
	  DOUBLE PRECISION a, b, x, y, alpha, mu, kappa, zeta, epsilon_sq, epsilon_b, delta 
	  DOUBLE PRECISION k1, k2, k3, k4, k5, k6, k7, k8, kn7, M(3,3)
	  INTEGER D
	  
	  !startpoint BCs
	  FB(1)= U0(1) - (PAR(12) + PAR(21)*PAR(22))
	  FB(2)= U0(2) - PAR(13)
	  FB(3)= U0(3) - (PAR(14) + PAR(21)*PAR(23))
	  FB(4)= U0(4) - (PAR(15) + PAR(21)*PAR(24))
	  
	  !endpoint BCs
	  FB(5)= U1(1) - PAR(16)
	  FB(6)= U1(2) - PAR(17)
	  FB(7)= U1(3) - PAR(18)
	  FB(8)= U1(4) - PAR(19)

      END SUBROUTINE BCND

!---------------------------------------------------------------------- 
      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
      DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

      END SUBROUTINE ICND

      SUBROUTINE PVLS
      END SUBROUTINE PVLS

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
!----------------------------------------------------------------------   

	  SUBROUTINE JACOB(NDIM,PAR,U,M)
!     ---------- ----	
	  IMPLICIT NONE
	  INTEGER, INTENT(IN):: NDIM
	  DOUBLE PRECISION, INTENT(IN):: PAR(NDIM), U(NDIM)
	  DOUBLE PRECISION, INTENT(OUT):: M(NDIM,NDIM)
	  
	 DOUBLE PRECISION a, b, x, y, Period, alpha, mu, kappa, zeta, epsilon_sq, epsilon_b, delta, k1, k2, k3, k4, k5, k6,k7,k8,kn7
	  
	  k1= PAR(1)
	  k2= PAR(2)
	  k3= PAR(3)
	  k4= PAR(4)
	  k5= PAR(5)
	  k6= PAR(6)
	  k7= PAR(7)
	  kn7= PAR(8)
	  k8= PAR(9)
      
	  mu= k7/k8
	  alpha= (k1*k5*kn7)/(k3*k8*sqrt(2*k2*k8))
	  epsilon_b= ((k1**2)*k5)/(2*k2*k3*k8)
	  kappa= (sqrt(2*k2*k8))/k5
	  epsilon_sq= (k3*k8)/(k1*k5)
	  zeta= k4/(sqrt(2*k2*k8))
	  delta= k6/k8
	  
	  a= U(1)
	  b= U(2)
	  x= U(3)
	  y= U(4)
	
	  M(1,1)= - alpha - b*y
	  M(1,2)= 0
	  M(1,3)= -a*b
	  
	  M(2,1)= (3*b*y)/epsilon_sq
	  M(2,2)= -(2*x - b + zeta)/epsilon_sq
	  M(2,3)= (3*a*b)/epsilon_sq
	  
	  M(3,1)= -(b*kappa*y)/epsilon_sq
	  M(3,2)= (2*kappa*x)/epsilon_sq
	  M(3,3)= -(kappa*(a*b + 1))/epsilon_sq
	  !WRITE(*,*) M(1,1)
  	END SUBROUTINE JACOB
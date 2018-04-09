!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   Stable manifold of slow manifold of upper branch of critical manifold.
!   We are using the same parameter values as in Desroches, Krauskopf, Osinga (2008).
!   The system has been rescaled to a system that is appropriate for capturing SAOs
!   in MMOs in Olsen model phase space.  This sytem can be found in Kuehn and 
!   Szmolyan (2015).
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

      DOUBLE PRECISION Period, alpha, mu, kappa, zeta, epsilon, epsilon_b

	  PAR(1:37)= (/ 0.34, 250.0, 0.035, 20.0, 5.35, 1E-4, 0.8, 0.1, 0.825,& !PAR(1:9) table 1 parameters
	  0.0, 0.0,& !PAR(10:11) PAR(11) is where AUTO stores the value of the period so we want to avoid it
	  10.436686738509971973100582955989, 0.75, 0.14087341527311976705996439770791, 0.0022481206783085090720331432231075,& !PAR(12:15) initial endpoint
	  10.436686738509971973100582955989, 0.75, 0.14087341527311976705996439770791, 0.0022481206783085090720331432231075,& !PAR(16:19) initial startpoint
	  -0.057393376079136078321021937571087, 0.99785001695446466542873519512164, 0.031644020702392266436058747410753,& !PAR(20:22) unstable eigenvector of initial startpoint
	  4.2228130176691198153663779476142, 0.0,&  !PAR(23:24) unstable eigenvalue and radius
	  -0.090816196238710541768445151887843, 0.0,& !PAR(25:26) weak stable eigenvalue and radius 
	  0.99999707302971025753674627752219, -0.0024046247390376816073830879728421, -0.00026779073310667623012098205397044,& !PAR(27:29) weak stable eigenvector of initial startpoint
	  -634.70400749900674383867893956826, 0.0,& !PAR(30:31) strong stable eigenvalue and radius
	  0.010123197508043369256178567694836, -0.57122343907965346994992714410738, 0.82073217526683244163977033521288,&  !PAR(32:34) strong stable eigenvector of initial startpoint
	  0.0, 0.0, 0.0 /) !PAR(35:37) variable used for releasing startpoint from unstable eigenvector in second step
	  
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
	  
	  DOUBLE PRECISION a1, b1, x1, y1, a2, b2, x2, y2, alpha, mu, kappa, zeta, epsilon_sq, epsilon_b, delta 
	  DOUBLE PRECISION k1, k2, k3, k4, k5, k6, k7, k8, kn7, M(3,3)
	  INTEGER D
	  
	  !endpoint BCs
	  FB(1)= U0(1) - (PAR(12) + PAR(26)*PAR(27) + PAR(31)*PAR(32))
	  FB(2)= U0(2) - PAR(13)
	  FB(3)= U0(3) - (PAR(14) + PAR(26)*PAR(28) + PAR(31)*PAR(33))
	  FB(4)= U0(4) - (PAR(15) + PAR(26)*PAR(29) + PAR(31)*PAR(34))
	  
	  !startpoint BCs
	  FB(5)= U1(1) - (PAR(16) + PAR(24)*PAR(20) + PAR(35))
	  FB(6)= U1(2) - PAR(17)
	  FB(7)= U1(3) - (PAR(18) + PAR(24)*PAR(21) + PAR(36))
	  FB(8)= U1(4) - (PAR(19) + PAR(24)*PAR(22) + PAR(37))
	  
	  IF(NBC==8) RETURN
	  
	  !calculate new fast subsystem steady state for endpoint that depends on B=PAR(13)
	  
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
	  
	  a1= PAR(12)
	  b1= PAR(13)
	  x1= PAR(14)
	  y1= PAR(15)
	  
      FB(9)= (mu - alpha*a1 - a1*b1*y1) 
      FB(10)= (b1*x1-x1**2 + 3*a1*b1*y1 - zeta*x1 + delta)/(epsilon_sq)
   	  FB(11)= kappa*(x1**2 - y1 - a1*b1*y1)/(epsilon_sq)
	  
	  !endpoint strong stable eigenvector
	  
	  D = NDIM - 1
	  CALL JACOB(D,PAR,PAR(12:15),M)
	  
	  FB(12:14)= MATMUL(M, PAR(32:34)) - PAR(30)*PAR(32:34)
	  FB(15)= DOT_PRODUCT(PAR(32:34),PAR(32:34))-1.0
	  
	  !endpoint weak stable eigenvector
	  
	  D = NDIM - 1
	  CALL JACOB(D,PAR,PAR(12:15),M)
	  
	  FB(16:18)= MATMUL(M, PAR(27:29)) - PAR(25)*PAR(27:29)
	  FB(19)= DOT_PRODUCT(PAR(27:29),PAR(27:29))-1.0
	  
	  !startpoint equilibrium
	  
	  a2= PAR(16)
	  b2= PAR(17)
	  x2= PAR(18)
	  y2= PAR(19)
	  
      FB(20)= (mu - alpha*a2 - a2*b2*y2) 
      FB(21)= (b2*x2-x2**2 + 3*a2*b2*y2 - zeta*x2 + delta)/(epsilon_sq)
   	  FB(22)= kappa*(x2**2 - y2 - a2*b2*y2)/(epsilon_sq)
	  
	  !startpoint weak eigenvector
	  CALL JACOB(D,PAR,PAR(16:19),M)
	  
	  FB(23:25)= MATMUL(M, PAR(20:22)) - PAR(23)*PAR(20:22)
	  FB(26)= DOT_PRODUCT(PAR(20:22),PAR(20:22))-1.0

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
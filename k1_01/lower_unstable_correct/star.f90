!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   one_slow :  slow manifold of system corresponding to third bullet
!	point in Kuhen and Szmolyan 2015
!   Use method from Farjami, Kirk, Osinga 2017 (?) to find stable and 
!	unstable manifolds of slow manifolds.
!   NOTE the endpoint and startpoint with respect to the flow are 
!	respectively the startpoint and endpoint with respect to the 
!	boundary conditions.  This is because we're continuing in backwards
!	time.
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

	  PAR(1:37)= (/ 0.1, 250.0, 0.035, 20.0, 5.35, 1E-5, 0.8, 0.1, 0.825,& !PAR(1:9) table 1 parameters
	  0.0, 0.0,& !PAR(10:11) PAR(11) is where AUTO stores the value of the period so we want to avoid it
	  1.9788805938744591214490646144441, 0.45, 1.2943670588295220112627801603527, 0.88621496482816810457094286066201,& !PAR(12:15) initial startpoint
	  1.9788805938744591214490646144441, 0.45, 1.2943670588295220112627801603527, 0.88621496482816810457094286066201,& !PAR(16:19) initial endpoint
	  -0.054113704099610944296001741255979, 0.99849139466813268778285899059149, 0.0093081578363839335758060423184044,& !PAR(20:22) unstable eigenvector of initial startpoint
	  1.5504068239584284429516954472173, 0.0,&  !PAR(23:24) unstable eigenvalue and radius
	  -0.091187050327067803627843733524516, 0.0,& !PAR(25:26) weak stable eigenvalue and radius 
	  -0.33346660064130737047731088355735, 0.52657049771917792536499935951078, 0.78199970408466783390660184501199,& !PAR(27:29) weak stable eigenvector of initial startpoint
	  -746.65728903445486119436933678594, 0.0,& !PAR(30:31) strong stable eigenvalue and radius
	  -0.13577820633870756410460892892982, 0.58428865940463036139575620763901, 0.80010689359271440375124350611793,&  !PAR(32:34) strong stable eigenvector of initial startpoint
	  10.606801438770783490352775268692, 0.047916398397789360680124248762032, 0.00021770859363502248771163454388334 /) !PAR(35:37) variables used for releasing startpoint from unstable eigenvector in second step
	  
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
	  
	  !endpoint BCs
	  FB(1)= U0(1) - (PAR(12) + PAR(26)*PAR(27) + PAR(31)*PAR(32))
	  FB(2)= U0(2) - PAR(13)
	  FB(3)= U0(3) - (PAR(14) + PAR(26)*PAR(28) + PAR(31)*PAR(33))
	  FB(4)= U0(4) - (PAR(15) + PAR(26)*PAR(29) + PAR(31)*PAR(34))
	  
	  !startpoint BCs
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

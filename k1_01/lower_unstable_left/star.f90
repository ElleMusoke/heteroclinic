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

	  !PAR(1:37)= (/ 0.1, 250.0, 0.035, 20.0, 5.35, 1E-4, 0.8, 0.1, 0.825,& !PAR(1:9) table 1 parameters
	  !0.0, 0.0,& !PAR(10:11) PAR(11) is where AUTO stores the value of the period so we want to avoid it
	  !0.81121984904465478078341498694954, 0.8000, 1.0217186418928954008370196563427, 0.6330650413608834359467877018843,& !PAR(12:15) initial endpoint
	  !0.81121984904465478078341498694954, 0.8000, 1.0217186418928954008370196563427, 0.6330650413608834359467877018843,& !PAR(16:19) initial startpoint
	  !-0.078125183105597024611054555384932, 0.61135628225381851355442321377554, 0.77461280007098070706464062296643,& !PAR(20:22) real part of complex conjugate eigenvectors
	  !0.0, 0.0,&  !PAR(23:24) strong unstable eigenvalue of startpoint and radius
	  !1.5999826631956767179885167014501, 0.0,& !PAR(25:26) weak unstable eigenvalue and radius 
	  !0.13339172174634352552632137932572, 0.048183092981327620190914515143416, 0.0,& !PAR(27:29) imaginary part of complex conjugate eigenvectors
	  !0.0, 0.0,& !PAR(30:31) strong stable eigenvalue and radius
	  !0.0, 0.0, 0.0,&  !PAR(32:34) strong stable eigenvector of initial startpoint is not used
	  !10.606801438770783490352775268692, 0.047916398397789360680124248762032, 0.00021770859363502248771163454388334 /) !PAR(35:37) variables used for releasing startpoint from unstable eigenvector in second step
	  
	  PAR(1:37)= (/ 0.1, 250.0, 0.035, 20.0, 5.35, 1E-4, 0.8, 0.1, 0.825,& !PAR(1:9) table 1 parameters
	  0.0, 0.0,& !PAR(10:11) PAR(11) is where AUTO stores the value of the period so we want to avoid it
	  4.5061896805176323384727756760086, 0.29, 0.99301715444958109099407422960427, 0.42746887603328361760208663377664,& !PAR(12:15) initial endpoint
	  4.5061896805176323384727756760086, 0.29, 0.99301715444958109099407422960427, 0.42746887603328361760208663377664,& !PAR(16:19) initial startpoint
	  -0.088087775136734942290496943739116, 0.76997426705126931209497411969253, 0.63196532495883081505689389538914,& !PAR(20:22) strong unstable eigenvector
	  9.160105657877419427176574427416, 0.0,&  !PAR(23:24) strong unstable eigenvalue of startpoint and radius
	  0.11962652483774959679780890038609, 0.0,& !PAR(25:26) weak unstable eigenvalue and radius 
	  0.94446596811088211746145172282912, -0.22232240034899249933270853365139, -0.24198509330832034330998279666784,& !PAR(27:29) weak unstable eigenvector
	  0.0, 0.0,& !PAR(30:31) strong stable eigenvalue and radius
	  0.0, 0.0, 0.0,&  !PAR(32:34) strong stable eigenvector of initial startpoint is not used
	  0.0, 0.0, 0.0/) !PAR(35:37) variables used for releasing startpoint from unstable eigenvector in second step
	  
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
	  FB(1)= U0(1) - (PAR(12) + PAR(24)*PAR(20) + PAR(26)*PAR(27))
	  FB(2)= U0(2) - PAR(13)
	  FB(3)= U0(3) - (PAR(14) + PAR(24)*PAR(21) + PAR(26)*PAR(28))
	  FB(4)= U0(4) - (PAR(15) + PAR(24)*PAR(22) + PAR(26)*PAR(29))
	  
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

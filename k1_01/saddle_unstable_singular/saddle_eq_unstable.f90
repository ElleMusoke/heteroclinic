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
	   x= U(2)
	   y= U(3)
	  
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
	   b=PAR(10)
	      
	   Period= PAR(11)
	   
	   mu= k7/k8
	   alpha= (k1*k5*kn7)/(k3*k8*sqrt(2*k2*k8))
	   epsilon_b= ((k1**2)*k5)/(2*k2*k3*k8)
	   kappa= (sqrt(2*k2*k8))/k5
	   epsilon_sq= (k3*k8)/(k1*k5)
	   zeta= k4/(sqrt(2*k2*k8))
	   delta= k6/k8
	  

       F(1)= Period*(mu - alpha*a - a*b*y) 
       F(2)= Period*(b*x-x**2 + 3*a*b*y - zeta*x + delta)/(epsilon_sq)
	   F(3)= Period*kappa*(x**2 - y - a*b*y)/(epsilon_sq)

      END SUBROUTINE FUNC
!---------------------------------------------------------------------- 
      SUBROUTINE STPNT(NDIM,U,PAR,T) 
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
	  DOUBLE PRECISION rad,h,theta,cs,sn,pi,s1,s2

	  PAR(1:25)= (/ 0.1, 250.0, 0.035, 20.0, 5.35, 1E-5, 0.8, 0.1, 0.825,& !PAR(1:9) table 1 parameters
	  0.32254818477151114657135458312817, 0.0,& !PAR(10) value for parameter B, PAR(11) is where AUTO stores the value of the period so we want to avoid it
	  3.5545054547341554928262530345153, 1.0992690167750540244021312593748, 0.56295959717795151286975553115473, & !PAR(12:14) equilibrium of full system
	  0.84546212014046716836108954944052, -0.33611863044884509943475501922561, -0.41499165012416690799563675007641, & !PAR(15:17) weak unstable eigenvector
	  -0.094008401683558401105574562876035, 0.71037835438813345045561415765963, 0.69751345078766059197251914234036, &  !PAR(18:20) strong unstable eigenvector
	  0.0, 0.0, 0.001,&  !PAR(21:23) arclength, fraction of 2*pi and radius
	  0.0, 0.0 /)
	  
	  PAR(27:28) = (/ 1.0, 1.0 /)
	  
      pi=4*ATAN(1.d0)
      rad=PAR(23)
      theta=2*pi*PAR(22)
      cs=COS(theta)
      sn=SIN(theta)
	  s1=1.0/PAR(27)
	  s2=1.0/PAR(28)
	  
	  PAR(41:43) = PAR(12:14) + rad*(s1*sn*PAR(15:17) + s2*cs*PAR(18:20)) !initialize PAR(41:44) AND U(1:4) with same thing, PAR(41:44) will keep track of the startpoint
	  
	  U(1:3)= PAR(41:43)
	  
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
	  DOUBLE PRECISION rad,h,theta,cs,sn,pi,s1,s2
	  
      pi=4*ATAN(1.d0)
      rad=PAR(23)
      theta=2*pi*PAR(22)
      cs=COS(theta)
      sn=SIN(theta)
	  s1=1.0/PAR(27)
	  s2=1.0/PAR(28)
	  
	  !startpoint BCs
	  FB(1:3)= U0(1:3) - (PAR(12:14) + rad*(s1*sn*PAR(15:17) + s2*cs*PAR(18:20))) !force startpoint to stay near saddle in the direction of the unstable eigenvector
	  
	  !endpoint BCs
	  FB(4:6)= U1(1:3) - PAR(41:43) !let endpoint go wherever
	  
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
	  
      DOUBLE PRECISION FF(3),DFDU(1),DFDP(1)

      CALL FUNC(NDIM,U,ICP,PAR,0,FF,DFDU,DFDP) 
      FI(1)=SQRT(FF(1)**2 + FF(2)**2 + FF(3)**2)*PAR(11)/2 - PAR(24) 
	  

      END SUBROUTINE ICND

      SUBROUTINE PVLS
      END SUBROUTINE PVLS

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
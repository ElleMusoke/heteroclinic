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
	  DOUBLE PRECISION rad,h,theta,cs,sn,pi,s1,s2

	  PAR(1:28)= (/ 0.1, 250.0, 0.035, 20.0, 5.35, 1E-5, 0.8, 0.1, 0.825,& !PAR(1:9) table 1 parameters
	  0.0, 0.0,& !PAR(10:11) PAR(11) is where AUTO stores the value of the period so we want to avoid it
	  3.5543523160658533377274211068911, 0.32252284625797180591497459978325, &
	  1.0993120635759691410092529826354, 0.56304026879390164145746976188598, & !PAR(12:15) equilibrium of full system
	  0.8216310390323403202136272807492, 0.0051178787024727747100739250767541, &
	  -0.35991370425517227805285100800578, -0.44199362948524950116946284058798, & !PAR(16:19) weak unstable eigenvector
	  -0.094368550455394742772297980659396, -0.00045691357544406957825645796388328, &
	  0.71001423820984576087018741045438, 0.69783533118797058479699059809073, &  !PAR(20:23) strong unstable eigenvector
	  0.0, 0.0, 0.001,&  !PAR(24:26) arclength, fraction of 2*pi and radius
	  1.0, 1.0/)
	  
      pi=4*ATAN(1.d0)
      rad=PAR(26)
      theta=2*pi*PAR(25)
      cs=COS(theta)
      sn=SIN(theta)
	  s1=1.0/PAR(27)
	  s2=1.0/PAR(28)
	  
	  PAR(41:44) = PAR(12:15) + rad*(s1*sn*PAR(16:19) + s2*cs*PAR(20:23)) !initialize PAR(41:44) AND U(1:4) with same thing, PAR(41:44) will keep track of the startpoint
	  
	  U(1:4)= PAR(41:44)
	  
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
      rad=PAR(26)
      h=PAR(25)
      theta=2*pi*h
      cs=COS(theta)
      sn=SIN(theta)
	  s1=1.0/PAR(27)
	  s2=1.0/PAR(28)
	  
	  !startpoint BCs
	  FB(1:4)= U0(1:4) - (PAR(12:15) + rad*(s1*sn*PAR(16:19) + s2*cs*PAR(20:23))) !force startpoint to stay near saddle in the direction of the unstable eigenvector
	  
	  !endpoint BCs
	  FB(5:8)= U1(1:4) - PAR(41:44) !let endpoint go wherever
	  
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
	  
      DOUBLE PRECISION FF(4),DFDU(1),DFDP(1)

       CALL FUNC(NDIM,U,ICP,PAR,0,FF,DFDU,DFDP) 
       FI(1)=SQRT(FF(1)**2 + FF(2)**2 + FF(3)**2 + FF(4)**2) - PAR(24) 
	  

      END SUBROUTINE ICND

      SUBROUTINE PVLS
      END SUBROUTINE PVLS

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
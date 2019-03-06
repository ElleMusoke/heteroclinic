      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION a, b, x, y, Period, alpha, mu, kappa, zeta, epsilon_sq, epsilon_b, delta, k1, k2, k3, k4, k5, k6,k7,k8,kn7
	  DOUBLE PRECISION a2, x2, y2, Period2
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
	   
	   a2= U(4)
	   x2= U(5)
	   y2= U(6)
	   
	   Period2=PAR(12)
	   
       F(4)= Period2*(mu - alpha*a2 - a2*b*y2) 
       F(5)= Period2*(b*x2-x2**2 + 3*a2*b*y2 - zeta*x2 + delta)/(epsilon_sq)
	   F(6)= Period2*kappa*(x2**2 - y2 - a2*b*y2)/(epsilon_sq)
	   
      END SUBROUTINE FUNC
!---------------------------------------------------------------------- 
      SUBROUTINE STPNT(NDIM,U,PAR,T) 
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION Period, alpha, mu, kappa, zeta, epsilon, epsilon_b
	  DOUBLE PRECISION rad,theta,cs,sn,rad2,theta2,cs2,sn2, pi

	  PAR(1:12)= (/ 0.1, 250.0, 0.035, 20.0, 5.35, 1E-5, 0.8, 0.1, 0.825, 0.7, 0.0, 0.0 /) !PAR(1:9) table parameters, PAR(10) slow variable b, PAR(11:12) Period, Period2
	  
	  !TOP ORBIT STUFF
	  PAR(40:42)= (/ 10.337919540910958466473353904001, 0.17402345878397065077820017696596, 0.0036768049061260305334243335255799 /) 
	  PAR(43:45)= (/ 10.337919540910958466473353904001, 0.17402345878397065077820017696596, 0.0036768049061260305334243335255799 /) !PAR(40:42), PAR(43:45) startpoint and endpoint for top orbit
	  PAR(48:50)= (/ 0.999995, -0.00323309, -0.000449168 /) !weak stable eigenvector of endpoint
	  PAR(51:53)= (/ 0.00999989, -0.566901, 0.823725 /) !strong stable eigenvector of endpoint
	  PAR(54:55)= (/ -0.0905496, -596.192 /) !weak and strong eigenvalues
	  
	  PAR(25:26)= (/ 0.0, 0.00109/) !theta and radius
	  
      pi=4*ATAN(1.d0)
      rad=PAR(26)
      theta=2*pi*PAR(25)
      cs=COS(theta)
      sn=SIN(theta)
	  
	  PAR(30:32) = PAR(43:45) + rad*(sn*PAR(48:50) + cs*PAR(51:53)) !initialize PAR(30:32) AND U(1:3) with same thing, PAR(30:32) will keep track of the startpoint
	  
	  U(1:3)= PAR(43:45) + rad*(sn*PAR(48:50) + cs*PAR(51:53))
	  
	  !BOTTOM ORBIT STUFF
	  PAR(56:61)= (/ 0.940272, 1.492271, 1.342954, 0.940272, 1.492271, 1.342954 /) !PAR(56:59), PAR(60:63) startpoint and endpoint for bottom orbit
	  PAR(64:66)= (/ - 0.0537135361608340861758282121244, 0.47335191629441716850755227165349, 0.87214305967004426414910112371712 /) !real part of unstable complex conjugate eigenvectors of saddle equilibrium at B=0.7
	  PAR(67:69)= (/ -0.099631805485786386764404684114047, -0.049928009959519587263507064053846, 0.0 /) !imaginary part of unstable complex conjugate eigenvectors of stable equilibrium at B=0.7
	  PAR(70:71)= (/ 0.0, 0.0 /) !real and imaginary radius
	  PAR(103:105)= (/ 0.940272, 1.492271, 1.342954 /) !point on critical manifold corresponding to bottom endpoint B-coordinate
	  PAR(106:107)= (/ 1.3753860940993691676233022299789, -4.4640865481402179732661936757437 /) !real (PAR(106)) and imaginary (PAR(107)) parts of eigenvalues
	  PAR(83)=0.0001
	  
      rad2=PAR(83)
      theta2=2*pi*PAR(84)
      cs2=COS(theta2)
      sn2=SIN(theta2)
	  
	  PAR(80:82) = PAR(56:58) + rad2*(sn2*PAR(64:66) + cs2*PAR(67:69)) !initialize PAR(80:82) AND U(4:6) with same thing, PAR(80:82) will keep track of the endpoint
	  U(4:6)= PAR(80:82)
	  
	  !Lin vector
	  PAR(74)= (U(2)-U(5))/SQRT((U(2)-U(5))**2+(U(3)-U(6))**2) !X
	  PAR(75)= (U(3)-U(6))/SQRT((U(2)-U(5))**2+(U(3)-U(6))**2) !Y
	  
	  !vector orthogonal to Lin vector
	  PAR(76)= -(U(3)-U(6))/SQRT((U(2)-U(5))**2+(U(3)-U(6))**2) !Lin vector -Y
	  PAR(77)= (U(2)-U(5))/SQRT((U(2)-U(5))**2+(U(3)-U(6))**2) !Lin vector X
	  
	  !Distance between endpt of bottom orbit and stpt of top orbit along Lin vector (squared)
	  PAR(100)= DOT_PRODUCT(PAR(74:75), (/ (U(2)-U(5)), (U(3)-U(6)) /)) 
	  
	  !used to keep A coordinates closed in fold steps
!	  PAR(102)= PAR(60)-PAR(44)
	  
      END SUBROUTINE STPNT
!---------------------------------------------------------------------- 
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
	  
	  DOUBLE PRECISION a, b, x, y, alpha, mu, kappa, zeta, epsilon_sq, epsilon_b, delta, k1, k2, k3, k4, k5, k6,k7,k8,kn7
	  DOUBLE PRECISION rad,theta,cs,sn,rad2,theta2,cs2,sn2,pi, M(3,3), K(3,3), a2, x2, y2, v1, v2, v3
	  INTEGER D
	  
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
   
   mu= k7/k8
   alpha= (k1*k5*kn7)/(k3*k8*sqrt(2*k2*k8))
   epsilon_b= ((k1**2)*k5)/(2*k2*k3*k8)
   kappa= (sqrt(2*k2*k8))/k5
   epsilon_sq= (k3*k8)/(k1*k5)
   zeta= k4/(sqrt(2*k2*k8))
   delta= k6/k8
	  
	  
	  !TOP ORBIT
	  
      pi=4*ATAN(1.d0)
      rad=PAR(26)
      theta=2*pi*PAR(25)
      cs=COS(theta)
      sn=SIN(theta)
	  
	  !endpoint BCs
	  FB(1:3)= U0(1:3) - (PAR(43:45) + rad*(sn*PAR(48:50) + cs*PAR(51:53)))
	  
	  !startpoint BCs
	  FB(4:6)= U1(1:3) - (PAR(30:32))
	  
	  
	  
	  !BOTTOM ORBIT
	  
      rad2=PAR(83)
      theta2=2*pi*PAR(84)
      cs2=COS(theta2)
      sn2=SIN(theta2)
	  
	  !startpoint BCs
	  FB(7:9)= U0(4:6) - (PAR(56:58) + rad2*(sn2*PAR(64:66) + cs2*PAR(67:69)))

	  !endpoint BCs
	  FB(10:12)= U1(4:6) - (PAR(80:82))
	  	    
	  !making sure that the endpoints stay on the span of the Lin vector by requiring that their difference is normal to the two normal vectors
	  FB(13)= DOT_PRODUCT( (/ U1(2)-U1(5), U1(3)-U1(6) /), PAR(76:77) )
!	  
!	  !distance between endpoints to be closed in last step
	  FB(14)= DOT_PRODUCT(PAR(74:75), (/ (U1(2)-U1(5)), (U1(3)-U1(6)) /)) - PAR(100)
	  
	  !calculating new stable and unstable eigenspaces for each new B
	  
	  !TOP ORBIT
	  
	  !new equilibrium
	  
      a= PAR(43)
      x= PAR(44)
      y= PAR(45)
	  
      FB(15)= (mu - alpha*a - a*b*y)
      FB(16)= (b*x-x**2 + 3*a*b*y - zeta*x + delta)/(epsilon_sq)
   	  FB(17)= kappa*(x**2 - y - a*b*y)/(epsilon_sq)
	  
	  !endpoint strong stable eigenvector
	  
	  D = 3
	  CALL JACOB(D,PAR,(/ PAR(43:45) /),M)
	  
	  FB(18:20)= MATMUL(M, PAR(51:53)) - PAR(55)*PAR(51:53)
	  FB(21)= DOT_PRODUCT(PAR(51:53),PAR(51:53))-1.0
	  
	  !endpoint weak stable eigenvector
	  
	  FB(22:24)= MATMUL(M, PAR(48:50)) - PAR(54)*PAR(48:50)
	  FB(25)= DOT_PRODUCT(PAR(48:50),PAR(48:50))-1.0
	  
	  !BOTTOM ORBIT
	  
	  !new equilibrium
	  
      a2= PAR(56)
      x2= PAR(57)
      y2= PAR(58)
	  
      FB(26)= (mu - alpha*a2 - a2*b*y2)
      FB(27)= (b*x2-x2**2 + 3*a2*b*y2 - zeta*x2 + delta)/(epsilon_sq)
   	  FB(28)= kappa*(x2**2 - y2 - a2*b*y2)/(epsilon_sq)
	  
	  CALL JACOB(D,PAR,(/ PAR(56:58) /),K)
	  
	  FB(29)= SQRT(PAR(64)**2 + PAR(67)**2 + PAR(65)**2 + PAR(68)**2 + PAR(66)**2 + PAR(69)**2) - 1.0
	  
	  FB(30:32)= MATMUL(K, PAR(64:66)) - (PAR(106)*PAR(64:66)-PAR(107)*PAR(67:69))
	  FB(33:35)= MATMUL(K, PAR(67:69)) - (PAR(107)*PAR(64:66)+PAR(106)*PAR(67:69))
	  
	  FB(36)= SQRT(PAR(64)**2 + PAR(65)**2 + PAR(66)**2) - 0.99377094803436699104907010676142
	  
	  IF(NBC==36) RETURN	  
	  !making sure the orthogonal vector is unit
	  FB(37)= DOT_PRODUCT(PAR(76:77), PAR(76:77))-1.0
	  
	  !Lin vector
	  FB(38)= (U1(2)-U1(5))/SQRT((U1(2)-U1(5))**2+(U1(3)-U1(6))**2)-PAR(74)
	  FB(39)= (U1(3)-U1(6))/SQRT((U1(2)-U1(5))**2+(U1(3)-U1(6))**2)-PAR(75)
	  
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
	  b= PAR(10)
	  x= U(2)
	  y= U(3)
	
	  M(1,1)= - alpha - b*y
	  M(1,2)= 0
	  M(1,3)= -a*b
	  
	  M(2,1)= (3*b*y)/epsilon_sq
	  M(2,2)= -(2*x - b + zeta)/epsilon_sq
	  M(2,3)= (3*a*b)/epsilon_sq
	  
	  M(3,1)= -(b*kappa*y)/epsilon_sq
	  M(3,2)= (2*kappa*x)/epsilon_sq
	  M(3,3)= -(kappa*(a*b + 1))/epsilon_sq
	  
  	END SUBROUTINE JACOB
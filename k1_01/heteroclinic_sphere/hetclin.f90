      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION a, b, x, y, Period, alpha, mu, kappa, zeta, epsilon_sq, epsilon_b, delta, k1, k2, k3, k4, k5, k6,k7,k8,kn7
	  DOUBLE PRECISION a2, b2, x2, y2, Period2
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
	   
	   a2= U(5)
	   b2= U(6)
	   x2= U(7)
	   y2= U(8)
	   
	   Period2=PAR(12)
	   
       F(5)= Period2*(mu - alpha*a2 - a2*b2*y2) 
       F(6)= Period2*epsilon_b*(1 - b2*x2 - a2*b2*y2)
       F(7)= Period2*(b2*x2-x2**2 + 3*a2*b2*y2 - zeta*x2 + delta)/(epsilon_sq)
	   F(8)= Period2*kappa*(x2**2 - y2 - a2*b2*y2)/(epsilon_sq)
	   
      END SUBROUTINE FUNC
!---------------------------------------------------------------------- 
      SUBROUTINE STPNT(NDIM,U,PAR,T) 
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION Period, alpha, mu, kappa, zeta, epsilon, epsilon_b

	  PAR(1:12)= (/ 0.1, 250.0, 0.035, 20.0, 5.35, 1E-5, 0.8, 0.1, 0.825, 0.0, 0.0, 0.0 /) !PAR(1:9) table parameters, PAR(10) placeholder, PAR(11:12) Period, Period2
	  
	  !top orbit specific stuff
	  PAR(40:47)= (/ 10.625702, 0.95, 0.0197257, 0.0000350720, 10.625702, 0.95, 0.0197257, 0.0000350720 /) !PAR(40:43), PAR(44:47) startpoint and endpoint for top orbit
	  PAR(48:50)= (/ 1.0, -0.000234108, -0.00000383610 /) !weak stable eigenvector of endpoint
	  PAR(51:53)= (/ 0.0104793, -0.583442, 0.812087 /) !strong stable eigenvector of endpoint
	  PAR(54:55)= (/ 0.0, 0.0 /) !weak and strong radius
	  
	  U(1)= PAR(40)
	  U(2)= PAR(41)
	  U(3)= PAR(42)
	  U(4)= PAR(43)
	  
	  !bottom orbit specific stuff
	  PAR(56:63)= (/ 0.940272, 0.7, 1.492271, 1.342954, 0.940272, 0.7, 1.492271, 1.342954 /) !PAR(56:59), PAR(60:63) startpoint and endpoint for top orbit
	  PAR(64:66)= (/ -0.0537135, 0.473352, 0.872143 /) !real part of unstable complex conjugate eigenvectors of saddle equilibrium at B=0.7
	  PAR(67:69)= (/ 0.0996318, 0.0499280, 0.0 /) !imaginary part of unstable complex conjugate eigenvectors of stable equilibrium at B=0.7
	  PAR(70:71)= (/ 0.0, 0.0 /) !real and imaginary radius
	  PAR(103:106)= (/ 0.940272, 0.7, 1.492271, 1.342954 /) !point on critical manifold corresponding to bottom endpoint B-coordinate
	  PAR(107)= 0.0 !distance of bottom orbit endpoint from critical manifold
	  
	  U(5)= PAR(56)
	  U(6)= PAR(57)
	  U(7)= PAR(58)
	  U(8)= PAR(59)
	  
	  !Lin vector
	  PAR(73)= (U(2)-U(6))/SQRT((U(2)-U(6))**2+(U(3)-U(7))**2+(U(4)-U(8))**2)
	  PAR(74)= (U(3)-U(7))/SQRT((U(2)-U(6))**2+(U(3)-U(7))**2+(U(4)-U(8))**2)
	  PAR(75)= (U(4)-U(8))/SQRT((U(2)-U(6))**2+(U(3)-U(7))**2+(U(4)-U(8))**2)
	  
	  !first vector orthogonal to Lin vector
	  PAR(77)= 0
	  PAR(78)= (U(4)-U(8))/SQRT((U(3)-U(7))**2+(U(4)-U(8))**2)
	  PAR(79)= -(U(3)-U(7))/SQRT((U(3)-U(7))**2+(U(4)-U(8))**2)
	  
	  !second vector orthogonal to Lin vector
	  PAR(80)= (U(4)-U(8))/SQRT((U(2)-U(6))**2+(U(4)-U(8))**2)
	  PAR(81)= 0
	  PAR(82)= -(U(2)-U(6))/SQRT((U(2)-U(6))**2+(U(4)-U(8))**2)
	  
	  !Distance between endpt of bottom orbit and stpt of top orbit along Lin vector
	  PAR(100)= DOT_PRODUCT(PAR(73:75), (/ (U(2)-U(6)), (U(3)-U(7)), (U(4)-U(8)) /)) 
	  
	  !top orbit startpoint distance from point on C^b corresponding to bottom endpoint
	  PAR(102)= SQRT((PAR(44)-PAR(103))**2 + (PAR(45)-PAR(104))**2 + (PAR(46)-PAR(105))**2 + (PAR(47)-PAR(106))**2)
	  PAR(108)= 0.0
	  
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
	  INTEGER D
	  
	  !top orbit
	  !endpoint BCs
	  FB(1)= U0(1) - (PAR(40) + PAR(54)*PAR(48) + PAR(55)*PAR(51))
	  FB(2)= U0(2) - PAR(41)
	  FB(3)= U0(3) - (PAR(42) + PAR(54)*PAR(49) + PAR(55)*PAR(52))
	  FB(4)= U0(4) - (PAR(43) + PAR(54)*PAR(50) + PAR(55)*PAR(53))
	  !startpoint BCs
	  FB(5)= U1(1) - (PAR(44))
	  FB(6)= U1(2) - (PAR(45))
	  FB(7)= U1(3) - (PAR(46))
	  FB(8)= U1(4) - (PAR(47))
	  
	  !bottom orbit
	  !startpoint BCs
	  FB(9)= U0(5) - (PAR(56) + PAR(70)*PAR(64) + PAR(71)*PAR(67))
	  FB(10)= U0(6) - PAR(57)
	  FB(11)= U0(7) - (PAR(58) + PAR(70)*PAR(65) + PAR(71)*PAR(68))
	  FB(12)= U0(8) - (PAR(59) + PAR(70)*PAR(66) + PAR(71)*PAR(69))
	  !endpoint BCs
	  FB(13)= U1(5) - (PAR(60))
	  FB(14)= U1(6) - (PAR(61))
	  FB(15)= U1(7) - (PAR(62))
	  FB(16)= U1(8) - (PAR(63))
	  	    
	  !making sure that the endpoints stay on the span of the Lin vector by requiring that their difference is normal to the two normal vectors
	  FB(17)= DOT_PRODUCT( (/ (U1(2)-U1(6)), (U1(3)-U1(7)), (U1(4)-U1(8)) /), PAR(77:79) )
	  FB(18)= DOT_PRODUCT( (/ (U1(2)-U1(6)), (U1(3)-U1(7)), (U1(4)-U1(8)) /), PAR(80:82) )
	  
	  !distance between endpoints to be closed in last step
	  FB(19)= DOT_PRODUCT(PAR(73:75), (/ (U1(2)-U1(6)), (U1(3)-U1(7)), (U1(4)-U1(8)) /)) - PAR(100)
	  
	  !distance of bottom endpoint from critical
	  FB(20)= SQRT((PAR(60)-PAR(103))**2 + (PAR(61)-PAR(104))**2 + (PAR(62)-PAR(105))**2 + (PAR(63)-PAR(106))**2) - PAR(107)
	  
	  !distance of top startpoint to point on critical manifold corresponding to bottom endpoint
	  FB(21)= SQRT((PAR(44)-PAR(103))**2 + (PAR(45)-PAR(104))**2 + (PAR(46)-PAR(105))**2 + (PAR(47)-PAR(106))**2) - PAR(102)
	  
	  IF(NBC==21) RETURN
	  
	  !Lin vector
	  FB(22)= (U1(2)-U1(6))/SQRT((U1(2)-U1(6))**2+(U1(3)-U1(7))**2+(U1(4)-U1(8))**2)-PAR(73)
	  FB(23)= (U1(3)-U1(7))/SQRT((U1(2)-U1(6))**2+(U1(3)-U1(7))**2+(U1(4)-U1(8))**2)-PAR(74)
	  FB(24)= (U1(4)-U1(8))/SQRT((U1(2)-U1(6))**2+(U1(3)-U1(7))**2+(U1(4)-U1(8))**2)-PAR(75)
	  
	  FB(25)= DOT_PRODUCT(PAR(77:79), PAR(77:79))-1.0
	  FB(26)= DOT_PRODUCT(PAR(80:82), PAR(80:82))-1.0
	  
	  IF(NBC==26) RETURN
	  
	  a=par(103)
	  b=par(104)
	  x=par(105)
	  y=par(106)
	  
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
	  
      FB(27)= (mu - alpha*a - a*b*y) 
      FB(28)= (b*x-x**2 + 3*a*b*y - zeta*x + delta)/(epsilon_sq)
      FB(29)= kappa*(x**2 - y - a*b*y)/(epsilon_sq)
	  
	  !used for freeing up the bottom B endpoint from the defining sphere of the 3D leaf
	  FB(30)= (PAR(61)-PAR(104))
	  
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
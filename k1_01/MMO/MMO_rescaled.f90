!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   olsen.f90 - compute unstable manifolds of branch A = 8
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!--------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the algebraic or ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   Equation or ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  !set variable type for parameters and state variables
  DOUBLE PRECISION a, b, x, y, Period, alpha, mu, kappa, zeta, epsilon_sq, epsilon_b, delta, k1, k2, k3, k4, k5, k6,k7,k8,kn7

  !define state variables
  
  a = U(1)
  b = U(2)
  x = U(3)
  y = U(4)
  
  !define parameters
  
  k1 = par(1)
  k2 = par(2)
  k3 = par(3)
  k4 = par(4)
  k5 = par(5)
  k6 = par(6)
  k7 = par(7)
  kn7 =par(8)
  k8 = par(9)
    
  Period= PAR(11)
  
  mu= k7/k8
  alpha= (k1*k5*kn7)/(k3*k8*sqrt(2*k2*k8))
  epsilon_b= ((k1**2)*k5)/(2*k2*k3*k8)
  kappa= (sqrt(2*k2*k8))/k5
  epsilon_sq= (k3*k8)/(k1*k5)
  zeta= k4/(sqrt(2*k2*k8))
  delta= k6/k8
  
  F(1)= (mu - alpha*a - a*b*y) 
  F(2)= epsilon_b*(1 - b*x - a*b*y)
  F(3)= (b*x-x**2 + 3*a*b*y - zeta*x + delta)/(epsilon_sq)
  F(4)= kappa*(x**2 - y - a*b*y)/(epsilon_sq)
  
  
  
END SUBROUTINE FUNC

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)
!--------- -----

! Input arguments :
!      NDIM   :   Dimension of the algebraic or ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values

! Note : For time- or space-dependent solutions this subroutine has
!        the scalar input parameter T contains the varying time or space
!        variable value.
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T
  
  DOUBLE PRECISION rad

! Initialize the equation parameters
  PAR(1:10) = (/ 0.1, 250.0, 0.035, 20.0, 0.1, 1e-5, 0.8, 0.1, 0.825, 1.0 /)
  

! Equilibrium of full system
  U(1) = 24.202251400006826732601648738031
  U(2) = 0.056439545712523865804148014425528
  U(3) = 1.268114854230984828407632827555
  U(4) = 0.67968710980537210679562952825643
  
END SUBROUTINE STPNT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!--------- ----

! Boundary Conditions

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      PAR    :   Equation parameters
!      ICP    :   Array indicating the free parameter(s)
!      NBC    :   Number of boundary conditions
!      U0     :   State variable values at the left boundary
!      U1     :   State variable values at the right boundary

! Values to be returned :
!      FB     :   The values of the boundary condition functions 

! Normally unused Jacobian arguments : IJAC, DBC (see manual)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
  DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

END SUBROUTINE BCND

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!--------- ----

! Integral Conditions

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      PAR    :   Equation parameters
!      ICP    :   Array indicating the free parameter(s)
!      NINT   :   Number of integral conditions
!      U      :   Value of the vector function U at `time' t

! The following input arguments, which are normally not needed,
! correspond to the preceding point on the solution branch
!      UOLD   :   The state vector at 'time' t
!      UDOT   :   Derivative of UOLD with respect to arclength
!      UPOLD  :   Derivative of UOLD with respect to `time'

! Normally unused Jacobian arguments : IJAC, DINT

! Values to be returned :
!      FI     :   The value of the vector integrand 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
  DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
  

END SUBROUTINE ICND

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
!--------- ----
!
! Defines the objective function for algebraic optimization problems
!
! Supplied variables :
!      NDIM   :   Dimension of the state equation
!      U      :   The state vector
!      ICP    :   Indices of the control parameters
!      PAR    :   The vector of control parameters
!
! Values to be returned :
!      FS      :   The value of the objective function
!
! Normally unused Jacobian argument : IJAC, DFDP

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: FS
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM),DFDP(*)

!X FS=

END SUBROUTINE FOPT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE PVLS(NDIM,U,PAR)
!--------- ----

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

!---------------------------------------------------------------------- 
! NOTE : 
! Parameters set in this subroutine should be considered as ``solution 
! measures'' and be used for output purposes only.
! 
! They should never be used as `true'' continuation parameters. 
!
! They may, however, be added as ``over-specified parameters'' in the 
! parameter list associated with the AUTO-Constant NICP, in order to 
! print their values on the screen and in the ``p.xxx file.
!
! They may also appear in the list associated with AUTO-Constant NUZR.
!
!---------------------------------------------------------------------- 
! For algebraic problems the argument U is, as usual, the state vector.
! For differential equations the argument U represents the approximate 
! solution on the entire interval [0,1]. In this case its values must 
! be accessed indirectly by calls to GETP, as illustrated below.
!---------------------------------------------------------------------- 
!
! Set PAR(2) equal to the L2-norm of U(1)
!X PAR(2)=GETP('NRM',1,U)
!
! Set PAR(3) equal to the minimum of U(2)
!X PAR(3)=GETP('MIN',2,U)
!
! Set PAR(4) equal to the value of U(2) at the left boundary.
!X PAR(4)=GETP('BV0',2,U)
!
! Set PAR(5) equal to the pseudo-arclength step size used.
!X PAR(5)=GETP('STP',1,U)
!
!---------------------------------------------------------------------- 
! The first argument of GETP may be one of the following:
!        'NRM' (L2-norm),     'MAX' (maximum),
!        'INT' (integral),    'BV0 (left boundary value),
!        'MIN' (minimum),     'BV1' (right boundary value).
!
! Also available are
!   'STP' (Pseudo-arclength step size used).
!   'FLD' (`Fold function', which vanishes at folds).
!   'BIF' (`Bifurcation function', which vanishes at singular points).
!   'HBF' (`Hopf function'; which vanishes at Hopf points).
!   'SPB' ( Function which vanishes at secondary periodic bifurcations).
!---------------------------------------------------------------------- 


END SUBROUTINE PVLS
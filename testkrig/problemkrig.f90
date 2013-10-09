program problemKriging

 use dimkrig

  implicit none
  !
  !     include the Ipopt return codes
  !
  include 'IpReturnCodes.inc'
  include 'mpif.h'
  !
  !     Size of the problem (number of variables and equality constraints)
  !
  integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
  parameter  (N = 3, M = 8, NELE_JAC = 24, NELE_HESS = 6)
  parameter  (IDX_STY = 1 )
  !
  !     Space for multipliers and constraints
  !
  double precision LAM(M)
  double precision G(M)
  !
  !     Vector of variables
  !
  double precision X(N)
  !
  !     Vector of lower and upper bounds
  !
  double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
  double precision G_L(M), G_U(M)
  !
  !     Private data for evaluation routines
  !     This could be used to pass double precision and integer arrays untouched
  !     to the evaluation subroutines EVAL_*
  !
  double precision DAT(2000)
  integer IDAT(2)
  !
  !     Place for storing the Ipopt Problem Handle
  !
  integer*8 IPROBLEM
  integer*8 IPCREATE
  !
  integer IERR
  integer IPSOLVE, IPADDSTROPTION
  integer IPADDNUMOPTION, IPADDINTOPTION
  integer IPOPENOUTPUTFILE
  !
  double precision F,Fs,sigmax(N)
  integer i,kprob

  double precision  infbound
  parameter        (infbound = 1.d+20)
  !
  !     The following are the Fortran routines for computing the model
  !     functions and their derivatives - their code can be found further
  !     down in this file.
  !
  external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS, ITER_CB

  integer::probtype

  call MPI_START

  probtype=1

  kprob=1

  sigmax(1)=0.05
  sigmax(2)=0.05
  sigmax(3)=0.05


  do i=i,n
     dat(i)=sigmax(i)
  end do

  IDAT(1)=kprob
  IDAT(2)=0
  IDAT(3)=probtype

  !
  !     Set initial point and bounds:
  !
  do i=1,N
     X(i) = 1.0
     X_L(i) = 0.010
     X_U(i) = infbound
  end do

  !
  !     Set bounds for the constraints
  !
  do i=1,M
     G_L(i)=-infbound
     G_U(i)=0.d0
  end do

  !
  !     First create a handle for the Ipopt problem (and read the options
  !     file)
  !

  IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
  if (IPROBLEM.eq.0) then
     write(*,*) 'Error creating an Ipopt Problem handle.'
     call stop_all
  endif
  !
  !     Open an output file
  !
  IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
  if (IERR.ne.0 ) then
     write(*,*) 'Error opening the Ipopt output file.'
     goto 9000
  endif
  !

  !!
  !!     Set a callback function to give you control once per iteration.
  !!     You can use it if you want to generate some output, or to stop
  !!     the optimization early.
  !!
  call IPSETCALLBACK(IPROBLEM, ITER_CB)

  !
  !     Call optimization routine
  !

  if (id_proc.eq.0) then
     IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 0)
     if (IERR.ne.0 ) goto 9990
  else
     IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 0)
     if (IERR.ne.0 ) goto 9990
  end if

  IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)

  !
  !     Output:
  !
  if (id_proc.eq.0) then

     if( IERR.eq.IP_SOLVE_SUCCEEDED .or. IERR.eq.5) then
        write(*,*)
        write(*,*) 'The solution was found.'
        write(*,*)
     else
        write(*,*)
        write(*,*) 'An error occoured.'
        write(*,*) 'The error code is ',IERR
        write(*,*)
     endif

     write(*,*) 'The final value of the objective function is ',F
     write(*,*)
     write(*,*) 'The optimal values of X are:'
     write(*,*)
     do i = 1, N
        write(*,*) 'X  (',i,') = ',X(i)
     enddo
     write(*,*)
     write(*,*) 'The multipliers for the equality constraints are:'
     write(*,*)
     do i = 1, M
        write(*,*) 'LAM(',i,') = ',LAM(i)
     enddo
     write(*,*)
     write(*,*) 'Weight and its variance:',DAT(N+2),DAT(N+3)

  end if
  !
9000 continue
  !
  !     Clean up
  !
  call IPFREE(IPROBLEM)

  call stop_all
  !
9990 continue
  write(*,*) 'Error setting an option'
  goto 9000

end program problemKriging
!
! =============================================================================
!
!                    Computation of objective function
!
! =============================================================================
!

subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
  use dimkrig

  implicit none
  integer N, NEW_X,I,probtype
  double precision F, X(N),sigmax(N),fmeantmp,fvartmp,fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  double precision DAT(*)
  integer IDAT(*),kprob,NMC
  integer IERR
  double precision fmin,fmax,gradmin(N-1),gradmax(N-1),gtol,low(N-1),up(N-1),Xsave(N)
  double precision  rho, L, sigmay, pi, p, E, Fs 

  !      if (id_proc.eq.0) print *,'Calculate Objective',X

  kprob=IDAT(1)
  probtype=IDAT(3)

  do i=1,N
     sigmax(i)=DAT(i)
  end do

  !---- MEAN and VARIANCE OF worst OBJECTIVE FUNCTION
!  call  Krigingestimate(ndimin,ndimint,xavgin,xstdin,fctin,fctindxin,nptsin,statin,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

  call  Krigingestimate(N,N,X,sigmax,11,0,60,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

  if (IDAT(2).eq.1) then ! Deterministic with PC
     fvartmp=0.0d0
     fvarprimetmp=0.0d0
  end if


  dat(N+2)=fmeantmp
  dat(N+3)=fvartmp

!  print*,fmeantmp,fvartmp

  !---- COMBINED OBJECTIVE FUNCTION

  F=fmeantmp+fvartmp

  IERR = 0
  return

end subroutine EV_F

!
! =============================================================================
!
!                     Computation of constraints
!
! =============================================================================
!
subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
  use dimkrig

  implicit none
  integer N, NEW_X, M,probtype
  double precision G(M), X(N), sigmax(N), cmean(M), cstd(M), fmeantmp, fvartmp
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n),dc(M,N)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  integer IDAT(*),kprob,NMC
  integer IERR, i, j, cnt
  double precision fmin,fmax,gradmin(N-1),gradmax(N-1),gtol,low(N-1),up(N-1),Xsave(N)
  double precision  rho, L, sigmay, pi, p, E, Fs 


  kprob=IDAT(1)
  probtype=IDAT(3)

  do i=1,N
     sigmax(i)=DAT(i)
  end do

  do i=1,M

     !---- MEAN OF INEQUALITY CONSTRAINT i

!  call  Krigingestimate(ndimin,ndimint,xavgin,xstdin,fctin,fctindxin,nptsin,statin,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

  call  Krigingestimate(N,N,X,sigmax,11,i,60,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

  if (IDAT(2).eq.1) then ! Deterministic with PC
     fvartmp=0.0d0
     fvarprimetmp=0.0d0
  end if

     G(i)=fmeantmp+fvartmp

  end do

  !Just printing

  if (id_proc.eq.0) then
     print*,''
     write(*,'(4x,a)') '>>Normalized Constraint Values:'
     do i=1,8
        write(*,'(E13.2)'),g(i)
     end do
     print*,''
  end if

  IERR = 0
  return
end subroutine EV_G

!
! =============================================================================
!
!                Computation of gradient of objective function
!
! =============================================================================
!
subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
  use dimkrig

  implicit none
  integer N, NEW_X,i,probtype
  double precision GRAD(N), X(N), sigmax(N), fmeantmp, fvartmp
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  integer IDAT(*),kprob,NMC
  integer IERR
  double precision  rho, L, sigmay, pi, p, E, Fs 

  kprob=IDAT(1)
  probtype=IDAT(3)

  do i=1,N
     sigmax(i)=DAT(i)
  end do

!  call  Krigingestimate(ndimin,ndimint,xavgin,xstdin,fctin,fctindxin,nptsin,statin,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

  call  Krigingestimate(N,N,X,sigmax,11,0,60,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

  !---- GRADIENT OF OBJECTIVE FUNCTION

  do i=1,n


     if (IDAT(2).eq.1) then ! Deterministic with PC
        fvartmp=0.0d0
        fvarprimetmp=0.0d0
     end if


     GRAD(i)=fmeanprimetmp(i)+fvarprimetmp(i)

  end do

  IERR = 0
  return
end subroutine EV_GRAD_F

!
! =============================================================================
!
!                Computation of Jacobian of constraints
!
! =============================================================================
!
subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,IDAT, DAT, IERR)
  use dimkrig

  implicit none
  integer TASK, N, NEW_X, M, NZ
  double precision X(N), A(NZ),dc(M,N), sigmax(N), fmeantmp, fvartmp
  integer ACON(NZ), AVAR(NZ), I, J, K, cnt, NMC
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  double precision  rho, L, sigmay, pi, p, E, Fs
  integer IDAT(*)
  integer IERR, kprob,probtype


  if( TASK.eq.0 ) then 
     !
     !     structure of Jacobian:
     !

     ACON(1) = 1
     AVAR(1) = 1

     ACON(2) = 1
     AVAR(2) = 2

     ACON(3) = 1
     AVAR(3) = 3



     ACON(4) = 2
     AVAR(4) = 1

     ACON(5) = 2
     AVAR(5) = 2

     ACON(6) = 2
     AVAR(6) = 3



     ACON(7) = 3
     AVAR(7) = 1

     ACON(8) = 3
     AVAR(8) = 2

     ACON(9) = 3
     AVAR(9) = 3



     ACON(10) = 4
     AVAR(10) = 1

     ACON(11) = 4
     AVAR(11) = 2

     ACON(12) = 4
     AVAR(12) = 3



     ACON(13) = 5
     AVAR(13) = 1

     ACON(14) = 5
     AVAR(14) = 2

     ACON(15) = 5
     AVAR(15) = 3



     ACON(16) = 6
     AVAR(16) = 1

     ACON(17) = 6
     AVAR(17) = 2

     ACON(18) = 6
     AVAR(18) = 3



     ACON(19) = 7
     AVAR(19) = 1

     ACON(20) = 7
     AVAR(20) = 2

     ACON(21) = 7
     AVAR(21) = 3

     ACON(22) = 8
     AVAR(22) = 1

     ACON(23) = 8
     AVAR(23) = 2

     ACON(24) = 8
     AVAR(24) = 3

  else


     !---- TOTAL GRADIENT OF CONSTRAINTS 

     kprob=IDAT(1)
     probtype=IDAT(3)

     do i=1,N
        sigmax(i)=DAT(i)
     end do


     dc(:,:)=0.0

     do i=1,M

        !---- MEAN OF INEQUALITY CONSTRAINT i

!  call  Krigingestimate(ndimin,ndimint,xavgin,xstdin,fctin,fctindxin,nptsin,statin,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

  call  Krigingestimate(N,N,X,sigmax,11,i,60,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)
        
        if (IDAT(2).eq.1) then ! Deterministic with PC
           fvartmp=0.0d0
           fvarprimetmp=0.0d0
        end if
        kprob=1.0
        do j=1,N
           dc(i,j)=fmeanprimetmp(j)
           if (fvartmp.ne.0.0) then
              dc(i,j)=dc(i,j)+kprob*fvarprimetmp(j)/(2.0*sqrt(fvartmp))
           endif
        end do

     end do

     A(1)=dc(1,1)
     A(2)=dc(1,2)
     A(3)=dc(1,3)

     A(4)=dc(2,1)
     A(5)=dc(2,2)
     A(6)=dc(2,3)

     A(7)=dc(3,1)
     A(8)=dc(3,2)
     A(9)=dc(3,3)

     A(10)=dc(4,1)
     A(11)=dc(4,2)
     A(12)=dc(4,3)

     A(13)=dc(5,1)
     A(14)=dc(5,2)
     A(15)=dc(5,3)

     A(16)=dc(6,1)
     A(17)=dc(6,2)
     A(18)=dc(6,3)

     A(19)=dc(7,1)
     A(20)=dc(7,2)
     A(21)=dc(7,3)

     A(22)=dc(8,1)
     A(23)=dc(8,2)
     A(24)=dc(8,3) 

     !if (id_proc.eq.0) print *,'Cons Gradients',jac(1:6)

  end if


  IERR = 0
  return
end subroutine EV_JAC_G
!
! =============================================================================
!
!                Computation of Hessian of Lagrangian
!
! =============================================================================
!
subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
  implicit none
  integer TASK, N, NEW_X, M, NEW_LAM, NNZH,  i,j,ii
  double precision X(N), OBJFACT, LAM(M), HESS(NNZH), sigmax(N)
  integer IRNH(NNZH), ICNH(NNZH)
  double precision::fmeantmp,fvartmp
  double precision OBJHESS(NNZH),CONHESS(M,NNZH)
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n)
  integer IDAT(*), kprob
  integer IERR,cnt
  double precision  rho, L, sigmay, pi, p, E, Fs 
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  double precision :: hesstmp


!!$

!!$        rho=0.2836
!!$        sigmay=36260.0
!!$        p=25000.0
!!$        L=5.0
!!$        E=30e6
!!$        pi=4.0*atan(1.0)

  Fs=DAT(1)

  kprob=IDAT(1)
  do i=1,N
     sigmax(i)=DAT(i+1)
  end do


  if( TASK.eq.0 ) then
     !
     !     structure of sparse Hessian (lower triangle):
     !
     IRNH(1) = 1
     ICNH(1) = 1

     IRNH(2) = 2
     ICNH(2) = 1

     IRNH(3) = 2
     ICNH(3) = 2

     IRNH(4) = 3
     ICNH(4) = 1

     IRNH(5) = 3
     ICNH(5) = 2

     IRNH(6) = 3
     ICNH(6) = 3

  else

!!$     
!!$     do ii=0,m
!!$
!!$        !      call PCestimate(dim,xavgin,xstdin,fctin,fctindxin,orderinitial,orderfinal,statin,fmeanout,fvarout,fmeanprimeout,fvarprimeout)
!!$        call  PCestimate(N,x,sigmax,11,ii,3,3,0,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)
!!$
!!$        if (ii.eq.0) then
!!$           
!!$           cnt=0
!!$           do i=1,N
!!$              do j=1,N
!!$                 if (i.le.j) then
!!$                    cnt=cnt+1
!!$                    objhess(cnt)=fmeandbleprimetmp(i,j)+kprob*fvardbleprimetmp(i,j)
!!$                 end if
!!$              end do
!!$           end do
!!$
!!$
!!$        else
!!$           
!!$           cnt=0
!!$           do i=1,N
!!$              do j=1,N
!!$                 if (i.le.j) then
!!$                    cnt=cnt+1
!!$                    conhess(ii,cnt)=fmeandbleprimetmp(i,j)+kprob*fvardbleprimetmp(i,j)
!!$                 end if
!!$              end do
!!$           end do
!!$
!!$        end if
!!$
!!$     end do
!!$
!!$     ! Assemble all into HESS
!!$     
!!$     HESS(:)=0.0
!!$     do i=1,NNZH
!!$        hesstmp=0.0
!!$        do j=1,m
!!$           hesstmp=hesstmp+lam(j)*conhess(j,i)
!!$        end do
!!$        hess(i)=hesstmp+objhess(i)
!!$     end do
     
     IERR = 0

  endif

  return
end subroutine EV_HESS











!
! =============================================================================
!
!                   Callback method called once per iteration
!
! =============================================================================
!
subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,DAT, ISTOP)
  use dimkrig

  implicit none
  integer ALG_MODE, ITER_COUNT, LS_TRIAL
  double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
  double precision ALPHA_DU, ALPHA_PR
  double precision DAT(*)
  integer IDAT(*)
  integer ISTOP
  !
  !     You can put some output here
  !
  if (id_proc.eq.0) then

     if (ITER_COUNT .eq.0) then
        write(*,*) 
        write(*,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
     end if

     write(*,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU

  end if
  !
  !     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
  !     simple example.
  !
  if (ITER_COUNT .gt. 1 .and. DNORM.le.1D-03) ISTOP = 1

  return
end subroutine ITER_CB

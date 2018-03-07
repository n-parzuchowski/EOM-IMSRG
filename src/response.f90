module response
  use EOM_IMSRG
  use basic_IMSRG
  use cross_coupled
  implicit none

contains
!==============================================================================================
!==============================================================================================
subroutine compute_response_function(jbas,HS,OP)
  
  integer :: N 
  integer :: nev
  type(spd) :: jbas
  type(sq_op) :: op,V1,Q1,Q2,w1,w2,PIVOT,HS
  type(cc_mat) :: OpCC,QCC,WCC
  type(ex_pandya_mat) :: QPP,WPP
  type(ex_cc_mat) :: OpPP
  real(8),allocatable,dimension(:) :: workl,D,eigs,resid,work,workD
  real(8),allocatable,dimension(:,:) :: V,Z,VX,XXX 
  integer :: i,j,k,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  real(8) ::  x,tol,y,sigma,t1,t2,strength,sm
  character(1) :: BMAT,HOWMNY 
  character(2) :: which
  logical :: rvec
  logical,allocatable,dimension(:) :: selct
!!! this is how many valid eigenvalues you seek, multiply by ten for
!!! number of pivots 
  nev = 100
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Q1%pphh_ph=.true.
  Q2%pphh_ph=.true.
  w1%pphh_ph=.false.
  w2%pphh_ph=.false.
  
  call duplicate_sq_op(OP,Q1,'y') !workspace
  call duplicate_sq_op(OP,Q2,'y') !workspace
  call duplicate_sq_op(Q1,w1,'w') !workspace
  call duplicate_sq_op(Q1,w2,'w') !workspace
  
  call init_ph_mat(Q1,OpPP,jbas) !cross coupled ME
  call init_ph_mat(Q1,QPP,jbas) !cross coupled ME
  call init_ph_wkspc(QPP,WPP) 
  
  h = OP%belowEF !holes
  p = OP%Nsp-h  !particles

  print*,
  print*,
  print*, "==========================================="
  print*, " COMPUTING "//OP%trans_label//" RESPONSE FUNCTION FOR "&
       //nucleus_name(HS%Aneut,HS%Aprot)
  print*, "==========================================="
  
  sps = 0 
  do ix = 1,p
     do jx = 1,h
        
        i = jbas%parts(ix)
        j = jbas%holes(jx) 
  
        if (triangle(jbas%jj(i),jbas%jj(j),OP%rank)) then  

           if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
           if (mod(jbas%ll(i) + jbas%ll(j) + OP%dpar/2,2) .ne. 0 ) cycle
                     
           sps = sps+1
        end if 
     end do
  end do
  

  ! tensor case
  tps = 0 
  do q = 1, OP%nblocks
     
     do II = 1,size( OP%tblck(q)%tgam(3)%X(:,1) ) 
        do JJ = 1, size( OP%tblck(q)%tgam(3)%X(1,:) )  
           
           if (mod(OP%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( OP%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                   OP%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
           end if
           
           if (mod(OP%tblck(q)%Jpair(2)/2,2) == 1) then 
              
              if ( OP%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                   OP%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
           end if
           
           tps = tps+ 1
        end do
     end do
     
     if (OP%tblck(q)%Jpair(1) == OP%tblck(q)%Jpair(2)) cycle
     
     do II = 1,size( OP%tblck(q)%tgam(7)%X(:,1) ) 
        do JJ = 1, size( OP%tblck(q)%tgam(7)%X(1,:) )  
           
           if (mod(OP%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( OP%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                   OP%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
           end if
           
           if (mod(OP%tblck(q)%Jpair(2)/2,2) == 1) then 
              
              if ( OP%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                   OP%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
           end if
           
           tps = tps+ 1
        end do
     end do
     
  end do
           
  print*, '1p1h Amplitudes: ', sps
  print*, '2p2h Amplitudes: ', tps
  N = sps + tps ! number of ph and pphh SDs 
  Q1%neq = N
  Q2%neq = N
  ido = 0  ! status integer is 0 at start
  BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
  which = 'SM' ! compute smallest eigenvalues in magnitude ('SA' is algebraic). 
  tol = 1.0E-10 ! error tolerance?
  info = 1 ! THIS MEANS USE THE PIVOT I TELL YOU TO. 
  ncv = 10*nev ! number of lanczos vectors I guess
  lworkl = ncv*(ncv+8) 
  allocate(V(N,NCV),VX(N,NCV),workl(lworkl))
  LDV = N  
  ishift = 1
  mxiter = 1   !!!! THIS TURNS OFF IMPLICIT RESTART
  mode = 1
   
  allocate(eigs(N),resid(N),work(10*N),workD(3*N)) 
  
  iparam(1) = ishift
  iparam(3) = mxiter
  iparam(7) = mode
  i = 0

  call rewrap_tensor( V(:,1), OP  ,N ,jbas) !!! GET PIVOT

  ! normalize pivot
  strength = sum(V(:,1)**2)
  V(:,1) = V(:,1)/sqrt(strength) 
  resid = V(:,1)

  do 
     ! so V is the krylov subspace matrix that is being diagonalized
     ! it does not need to be initialized, so long as you have the other 
     ! stuff declared right, the code should know this. 
     call dsaupd ( ido, bmat, N, which, nev, tol, resid, &
          ncv, v, ldv, iparam, ipntr, workd, workl, &
          lworkl, info )
     ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
     call progress_bar( i )
     i=i+1 

     if ( ido /= -1 .and. ido /= 1 ) then
        exit
     end if
     if ( i > NCV ) STOP 'response failed to converge' 
     call matvec_nonzeroX_prod(N,HS,Q1,Q2,w1,w2,OpPP,QPP,WPP,jbas, workd(ipntr(1)), workd(ipntr(2)) ) 
     
  end do
  
  rvec= .true. 
  howmny = 'A'

  VX = V 

  allocate(selct(NCV)) 
  allocate(D(NEV)) 
  allocate(Z(N,NEV)) 
  ldz = N  
  call dseupd( rvec, howmny, selct, d, Z, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, info )  


  open (unit=82, file = trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_gaussianfold_'//&
       OP%trans_label//'_response.dat')
  
  open (unit=83, file = trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_discrete_'//&
       OP%trans_label//'_response.dat')
    
  x = 0.d0 
  do i = 1, 10000
     sm = 0.d0 
     do j = 1,NEV
        sm = sm + strength * sum(Z(:,j)*VX(:,1))**2 * gaussian(x,D(j),0.5d0)        
     end do
     
     write(82,'(2(f25.14))') x ,  sm
     x = x + .01d0
  end do
  
  do j = 1, NEV
     write(83,'(2(f25.14))') D(j)-0.00001, 0.0 
     write(83,'(2(f25.14))') D(j),  strength * sum(Z(:,j)*VX(:,1))**2
     write(83,'(2(f25.14))') D(j)+0.00001, 0.0
  end do


  close(82)
  close(83)
  
end subroutine COMPUTE_RESPONSE_FUNCTION

real(8) function gaussian(x,E,sig)
  implicit none

  real(8),intent(in) :: E,sig,x
  real(8) :: div
  
  div = sqrt(2*PI_const)*sig 
  gaussian = exp( -(x-E)**2 /(2.d0*sig**2) )/div 
end function gaussian

real(8) function dirac_delta(x,E,sig)
  implicit none

  real(8),intent(in) :: E,sig,x
  
  if ( abs(x-E) >1e-3) then 
     dirac_delta = 0.0
  else
     dirac_delta = 1.0
  end if
  
end function dirac_delta
   

end module response



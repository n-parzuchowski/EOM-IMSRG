module EOM_IMSRG
  use basic_IMSRG
  use EOM_scalar_commutators
  use EOM_TS_commutators
  use operators
  use isospin_operators
  use EOM_dTZ_commutators
  implicit none
  
contains 


subroutine calculate_excited_states(J,PAR,Numstates,HS,jbas,ladder_ops) 
  ! this is the control routine for the EOM-IMSRG calculation
  implicit none
  
  real(8) :: BE,Mfi ,SD_shell_content,dEtrips,dcgi,dcgi00
  real(8) :: t1,t2,t0,omp_get_wtime,XX,QQ,sm,sm2 
  type(obsv_mgr) :: transitions, moments 
  type(spd) :: jbas
  type(sq_op) :: HS ,newladder
  type(sq_op),dimension(Numstates) :: ladder_ops
  integer :: J,PAR,Numstates,i,q,aa,jj,istart,ist,prots,neuts
  character(2) :: Jlabel,Plabel,betalabel  
  character(2) :: statelabel
  REAL(8),dimension(Numstates) :: Es,BEs ,moms ,trips

  ladder_ops%herm = 1
  ladder_ops%rank = J 
  ladder_ops%dpar = PAR
  ladder_ops%pphh_ph = .true. 
  
  prots = 0 
  neuts = 0 
  do i = 1, jbas%total_orbits,2 
     prots = prots + (jbas%jj(i)+1)*jbas%con(i) 
  end do 
  do i = 2, jbas%total_orbits,2 
     neuts = neuts + (jbas%jj(i)+1)*jbas%con(i) 
  end do   
  
  if ( allocated(phase_hh) ) then
     deallocate(phase_hh,phase_pp)
  end if
  call allocate_tensor(jbas,ladder_ops(1),HS)   
  do q = 1,ladder_ops(1)%nblocks
     ladder_ops(1)%tblck(q)%lam(1) = 1 
  end do
    
  do i = 2, Numstates
     call duplicate_sq_op(ladder_ops(1),ladder_ops(i),'y')
  end do
  
  print* 
  write(*,'((A55),(I1),(A3),(I1))') 'EXECUTING EOM CALCULATION'// &
       ' FOR EXCITED STATES: J=',J/2,' P=',PAR  
  print*

  ! heavy lifting done here 
  call lanczos_diagonalize(jbas,HS,ladder_ops,Numstates)

  print*
  write(*,'((A21),(f16.9))') 'Ground State Energy: ',HS%E0 
  print*
  print*, 'EXCITED STATE ENERGIES:'
  print*, '==========================================='
  print*, '      dE         E_0 + dE         n(1p1h)  ' 
  print*, '==========================================='

  do i = 1, Numstates
     
     Es(i) = ladder_ops(i)%E0          
     write(*,'(6(f16.9))') ladder_ops(i)%E0 , ladder_ops(i)%E0+HS%E0, sum(ladder_ops(i)%fph**2)
       
  end do
  
  ! WRITE STUFF TO FILES. 
  write( Jlabel ,'(I2)') J
  write( betalabel ,'(I2)') nint(HS%lawson_beta)
  if (PAR == 0 ) then 
     Plabel ='+'
  else 
     Plabel ='-'
  end if
  
  Jlabel = adjustl(Jlabel)
  Plabel = adjustl(Plabel)
  betalabel = adjustl(betalabel)
  open(unit=72,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//&
       '_EOM_spec_law'//trim(betalabel)//'.dat')  
  do i = 1, Numstates
     write(72,'(3(f16.9))') ladder_ops(i)%E0,ladder_ops(i)%E0+HS%E0,sum(ladder_ops(i)%fph**2)
  end do

  close(72)
  
end subroutine calculate_excited_states
!==================================================================================
!==================================================================================
subroutine calculate_isospin_states(J,PAR,dTZ,Numstates,HS,jbas,ladder_ops) 
  !control routine for isospin projection changing EOM
  implicit none
  
  real(8) :: BE,Mfi ,SD_shell_content,dEtrips,dcgi,dcgi00
  real(8) :: t1,t2,t0,omp_get_wtime,XX,QQ,sm,sm2 
  type(obsv_mgr) :: transitions, moments 
  type(spd) :: jbas
  type(sq_op) :: HS ,newladder
  type(iso_ladder),dimension(Numstates) :: ladder_ops
  integer :: Numstates,i,q,aa,jj,istart,ist,prots,neuts
  character(2) :: Jlabel,Plabel,betalabel
  character(3) :: tzlab 
  character(2) :: statelabel
  REAL(8),dimension(Numstates) :: Es,BEs ,moms ,trips
  integer,intent(in) :: J,PAR,DTz

  ladder_ops%herm = 1
  ladder_ops%rank = J
  ladder_ops%dTZ = dTZ 
  ladder_ops%dpar = PAR

  prots = 0 
  neuts = 0 
  do i = 1, jbas%total_orbits,2 
     prots = prots + (jbas%jj(i)+1)*jbas%con(i) 
  end do 
  do i = 2, jbas%total_orbits,2 
     neuts = neuts + (jbas%jj(i)+1)*jbas%con(i) 
  end do   

  call allocate_isospin_ladder(jbas,ladder_ops(1),HS)   
  
  do i = 2, Numstates
     call duplicate_isospin_ladder(ladder_ops(1),ladder_ops(i))
  end do
  
  print* 
  write(*,'((A55),(I1),(A3),(I1),(A5),(I1))') 'EXECUTING EOM CALCULATION'// &
       ' FOR EXCITED STATES: J=',J/2,' P=',PAR,' dTZ=',dTZ   
  print*

  t1 = omp_get_wtime() 
  call lanczos_isospin_changer(jbas,HS,ladder_ops,Numstates)
  t2 = omp_get_wtime()
  print*, 'LANCZOS RUN TIME:', t2-t1

     
  print*
  write(*,'((A21),(f16.9))') 'Ground State Energy: ',HS%E0 
  print*
  print*, 'EXCITED STATE ENERGIES:'
  print*, '================================================='
  print*, '      dE              E_0 + dE         n(1p1h)  '
  print*, '================================================='
  do i = 1, Numstates
     Es(i) = ladder_ops(i)%E0     
     write(*,'(3(f16.9))') ladder_ops(i)%E0 , ladder_ops(i)%E0+HS%E0,sum(ladder_ops(i)%fph**2)
       
  end do

end subroutine calculate_isospin_states
!==============================================================================================
!==============================================================================================
subroutine LANCZOS_DIAGONALIZE(jbas,OP,Vecs,nev)    
  
  integer :: N 
  integer,intent(in) :: nev
  type(spd) :: jbas
  type(sq_op) :: op,V1,Q1,Q2,w1,w2
  type(sq_op),dimension(nev) :: Vecs
  type(cc_mat) :: OpCC,QCC,WCC
  type(ex_pandya_mat) :: QPP,WPP
  type(ex_cc_mat) :: OpPP
  real(8),allocatable,dimension(:) :: workl,D,eigs,resid,work,workD
  real(8),allocatable,dimension(:,:) :: V,Z
  integer :: i,j,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  real(8) ::  x,tol,y,sigma,t1,t2
  character(1) :: BMAT,HOWMNY 
  character(2) :: which
  logical :: rvec
  logical,allocatable,dimension(:) :: selct

  Q1%pphh_ph=.true.
  Q2%pphh_ph=.true.
  w1%pphh_ph=.false.
  w2%pphh_ph=.false.

  call duplicate_sq_op(vecs(1),w1,'w') !workspace
  call duplicate_sq_op(vecs(1),w2,'w') !workspace
  call duplicate_sq_op(vecs(1),Q1,'y') !workspace
  call duplicate_sq_op(vecs(1),Q2,'y') !workspace
 
  call init_ph_mat(Op,OpPP,jbas) !cross coupled ME
  call init_ph_mat(vecs(1),QPP,jbas) !cross coupled ME
  call init_ph_wkspc(QPP,WPP) 
    
  h = OP%belowEF !holes
  p = OP%Nsp-h  !particles
 
  sps = 0 
  do ix = 1,p
     do jx = 1,h
        
        i = jbas%parts(ix)
        j = jbas%holes(jx) 
  
        if (triangle(jbas%jj(i),jbas%jj(j),vecs(1)%rank)) then  

           if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
           if (mod(jbas%ll(i) + jbas%ll(j) + vecs(1)%dpar/2,2) .ne. 0 ) cycle
        
           sps = sps+1
        end if 
     end do
  end do
  
  tps = 0 
  do q = 1, vecs(1)%nblocks

     do II = 1,size( vecs(1)%tblck(q)%tgam(3)%X(:,1) ) 
        do JJ = 1, size( vecs(1)%tblck(q)%tgam(3)%X(1,:) )  

           if (mod(vecs(1)%tblck(q)%Jpair(1)/2,2) == 1) then 

              if ( vecs(1)%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                   vecs(1)%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
           end if

           if (mod(vecs(1)%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( vecs(1)%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                   vecs(1)%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
           end if

           tps = tps+ 1
        end do
     end do

     if (vecs(1)%tblck(q)%Jpair(1) == vecs(1)%tblck(q)%Jpair(2)) cycle

     do II = 1,size( vecs(1)%tblck(q)%tgam(7)%X(:,1) ) 
        do JJ = 1, size( vecs(1)%tblck(q)%tgam(7)%X(1,:) )  

           if (mod(vecs(1)%tblck(q)%Jpair(1)/2,2) == 1) then 

              if ( vecs(1)%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                   vecs(1)%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
           end if

           if (mod(vecs(1)%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( vecs(1)%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                   vecs(1)%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
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
  Vecs%neq = N
  ido = 0  ! status integer is 0 at start
  BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
  which = 'SM' ! compute smallest eigenvalues in magnitude ('SA') is algebraic. 
  tol = 0.0E+00 ! error tolerance, machine precision. 
  info = 0
  ncv = 5*nev ! number of lanczos vectors I guess
  lworkl = ncv*(ncv+8) 
  allocate(V(N,NCV),workl(lworkl))
  LDV = N  
  ishift = 1
  mxiter = 500 
  mode = 1
   
  allocate(eigs(N),resid(N),work(10*N),workD(3*N)) 
  
  iparam(1) = ishift
  iparam(3) = mxiter
  iparam(7) = mode
  i = 0

  do 
     ! so V is the krylov subspace matrix that is being diagonalized
     ! it does not need to be initialized, so long as you have the other 
     ! stuff declared right, the code should know this. 
     call dsaupd ( ido, bmat, N, which, nev, tol, resid, &
          ncv, v, ldv, iparam, ipntr, workd, workl, &
          lworkl, info ) ! ARPACK IMPLICITLY RESTARTED LANCZOS ROUTINE
     
     ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
     call progress_bar( i )

     i=i+1 

     if ( ido /= -1 .and. ido /= 1 ) then
        exit
     end if
     
     call matvec_nonzeroX_prod(N,OP,Q1,Q2,w1,w2,OpPP,QPP,WPP,jbas, workd(ipntr(1)), workd(ipntr(2)) ) 
     
  end do
  write(6,*) 
  ! the ritz values are out of order right now. Need to do post
  ! processing to fix this, and get the eigenvectors
  rvec= .true. 
  howmny = 'A'
  
  allocate(selct(NCV)) 
  allocate(D(NEV)) 
  allocate(Z(N,NEV)) 
  ldz = N  
  call dseupd( rvec, howmny, selct, d, Z, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, info ) ! ARPACK POST PROCESSING FOR DSAUPD 
  
  ! right now Z contains the eigenvectors in the columns
  ! d contains the eigenvalues in the same order. 

  do i = 1, nev
     call unwrap_tensor(Z(:,i),Vecs(i),N,jbas) 
     Vecs(i)%E0 = d(i)
  end do 
      
end subroutine LANCZOS_DIAGONALIZE
!===================================================================
!===================================================================
subroutine LANCZOS_ISOSPIN_CHANGER(jbas,OP,Vecs,nev)    
  
  integer :: N 
  integer,intent(in) :: nev
  type(spd) :: jbas
  type(sq_op) :: op
  type(iso_ladder) :: V1,Q1,Q2
  type(iso_ladder),dimension(nev) :: Vecs
  type(cc_mat) :: OpCC
  real(8),allocatable,dimension(:) :: workl,D,eigs,resid,work,workD
  real(8),allocatable,dimension(:,:) :: V,Z
  integer :: i,j,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  real(8) ::  x,tol,y,sigma,t1,t2
  character(1) :: BMAT,HOWMNY 
  character(2) :: which
  logical :: rvec
  logical,allocatable,dimension(:) :: selct

  Q1%pphh_ph=.true.
  Q2%pphh_ph=.true.

  call duplicate_isospin_ladder(vecs(1),Q1) !workspace
  call duplicate_isospin_ladder(vecs(1),Q2) !workspace
   
  h = OP%belowEF !holes
  p = OP%Nsp-h  !particles
 
  sps = 0 
  do ix = 1,p
     do jx = 1,h
        
        i = jbas%parts(ix)
        j = jbas%holes(jx) 
  
        if (triangle(jbas%jj(i),jbas%jj(j),vecs(1)%rank)) then  

           if (jbas%itzp(i)-vecs(1)%dTz*2 .ne. jbas%itzp(j) ) cycle
           if (mod(jbas%ll(i) + jbas%ll(j) + vecs(1)%dpar/2,2) .ne. 0 ) cycle
        
           sps = sps+1
        end if 
     end do
  end do
  

  ! tensor case
  tps = 0 
  do q = 1, vecs(1)%nblocks
     
     do II = 1,vecs(1)%tblck(q)%npp1 
        do JJ = 1, vecs(1)%tblck(q)%nhh2 

           if (mod(vecs(1)%tblck(q)%Jpair(1)/2,2) == 1) then 

              if ( vecs(1)%tblck(q)%qn(1)%Y(II,1) == &
                   vecs(1)%tblck(q)%qn(1)%Y(II,2) ) cycle
           end if

           if (mod(vecs(1)%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( vecs(1)%tblck(q)%qn(2)%Y(JJ,1) == &
                   vecs(1)%tblck(q)%qn(2)%Y(JJ,2) ) cycle
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
  Vecs%neq = N
  ido = 0  ! status integer is 0 at start
  BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
  which = 'SM' ! compute smallest eigenvalues in magnitude ('SA') is algebraic. 
  tol = 0.0E+00 ! error tolerance? (wtf zero?) 
  info = 0
  ncv = 5*nev ! number of lanczos vectors I guess
  lworkl = ncv*(ncv+8) 
  allocate(V(N,NCV),workl(lworkl))
  LDV = N  
  ishift = 1
  mxiter = 500 
  mode = 1
   
  allocate(eigs(N),resid(N),work(10*N),workD(3*N)) 

  iparam(1) = ishift
  iparam(3) = mxiter
  iparam(7) = mode
  i = 0

  do 
     ! so V is the krylov subspace matrix that is being diagonalized
     ! it does not need to be initialized, so long as you have the other 
     ! stuff declare right, the code should know this. 
     call dsaupd ( ido, bmat, N, which, nev, tol, resid, &
          ncv, v, ldv, iparam, ipntr, workd, workl, &
          lworkl, info )
     ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
     call progress_bar( i )
     i=i+1 
     
     if ( ido /= -1 .and. ido /= 1 ) then
        exit
     end if
     
     call matvec_isoX_prod(N,OP,Q1,Q2,jbas, workd(ipntr(1)), workd(ipntr(2)) ) 
          
  end do
  write(6,*) 
  ! the ritz values are out of order right now. Need to do post
  ! processing to fix this, and get the eigenvectors
  rvec= .true. 
  howmny = 'A'
  
  allocate(selct(NCV)) 
  allocate(D(NEV)) 
  allocate(Z(N,NEV)) 
  ldz = N  
  call dseupd( rvec, howmny, selct, d, Z, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, info )
  
  ! right now Z contains the eigenvectors in the columns
  ! d contains the eigenvalues in the same order. 
  do i = 1, nev
     call unwrap_iso_ladder(Z(:,i),Vecs(i),N,jbas) 
     Vecs(i)%E0 = d(i)
  end do 
      
end subroutine LANCZOS_ISOSPIN_CHANGER
!=========================================================================
!=========================================================================
subroutine matvec_nonzeroX_prod(N,OP,Q_op,Qout,w1,w2,OpCC,QCC,WCC,jbas,v,w) 
  implicit none 
  
  integer :: N ,q,a,b,c,d,i,j,k,l,Jtot
  real(8) :: sm
  type(sq_op) :: OP, Q_op ,Qout,w1,w2
  type(ex_pandya_mat) :: QCC,WCC
  type(ex_cc_mat) :: OpCC
  type(spd) :: jbas
  real(8),dimension(N) :: v,w 
  real(8) :: coef9,dfact0

  ! the name says mat-vec product, but really this
  ! is a commutator. 
  
  ! the commutator here is equivalent to the matrix-vector product 
  ! with only connected diagrams retained. 

  
  ! FIRST WE NEED TO CONVERT v TO a (SQ_OP) variable

  call unwrap_tensor(v,Q_op,N,jbas)
  ! now we have two sq_op operators which can be used with my commutator expressions. Noice. 
  
  call EOM_generalized_pandya(Q_op,QCC,jbas)
  call calculate_cross_coupled_pphh(Op,OpCC,jbas) 
  
  call EOM_TS_commutator_111(Op,Q_op,Qout,jbas) 
  call EOM_TS_commutator_121(Op,Q_op,Qout,jbas)
  call EOM_TS_commutator_211(OpCC,Q_op,Qout,jbas) 
  call EOM_TS_commutator_122(Op,Q_op,Qout,jbas)
  call EOM_TS_commutator_212(Op,Q_op,Qout,jbas)

  call EOM_TS_commutator_222_pp_hh(Op,Q_op,Qout,w1,w2,jbas)  
  call EOM_TS_commutator_221(w1,w2,Op%herm*Q_op%herm,Qout,jbas)
  call EOM_TS_commutator_222_ph(OpCC,QCC,Qout,WCC,jbas)
  
  ! Okay, now we have the "matrix vector product" So lets put it back in vector form:
  call rewrap_tensor(w,Qout,N,jbas) 

end subroutine matvec_nonzeroX_prod
!======================================================================================
!======================================================================================
subroutine matvec_ISOX_prod(N,OP,Q_op,Qout,jbas,v,w) 
  implicit none 
  
  integer :: N ,q,a,b,c,d,i,j,k,l,Jtot
  real(8) :: sm
  type(sq_op) :: OP
  type(iso_ladder) :: Q_op ,Qout
  type(spd) :: jbas
  real(8),dimension(N) :: v,w 
  real(8) :: coef9,dfact0

  ! the name says mat-vec product, but really this
  ! is a commutator. 
  
  ! the commutator here is equivalent to the matrix-vector product 
  ! with only connected diagrams retained. 

  
  ! FIRST WE NEED TO CONVERT v TO a (SQ_OP) variable

  call unwrap_iso_ladder(v,Q_op,N,jbas)
  ! now we have two sq_op operators which can be used with my commutator expressions. Noice. 
  call EOM_dTZ_commutator_111(Op,Q_op,Qout,jbas) 
  call EOM_dTZ_commutator_121(Op,Q_op,Qout,jbas)
  call EOM_dTZ_commutator_211(Op,Q_op,Qout,jbas) 
  call EOM_dTZ_commutator_122(Op,Q_op,Qout,jbas)
  call EOM_dTZ_commutator_212(Op,Q_op,Qout,jbas)
  
  call EOM_dTZ_commutator_222_pp_hh(Op,Q_op,Qout,jbas)
  
  call EOM_dTZ_commutator_222_ph(Op,Q_op,Qout,jbas)
  
  ! Okay, now we have the "matrix vector product" So lets put it back in vector form:
  call rewrap_iso_ladder(w,Qout,N,jbas) 

end subroutine
!======================================================================================
!======================================================================================
subroutine progress_bar( step  )  
  implicit none
  
  integer :: step
 
  if ( step .ne. 0) then 
     ! hold the backspace key down
     write(6,'(A)',advance='no') char(8)//char(8)//char(8)//char(8)//char(8) &
       //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)&
       //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8) &
       //char(8)//char(8)//char(8)//char(8)
     flush 6
  end if 

  write(6,'((A19),(I5))',advance='no') 'Lanczos iteration: ',step
  flush 6 

end subroutine 
!=====================================================
!=====================================================
real(8) function EOM_triples(H,Xdag,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H,Xdag
  
  EOM_triples = tensor_Triples(H,Xdag,jbas)

end function EOM_triples
!=====================================================
!=====================================================
real(8) function tensor_triples(H,Xdag,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(tpd),allocatable,dimension(:) :: threebas
  type(sq_op) :: H,Xdag,Xrag
  integer :: a,b,c,i,j,k,jtot1,jtot2,Jab,Jij,g1
  integer :: ja,jb,jc,ji,jj,jk,AAA,q,q2,TZ,PAR
  integer :: ax,bx,cx,ix,jx,kx,III,rank,fails
  integer :: jab_min,jab_max,jij_min,jij_max
  integer :: J_min, J_max,x,total_threads,thread
  real(8) :: faa,fbb,fcc,fii,fjj,fkk,Gabab,Gkbkb,Gkckc
  real(8) :: Gacac,Gbcbc,Gijij,Gikik,Gjkjk,Giaia
  real(8) :: Gibib,Gicic,Gjaja,Gjbjb,Gjcjc,Gkaka  
  real(8) :: sm,denom,dlow,w,w_test,min_Denom,pre1,pre2
  
  sm = 0.d0   
  call enumerate_three_body(threebas,jbas)  
  total_threads = size(threebas(1)%direct_omp) - 1
  rank = Xdag%rank
  fails = 0
!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(threebas,jbas,H,Xdag) & 
!$OMP& REDUCTION(+:sm) REDUCTION(+:global_counter1) REDUCTION(+:global_counter2)
  
  do thread = 1, total_threads
  do q = 1+threebas(1)%direct_omp(thread),&
       threebas(1)%direct_omp(thread+1)
  !!do q = 1, size(threebas)

     jtot1 = threebas(q)%chan(1)      
     TZ = threebas(q)%chan(2)      
     PAR = threebas(q)%chan(3)      

     do AAA = 1, size(threebas(q)%ppp(:,1)) 
        pre1 = 1.d0
        a = threebas(q)%ppp(AAA,1)
        b = threebas(q)%ppp(AAA,2)
        c = threebas(q)%ppp(AAA,3)
        
        if (a==b) then 
           if (a==c) then 
              pre1 = 6.d0
           else
              pre1 = 2.d0
           end if
        else if (a==c) then 
           pre1 = 2.d0
        else if (b==c) then 
           pre1 = 2.d0 
        else
           pre1 = 1.d0 
        end if

        ja = jbas%jj(a)      
        jb = jbas%jj(b) 
        jc = jbas%jj(c) 

        jab_min = abs(ja-jb) 
        jab_max = ja+jb

        faa = f_elem(a,a,H,jbas)
        fbb = f_elem(b,b,H,jbas)
        fcc = f_elem(c,c,H,jbas)
        Gabab = twobody_monopole(a,b,ja,jb,H,jbas) 
        Gacac = twobody_monopole(a,c,ja,jc,H,jbas) 
        Gbcbc = twobody_monopole(b,c,jb,jc,H,jbas) 
        
        do jab = jab_min,jab_max,2
           
           if ( .not. (triangle(jtot1,jc,Jab))) cycle
           if ((a==b).and.(mod(Jab/2,2)==1)) cycle
                 
           do q2 = 1, size(threebas)

              jtot2 = threebas(q2)%chan(1)

              if (threebas(q2)%chan(2) .ne. TZ ) cycle
              if (mod(threebas(q2)%chan(3)+Xdag%dpar/2,2) .ne. PAR) cycle
              if (.not. triangle(jtot1,jtot2,rank) ) cycle 

              do III = 1, size(threebas(q2)%hhh(:,1)) 
                 
                 pre2=1.d0
                 i = threebas(q2)%hhh(III,1)
                 j = threebas(q2)%hhh(III,2)
                 k = threebas(q2)%hhh(III,3)
                 
                 if (i==j) then 
                    if (i==k) then 
                       pre2 = 6.d0
                    else
                       pre2 = 2.d0
                    end if
                 else if (i==k) then 
                    pre2 = 2.d0
                 else if (j==k) then 
                    pre2 = 2.d0 
                 else
                    pre2 = 1.d0 
                 end if

                 ji = jbas%jj(i)       
                 jj = jbas%jj(j) 
                 jk = jbas%jj(k)  

                 jij_min = abs(ji-jj) 
                 jij_max = ji+jj

                 fii = f_elem(i,i,H,jbas)
                 fjj = f_elem(j,j,H,jbas)
                 fkk = f_elem(k,k,H,jbas)

                 Gijij = twobody_monopole(i,j,ji,jj,H,jbas) 
                 Gikik = twobody_monopole(i,k,ji,jk,H,jbas) 
                 Gjkjk = twobody_monopole(j,k,jj,jk,H,jbas) 

                 Giaia = twobody_monopole(i,a,ji,ja,H,jbas) 
                 Gibib = twobody_monopole(i,b,ji,jb,H,jbas) 
                 Gicic = twobody_monopole(i,c,ji,jc,H,jbas) 

                 Gjaja = twobody_monopole(j,a,jj,ja,H,jbas) 
                 Gjbjb = twobody_monopole(j,b,jj,jb,H,jbas) 
                 Gjcjc = twobody_monopole(j,c,jj,jc,H,jbas) 

                 Gkaka = twobody_monopole(k,a,jk,ja,H,jbas) 
                 Gkbkb = twobody_monopole(k,b,jk,jb,H,jbas) 
                 Gkckc = twobody_monopole(k,c,jk,jc,H,jbas) 

                 denom = (Xdag%E0-(faa+fbb+fcc-fii-fjj-fkk+Gabab+&
                      Gacac+Gbcbc+Gijij+Gikik+Gjkjk-Giaia&
                      -Gibib-Gicic-Gjaja-Gjbjb-Gjcjc-Gkaka-&
                      Gkbkb-Gkckc))*pre1*pre2
              
              
                 do jij = jij_min, jij_max,2
                    
                    if ( .not. (triangle(jtot2,jk,Jij))) cycle
                    if ((i==j).and.(mod(Jij/2,2)==1)) cycle
                    w = EOM_TS_commutator_223_single(H,Xdag,a,b,c,i,j,k,jtot1,jtot2,jab,jij,jbas)
                    sm = sm + w*w/denom*(jtot1+1.d0)/(rank+1.d0)

                 end do
              end do

           end do
        end do
     end do
  end do
  end do
 !$OMP END PARALLEL DO 
  tensor_triples = sm
  
end function tensor_triples
!======================================================================================
!======================================================================================
subroutine unwrap_tensor( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX,qx
  real(8),dimension(N) :: v
  type(sq_op) :: AX 
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  

  do ix = 1,parts
     do jx = 1,holes
        
        ii = jbas%parts(ix)
        JJ = jbas%holes(jx) 
  
        if (triangle(jbas%jj(II),jbas%jj(JJ),AX%rank)) then  

           if (jbas%itzp(II) .ne. jbas%itzp(JJ) ) cycle
           if (mod(jbas%ll(II) + jbas%ll(JJ) + AX%dpar/2,2) .ne. 0 ) cycle
        
           AX%fph(IX,JX) = v(i) 
           i = i+ 1
        end if 
     end do
  end do
  
  ! tensor case
  
  do q = 1, AX%nblocks
     
     do II = 1,size( AX%tblck(q)%tgam(3)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(3)%X(1,:) )  
           
           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( AX%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
           end if

           AX%tblck(q)%tgam(3)%X(II,JJ) = v(i)
           i = i + 1
        end do
     end do

     ! IF THE Js ARE THE SAME THEN WE NEED TO MAKE SURE THAT IS REFLECTED IN 
     ! THE TRANSPOSE
     if (AX%tblck(q)%Jpair(1) == AX%tblck(q)%Jpair(2)) then 
        if (AX%dpar == 2 ) then 
           qx = tensor_block_index(AX%tblck(q)%Jpair(1),AX%tblck(q)%Jpair(2)&
                ,AX%rank,AX%tblck(q)%lam(3),mod(AX%tblck(q)%lam(2)+1,2))
           AX%tblck(qx)%tgam(7)%X = Transpose(AX%tblck(q)%tgam(3)%X)  
        else       
           AX%tblck(q)%tgam(7)%X = Transpose(AX%tblck(q)%tgam(3)%X)  
        end if
        
        cycle 
     
     end if
     
     ! OTHERWISE BUSINESS AS USUAL.
     do II = 1,size( AX%tblck(q)%tgam(7)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(7)%X(1,:) )  

           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
           end if

           AX%tblck(q)%tgam(7)%X(II,JJ) = v(i)
           i = i + 1
        end do
     end do

  end do

end subroutine unwrap_tensor
!======================================================================================
!======================================================================================
subroutine unwrap_iso_ladder( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX,qx
  real(8),dimension(N) :: v
  type(iso_ladder) :: AX 
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  

  do ix = 1,parts
     do jx = 1,holes
        
        ii = jbas%parts(ix)
        JJ = jbas%holes(jx) 
  
        if (triangle(jbas%jj(II),jbas%jj(JJ),AX%rank)) then  

           if (jbas%itzp(II)-AX%dTz*2 .ne. jbas%itzp(JJ) ) cycle
           if (mod(jbas%ll(II) + jbas%ll(JJ) + AX%dpar/2,2) .ne. 0 ) cycle
        
           AX%fph(IX,JX) = v(i) 
           i = i+ 1
        end if 
     end do
  end do
  
  ! tensor case
  
  do q = 1, AX%nblocks
     
     do II = 1,AX%tblck(q)%npp1
        do JJ = 1, AX%tblck(q)%nhh2
           
           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( AX%tblck(q)%qn(1)%Y(II,1) == &
                   AX%tblck(q)%qn(1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%qn(2)%Y(JJ,1) == &
                   AX%tblck(q)%qn(2)%Y(JJ,2) ) cycle
           end if

           AX%tblck(q)%Xpphh(II,JJ) = v(i)
           i = i + 1
        end do
     end do

  end do

end subroutine unwrap_iso_ladder
!======================================================================================
!======================================================================================
subroutine rewrap_iso_ladder( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX,qx
  real(8),dimension(N) :: v
  type(iso_ladder) :: AX 
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  

  do ix = 1,parts
     do jx = 1,holes
        
        ii = jbas%parts(ix)
        JJ = jbas%holes(jx) 
  
        if (triangle(jbas%jj(II),jbas%jj(JJ),AX%rank)) then  

           if (jbas%itzp(II)-AX%dTz*2 .ne. jbas%itzp(JJ) ) cycle
           if (mod(jbas%ll(II) + jbas%ll(JJ) + AX%dpar/2,2) .ne. 0 ) cycle
        
           v(i) = AX%fph(IX,JX)
           i = i+ 1
        end if 
     end do
  end do
  
  ! tensor case
  
  do q = 1, AX%nblocks
     
     do II = 1,AX%tblck(q)%npp1
        do JJ = 1, AX%tblck(q)%nhh2
           
           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( AX%tblck(q)%qn(1)%Y(II,1) == &
                   AX%tblck(q)%qn(1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%qn(2)%Y(JJ,1) == &
                   AX%tblck(q)%qn(2)%Y(JJ,2) ) cycle
           end if

           v(i) = AX%tblck(q)%Xpphh(II,JJ)
           i = i + 1
        end do
     end do

  end do

end subroutine rewrap_iso_ladder
!======================================================================================
!======================================================================================
subroutine rewrap_tensor( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX,qx
  real(8),dimension(N) :: v
  type(sq_op) :: AX 
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  
  do ix = 1,parts
     do jx = 1,holes
        
        ii = jbas%parts(ix)
        JJ = jbas%holes(jx) 
  
        if (triangle(jbas%jj(II),jbas%jj(JJ),AX%rank)) then  

           if (jbas%itzp(II) .ne. jbas%itzp(JJ) ) cycle
           if (mod(jbas%ll(II) + jbas%ll(JJ) + AX%dpar/2,2) .ne. 0 ) cycle
        
           v(i) = AX%fph(IX,JX) 
           i = i + 1
        end if 
     end do
  end do
  
  ! tensor case
  
  do q = 1, AX%nblocks
     
     do II = 1,size( AX%tblck(q)%tgam(3)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(3)%X(1,:) )  
           
           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( AX%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
           end if

           v(i) = AX%tblck(q)%tgam(3)%X(II,JJ) 
           i = i + 1
        end do
     end do

     ! IF THE Js ARE THE SAME THEN WE NEED TO MAKE SURE THAT IS REFLECTED IN 
     ! THE TRANSPOSE
     if (AX%tblck(q)%Jpair(1) == AX%tblck(q)%Jpair(2)) then 
        
        cycle 
     
     end if
     
     ! OTHERWISE BUSINESS AS USUAL.
     do II = 1,size( AX%tblck(q)%tgam(7)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(7)%X(1,:) )  

           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
           end if

           v(i) = AX%tblck(q)%tgam(7)%X(II,JJ) 
           i = i + 1
        end do
     end do

  end do

end subroutine rewrap_tensor
!======================================================
!======================================================
subroutine write_ladder_operators(AX,jbas) 
  ! write the ladder operators to file in the playplace directory
  type(sq_op),dimension(:) :: AX 
  type(spd) :: jbas
  integer Atot,Ntot,q,i,numstates,neq
  real(8),allocatable,dimension(:):: outvec 
  character(200) :: prefix2,stringout
  integer(c_int) :: rx,filehandle
  logical :: isthere
  
  numstates = size(AX)

  if (prefix(1:8) == 'testcase') return  

  do i = 1,200
     if (prefix(i:i+1) == 'hw') exit
  end do
  
  prefix2(1:i+6)=prefix(1:i+6) 
  print*, 'Writing normal ordered ladder operators to ',&
       playplace//trim(adjustl(prefix2(1:i+6)))//&
       '_ladder.gz'
  neq = AX(1)%neq
  
  allocate(outvec(numstates*neq)) 

  do q = 1, numstates
     call rewrap_tensor(outvec((q-1)*neq+1:q*neq),AX(q),neq,jbas) 
  end do
  

  filehandle = gzOpen(playplace//trim(adjustl(prefix2(1:i+6)))//&
       '_ladder.gz'//achar(0),'w'//achar(0)) 


  write(stringout(1:20),'(d20.14)') dfloat(neq) 
  stringout = adjustl(trim(adjustl(stringout(1:20)))//' XXXX')
  call write_gz(filehandle,stringout) 

  do q = 1,numstates
       write(stringout(1:20),'(d20.14)') AX(q)%E0 
       stringout = adjustl(trim(adjustl(stringout(1:20)))//' XXXX')
       call write_gz(filehandle,stringout) 
  end do 
  
  do q =1,neq*numstates
     write(stringout(1:20),'(d20.14)') outvec(q)
     stringout = adjustl(trim(adjustl(stringout(1:20)))//' XXXX')
     call write_gz(filehandle,stringout)    
  end do
   
  rx = gzClose(filehandle) 
end subroutine write_ladder_operators
!======================================================
!======================================================
logical function read_ladder_operators(AX,jbas) 
  ! read the ladder operators from file in the playplace directory
  type(sq_op),dimension(:) :: AX 
  type(spd) :: jbas
  integer Atot,Ntot,q,i,numstates,neq
  real(8),allocatable,dimension(:):: outvec 
  character(200) :: prefix2
  character(20) :: instring
  integer(c_int) :: rx,filehandle
  logical :: isthere
  real(8) :: neq_float
  
  read_ladder_operators = .true. 
  if (prefix(1:8) == 'testcase') return  

  do i = 1,200
     if (prefix(i:i+1) == 'hw') exit
  end do
  
  prefix2(1:i+6)=prefix(1:i+6) 
    
  inquire(file=playplace//trim(adjustl(prefix2(1:i+6)))//&
       '_ladder.gz',exist=isthere)
  
  if ( .not. isthere ) then 
     return
  end if 
  
  numstates = size(AX)

  print*, 'Reading normal ordered ladder operators from ',&
       playplace//trim(adjustl(prefix2(1:i+6)))//&
       '_ladder.gz'

  
  filehandle = gzOpen(playplace//trim(adjustl(prefix2(1:i+6)))//&
       '_ladder.gz'//achar(0),'r'//achar(0)) 
  
  instring = read_normal_gz(filehandle) 

  read(instring,'(d20.14)') neq_float 
  
  neq = nint(neq_float) 
  allocate(outvec(numstates*neq)) 
  
  do q =1,numstates
     instring = read_normal_gz(filehandle) 
     read(instring,'(d20.14)') AX(q)%E0
  end do

  do q =1,neq*numstates
     instring = read_normal_gz(filehandle) 
     read(instring,'(d20.14)') outvec(q)
  end do

  do q = 1, numstates
     call unwrap_tensor( outvec((q-1)*neq+1:q*neq), AX(q) ,neq ,jbas) 
  end do

  rx = gzClose(filehandle) 
  read_ladder_operators = .false. 

end function read_ladder_operators
!===============================================================================
!===============================================================================
real(8) function W_mscheme(p,mp,q,mq,r,mr,s,ms,t,mt,u,mu,MU_B,AA,BB,jbas)
  ! for testing
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB
  integer :: p,mp,q,mq,r,mr,s,ms,t,mt,u,mu,mu_B
  integer :: ia,a,ma,j1,m1,j2,m2,Jpq,Mpq,Jst,Mst
  integer :: ja,jp,jq,jr,js,jt,ju
  real(8) :: sm 
  
  jp = jbas%jj(p)
  jq = jbas%jj(q)
  jr = jbas%jj(r)
  js = jbas%jj(s)
  jt = jbas%jj(t)
  ju = jbas%jj(u) 
  sm = 0.d0 

  do ia = 1,AA%belowEF
     a = jbas%holes(ia) 
     ja = jbas%jj(a) 
     
     do ma = -1*ja,ja,2 
        
        sm = sm - v_mscheme(p,mp,a,ma,s,ms,t,mt,AA,jbas)* &
             tensor_mscheme(q,mq,r,mr,a,ma,u,mu,mu_B,BB,jbas)

        sm = sm - v_mscheme(p,mp,a,ma,t,mt,u,mu,AA,jbas)* &
             tensor_mscheme(q,mq,r,mr,a,ma,s,ms,mu_B,BB,jbas)

        sm = sm - v_mscheme(p,mp,a,ma,u,mu,s,ms,AA,jbas)* &
             tensor_mscheme(q,mq,r,mr,a,ma,t,mt,mu_B,BB,jbas)

        sm = sm - v_mscheme(q,mq,a,ma,s,ms,t,mt,AA,jbas)* &
             tensor_mscheme(r,mr,p,mp,a,ma,u,mu,mu_B,BB,jbas)

        sm = sm - v_mscheme(q,mq,a,ma,t,mt,u,mu,AA,jbas)* &
             tensor_mscheme(r,mr,p,mp,a,ma,s,ms,mu_B,BB,jbas)

        sm = sm - v_mscheme(q,mq,a,ma,u,mu,s,ms,AA,jbas)* &
             tensor_mscheme(r,mr,p,mp,a,ma,t,mt,mu_B,BB,jbas)

        sm = sm - v_mscheme(r,mr,a,ma,s,ms,t,mt,AA,jbas)* &
             tensor_mscheme(p,mp,q,mq,a,ma,u,mu,mu_B,BB,jbas)

        sm = sm - v_mscheme(r,mr,a,ma,t,mt,u,mu,AA,jbas)* &
             tensor_mscheme(p,mp,q,mq,a,ma,s,ms,mu_B,BB,jbas)

        sm = sm - v_mscheme(r,mr,a,ma,u,mu,s,ms,AA,jbas)* &
             tensor_mscheme(p,mp,q,mq,a,ma,t,mt,mu_B,BB,jbas)
     end do
  end do

  do ia = 1,AA%Nsp-AA%belowEF

     a = jbas%parts(ia) 
     ja = jbas%jj(a) 
     
     do ma = -1*ja,ja,2 

        sm = sm + tensor_mscheme(p,mp,a,ma,s,ms,t,mt,mu_B,BB,jbas)* &
             v_mscheme(q,mq,r,mr,a,ma,u,mu,AA,jbas)

        sm = sm + tensor_mscheme(p,mp,a,ma,t,mt,u,mu,mu_B,BB,jbas)* &
             v_mscheme(q,mq,r,mr,a,ma,s,ms,AA,jbas)

        sm = sm + tensor_mscheme(p,mp,a,ma,u,mu,s,ms,mu_B,BB,jbas)* &
             v_mscheme(q,mq,r,mr,a,ma,t,mt,AA,jbas)

        sm = sm + tensor_mscheme(q,mq,a,ma,s,ms,t,mt,mu_B,BB,jbas)* &
             v_mscheme(r,mr,p,mp,a,ma,u,mu,AA,jbas)

        sm = sm + tensor_mscheme(q,mq,a,ma,t,mt,u,mu,mu_B,BB,jbas)* &
             v_mscheme(r,mr,p,mp,a,ma,s,ms,AA,jbas)

        sm = sm + tensor_mscheme(q,mq,a,ma,u,mu,s,ms,mu_B,BB,jbas)* &
             v_mscheme(r,mr,p,mp,a,ma,t,mt,AA,jbas)

        sm = sm + tensor_mscheme(r,mr,a,ma,s,ms,t,mt,mu_B,BB,jbas)* &
             v_mscheme(p,mp,q,mq,a,ma,u,mu,AA,jbas)

        sm = sm + tensor_mscheme(r,mr,a,ma,t,mt,u,mu,mu_B,BB,jbas)* &
             v_mscheme(p,mp,q,mq,a,ma,s,ms,AA,jbas)

        sm = sm + tensor_mscheme(r,mr,a,ma,u,mu,s,ms,mu_B,BB,jbas)* &
             v_mscheme(p,mp,q,mq,a,ma,t,mt,AA,jbas)
        
     end do
  end do
  W_mscheme = sm 

end function W_mscheme
!=====================================================================================
!=====================================================================================
real(8)  function W_via_mscheme(p,q,r,s,t,u,j1,j2,Jpq,Jst,AA,BB,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  integer :: p,q,r,s,t,u,j1,j2,Jpq,Jst
  integer :: mp,mq,mr,ms,mt,mu,m1,m2,Mpq,Mst
  real(8) :: sm,dcgi00,dcgi
  integer :: jp,jq,jr,js,jt,ju,rank,Mew
  
  sm = dcgi00()
  sm = 0.d0 

  jp = jbas%jj(p)
  jq = jbas%jj(q)
  jr = jbas%jj(r)
  js = jbas%jj(s)
  jt = jbas%jj(t)
  ju = jbas%jj(u)
  
  rank = BB%rank

  do mp = -1*jp,jp,2
     do mq = -1*jq,jq,2
        Mpq = mp+mq
        do mr = -1*jr,jr,2
           do ms = -1*js,js,2
              do mt = -1*jt,jt,2
                 Mst = ms+mt
                 do mu = -1*ju,ju,2
                    do m1 = -1*j1,j1,2
                       do m2 = -1*j2,j2,2
                          do Mew = -1*Rank,Rank,2 
                             
                             sm = sm + W_mscheme(p,mp,q,mq,r,mr,s,ms,t,mt,u,mu,Mew,AA,BB,jbas)&
                                  * dcgi(jp,mp,jq,mq,Jpq,Mpq) * dcgi(Jpq,Mpq,jr,mr,j1,m1) &
                                  * dcgi(js,ms,jt,mt,Jst,Mst) * dcgi(Jst,Mst,ju,mu,j2,m2) &
                                  * dcgi(j2,m2,RANK,Mew,j1,m1)/sqrt(j1+1.d0) 
                             
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do 
 end do
  W_via_mscheme = sm
  
end function W_via_mscheme
!================================================
!================================================
integer function read_eom_file(trs,mom,eom_states,jbas)
  ! read the inputs for the EOM calculation 
  implicit none

  type(spd) :: jbas
  type(eom_mgr) :: eom_states
  type(obsv_mgr) :: trs,mom
  character(2) :: op,init,fin
  integer :: ist,num_trans,num_mom,i,num_jpi,N ,uniq
  integer :: totstates
  logical :: isthere
  
  N =jbas%total_orbits
  if (trim(INI_DIR) == './' ) then
     inquire(file=trim(INI_DIR)//trim(eomfile),exist=isthere)
     if (isthere) then 
        open(unit=44,file=trim(INI_DIR)//trim(eomfile))
     else
        open(unit=44,file='../../inifiles/'//trim(eomfile))        
     end if
  else
     open(unit=44,file=trim(INI_DIR)//trim(eomfile))
  end if
  read(44,*);read(44,*);read(44,*)
  ! read transition types
  read(44,*) num_jpi

  eom_states%num = num_jpi
  allocate(eom_states%name(num_jpi))
  allocate(eom_states%dTz(num_jpi)) 
  allocate(eom_states%ang_mom(num_jpi))
  allocate(eom_states%par(num_jpi))
  allocate(eom_states%number_requested(num_jpi))
  totstates = 0 
  read(44,*)
  do i = 1, num_jpi
     read(44,*) eom_states%name(i),eom_states%dTz(i),eom_states%number_requested(i)
     read(eom_states%name(i)(1:1),'(I1)') eom_states%ang_mom(i)
     eom_states%ang_mom(i) =  eom_states%ang_mom(i) *2
     if (eom_states%name(i)(2:2) == '+') then
        eom_states%par(i) = 0
     else
        eom_states%par(i) = 2
     end if
     totstates =totstates + eom_states%number_requested(i) 
     if (eom_states%dTz(i).ne.0) then
        eom_states%total_dtz = eom_states%total_dtz +  eom_states%number_requested(i) 
     end if
  end do

  uniq = num_jpi + 2 ! plus one for operator, one for Hcm.
  ! too lazy to check if the operator has the same structure as other stuff

  read(44,*)
  read(44,*) op
  trs%oper= op
  mom%oper= op
  read(44,*)

  read(44,*) num_trans
  trs%num = num_trans

  allocate(trs%Jpi1(num_trans))
  allocate(trs%dTz(num_trans)) 
  allocate(trs%Jpi2(num_trans))

  read(44,*)
  do i=1,num_trans
     read(44,*) trs%Jpi1(i),trs%Jpi2(i),trs%dtz(i)
  end do 

  read(44,*)
  read(44,*) num_mom
  mom%num = num_mom

  allocate(mom%Jpi1(num_mom))
  allocate(mom%dTz(num_trans)) 
  read(44,*)
  do i=1,num_mom
     read(44,*) mom%Jpi1(i),mom%dtz(i)
  end do 

  allocate(jbas%xmap_tensor(uniq,N*(N+1)/2))
  allocate(half6j(uniq)) 
  read_eom_file = totstates

end function read_eom_file

end module EOM_IMSRG

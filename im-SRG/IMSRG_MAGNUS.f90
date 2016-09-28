module IMSRG_MAGNUS
  use commutators
  use TS_commutators
  use operators
  use HF_mod 
  implicit none 

  !interface such that dTZ operators and sq_op use the same function
  interface transform_observable_BCH
     module procedure BCH_ISOTENSOR,transform_sq_op_BCH
  end interface transform_observable_BCH
  
contains

subroutine magnus_decouple(HS,G,jbas,quads,trips,build_generator)
  ! runs IMSRG using magnus expansion method
  implicit none 
  
  real(8),parameter :: conv_crit = 1d-6
  type(spd) :: jbas
  type(sq_op) :: H,G,HS,AD
  type(sq_op) :: DG,CR
  type(cc_mat) :: GCC,ADCC,WCC 
  real(8) :: ds,s,sx,Eold,E_mbpt2,crit,nrm1,nrm2,dcgi00
  real(8) :: omp_get_wtime,t1,t2
  character(1) :: quads,trips
  integer :: Atot,Ntot,nh,np,nb,q,steps,i,j,chk
  logical :: trip_calc,xxCR,not_chkpoint_restart,first

  ! dummy variable for function to build generator
  external :: build_generator 

  ! need to make some new objects to hold auxillary operators 
  HS%neq = 1
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  if (.not. allocated(G%mat)) call duplicate_sq_op(HS,G) !magnus operator
  call duplicate_sq_op(HS,AD) !workspace

  call duplicate_sq_op(HS,DG) !magnus operator
  G%herm = -1 
  DG%herm = -1 

  ! construct generator of the transformation
  call build_generator(HS,DG,jbas)  
  call copy_sq_op(HS,H) 

  ! start parameters
  s = 0.d0 
  ds = 0.5d0
  steps = 0

  ! relevant for checkpointing
  ! if we are checkpointing we write
  ! omega to file every "chk" steps and then resubmit the job. 
  if (HS%eMax.ge.14) then
     ! eMax=14 is the largest I've done, and it's pretty demanding.
     ! this takes 10 to 15 hours depending on a variety of factors
     chk = 8
  else
     chk = 24
  end if 
  
  if (checkpointing) then 
     not_chkpoint_restart = read_omega_checkpoint(G,sx) 
  else
     not_chkpoint_restart = .true. 
  end if 
  
  if (not_chkpoint_restart) then 
  
     ! start from the beginning.
     open(unit=36,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_0b_magnus_flow.dat')
     
     write(36,'(I6,4(e15.7))') steps,s,H%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     Eold=0.
     E_mbpt2 = mbpt2(HS,jbas) ! many body perturbation theory at second order
     crit=abs(E_mbpt2) !convergence is based on the magnitude of MBPT(2), which is resummed into E_HF,
     first = .false.  
  else 

     ! CHECKPOINT RESTART
     open(unit=36,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_0b_magnus_flow.dat',position='append')

     ! trnasform to where we left off 
     call BCH_EXPAND(HS,G,H,jbas,quads) 
     call build_generator(HS,DG,jbas)  

     ! advance to the current value of s 
     s= sx
  
     steps = nint(sx/ds)
     first = .true. 
  end if 
  

  do while (crit > conv_crit) 
 
     steps = steps + 1

     ! if checkpointing and you reach a checkpoint, write and quit. 
     if (checkpointing) then 
        if (mod(steps,chk)==0) then 
           if (first) then
              ! don't write and quit at the first iteration.               
              first = .false.
           else
              call write_omega_checkpoint(G,s)
           end if
        end if
     end if 

     ! magnus expansion to calculate transformation derivative 
     call MAGNUS_EXPAND(DG,G,AD,jbas)

     ! euler solve the magnus differential equation
     call euler_step(G,DG,s,ds) 

     ! apply the current transformation to H 
     call BCH_EXPAND(HS,G,H,jbas,quads) 

     ! update generator 
     call build_generator(HS,DG,jbas)   

     
     E_mbpt2 = mbpt2(HS,jbas) 
     crit = abs(E_mbpt2)  

     write(36,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     
  end do

  ! HS is now the transformed hamiltonian
  ! exp_omega now holds the IM-SRG unitary transformation. 
end subroutine magnus_decouple
!=========================================================================
!=========================================================================
subroutine transform_sq_op_BCH(Op,G,jbas,quads)
  ! BCH expansion for some operator
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: Op,G,Oevolved
  character(1) :: quads
  
  if (Op%rank > 0) then      
     call BCH_TENSOR(G,Op,jbas,quads)    
  else
     call duplicate_sq_op(Op,Oevolved) 
     call BCH_EXPAND(Oevolved,G,Op,jbas,quads)
     call copy_sq_op(Oevolved,Op) 
  end if
  
end subroutine transform_sq_op_BCH
!=========================================================================
!=========================================================================
subroutine BCH_EXPAND(HS,G,H,jbas,quads) 
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) :: adnorm,fullnorm,s,advals(20),sm,sm2,coef
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 


  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  call init_ph_mat(AD,ADCC,jbas) !cross coupled ME
  call duplicate_ph_mat(ADCC,GCC) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME

  advals = 0.d0 
  coef = 1.d0 

  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.
  if (quads=='y') call duplicate_sq_op(H,INT3)

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(H%E0)   

  do iw = 2 ,20

     coef = coef/(iw-1.d0) 
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)
     !now: INT2 = [ G , AD ]  
        
! zero body commutator 
     call calculate_cross_coupled(AD,ADCC,jbas)
     call calculate_cross_coupled(G,GCC,jbas) 
  
     INT2%E0 = commutator_110(G,AD,jbas) + commutator_220(G,AD,jbas)
     
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    
          
     call commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)

     call commutator_221(G,AD,INT2,w1,w2,jbas)
     call commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)
     
     ! Here is where you add perturbative quadrupoles if desired 
     if (quads=='y') then 
        call add_sq_op(INT3, 1.d0 , INT2, 1.d0, INT2)          
        call restore_quadrupoles(AD,G,w1,w2,INT3,jbas) 
     end if 

     ! so now just add INT1 + c_n * INT2 to get current value of HS
     call add_sq_op(INT1 ,1.d0 , INT2 , coef , HS )   !basic_IMSRG
          
     advals(iw) = mat_frob_norm(INT2)*coef
     if (advals(iw) < conv) exit
  end do 
  if (iw > 20) STOP 'clipping BCH'
  
end subroutine BCH_EXPAND
!=========================================================================
!=========================================================================
subroutine BCH_EXPAND_1b(HS,G,H,jbas,quads) 
  ! for 1-body transformations
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) :: adnorm,fullnorm,s,advals(15),sm,sm2,coef
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 

  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1

  advals = 0.d0 
  coef = 1.d0 

  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(H%E0)   

  do iw = 2 ,30
     coef = coef/(iw-1.d0) 
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)    
     !now: INT2 = [ G , AD ]  
        
! zero body commutator 
     INT2%E0 = commutator_110(G,AD,jbas)
     
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    
     
     ! so now just add INT1 + c_n * INT2 to get current value of HS
     call add_sq_op(INT1 ,1.d0 , INT2 , coef , HS )   !basic_IMSRG
          
     advals(iw) = mat_frob_norm(INT2)*coef
     if (advals(iw) < conv) exit
     
  end do 
  
end subroutine BCH_EXPAND_1b
!=========================================================================
!=========================================================================
subroutine BCH_TENSOR(G,HS,jbas,quads) 
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw,omp_get_num_threads
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT,threads
  type(spd) :: jbas
  type(sq_op) :: G, ETA, INT2, INT3,HS, AD,w1,w2
  type(pandya_mat) :: WCC,ADCC
  type(cc_mat) :: GCC 
  real(8) ::  coef,adnorm,fullnorm,s,advals(30),sm,sm2,dcgi,dcgi00
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 
  
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  INT2%herm = 1
  AD%herm = 1
  call allocate_small_tensor_CCMAT(AD,ADCC,jbas)
  call init_ph_mat(G,GCC,jbas) !cross coupled ME
  
  advals = 0.d0 
  coef = 1.d0
  ! intermediates must be HERMITIAN
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.

  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(HS%E0)   
  do iw = 2 ,30
     
     coef = coef/(iw-1.d0)
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)    
!     call clear_sq_op(INT3)    
     !now: INT2 = [ G , AD ]  
        
     call calculate_cross_coupled(G,GCC,jbas)      

     call TS_commutator_111(G,AD,INT2,jbas) 
     call TS_commutator_121(G,AD,INT2,jbas)      

     call TS_commutator_122(G,AD,INT2,jbas)   
     call TS_commutator_212(G,AD,INT2,jbas)  
     call TS_commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)
     call TS_commutator_221(w1,w2,G%herm*AD%herm,INT2,jbas)

     call TS_commutator_211(GCC,AD,INT2,jbas)
     call TS_commutator_222_ph(GCC,ADCC,AD,INT2,jbas)       

     ! so now just add HS + c_n * INT2 to get current value of HS
  !   call append_operator( INT3 , 1.d0 , INT2 )   !basic_IMSRG
     call append_operator( INT2 , coef , HS )   !basic_IMSRG
     
         
     advals(iw) = mat_frob_norm(INT2)*coef
    if (advals(iw) < conv) exit
     
  end do
 
end subroutine BCH_TENSOR
!=========================================================================
!=========================================================================
subroutine BCH_ISOTENSOR(HS,G,jbas,quads) 
  use operator_commutators
  use isospin_operators
  implicit none 
  
  real(8), parameter :: conv = 1e-4
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw,omp_get_num_threads
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT,threads
  type(spd) :: jbas
  type(sq_op) :: G
  type(iso_operator) :: INT2, INT3,HS, AD 
  type(pandya_mat) :: WCC,ADCC
  type(cc_mat) :: GCC 
  real(8) ::  coef,adnorm,fullnorm,s,advals(30),sm,sm2,dcgi,dcgi00
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 

  call duplicate_isospin_operator(HS,INT2) !workspace
 ! call duplicate_sq_op(HS,INT3) !workspace
  call duplicate_isospin_operator(HS,AD) !workspace
  INT2%herm = 1
  AD%herm = 1
  
  advals = 0.d0 
  coef = 1.d0
  ! intermediates must be HERMITIAN
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.

  call copy_isospin_operator( HS , INT2 )
 
  advals(1) = abs(HS%E0)   
  do iw = 2 ,30
     
     coef = coef/(iw-1.d0)
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_isospin_operator( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_isospin_operator(INT2)    
!     call clear_sq_op(INT3)    
     !now: INT2 = [ G , AD ]  

     call operator_commutator_111(G,AD,INT2,jbas) 
     call operator_commutator_211(G,AD,INT2,jbas)
     call operator_commutator_121(G,AD,INT2,jbas)      
     call operator_commutator_122(G,AD,INT2,jbas)   
     call operator_commutator_212(G,AD,INT2,jbas)  
     call operator_commutator_222_pp_hh(G,AD,INT2,jbas)    
     call operator_commutator_222_ph(G,AD,INT2,jbas)       
     ! so now just add HS + c_n * INT2 to get current value of HS
     !   call append_operator( INT3 , 1.d0 , INT2 )   !basic_IMSRG
    
     call append_isospin_operator( INT2 , coef , HS )   !basic_IMSRG
              
     advals(iw) = iso_frob_norm(INT2)*coef
     if (advals(iw) < conv) exit
     
  end do
 
end subroutine BCH_ISOTENSOR
!=========================================================================
!=========================================================================
subroutine CR_EXPAND(HS,G,H,jbas,quads) 
  ! this is the BCH expansion less one term, where the number of
  ! nested commutators is reducted by 1, This allows us to construct
  ! the induced three-body force by taking one last commutator [exp_omega,HS]_{223} 
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) ::  coef,adnorm,fullnorm,s,advals(14),sm,sm2,dcgi,dcgi00
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 

  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  call init_ph_mat(AD,ADCC,jbas) !cross coupled ME
  call duplicate_ph_mat(ADCC,GCC) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME


  advals = 0.d0 
  
  coef = 1.d0
  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.
  if (quads=='y') call duplicate_sq_op(H,INT3)

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(H%E0)   

  do iw = 2 ,14
     coef= coef/float(iw)
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)    
     !now: INT2 = [ G , AD ]  
        
! zero body commutator
 
     call calculate_cross_coupled(AD,ADCC,jbas)
     call calculate_cross_coupled(G,GCC,jbas) 
 
     INT2%E0 = commutator_110(G,AD,jbas) + commutator_220(G,AD,jbas)

     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    

     call commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)
  
     call commutator_221(G,AD,INT2,w1,w2,jbas)
     call commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)

     ! so now just add INT1 + c_n * INT2 to get current value of HS
     if (quads=='y') then
        call add_sq_op(INT3, 1.d0 , INT2, 1.d0, INT2)        
        call restore_quadrupoles(AD,G,w1,w2,INT3,jbas) 
     end if           
     
     call add_sq_op(INT1 ,1.d0 , INT2 , coef , HS )   !basic_IMSRG
         
     advals(iw) = mat_frob_norm(INT2)*coef
     if (advals(iw) < conv) exit
     
  end do 
 
end subroutine 
!===============================================================
!===============================================================
subroutine MAGNUS_EXPAND(DG,G,AD,jbas)
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,q,j,k,l,ry
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2,DG
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(7),adnorm,fullnorm,s,advals(7) 
  character(3) :: args
  
  call duplicate_sq_op(DG,w1) !workspace
  call duplicate_sq_op(DG,w2) !workspace
  call duplicate_sq_op(DG,INT1) !workspace
  call duplicate_sq_op(DG,INT2) !workspace
  INT2%herm = -1
  INT1%herm = -1 
  AD%herm = -1
  
  call init_ph_mat(AD,ADCC,jbas) !cross coupled ME
  call duplicate_ph_mat(ADCC,GCC) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME

  advals = 0.d0 
  ! Intermediates are ANTI-HERMITIAN 

  ! Bernoulli numbers
  cof = (/1.d0,-0.5d0,.0833333333d0,0.d0,-0.00138888888d0,0.d0,3.306878d-5/) 
  ! nested derivatives not so important
     
  !! same deal as BCH expansion, which is explained ad nauseam above. 
   
  call copy_sq_op( DG , INT2 )
  advals(1) = mat_frob_norm(INT2)  
  
  fullnorm = mat_frob_norm(G)   
  if (fullnorm < 1e-9) return
  
  q = 1
 ! return
  do i = 2 , 7

     call copy_sq_op( DG , INT1) 
     call copy_sq_op( INT2 , AD ) 
  
     adnorm = advals(i-q) 
  
     if  (abs(cof(i)) > 1e-6) then  
        q = 1 
     else 
        q = 2
     end if 
     
     if ( abs(adnorm/fullnorm) < conv ) exit
          
     call calculate_cross_coupled(AD,ADCC,jbas)
     call calculate_cross_coupled(G,GCC,jbas) 
 
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    
  
     call commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)   
     call commutator_221(G,AD,INT2,w1,w2,jbas)     
     call commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)
                 
     call add_sq_op(INT1 , 1.d0 , INT2 , cof(i) , DG ) !ME_general
     
     advals(i) = mat_frob_norm(INT2)*abs(cof(i))
    
  end do   
     
end subroutine 
!=====================================================
!=====================================================
subroutine MAGNUS_EXPAND_1b(DG,G,AD,jbas)
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,q,j,k,l,ry
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2,DG
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(7),adnorm,fullnorm,s,advals(7) 
  character(3) :: args
  
  call duplicate_sq_op(DG,INT1) !workspace
  call duplicate_sq_op(DG,INT2) !workspace
  INT2%herm = -1
  INT1%herm = -1 
  AD%herm = -1
  
  advals = 0.d0 
  ! Intermediates are ANTI-HERMITIAN 
  
  cof = (/1.d0,-0.5d0,.0833333333d0,0.d0,-0.00138888888d0,0.d0,3.306878d-5/) 
  ! nested derivatives not so important
     
  !! same deal as BCH expansion, which is explained ad nauseam above. 
   
  call copy_sq_op( DG , INT2 )
  advals(1) = mat_frob_norm(INT2)  
  
  fullnorm = mat_frob_norm(G)   
  if (fullnorm < 1e-9) return
  
  q = 1
 ! return
  do i = 2 , 7
     call copy_sq_op( DG , INT1) 
     call copy_sq_op( INT2 , AD ) 
  
     adnorm = advals(i-q) 
  
     if  (abs(cof(i)) > 1e-6) then  
        q = 1 
     else 
        q = 2
     end if 
     
     if ( abs(adnorm/fullnorm) < conv ) exit
  
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    
                 
     call add_sq_op(INT1 , 1.d0 , INT2 , cof(i) , DG ) !ME_general
     
     advals(i) = mat_frob_norm(INT2)*abs(cof(i))
    
  end do   
     
end subroutine MAGNUS_EXPAND_1b
!=====================================================
!=====================================================
subroutine euler_step(G,DG,s,stp)
  ! forward euler step
  implicit none 

  integer :: i
  real(8) :: s , stp
  type(sq_op) :: G , DG

  call add_sq_op(G,1.d0,DG,stp,G) !durr probably wrong. 
  s = s + stp 

end subroutine 
!=====================================================
!=====================================================
real(8) function restore_triples(H,OM,jbas) 
  ! this is the MBPT(2) correction for a three-body force.
  ! the actual correction comes in at fourth order, as the
  ! induced three-body force is 2nd order in MBPT. 
  implicit none 
  
  type(spd) :: jbas
  type(tpd),allocatable,dimension(:) :: threebas
  type(sq_op) :: H,OM
  integer :: a,b,c,i,j,k,Jtot,Jab,Jij,g1
  integer :: ja,jb,jc,ji,jj,jk,AAA,q
  integer :: ax,bx,cx,ix,jx,kx,III
  integer :: jab_min,jab_max,jij_min,jij_max
  integer :: J_min, J_max,x,total_threads,thread
  real(8) :: faa,fbb,fcc,fii,fjj,fkk,Gabab,Gkbkb,Gkckc
  real(8) :: Gacac,Gbcbc,Gijij,Gikik,Gjkjk,Giaia
  real(8) :: Gibib,Gicic,Gjaja,Gjbjb,Gjcjc,Gkaka    
  real(8) :: sm,denom,dlow,w,pre1,pre2
  
  
  sm = 0.d0   
  call enumerate_three_body(threebas,jbas)  
  total_threads = size(threebas(1)%direct_omp) - 1

  !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(threebas,jbas,H,OM) & 
  !$OMP& REDUCTION(+:sm)  
  
  do thread = 1, total_threads
  do q = 1+threebas(1)%direct_omp(thread),&
       threebas(1)%direct_omp(thread+1)
  
     Jtot = threebas(q)%chan(1) 
     do AAA = 1, size(threebas(q)%ppp(:,1)) 

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
        
        do III = 1, size(threebas(q)%hhh(:,1)) 
             
           i = threebas(q)%hhh(III,1)
           j = threebas(q)%hhh(III,2)
           k = threebas(q)%hhh(III,3)

           ! prefactors to to account for
           ! double counting. 
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
  
           ! THIS IS THE EPSTEIN NESBET DENOMINATOR     
           denom = -1*(faa+fbb+fcc-fii-fjj-fkk+Gabab+&
                Gacac+Gbcbc+Gijij+Gikik+Gjkjk-Giaia&
                -Gibib-Gicic-Gjaja-Gjbjb-Gjcjc-Gkaka-&
                Gkbkb-Gkckc)*pre1*pre2 

           do jab = jab_min,jab_max,2
               
              if ( .not. (triangle(Jtot,jc,Jab))) cycle
              if ((a==b) .and. (mod(Jab/2,2)==1)) cycle               
              do jij = jij_min, jij_max,2
                 
                 if ( .not. (triangle(Jtot,jk,Jij))) cycle
                 if ((i==j) .and. (mod(Jij/2,2)==1)) cycle
                 w = commutator_223_single(OM,H,a,b,c,i,j,k,Jtot,jab,jij,jbas)
                 sm = sm + w*w/denom*(Jtot+1.d0)
              end do
           end do

        end do
     end do
  end do
  end do 
 !$OMP END PARALLEL DO 
  restore_triples = sm

end function restore_triples
!================================================
!================================================
end module
!================================================
!================================================
subroutine restore_quadrupoles( X , OM, w1,w2, RES,jbas ) 
  ! perturbative quadrupoles restoration, these are best included
  ! on the fly in the BCH expansion, and this routine costs NOTHING
  ! compared with the triples. In fact the scaling is the same as the IM-SRG calculation. 
  use basic_IMSRG
  implicit none
  
  type(sq_op) :: X,OM,INT1,INT2,RES
  type(spd) :: jbas
  type(sq_op) :: w1,w2
  integer :: p1x,p2x,h1x,h2x
  integer :: i,j,b,q,Abody,Ntot,nh,np,nb,a,c,p1,p2,h1,h2,cx
  integer :: AA,II,Jtot,ik,jk,ck,ji,jj,ti,tj,li,lj,jc,JT,pm,jp1,jp2,ja,jb
  integer :: kx,dx,k,d
  real(8) :: sm,sm2
  
  Abody = X%belowEF
  Ntot = X%Nsp
  
  ! I think this is all I need here (these are the one -body insertions)
  allocate(INT1%fpp(Ntot-Abody,Ntot-Abody),INT1%fhh(Abody,Abody)) 
  INT1%fpp = 0.d0
  INT1%fhh = 0.d0 
 
  allocate(INT2%fpp(Ntot-Abody,Ntot-Abody),INT2%fhh(Abody,Abody)) 
  INT2%fpp = 0.d0
  INT2%fhh = 0.d0
 
  
  pm = X%herm*OM%herm

  ! w1=X.OM w2=OM.OM matrices are made 
  call build_intermediates_For_intermediates(X,OM,w1,w2,jbas)

! fhh
  do i = 1 , Abody
     ik = jbas%holes(i) 
     ji = jbas%jj(ik) 
     li = jbas%ll(ik) 
     ti = jbas%itzp(ik) 
     
     do j = 1 , Abody
        
        jk = jbas%holes(j) 
        jj = jbas%jj(jk) 
        if (ji .ne. jj)  cycle
        lj = jbas%ll(jk) 
        if (li .ne. lj)  cycle
        tj = jbas%itzp(jk)
        if (tj .ne. ti) cycle 
                
        sm = 0.d0 
        sm2 = 0.d0 
        do c = 1, Abody
           ck = jbas%holes(c) 
           jc = jbas%jj(ck)
           do JT = abs(jc - ji),jc+ji,2
              sm = sm + v_elem(ck,ik,ck,jk,JT,w1,jbas)*(JT + 1)
              sm2 = sm2 + v_elem(ck,ik,ck,jk,JT,w2,jbas)*(JT + 1) 
              ! w2 is subtracted, because the commutator in this case has a minus sign. 
           end do
        end do 
           
        INT1%fhh(i,j) = INT1%fhh(i,j) + sm / (ji + 1.d0 )
        INT2%fhh(i,j) = INT2%fhh(i,j) + sm2 / (ji + 1.d0 )
       ! nothing is hermitian or anti-hermitian here
     end do 
  end do       

! fpp
  do i = 1 , Ntot - Abody
     ik = jbas%parts(i) 
     ji = jbas%jj(ik) 
     li = jbas%ll(ik) 
     ti = jbas%itzp(ik) 
     
     do j = i , Ntot - Abody
        
        jk = jbas%parts(j) 
        jj = jbas%jj(jk) 
        if (ji .ne. jj)  cycle
        lj = jbas%ll(jk) 
        if (li .ne. lj)  cycle
        tj = jbas%itzp(jk)
        if (tj .ne. ti) cycle 
                
       
        sm = 0.d0 
        sm2 = 0.d0
        do c = 1, Ntot - Abody
           ck = jbas%parts(c) 
           jc = jbas%jj(ck)
           ! w2 matrix results from multiplying the hh channel
           do JT = abs(jc - ji),jc+ji,2
              ! the hermiticity of X and OM are exploited here. 
              sm = sm + X%herm*OM%herm*v_elem(ck,jk,ck,ik,JT,w1,jbas)*(JT + 1) 
              sm2 = sm2 + v_elem(ck,jk,ck,ik,JT,w2,jbas)*(JT + 1)
           end do 
        end do 
     
        INT1%fpp(i,j) = INT1%fpp(i,j) + sm / (ji + 1.d0) 
        INT2%fpp(i,j) = INT2%fpp(i,j) + sm2 / (ji + 1.d0) 
        
     end do 
  end do       

  !!! now add the new stuff to RES
  
  do q = 1, RES%nblocks
     
     JTot = RES%mat(q)%lam(1) 
     
     np = RES%mat(q)%npp
     nh = RES%mat(q)%nhh
     
     do AA = 1, np
        a = RES%mat(q)%qn(1)%Y(AA,1)
        b = RES%mat(q)%qn(1)%Y(AA,2)
        ja = jbas%jj(a)
        jb = jbas%jj(b)
      
        do II = 1, nh    
           i = RES%mat(q)%qn(3)%Y(II,1)
           j = RES%mat(q)%qn(3)%Y(II,2)
           ji = jbas%jj(i)
           jj = jbas%jj(j)
           
           RES%mat(q)%gam(3)%X(AA,II) = 0.d0
           do c = 1, Abody
              
              cx = jbas%holes(c) 
              h1 = hb4(i)+1 
              h2 = hb4(j)+1
              
              
              ! i've already built the INT2 minus signs into the one-body
              ! insertion
              RES%mat(q)%gam(3)%X(AA,II) = RES%mat(q)%gam(3)%X(AA,II) &
                   - INT1%fhh(h1,c) * v_elem(a,b,cx,j,Jtot,OM,jbas)  &
                   + (-1)**((ji+jj-Jtot)/2)*INT1%fhh(h2,c) * v_elem(a,b,cx,i,Jtot,OM,jbas) &
                   + INT2%fhh(h1,c) * v_elem(a,b,cx,j,Jtot,X,jbas) &
                   - (-1)**((ji+jj-Jtot)/2)*INT2%fhh(h2,c) * v_elem(a,b,cx,i,Jtot,X,jbas)
                   
          end do 
       
          do c = 1, Ntot-Abody
              
              cx = jbas%parts(c) 
              p1 = pb4(a)+1 
              p2 = pb4(b)+1
              
              
              ! i've already built the INT2 minus signs into the one-body
              ! insertion
               RES%mat(q)%gam(3)%X(AA,II) = RES%mat(q)%gam(3)%X(AA,II) &
                   - INT1%fpp(p1,c) * v_elem(cx,b,i,j,Jtot,OM,jbas) &
                    + (-1)**((ja+jb-Jtot)/2)* INT1%fpp(p2,c) * v_elem(cx,a,i,j,Jtot,OM,jbas) &
                    + INT2%fpp(p1,c) * v_elem(cx,b,i,j,Jtot,X,jbas) &
                    -(-1)**((ja+jb-Jtot)/2)*INT2%fpp(p2,c) * v_elem(cx,a,i,j,Jtot,X,jbas)
                   
          end do 
       end do 
    end do 
 end do 

end subroutine restore_quadrupoles
!========================================================================
!========================================================================
subroutine build_intermediates_For_intermediates(L,R,w1,w2,jbas) 
  !VERIFIED
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediate matrices for the calculation of 
  ! one body insertions in restore_quadrupoles
  use basic_IMSRG
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  L,R,w1,w2
  integer :: q
  integer :: np,nb,nh,pm
  real(8) :: bet_off,al_off
  
  pm = R%herm*L%herm
   !construct temporary matrices
  do q = 1, L%nblocks
     
     nh = L%mat(q)%nhh
     np = L%mat(q)%npp
     nb = L%mat(q)%nph
        
     if (nh + np == 0 ) cycle 
     if (np + nb == 0 ) cycle 
     if (nh + nb == 0 ) cycle
     
  
             
     if (np*nh .ne. 0) then 
     !L_hhpp . R_pphh = W1_hhhh
     call dgemm('T','N',nh,nh,np,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(5)%X,nh)
    
     w1%mat(q)%gam(5)%X = w1%mat(q)%gam(5)%X * L%herm 
     
     !L_pphh . R_hhpp = W1_pppp
     call dgemm('N','T',np,np,nh,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(1)%X,np)
     
     w1%mat(q)%gam(1)%X = w1%mat(q)%gam(1)%X * R%herm
     
     
     !R_hhpp . R_pphh = W2_hhhh
     call dgemm('T','N',nh,nh,np,al,R%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w2%mat(q)%gam(5)%X,nh)
    
     w2%mat(q)%gam(5)%X = w2%mat(q)%gam(5)%X * R%herm 
     
     !R_pphh . R_hhpp = W2_pppp
     call dgemm('N','T',np,np,nh,al,R%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w2%mat(q)%gam(1)%X,np)
     
     w2%mat(q)%gam(1)%X = w2%mat(q)%gam(1)%X * R%herm
     
     end if
    
    
  end do

end subroutine build_intermediates_For_intermediates

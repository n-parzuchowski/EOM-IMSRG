module IMSRG_CANONICAL
  use commutators
  use operators
  use HF_mod 
  implicit none 
  
contains

subroutine discrete_decouple(HS,jbas,O1,O2) 
  ! runs discrete canonical transformation 
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps,i,j
  type(spd) :: jbas
  type(sq_op),optional :: O1,O2
  type(sq_op) :: H , G ,ETA, HS,INT1,INT2,AD,w1,w2,DG,Oevolv
  type(cc_mat) :: GCC,ADCC,WCC 
  real(8) :: ds,s,E_old,E_mbpt2,crit,nrm1,nrm2,wTs(2),Ecm(3),oldcrit
  
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  
  call init_ph_mat(HS,ADCC,jbas) !cross coupled ME
  call duplicate_ph_mat(ADCC,GCC) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME
  
  call copy_sq_op(HS,H) 
  
  s = 0.d0 
  ds = 1.d0
  
  crit = 10.
  steps = 0

  open(unit=36,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_0b_discrete_flow.dat')

if (present(O1)) then 
   if (present(O2)) then 
! two operators   
  do while (crit > 1e-6) 
 
    ! call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2) 
     call build_gs_white(HS,ETA,jbas) 
     call scale_sq_op(ETA,ds) 

     call BCH_EXPAND(HS,ETA,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
    
     call BCH_EXPAND(H,ETA,O2,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     call copy_sq_op(H,O2)

     call BCH_EXPAND(H,ETA,O1,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s)
     call copy_sq_op(H,O1)
    
     call copy_sq_op(HS,H)
     E_mbpt2 = mbpt2(HS,jbas) 

     crit = abs(E_mbpt2) 
     
     s = s + ds 
     steps = steps + 1
     
     write(36,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     
  end do   
else
! one operator
  do while (crit > 1e-6) 
 
    ! call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2) 
     call build_gs_white(HS,ETA,jbas) 
     call scale_sq_op(ETA,ds) 

     call BCH_EXPAND(HS,ETA,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
    
     call BCH_EXPAND(H,ETA,O1,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     call copy_sq_op(H,O1)
    
     call copy_sq_op(HS,H)
     E_mbpt2 = mbpt2(HS,jbas) 

     crit = abs(E_mbpt2) 
     
     s = s + ds 
     steps = steps + 1
     
     write(36,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     
  end do 
end if 
else
! no operators
  do while (crit > 1e-6) 
 
    ! call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2) 
     call build_gs_white(HS,ETA,jbas) 
     call scale_sq_op(ETA,ds) 

     call BCH_EXPAND(HS,ETA,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
      
     call copy_sq_op(HS,H)
     E_mbpt2 = mbpt2(HS,jbas) 

     oldcrit = crit
     crit = abs(E_mbpt2) 
     
     s = s + ds 
     steps = steps + 1

     write(36,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     
  end do 
  
end if 

  close(36)
end subroutine  


subroutine discrete_TDA( HS , TDA, jbas,O1,O1TDA,O2,O2TDA ) 
  ! runs IMSRG TDA decoupling using discrete transformation
  implicit none 

  ! SRG convergence / failsafe / error tolerances
  integer,parameter :: max_steps = 50
  real(8),parameter :: conv_crit = 1.d-5
  real(8),parameter :: relerr = 1.d-5, abserr = 1.d-5

  type(spd) :: jbas
  type(sq_op) :: HS,H,ETA,w1,w2,AD,INT1,INT2
  type(sq_op),optional :: O1,O2 
  type(cc_mat) :: HCC,OeCC,ADCC,GCC,WCC
  type(full_sp_block_mat) :: TDA
  type(full_sp_block_mat),optional :: O1TDA,O2TDA
  integer,dimension(5) :: iwork
  real(8),allocatable,dimension(:) :: cur_vec,work,E_old
  integer :: neq,iflag,Atot,Ntot,nh,np,nb,q,steps ,i 
  real(8) :: ds,s,crit,min_crit  
  
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  
  call init_ph_mat(HS,ADCC,jbas) !cross coupled ME
  call duplicate_ph_mat(ADCC,GCC) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME
  
  call copy_sq_op(HS,H) 
  
  s = 0.d0 
  ds = 0.001d0
  
  crit = 10.
  steps = 0

  call duplicate_ph_mat(ADCC,HCC) 
  HCC%herm = H%herm
  allocate(E_old(TDA%map(1)))
  call calculate_cross_coupled(H,HCC,jbas)
  call calc_TDA(TDA,H,HCC,jbas) 
  call diagonalize_blocks(TDA)

  E_old = TDA%blkM(1)%eigval
  open(unit=37,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_excited.dat')

  call write_excited_states(steps,s,TDA,H%E0,37) 
  
min_crit = 10000.d0
if (present(O1)) then 
   if (present(O2)) then 
! two operators   
  do while (crit > 1e-6) 
 
    ! call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2) 
     call build_specific_space(HS,ETA,jbas) 
     call scale_sq_op(ETA,ds) 

     call BCH_EXPAND(HS,ETA,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
    
     call BCH_EXPAND(H,ETA,O2,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     call copy_sq_op(H,O2)

     call BCH_EXPAND(H,ETA,O1,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s)
     call copy_sq_op(H,O1)
    
     call copy_sq_op(HS,H)

     s = s + ds 
     steps = steps + 1
  
     call calculate_cross_coupled(H,HCC,jbas)
     call calc_TDA(TDA,H,HCC,jbas) 
     call diagonalize_blocks(TDA)
  
     call write_excited_states(steps,s,TDA,H%E0,37) 
     
     ! convergence criteria
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(I4,7(e14.5))') steps,crit,E_old(1:6) 
     E_old = TDA%blkM(1)%eigval
 

     if (crit > 100*min_crit) then
        print*, 'convergence failed' 
        exit
     end if 
     min_crit = min(min_crit,crit) 
     if (crit < conv_crit) exit
     
  end do   
else
! one operator
  do while (crit > 1e-6) 
 
    ! call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2) 
     call build_specific_space(HS,ETA,jbas) 
     call scale_sq_op(ETA,ds) 

     call BCH_EXPAND(HS,ETA,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
    
     call BCH_EXPAND(H,ETA,O1,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     call copy_sq_op(H,O1)
    
     call copy_sq_op(HS,H)
     s = s + ds 
     steps = steps + 1
  
     call calculate_cross_coupled(H,HCC,jbas)
     call calc_TDA(TDA,H,HCC,jbas) 
     call diagonalize_blocks(TDA)
  
     call write_excited_states(steps,s,TDA,H%E0,37) 
     
     ! convergence criteria
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(I4,7(e14.5))') steps,crit,E_old(1:6) 
     E_old = TDA%blkM(1)%eigval
 

     if (crit > 100*min_crit) then
        print*, 'convergence failed' 
        exit
     end if 
     min_crit = min(min_crit,crit) 
     if (crit < conv_crit) exit
      
  end do 
end if 
else
! no operators
  do while (crit > 1e-6) 
 
    ! call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2) 
     call build_specific_space(HS,ETA,jbas) 
     call scale_sq_op(ETA,ds) 

     call BCH_EXPAND(HS,ETA,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
      
     call copy_sq_op(HS,H)
     s = s + ds 
     steps = steps + 1
  
     call calculate_cross_coupled(H,HCC,jbas)
     call calc_TDA(TDA,H,HCC,jbas) 
     call diagonalize_blocks(TDA)
  
     call write_excited_states(steps,s,TDA,H%E0,37) 
     
     ! convergence criteria
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(I4,7(e14.5))') steps,crit,E_old(1:6) 
     E_old = TDA%blkM(1)%eigval
 

     if (crit > 100*min_crit) then
        print*, 'convergence failed' 
        exit
     end if 
     min_crit = min(min_crit,crit) 
     if (crit < conv_crit) exit
      
  end do 
  
end if 

! calculate TDA matrices for operators
  if (present(O1)) then 
     call duplicate_ph_mat(HCC,OeCC)

     if (present(O2)) then 
        call duplicate_sp_mat(TDA,O2TDA)
        allocate(O2TDA%blkM(1)%labels(TDA%map(1),2)) 
        O2TDA%blkM(1)%labels = TDA%blkM(1)%labels      
        call calculate_cross_coupled(O2,OeCC,jbas) 
        call calc_TDA(O2TDA,O2,OeCC,jbas)
     end if 
        call duplicate_sp_mat(TDA,O1TDA)     
        allocate(O1TDA%blkM(1)%labels(TDA%map(1),2)) 
        O1TDA%blkM(1)%labels = TDA%blkM(1)%labels      
        call calculate_cross_coupled(O1,OeCC,jbas) 
        call calc_TDA(O1TDA,O1,OeCC,jbas)
        
  end if 

  close(37)
end subroutine  
  

!=========================================================================
!=========================================================================
subroutine BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
  implicit none 
  
  real(8), parameter :: conv = 1e-8
  integer :: trunc,i,m,n,q,j,k,l
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(15),adnorm,fullnorm,s,advals(15)
  character(3) :: args
  
  advals = 0.d0 
  
  cof = (/1.d0,1.d0,0.5d0,0.166666666666666666d0, &
       0.04166666666666666d0,0.0083333333333333333d0,&
       .001388888888888d0,1.984126984d-4,2.48015873016d-5,&
       2.75573192239d-6,2.75573192239d-7,2.505210839d-8, &
       2.087675698d-9,1.6059043837d-10,1.1470745598d-11/) 

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
 
  do i = 2 ,15

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
     
     call add_sq_op(INT1 , 1.d0 , INT2 , cof(i) , HS )   !basic_IMSRG
    
     advals(i) = abs(INT2%E0*cof(i))
   !  if (advals(i) < conv) exit
     
  end do 
 
end subroutine 
!===============================================================
!===============================================================
end module

program main_IMSRG
  use isospin_operators
  use HF_mod
  use IMSRG_ODE
  use IMSRG_MAGNUS
  use IMSRG_CANONICAL
  use operators
  use interaction_IO
  use EOM_IMSRG
  use brute_force_testing
  use three_body_routines
  ! ground state IMSRG calculation for nuclear system 
  implicit none

  type(spd) :: jbas,jbx
  type(sq_op) :: HS,pipj,rirj,r2_rms,Otrans,exp_omega,cr,H0,Hcm
  type(sq_op),allocatable,dimension(:) :: ladder_ops 
  type(iso_ladder),allocatable,dimension(:) :: isoladder_ops
  type(iso_operator) :: GT_trans
  type(full_sp_block_mat) :: coefs,TDA,ppTDA,rrTDA
  type(three_body_force) :: threebod
  type(obsv_mgr) :: trans,moments
  type(eom_mgr) :: eom_states
  character(200) :: inputs_from_command
  character(1) :: quads,trips,trans_type
  integer :: q,ham_type,PAR,Tz,trans_rank,method_int,ex_Calc_int
  integer :: totstates,numstates,oldnum,dTZ,oldnum_dTz,numstates_dTz
  real(8) :: hw ,omp_get_wtime,t1,t2,t3,t4,corr,de_trips
  logical :: hartree_fock,COM_calc,r2rms_calc,me2j,me2b,trans_calc
  logical :: skip_setup,skip_gs,do_HF,TEST_commutators,mortbin,decouple

  ! generator options 
  external :: build_gs_white,build_specific_space,build_gs_atan,build_gs_w2
  external :: build_ex_imtime,build_sd_shellmodel
!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  t1 = omp_get_wtime() ! record start time 

! arguments from command line (no arguments runs test case) 
  call getarg(1,inputs_from_command) 
  call getarg(2,resubmitter) 

  if (trim(inputs_from_command) == 'X') then 
     ! this is just a test run using the parameters of testcase.ini 
     test_commutators = .true.
     inputs_from_command = ''
  else
     test_commutators = .false.
  end if

  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,method_int,ex_calc_int,COM_calc,r2rms_calc,me2j,&
       me2b,mortbin,hw,skip_setup,skip_gs,quads,trips,threebod%e3max)

  call read_sp_basis(jbas,HS%Aprot,HS%Aneut,HS%eMax,HS%lmax,trips,jbx)

  call print_system(jbas) 
  
  if (TEST_COMMUTATORS) then 
     ! run this by typing ' X' after the input file in the command line
     call test
     stop 
  end if 

  ! allocate hamiltonian storage
  call allocate_blocks(jbas,HS)

  HS%herm = 1   ! hermiticity 
  HS%hospace = hw  

!=================================================================
! SET UP OPERATORS
!=================================================================
  if (COM_calc) then   ! THIS IS FOR CALCULATION OF CENTER OF MASS ENERGY
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,pipj)
     call calculate_pipj(pipj,jbas)
     call calculate_rirj(rirj,jbas)
  end if
   
  if (r2rms_calc) then  ! Point-Nucleon charge radius 
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,r2_rms) 
     call initialize_rms_radius(r2_rms,rirj,jbas) 
  end if 
  
!============================================================
!  CAN WE SKIP STUFF?  
!============================================================
  
  do_hf = .true. 
  IF (reading_decoupled) then 
     !Check if the entire IM-SRG calculation has been performed, and stored somewhere. 
     do_hf = read_twobody_operator(HS,'decoupled') 
     ! SEARCH "t2" or "91"
     if (.not. do_hf) goto 91  ! I use goto VERY sparingly, with the counter-intutive
                               ! goal of making things more readable.  
     ! depending on the boolean, we may skip the entire IM-SRG calculation. 
  end if
  
  if (reading_bare) then
     ! This is very useful if working with three-body operators,
     ! you only need to do Hartree Fock once, and then you can
     ! read in a normal ordered two-body force for all subsequent uses. 
     do_hf = read_twobody_operator(HS,'bare')
     
     if (.not. do_hf) then
        call read_umat(coefs,jbas)
        ! SEARCH "hartree fock" or "91"
        goto 90
        ! skipping the hartree fock/ normal ordering portion. 
     end if
  end if
  ! yes, goto 
!=============================================================
! READ INTERACTION FROM FILE
!=============================================================
  print*, 'reading 2-body interaction'
  call read_me2j_interaction(HS,jbas,jbx,ham_type)  ! me2j format 
    
!============================================================
! BUILD HAMILTONIAN
!============================================================

  ! set up the one-body hamiltonian, specified by ham_type, 1: intrinsic, 2: trapped, 3: Full T+V
  call calculate_h0_harm_osc(hw,jbas,HS,ham_type) 
  
  if (threebod%e3Max.ne.0) then
     print*, 'Reading Three Body Force From file'
     call allocate_three_body_storage(jbas,jbx,threebod,HS%eMax,HS%lmax)
     call read_me3j(threebod,jbas,jbx,HS%eMax,HS%lmax)
     ! 3-body forces can only be read from me3j in either .gz or .bin forms 
  end if 
    
  if (hartree_fock) then 
     call calc_HF(HS,threebod,jbas,coefs)   
  else 
     ! hartree fock implicitly normal orders, if we don't do HF, we still normal order. 
     call normal_order(HS,jbas) 
  end if

  ! at this point we are done with 3-body forces, thank god, they are very cumbersome.
  ! they should be mostly encoded in the normal ordered Hamiltonian: NO2B approx. 
  call deallocate_3b(threebod)

  ! lawson 0b term
12 HS%E0 = HS%E0 - HS%lawson_beta * 1.5d0* HS%com_hw
  
  if (writing_bare) then 
     ! write to file so you don't have to do HF in the future.  
     call write_twobody_operator(HS,'bare')
     ! we can store the transformation instead of storing every operator that one might want
     ! to investigate. 
     call write_umat(coefs)
  end if

  ! Normal order scalar observables 
90 if (hartree_fock) then    
     call observable_to_HF(pipj,coefs,jbas)
     call observable_to_HF(rirj,coefs,jbas)
     call observable_to_HF(r2_rms,coefs,jbas)
  else 
     call normal_order(pipj,jbas) 
     call normal_order(rirj,jbas)
     call normal_order(r2_rms,jbas)
  end if

print*, 'SCALAR OPERATOR SETUP COMPLETE' 
!============================================================
! IM-SRG CALCULATION 
!============================================================ 
 
! just a large series of IF statements deciding which type of
! calculation to run. (The call statements run the full calculation) 
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ground state decoupling
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     
  call print_header
  select case (method_int) 

          
  case (1) ! magnus 
     ! Magnus method is recommended due to explicit construction of the transformation. 
     call duplicate_sq_op(HS,exp_omega)
     !exp_omega is the anti-hermitian operator in the exponential U = e^{exp_omega} 
     exp_omega%herm = -1
     decouple=.true. ! this boolean determines if we need to run magnus,
                     ! or if it's already been done and stored 

     if ( reading_omega ) then
        !check if the calculation has been done before
        decouple = read_twobody_operator( exp_omega ,'omega' )
     end if 

     if (trips .ne. 'n') then
        ! For perturbative triples corrections, must keep bare Hamiltonian
        call duplicate_sq_op(HS,H0)
        call copy_sq_op(HS,H0)
     end if
    
     if ( decouple ) then  ! actual IM-SRG calculation 
        call magnus_decouple(HS,exp_omega,jbas,quads,trips,build_gs_w2)    
        if (writing_omega) call write_twobody_operator(exp_omega,'omega')
     else
        print*, 'READ TRANSFORMATION FROM FILE, SKIPPING IMSRG...' 
        call transform_observable_BCH(HS,exp_omega,jbas,quads) 
     end if 

     ! transform scalar observables.
     if (COM_calc) then 
        print*, 'TRANSFORMING Hcm'
        call transform_observable_BCH(pipj,exp_omega,jbas,quads)
        call transform_observable_BCH(rirj,exp_omega,jbas,quads)
        print*, 'CALCULTING Ecm' 
        call calculate_CM_energy(pipj,rirj,hw) 
     end if 
     
     if (r2rms_calc) then
        print*, 'TRANSFORMING RADIUS'
        call transform_observable_BCH(r2_rms,exp_omega,jbas,quads)
        call write_tilde_from_Rcm(r2_rms) 
     end if
        
  case (2) ! traditional flow equation, observables must be transformed simultaneously
     if (COM_calc) then 
        call decouple_hamiltonian(HS,jbas,build_gs_white,pipj,rirj)
        call calculate_CM_energy(pipj,rirj,hw)  ! this writes to file
     else if (r2rms_calc) then 
        call decouple_hamiltonian(HS,jbas,build_gs_white,r2_rms)
!        call write_tilde_from_Rcm(r2_rms)
        print*, sqrt(r2_rms%E0)
    else 
        call decouple_hamiltonian(HS,jbas,build_gs_white) 
    !    call discrete_decouple(HS,jbas) 
     end if
      
  case (3) !discrete product of exponential operators.
     ! observables must be transformed simultaneously. 
     if (COM_calc) then 
        call discrete_decouple(HS,jbas,pipj,rirj)
        call calculate_CM_energy(pipj,rirj,hw)  ! this writes to file
     else if (r2rms_calc) then 
        call discrete_decouple(HS,jbas,r2_rms)
        call write_tilde_from_Rcm(r2_rms)     
     else 
        call discrete_decouple(HS,jbas) 
     end if

  end select
!============================================================
! store hamiltonian in easiest format for quick reading
!============================================================
  if (writing_decoupled) then 
     call write_twobody_operator(HS,'decoupled')
  end if


!============TRIPLES MAGNUS=============================================
  if (trips == 'y') then ! this will be false if magnus is not chosen. 
     ! perturbative triples correction to energy
     print*, 'computing triples'
     t3 = omp_get_wtime() 
     corr =  restore_triples(H0,exp_omega,jbas) 
     t4 = omp_get_wtime()
     print*, 'FINAL ENERGY:', corr + HS%E0,t4-t3
     ! write result to file. 
     open(unit=39,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_magnus_triples.dat')
     write(39,'(2(e15.7))') HS%E0,HS%E0+corr
     close(39)
     
  else if (trips == 'C') then
     ! same as above, but with a more extensive approximate
     ! three-body force. This is the genius of Titus Morris.  
     t3 = omp_get_wtime()
     call duplicate_sq_op(H0,CR)         

     ! "completely renormalized" bit. The goal is to make this
     ! method look like CR-CC(2,3), requires a modified BCH expansion. 
     call CR_EXPAND(CR,exp_omega,H0,jbas,quads) 
     corr =  restore_triples(CR,exp_omega,jbas)
     t4 = omp_get_wtime()
     print*, 'FINAL ENERGY:', corr + HS%E0,t4-t3
     open(unit=39,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_magnus_cr_triples.dat')
     write(39,'(2(e15.7))') HS%E0,HS%E0+corr
     close(39)
  end if
!=======================================================================

91 t2 = omp_get_wtime() 
  write(*,'(A5,f12.7)') 'TIME:', t2-t1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  equation of motion calculation 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  ! Equations of Motion part. 
  if (ex_calc_int==1) then

     ! read_eom_file reads a secondary *.eom file that specifies all
     ! excited states and transitions/moments one might wish to calculate
     
     totstates=read_eom_file(trans,moments,eom_states,jbas)! total number of states

     ! ladder_ops stores eigenfunctions for the current nucleus
     allocate(ladder_ops(totstates-eom_states%total_dTz))
     ! isoladder_ops stores eigenfunctions for isobaric neighbors.  
     allocate(isoladder_ops(eom_states%total_dTz))     

     oldnum = 0
     oldnum_dTz = 0
     numstates = 0
     numstates_dTz = 0
     do q = 1,eom_states%num
        
        ! dTZ specifies which type of EOM we are performing.
        ! dTz == 0 is traditional excited state.
        ! else we are doing isospin projection changing EOM for a different nucleus. 
        
        if (eom_states%dTz(q) == 0 ) then 
           oldnum = oldnum + Numstates
           Numstates = eom_states%number_requested(q)        
           ladder_ops(1+oldnum:Numstates+oldnum)%xindx = q
           call calculate_excited_states(eom_states%ang_mom(q),eom_states%par(q),numstates,HS,&
                jbas,ladder_ops(1+oldnum:Numstates+oldnum))
        else
           oldnum_dTz = oldnum_dTz + Numstates_dTz
           Numstates_dTz = eom_states%number_requested(q)        
           isoladder_ops(1+oldnum_dTz:Numstates_dTz+oldnum_dTz)%xindx = q
           call calculate_isospin_states(eom_states%ang_mom(q),eom_states%par(q),eom_states%dTz(q),&
                numstates_dTZ,HS,jbas,isoladder_ops(1+oldnum_dTz:Numstates_dTz+oldnum_dTz))

        end if

     end do

     t2 = omp_get_wtime() ! we are done with EOM, check time. 
     write(*,'(A5,f12.7)') 'TIME:', t2-t1

     ! xindx specifies which hash tables we should be using to index the various tensor operators 
     Otrans%xindx = eom_states%num+1
     GT_Trans%xindx = Otrans%xindx
     
     trans_type = trans%oper(1:1)
     ! can be 'E','M' or 'G' 
     if (trans_type == 'G') then
        ! GAMOW TELLER: RANK=1, PAR=+
        call initialize_transition_operator('G',trans%oper(2:2),GT_trans,HS,jbas,trans_calc)
     else        
        ! ELECTROMAGNETIC TRANSITION 
        read(trans%oper(2:2),'(I1)') trans_rank        
        call initialize_transition_operator(trans_type,trans_rank,Otrans,HS,jbas,trans_calc)
     end if 
     
     if (trans_calc) then

        ! transform to hartree fock basis 
        if (hartree_fock) then 
           call observable_to_HF(Otrans,coefs,jbas)
        end if

        print*, 'Transforming transition operator...'    
        if (trans_type == 'G') then 
           if ( trans%num + moments%num > 0 ) call transform_observable_BCH(GT_trans,exp_omega,jbas,quads)
        else
           if ( trans%num + moments%num > 0 ) call transform_observable_BCH(Otrans,exp_omega,jbas,quads)
        end if

        ! Center of Mass calculation can be done for each state. 
        if (com_calc) then 
           Hcm%rank = 0
           Hcm%dpar = 0
           Hcm%xindx = Otrans%xindx + 1 
           call build_Hcm(pipj,rirj,Hcm,jbas)
        end if

        ! these routines will actually compute all of the observables you've requested
        if (trans_type == 'G') then 
           call EOM_beta_observables( ladder_ops, isoladder_ops, GT_trans, HS, Hcm,trans, moments,eom_states,jbas)
        else           
           call EOM_observables( ladder_ops, isoladder_ops, Otrans, HS, Hcm,trans, moments,eom_states,jbas)
        end if 
        
     end if

     deallocate(isoladder_ops,ladder_ops)
  end if

contains

subroutine test
  ! it's not recommended to run all of these. Comment them out
  ! and use for unit-testing. 
  call compare_tensor_scalar_commutator(jbas,-1,1) 
  call test_scalar_scalar_commutator(jbas,-1,1) 
  call test_EOM_scalar_scalar_commutator(jbas,1,1)
  call test_EOM_scalar_tensor_commutator(jbas,1,1,6,2)  
  deallocate(phase_hh,phase_pp)
  call test_scalar_tensor_commutator(jbas,-1,1,6,2) 
  call test_tensor_product(jbas,1,1,2,4,6,2,2,0) 
  call test_EOM_iso_commutator(jbas,1,1,4,0,0)
  call test_scalar_iso_commutator(jbas,-1,1,6,2,1) 
  call test_tensor_dTZ_product(jbas,1,1,2,0,2,0,0,0,-1)
end subroutine test

end program main_IMSRG
!=========================================================================
subroutine print_header
  implicit none 
  
  print* 
  print*, 'Constructing Basis...' 
  print*
  print*, '================================'//&
       '==================================='
  print*, '  iter        s            E0      '//&
       '    E0+MBPT(2)      |MBPT(2)|  '  
  print*, '================================'//&
       '==================================='

end subroutine   



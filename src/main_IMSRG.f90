program main_IMSRG
!!!===================================================================
!!!     Equations of Motion IMSRG for closed-shell nuclei 
!!!     Written by Nathan M. Parzuchowski       
!!!====================================================================
  use isospin_operators
  use HF_mod
  use response
  use IMSRG_ODE
  use IMSRG_MAGNUS
  use IMSRG_CANONICAL
  use operators
  use interaction_IO
  use EOM_IMSRG
  use brute_force_testing
  use three_body_routines
  use deuteron
  implicit none
  
  type(spd) :: jbas,jbx
  type(sq_op) :: HS,ETA,DH,w1,w2,rirj,pipj,r2_rms,Otrans,exp_omega
  type(sq_op) :: num,cr,H0,Hcm,Oscal 
  type(sq_op),allocatable,dimension(:) :: ladder_ops 
  type(iso_ladder),allocatable,dimension(:) :: isoladder_ops
  type(iso_operator) :: GT_trans
  type(cc_mat) :: CCHS,CCETA,WCC
  type(full_sp_block_mat) :: coefs,coefsT,TDA,ppTDA,rrTDA
  type(three_body_force) :: threebod
  type(obsv_mgr) :: trans,moments
  type(eom_mgr) :: eom_states
  character(200) :: inputs_from_command
  character(10) :: other_obs
  character(1) :: quads,trips,trans_type
  integer :: i,j,T,JTot,a,b,c,d,g,q,ham_type,j3,ix,jx,kx,lx,PAR,Tz,trans_rank
  integer :: np,nh,nb,k,l,m,n,method_int,mi,mj,ma,mb,j_min,ex_Calc_int,J1,J2
  integer :: na,la,lb,totstates,numstates,oldnum,qx,dTZ,oldnum_dTz,numstates_dTz
  real(8) :: sm,omp_get_wtime,t1,t2,bet_off,d6ji,gx,dcgi,dcgi00
  real(8) :: corr_op,pre,x,corr,de_trips
  logical :: hartree_fock,COM_calc,r2rms_calc,me2j,me2b,trans_calc
  logical :: skip_setup,skip_gs,do_HF,TEST_commutators,mortbin,decouple
  external :: build_gs_white,build_specific_space,build_gs_atan,build_gs_w2
  external :: build_ex_imtime,build_sd_shellmodel
!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  time0 = omp_get_wtime()

  call getarg(1,inputs_from_command) 
  call getarg(2,resubmitter) 
  
  if (trim(inputs_from_command) == 'X') then 
     test_commutators = .true.
     inputs_from_command = ''
  else
     test_commutators = .false.
  end if

  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,method_int,ex_calc_int,COM_calc,r2rms_calc,other_obs,me2j,&
       me2b,mortbin,jbas%hw,skip_setup,skip_gs,quads,trips,threebod%e3max)

  call read_sp_basis(jbas,HS%Aprot,HS%Aneut,HS%eMax,HS%lmax,trips,jbx)

  call print_system(jbas) 
  
  if (TEST_COMMUTATORS) then 
     ! run this by typing ' X' after the input file in the command line
     ! This takes forever, you might want to comment some of this out. 
     call test
     stop
  end if 

  call allocate_blocks(jbas,HS)

  HS%herm = 1
  HS%hospace = jbas%hw

!=================================================================
! SET UP OPERATORS
!=================================================================
  if (COM_calc) then 
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,pipj)
     call calculate_pipj(pipj,jbas)
     call calculate_rirj(rirj,jbas)
  end if
  
  if (r2rms_calc) then 
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,r2_rms) 
     call initialize_rms_radius(r2_rms,rirj,jbas) 
  end if 

 if (.not. (trim(other_obs)=="none")) then 
    call  duplicate_sq_op(HS,Oscal)
    call initialize_scalar_operator(other_obs,Oscal,jbas)   
 end if 
  
!============================================================
!  CAN WE SKIP STUFF?  
!============================================================
  
  do_hf = .true. 
  IF (reading_decoupled) then 
     do_hf = read_twobody_operator(HS,'decoupled') 
     if (.not. do_hf) goto 91 
  end if  
  
  if (reading_bare) then 
     do_hf = read_twobody_operator(HS,'bare')

     if (.not. do_hf) then       
        if (hartree_fock) call read_umat(coefs,jbas)
        goto 90   ! goto is used to skip things in the main program file 
     end if
  end if 

!=============================================================
! READ INTERACTION 
!=============================================================
  print*, 'reading 2-body interaction'
  ! for calculating COM expectation value
  
  if (me2j) then 
     call read_me2j_interaction(HS,jbas,jbx,ham_type) 
  else if (me2b) then
     ! pre normal ordered interaction with three body included at No2b       
     print*, 'READING PRE-NORMAL ORDERED INTERACTION FROM HEIKO' 
     call read_me2b_interaction(HS,jbas,ham_type) 
     goto 12 ! skip the normal ordering. 
     ! it's already done.  line 128 or search "bare" 
  else if (mortbin) then 
     call read_gz(HS,jbas,ham_type) 
  else
     call read_interaction(HS,jbas,ham_type)
  end if
  call print_time
!============================================================
! BUILD BASIS
!============================================================
  call calculate_h0_harm_osc(jbas,HS,ham_type) 

  if (HS%Aprot+HS%Aneut == 2) then
     call compute_deuteron_ground_state(HS,jbas,Oscal)    
     stop
  end if
  
  !! threebody machinery if needed
  if (threebod%e3Max.ne.0) then 
     print*, 'Reading Three Body Force From file'
     call allocate_three_body_storage(jbas,jbx,threebod,HS%eMax,HS%lmax)
     call read_me3j(threebod,jbas,jbx,HS%eMax,HS%lmax)
     call print_time
  end if

  if (hartree_fock) then 
     call calc_HF(HS,threebod,jbas,coefs)   
  else
     call calc_HO(HS,threebod,jbas,coefs)   
  end if

  call print_time

  call deallocate_3b(threebod)

  ! lawson 0b term
12 HS%E0 = HS%E0 - HS%lawson_beta * 1.5d0* HS%com_hw

  !! write bare operator to gzipped file (so you don't need to do HF again)   
  if (writing_bare) then 
     call write_twobody_operator(HS,'bare')    
     call write_umat(coefs)   ! store HF transformation 
  end if

  if (writing_human) then
     call write_onebody_human(HS,jbas,"ham0")   
     call write_twobody_human(HS,jbas,"ham0")
  end if

  ! Normal order scalar observables
  print*, "Normal ordering scalar observables..." 
90 if (hartree_fock) then    
     call observable_to_HF(Oscal,coefs,jbas) 
     call observable_to_HF(pipj,coefs,jbas)
     call observable_to_HF(rirj,coefs,jbas)
     call observable_to_HF(r2_rms,coefs,jbas)
  else 
     call normal_order(Oscal,jbas) 
     call normal_order(pipj,jbas) 
     call normal_order(rirj,jbas)
     call normal_order(r2_rms,jbas)
  end if

print*, 'BASIS SETUP COMPLETE' 
!============================================================
! IM-SRG CALCULATION 
!============================================================ 
 
! just a large series of IF statements deciding exactly which type of
! calculation to run. (The call statements run the full calculation) 
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ground state decoupling
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call print_time
     
call print_header

  select case (method_int) 

          
  case (1) ! magnus 
     
     call duplicate_sq_op(HS,exp_omega)
     exp_omega%herm = -1
     decouple=.true.

     !!! if we sucessfully read an omega, "decouple" is set to .false.
     if ( reading_omega ) then
        decouple = read_twobody_operator( exp_omega ,'omega' )
     end if 

     if (trips .ne. 'n') then 
        call duplicate_sq_op(HS,H0)
        call copy_sq_op(HS,H0)
     end if

     if ( decouple ) then 

        call magnus_decouple(HS,exp_omega,jbas,quads,trips,build_gs_w2)    
        open(unit=42,file=trim(OUTPUT_DIR)//&
             trim(adjustl(prefix))//'_Egs.dat')
        write(42,'(2(I4),e14.6)')  nint(HS%hospace), HS%eMax, HS%E0
        close(42)
        if (writing_omega) call write_twobody_operator(exp_omega,'omega')
        
 
     else
        print*, 'READ TRANSFORMATION FROM FILE, SKIPPING IMSRG...' 
        call transform_observable_BCH(HS,exp_omega,jbas,quads)
     end if 
     
     if (.not. (trim(other_obs)=="none")) then 
        call transform_observable_BCH(Oscal,exp_omega,jbas,quads)
        call write_scalar_operator(other_obs, Oscal,jbas)
     end if
     

     if (COM_calc) then 
        print*, 'TRANSFORMING Hcm'
        call transform_observable_BCH(pipj,exp_omega,jbas,quads)
        call transform_observable_BCH(rirj,exp_omega,jbas,quads)
        print*, 'CALCULTING Ecm' 
        call calculate_CM_energy(pipj,rirj) 
     end if 
     
     if (r2rms_calc) then
        print*, 'TRANSFORMING RADIUS'
        call transform_observable_BCH(r2_rms,exp_omega,jbas,quads)
        call write_tilde_from_Rcm(r2_rms) 
     end if
        
  case (2) ! traditional
     if (COM_calc) then 
        call decouple_hamiltonian(HS,jbas,build_gs_white,pipj,rirj)
        call calculate_CM_energy(pipj,rirj)  ! this writes to file
     else if (r2rms_calc) then 
        call decouple_hamiltonian(HS,jbas,build_gs_white,r2_rms)
        print*, sqrt(r2_rms%E0)
    else 
        call decouple_hamiltonian(HS,jbas,build_gs_white) 
     end if
     
  case (3)    ! product of discrete infinitesimal transformations
     if (COM_calc) then 
        call discrete_decouple(HS,jbas,pipj,rirj)
        call calculate_CM_energy(pipj,rirj)  ! this writes to file
     else if (r2rms_calc) then 
        call discrete_decouple(HS,jbas,r2_rms)
        print*, sqrt(r2_rms%E0)
     else 
        call discrete_decouple(HS,jbas) 
     end if

  end select
  call print_time
!============================================================
! store hamiltonian in easiest format for quick reading
!============================================================
  if (writing_decoupled) then 
     call write_twobody_operator(HS,'decoupled')
  end if

!============TRIPLES MAGNUS=============================================
  if (trips == 'y') then 
     print*, 'computing triples'
     t1 = omp_get_wtime()
     call restore_triples(H0,exp_omega,jbas,corr,corr_op) 
     t2 = omp_get_wtime()
     print*, 'FINAL ENERGY:', corr + HS%E0,t2-t1,HS%E0
     open(unit=39,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_magnus_triples.dat')
     write(39,'(2(e15.7))') HS%E0,HS%E0+corr
     close(39)
     call print_time
  else if (trips == 'C') then
     t1 = omp_get_wtime()
     call duplicate_sq_op(H0,CR)         
     ! completely renormalized bit.
     call CR_EXPAND(CR,exp_omega,H0,jbas,quads) 
     call restore_triples(CR,exp_omega,jbas,corr,corr_op)
     t2 = omp_get_wtime()
     print*, 'FINAL ENERGY:', corr + HS%E0,t2-t1
     open(unit=39,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_magnus_cr_triples.dat')
     write(39,'(2(e15.7))') HS%E0,HS%E0+corr
     close(39)
     call print_time
  end if
!=======================================================================
!  equations-of-motion calculation 
!=======================================================================
91 if (ex_calc_int==1) then
     print*, 'STARTING EXCITED STATES CALCULATION'
     totstates=read_eom_file(trans,moments,eom_states,jbas)! total number of states
     print*, totstates
     allocate(ladder_ops(totstates-eom_states%total_dTz))
     allocate(isoladder_ops(eom_states%total_dTz))     
     
     oldnum = 0
     oldnum_dTz = 0
     numstates = 0
     numstates_dTz = 0
     do q = 1,eom_states%num

        if (eom_states%dTz(q) == 0 ) then 
           oldnum = oldnum + Numstates
           Numstates = eom_states%number_requested(q)        
           ladder_ops(1+oldnum:Numstates+oldnum)%xindx = q
           call calculate_excited_states(eom_states%ang_mom(q),eom_states%par(q),numstates,HS,&
                jbas,ladder_ops(1+oldnum:Numstates+oldnum))
        
           call print_time
           if (eom_states%trips) then
              call print_triples_header
              do qx = 1+oldnum,Numstates+oldnum
                 
                 t1= omp_get_wtime()
                 dE_trips=EOM_triples(HS,ladder_ops(qx),jbas)  
                 t2= omp_get_wtime()
                                  
                 write(*,'(A2,4(f20.10))') eom_states%name(q),&
                      ladder_ops(qx)%E0,ladder_ops(qx)%E0+dE_trips,dE_trips,t2-t1
              end do
           end if

        else
   
           oldnum_dTz = oldnum_dTz + Numstates_dTz
           Numstates_dTz = eom_states%number_requested(q)        
           isoladder_ops(1+oldnum_dTz:Numstates_dTz+oldnum_dTz)%xindx = q
           call calculate_isospin_states(eom_states%ang_mom(q),eom_states%par(q),eom_states%dTz(q),&
                numstates_dTZ,HS,jbas,isoladder_ops(1+oldnum_dTz:Numstates_dTz+oldnum_dTz))
        
        end if        
               
    end do

    if ((trans%num + moments%num)*(Numstates+Numstates_dTz) == 0) then 
       trans_Calc=.false.
    else
       trans_Calc= .true.
    end if
    
    call print_time
    
     Otrans%xindx = eom_states%num+1
     GT_Trans%xindx = Otrans%xindx

     trans_type = trans%oper(1:1)

     if (trans_type == 'G') then 
        call initialize_transition_operator('G',trans%oper(2:2),GT_trans,HS,jbas,trans_calc)
     else        
        read(trans%oper(2:2),'(I1)') trans_rank        
        call initialize_transition_operator(trans_type,trans_rank,Otrans,HS,jbas,trans_calc)
     end if 

     if (writing_human) then
        call write_onebody_tensor_human(Otrans,jbas,"E2_0")   
        call write_twobody_tensor_human(Otrans,jbas,"E2_0")
     end if

     
     if (trans_calc) then 
        print*, 'Transforming transition operator...' 

        if (trans_type == 'G') then 
           if (hartree_fock) then 
              call observable_to_HF(GT_trans,coefs,jbas)
           end if
           if ( trans%num + moments%num > 0 ) call transform_observable_BCH(GT_trans,exp_omega,jbas,quads)
                                            
        else
           if (hartree_fock) then 
              call observable_to_HF(Otrans,coefs,jbas)
           end if

           if ( trans%num + moments%num > 0 ) call transform_observable_BCH(Otrans,exp_omega,jbas,quads)
        end if
   
     end if
     call print_time

     if (writing_human) then
        ! make sure the string is 4 characters
        call write_onebody_human(HS,jbas,"ham ")   
        call write_twobody_human(HS,jbas,"ham ")
        call write_onebody_tensor_human(Otrans,jbas,"E2  ")   
        call write_twobody_tensor_human(Otrans,jbas,"E2  ")
     end if

     if (com_calc) then 
        Hcm%rank = 0
        Hcm%dpar = 0
        Hcm%xindx = Otrans%xindx + 1 
        call build_Hcm(pipj,rirj,Hcm,jbas)
     end if
     
     if (trans_type == 'G') then 
        call EOM_beta_observables( ladder_ops, isoladder_ops, GT_trans, HS, Hcm,trans, moments,eom_states,jbas)
     else           
        call EOM_observables( ladder_ops, isoladder_ops, Otrans, HS, Hcm,trans, moments,eom_states,jbas)

        if (eom_states%response) then 
           call compute_response_function(jbas,HS,Otrans) 
        end if
        
     end if
     deallocate(isoladder_ops,ladder_ops)
  end if

  call print_time
contains

subroutine test
  ! call compare_tensor_scalar_commutator(jbas,-1,1) 
  ! stop
  ! deallocate(jbas%xmap) 
  ! call test_scalar_scalar_commutator(jbas,-1,1) 
  ! deallocate(jbas%xmap)
  ! call test_EOM_scalar_scalar_commutator(jbas,1,1)
  ! deallocate(jbas%xmap)
!   call test_EOM_scalar_tensor_commutator(jbas,1,1,6,2)  
!   deallocate(jbas%xmap,jbas%xmap_tensor,phase_hh,phase_pp)
!   deallocate(half6j%tp_mat)
!  call test_scalar_tensor_commutator(jbas,-1,1,6,2) 
 ! call test_tensor_product(jbas,1,1,2,2,2,2,0,2) 
!  call test_EOM_iso_commutator(jbas,1,1,4,0,0)
!  call test_scalar_iso_commutator(jbas,-1,1,6,2,1) 
  call test_tensor_dTZ_product(jbas,1,1,2,2,2,2,0,2,-1) 

end subroutine test
end program main_IMSRG
!=========================================================================
subroutine print_header
  implicit none 
  
  print*, '================================'//&
       '==================================='
  print*, '  iter        s            E0      '//&
       '    E0+MBPT(2)      |MBPT(2)|  '  
  print*, '================================'//&
       '==================================='

end subroutine print_header

subroutine print_triples_header
  implicit none
  
  print*
  print*, '====================================='//&
       '=============================================='
  print*, 'J^Pi      dE(1+2)            dE(1+2+3)'//&
       '             dE(3)               time    '
  print*, '====================================='//&
       '=============================================='
  
end subroutine print_triples_header

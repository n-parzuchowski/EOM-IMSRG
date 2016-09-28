module IMSRG_ODE
  use commutators
  use adams_ode 
  use operators
  implicit none 
  
contains
!==============================================================
subroutine decouple_hamiltonian(H,jbas,build_generator,O1,O2) 
  ! runs IMSRG using the specified derivative calculator
  implicit none 

  ! SRG convergence / failsafe / error tolerances
  integer,parameter :: max_steps = 10000
  real(8),parameter :: conv_crit = 1.d-6
  real(8),parameter :: relerr = 1.d-6, abserr = 1.d-6

  type(spd) :: jbas
  type(sq_op) :: H ,HOD
  type(sq_op),optional :: O1,O2
  type(cc_mat) :: HCC
  type(full_sp_block_mat) :: TDA
  integer,dimension(5) :: iwork
  real(8),allocatable,dimension(:) :: cur_vec,work
  integer :: neq,iflag,Atot,Ntot,nh,np,nb,q,steps  
  real(8) :: ds,s,E_old,crit,E_mbpt2
  logical :: com_calc 
  external :: build_generator,dHds,dHds_1op,dHds_2op

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! first figure out how many equations there are:
  Atot = H%belowEF
  Ntot = H%Nsp


  neq = 1 + Atot*Atot + Atot*(Ntot-Atot) + (Ntot - Atot)**2 

  do q = 1, H%nblocks

     nh = H%mat(q)%nhh
     np = H%mat(q)%npp
     nb = H%mat(q)%nph 

     neq = neq + (nh*nh+nh +  nb*nb+nb + np*np+np)/2 + nb*np + nh*np + nh*nb 
  end do

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  H%neq = neq 

  ! add more room if there are other operators
  if (present(O1)) then 
     O1%neq=H%neq
     if (present(O2)) then
        O2%neq = H%neq
        neq = 3*neq
     else 
        neq = 2*neq
     end if
  end if

  allocate(cur_vec(neq)) ! carries the system into SG solver
  allocate(work(100+21*neq))  ! memory eater

  ! parameters for solver
  work = 0.d0
  iwork = 0 
  iflag = 1

  ! flow equation variables
  ds = 0.1d0
  s = 0.d0 

  steps = 0 

  open(unit=36,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_0bflow.dat')

  E_mbpt2 = mbpt2(H,jbas) 
  crit = abs(E_mbpt2)

  write(36,'(I6,4(e17.9))') steps,s,H%E0,H%E0+E_mbpt2,crit
  write(*,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit

  if (present(O1)) then 
     if (present(O2)) then 

        ! main loop   (two operators) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        do while (steps < max_steps) 

           E_old = H%E0 
           ! send info to SG solver
           call vectorize(H,cur_vec(1:H%neq))
           call vectorize(O1,cur_vec(H%neq+1:2*H%neq))
           call vectorize(O2,cur_vec(2*H%neq+1:3*H%neq)) 

           call ode(dHds_2op,build_generator,neq,cur_vec,H,jbas,&
                s,s+ds,relerr,abserr,iflag,work,iwork) 

           call repackage(H,cur_vec(1:H%neq)) 
           call repackage(O1,cur_vec(H%neq+1:2*H%neq))
           call repackage(O2,cur_vec(2*H%neq+1:3*H%neq)) 

           steps = steps + 1

           ! weak convergence criteria, but it works

           E_mbpt2 = mbpt2(H,jbas) 
           crit = abs(E_mbpt2)

           write(36,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
           write(*,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit

           if (crit < conv_crit) exit

        end do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     else 

        ! main loop (one operator) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        do while (steps < max_steps) 

           E_old = H%E0 
           ! send info to SG solver
           call vectorize(H,cur_vec(1:H%neq))
           call vectorize(O1,cur_vec(H%neq+1:2*H%neq))

           call ode(dHds_1op,build_generator,neq,cur_vec,H,jbas,&
                s,s+ds,relerr,abserr,iflag,work,iwork) 

           call repackage(H,cur_vec(1:H%neq)) 
           call repackage(O1,cur_vec(H%neq+1:2*H%neq))   
           steps = steps + 1

           ! weak convergence criteria, but it works

           E_mbpt2 = mbpt2(H,jbas) 
           crit = abs(E_mbpt2)

           write(36,'(I6,5(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
           write(*,'(I6,5(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit

           if (crit < conv_crit) exit

        end do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     end if
  else
     ! main loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
     do while (steps < max_steps) 

        E_old = H%E0 
        ! send info to SG solver
        call vectorize(H,cur_vec)

        call ode(dHds,build_generator,neq,cur_vec,H,jbas,&
             s,s+ds,relerr,abserr,iflag,work,iwork) 

        call repackage(H,cur_vec) 

        steps = steps + 1

        ! weak convergence criteria, but it works

        E_mbpt2 = mbpt2(H,jbas) 
        crit = abs(E_mbpt2)

        write(36,'(I6,4(e17.9))') steps,s,H%E0,H%E0+E_mbpt2,crit
        write(*,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit

        if (crit < conv_crit) exit

     end do
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  end if

  close(36)
end subroutine decouple_hamiltonian
!================================================
!================================================
end module IMSRG_ODE
!These are external so that they can be passed to ODE   
!================================================
!================================================
subroutine dHds(t,yy,yp,HS,jbas,build_gen) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, (adams_ode.f90) modified to include HS,jbas
  use basic_IMSRG
  use commutators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,k,l,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,IX,JX,a,b,c,d
  integer :: ji,jj,jk,jl,ja,jb,jc,jd,J1,J2,q1,q2,ril,rkj,J3,J4,J5,J6,J7
  real(8) :: yp(*),yy(*),sm,x,d6ji,coef9,smr1,smr2,smr3,smr4
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2,w3,w4
  type(cc_mat) :: WCC,ETACC,HSCC,KCC 
  external :: build_gen

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call init_ph_mat(HS,HSCC,jbas) ! cross coupled ME
  call duplicate_ph_mat(HSCC,ETACC) !cross coupled ME
  call init_ph_wkspc(HSCC,WCC) ! workspace for CCME

  call build_gen(HS,ETA,jbas) ! constructs generator
  
  call calculate_cross_coupled(HS,HSCC,jbas)
  call calculate_cross_coupled(ETA,ETACC,jbas) 
   
  DH%E0 = commutator_110(ETA,HS,jbas) + commutator_220(ETA,HS,jbas)

  call commutator_111(ETA,HS,DH,jbas) 
  call commutator_121(ETA,HS,DH,jbas)
  call commutator_122(ETA,HS,DH,jbas)

  call commutator_222_pp_hh(ETA,HS,DH,w1,w2,jbas)

  call commutator_221(ETA,HS,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp)
  
end subroutine dHds
!==============================================================
!==============================================================
subroutine dHds_1op(t,yy,yp,HS,jbas,build_gen) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l,num_ops 
  real(8) :: yp(*),yy(*)
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2,O1,O2
  type(cc_mat) :: WCC,ETACC,HSCC
  external :: build_gen

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

  neq = 2*HS%neq
 
! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
 
  call duplicate_sq_op(HS,O1)   
  call repackage(O1,yy(HS%neq+1:2*HS%neq)) 
  
  call init_ph_mat(HS,HSCC,jbas) ! cross coupled ME
  call duplicate_ph_mat(HSCC,ETACC) !cross coupled ME
  call init_ph_wkspc(HSCC,WCC) ! workspace for CCME

  call build_gen(HS,ETA,jbas) ! constructs generator
  
  call calculate_cross_coupled(HS,HSCC,jbas)
  call calculate_cross_coupled(ETA,ETACC,jbas) 
   
  DH%E0 = commutator_110(ETA,HS,jbas) + commutator_220(ETA,HS,jbas)
 
  call commutator_111(ETA,HS,DH,jbas) 
  call commutator_121(ETA,HS,DH,jbas)
  call commutator_122(ETA,HS,DH,jbas)
  
  call commutator_222_pp_hh(ETA,HS,DH,w1,w2,jbas)
  
  call commutator_221(ETA,HS,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(1:HS%neq))
  
!===================================================================
  
  call calculate_cross_coupled(O1,HSCC,jbas)
  
  DH%E0 = commutator_110(ETA,O1,jbas) + commutator_220(ETA,O1,jbas)
 
  call commutator_111(ETA,O1,DH,jbas) 
  call commutator_121(ETA,O1,DH,jbas)
  call commutator_122(ETA,O1,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O1,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O1,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(HS%neq+1:2*HS%neq))

end subroutine dHds_1op
!==================================================================
!==================================================================
subroutine dHds_2op(t,yy,yp,HS,jbas,build_gen) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l,num_ops 
  real(8) :: yp(*),yy(*)
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2,O1,O2
  type(cc_mat) :: WCC,ETACC,HSCC 
  external :: build_gen
  
!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

  neq = 3*HS%neq
 
! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  
  call duplicate_sq_op(HS,O1) 
  call repackage(O1,yy(HS%neq+1:2*HS%neq)) 

  call duplicate_sq_op(HS,O2)
  call repackage(O2,yy(2*HS%neq+1:3*HS%neq))
  
  call init_ph_mat(HS,HSCC,jbas) ! cross coupled ME
  call duplicate_ph_mat(HSCC,ETACC) !cross coupled ME
  call init_ph_wkspc(HSCC,WCC) ! workspace for CCME

  call build_gen(HS,ETA,jbas) ! constructs generator
  
  call calculate_cross_coupled(HS,HSCC,jbas)
  call calculate_cross_coupled(ETA,ETACC,jbas) 
   
  DH%E0 = commutator_110(ETA,HS,jbas) + commutator_220(ETA,HS,jbas)
 
  call commutator_111(ETA,HS,DH,jbas) 
  call commutator_121(ETA,HS,DH,jbas)
  call commutator_122(ETA,HS,DH,jbas)
  
  call commutator_222_pp_hh(ETA,HS,DH,w1,w2,jbas)
  
  call commutator_221(ETA,HS,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(1:HS%neq))
  
!===================================================================
  
  call calculate_cross_coupled(O1,HSCC,jbas)
  
  DH%E0 = commutator_110(ETA,O1,jbas) + commutator_220(ETA,O1,jbas)
 
  call commutator_111(ETA,O1,DH,jbas) 
  call commutator_121(ETA,O1,DH,jbas)
  call commutator_122(ETA,O1,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O1,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O1,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(HS%neq+1:2*HS%neq))

!===================================================================
  call calculate_cross_coupled(O2,HSCC,jbas)
   
  DH%E0 = commutator_110(ETA,O2,jbas) + commutator_220(ETA,O2,jbas)

  call commutator_111(ETA,O2,DH,jbas) 
  call commutator_121(ETA,O2,DH,jbas)
  call commutator_122(ETA,O2,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O2,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O2,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(2*HS%neq+1:3*HS%neq))
end subroutine dHds_2op

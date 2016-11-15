module HF_mod
  use basic_IMSRG 
  use three_body_routines
  ! J-scheme Hartree-Fockery. 
  implicit none 
  
contains
!====================================================
subroutine calc_HF(H,THREEBOD,jbas,D)
  ! returns H in the normal orderd Hartree Fock basis
  implicit none 
  
  type(three_body_force) :: THREEBOD
  type(mono_3b) :: THREEBOD_MONO
  type(spd) :: jbas
  type(sq_op) :: H 
  type(full_sp_block_mat) :: T,F,Vgam,V3gam,rho,D,Dx
  integer :: q,r,i,j,k,l
  real(8) :: crit,sm
  logical :: tbforce=.false.
  
  if ( allocated( THREEBOD%mat) ) tbforce = .true.  
     
  ! allocate the workspace   
  call allocate_sp_mat(jbas,T) 
  call duplicate_sp_mat(T,F) 
  call duplicate_sp_mat(T,Vgam) 
  call duplicate_sp_mat(T,rho)
  call duplicate_sp_mat(T,D)
  call duplicate_sp_mat(T,Dx)
  call duplicate_sp_mat(T,V3gam)
  
  if ( tbforce ) then
     ! monopole matrix elements make life easier
     call allocate_mono(THREEBOD_MONO,jbas)      
  end if 
  
  call write_kin_matrix(T,H,jbas) 

  !initial eigenvectors
  do q = 1, D%blocks
     D%blkM(q)%eigval = 0.d0
     D%blkM(q)%matrix = 0.d0 
     do i =1,D%map(q)
        D%blkM(q)%matrix(i,i) = 1.d0
     end do 
  end do 

  crit = 10.d0 
  r = 0
     
  !!! HARTREE FOCK MAIN LOOP 
  print*, 'Computing Hartree Fock Basis...'
  do while (crit > 1e-6) 
     
     call density_matrix(rho,D,DX,jbas)      
     call gamma_matrix(Vgam,H,rho,jbas)
     if (tbforce) call gamma_matrix_three_body(V3gam,rho,THREEBOD,THREEBOD_MONO,jbas)
   
     ! fock matrix
     do q = 1,T%blocks
        F%blkM(q)%matrix = T%blkM(q)%matrix + Vgam%blkM(Q)%matrix&
             + V3gam%blkM(q)%matrix
     end do


     call diagonalize_blocks(F) 
     ! new eigenvectors
     ! calculate conv. criteria
     ! store eigenvalues
 
     crit = 0.d0 
     do q = 1,T%blocks
        D%blkM(q)%matrix = F%blkM(q)%matrix
        crit = crit + sqrt(sum((D%blkM(q)%eigval-F%blkM(q)%eigval)**2))
        D%blkM(q)%eigval = F%blkM(q)%eigval       
     end do

     r = r +1
 end do 

 do q = 1,T%blocks
   do i = 1, D%map(q) 
      if (D%blkM(q)%matrix(i,i) < 0.d0) then 
         ! make phase consistent for all eigenstates.
         ! helps for comparing matrix elements, does not affect observables
         D%blkM(q)%matrix(:,i) = D%blkM(q)%matrix(:,i)*(-1)
      end if 
   end do
   
   F%blkM(q)%matrix = T%blkM(q)%matrix + Vgam%blkM(Q)%matrix &
        + V3gam%blkM(q)%matrix

   T%blkM(q)%eigval = F%blkM(q)%eigval    
 end do
 
 call transform_1b_to_HF(D,Dx,F,H,jbas,T,Vgam,V3gam) 
  
 ! this needs to come after the transformation
 ! e_HF is calculated in the hartree fock basis
 H%E0 = e_HF(T,Vgam,V3gam,jbas)
 write(*,'(A22,f12.7)') ' Hartree Fock Energy: ',H%E0
 if ( tbforce ) call meanfield_2b(rho,H,THREEBOD,jbas) !horrifying normal ordering of three-body force
 call transform_2b_to_HF(D,H,jbas) 

end subroutine calc_HF
!====================================================
subroutine observable_to_HF(Op,coefs,jbas) 
  ! transforms observable to HF basis
  ! works for scalar observables 
  ! and ONE BODY TENSORS ONLY 
  implicit none 
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: coefs,T,Dx
  type(sq_op) :: Op 
  
  if (allocated(Op%tblck)) then 
     call transform_1b_to_HF_tensor(coefs,Op,jbas) 
     ! NOTE THAT SINCE THIS IS A 1-BODY OBSERVABLE,
     ! NORMAL ORDERING IS NOT NECESSARY
     ! THERE SHOULDN'T BE A ZERO BODY PIECE EITHER
     ! THAT WOULD BE A 0+ -> 0+ TRANSITION, WHICH IS 
     ! FORBIDDEN FOR ALL RANK> 0 TENSORS 
  else if (allocated(Op%mat)) then 
     call duplicate_sp_mat(coefs,T)
     call duplicate_sp_mat(coefs,Dx)   
     call write_kin_matrix(T,Op,jbas) 
     call transform_1b_to_HF(coefs,Dx,T,Op,jbas)
     call transform_2b_to_HF(coefs,Op,jbas) 
     call normal_order(Op,jbas)
  else 
     return
  end if
 
end subroutine observable_to_HF
!====================================================
subroutine write_kin_matrix(T,H,jbas)
  ! the kinetic energy matrix has already been calculated
  ! in the setup of the problem, here we are just writing it
  ! to a more convenient block matrix for the HF calculation
  implicit none 
  
  type(full_sp_block_mat) :: T
  type(sq_op) :: H
  type(spd) :: jbas
  integer :: q,i,j,n1,n2,c1,c2,cx,AX
  
  AX = H%belowEF
  
  do q = 1, T%blocks
    
     do i = 1, T%map(q)
        do j = i,T%map(q) 
           
           n1 = T%blkM(q)%states(i) 
           n2 = T%blkM(q)%states(j) 
           
           T%blkM(q)%matrix(i,j) = f_elem(n1,n2,H,jbas) 
           T%blkM(q)%matrix(j,i) = T%blkM(q)%matrix(i,j) 
           
           
           end do 
        end do 
     end do 
     
end subroutine write_kin_matrix
!==================================================== 
subroutine gamma_matrix(gam,int,rho,jbas) 
  ! hartree fock potential matrix
  implicit none
  
  type(full_sp_block_mat) :: gam,rho
  type(sq_op) :: int
  type(spd) :: jbas
  integer :: q,r,i,jmax,lmax,n1,n2,j,JJ,n3,n4,tzrho,tzfoc
  integer :: grho,hrho,qrho,jrho,lrho,PAR,jfoc,lfoc,TZ
  real(8) :: sm

  do q = 1, gam%blocks

     gam%blkM(q)%matrix = 0.d0 
     jfoc = rho%blkM(q)%lmda(2)
     
     ! loop over states in this block
     do i = 1, gam%map(q)
        do j = i,gam%map(q) 
        
           n1 = rho%blkM(q)%states(i) 
           n3 = rho%blkM(q)%states(j) 
           
           ! sum over blocks of the density matrix
           sm = 0.d0 
           do qrho =  1, rho%blocks
           
              jrho = rho%blkM(qrho)%lmda(2) 
                  
              ! sum over elements of the block
              do grho = 1,rho%map(qrho) 
                 do hrho = 1,rho%map(qrho) 
                    
                    n2 = rho%blkM(qrho)%states(grho) 
                    n4 = rho%blkM(qrho)%states(hrho) 
                       
                    ! sum over allowed JJ values
                    
                    do JJ = abs(jrho - jfoc),jrho+jfoc,2
                       
                       sm = sm + rho%blkM(qrho)%matrix(hrho,grho) &
                            * v_elem(n1,n2,n3,n4,JJ,int,jbas) &
                            * (JJ + 1.d0)/(jrho + 1.d0)
                             ! this (jrho+1.d0) comes from the density matrix
                       ! because the m-scheme density matrix has 1's instead
                       ! of 2j+1's 
                    end do 
                    
                 end do 
              end do 
           end do 
                      
           ! the matrix elements are multiplied by (2J+1)/(2j+1)/(2j'+1) 
           ! which amounts to averaging over J 
          
           gam%blkM(q)%matrix(i,j) = sm/(jfoc + 1.d0)
           gam%blkM(q)%matrix(j,i) = sm/(jfoc + 1.d0)

        end do 
     end do 

  end do 
  
end subroutine gamma_matrix
!==================================================== 
!====================================================
subroutine gamma_matrix_three_body(gam,rho,THREEBOD,TB_MONO,jbas) 
  ! hartree fock potential matrix
  implicit none
  
  type(three_body_force) :: THREEBOD
  type(mono_3b) :: TB_MONO
  type(full_sp_block_mat) :: gam,rho
  type(sq_op) :: int
  type(spd) :: jbas
  integer :: q,r,i,jmax,lmax,n1,n2,j,JJ,n3,n4,tzrho,tzfoc
  integer :: grho,hrho,qrho,jrho,IImono,JJmono,qmono,x1,x2,N
  integer :: g1rho,h1rho,q1rho,j1rho,jtot,n22,n44,aux
  integer :: g2rho,h2rho,q2rho,j2rho,lrho,PAR,jfoc,lfoc,TZ
  real(8) :: sm,sm_x,den1,den2
  
  N = size(jbas%con)
!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(rho,int,threebod,gam,jbas)
  do q = 1, gam%blocks

     gam%blkM(q)%matrix = 0.d0 
     jfoc = rho%blkM(q)%lmda(2)
     
     ! loop over states in this block
     do i = 1, gam%map(q)
        do j = i,gam%map(q) 
        
           n1 = rho%blkM(q)%states(i) 
           n3 = rho%blkM(q)%states(j) 
           
           sm = 0.d0
           ! sum over density matrices
           do q1rho =  1, rho%blocks              
              j1rho = rho%blkM(q1rho)%lmda(2) 
              do q2rho = 1,rho%blocks
                 j2rho = rho%blkM(q2rho)%lmda(2) 
                 
                 ! sum over elements of the block
                 do g1rho = 1,rho%map(q1rho) 
                    do h1rho = 1,rho%map(q1rho) 
                       
                       den1 = rho%blkM(q1rho)%matrix(g1rho,h1rho)
                       if ( abs(den1) <1e-6) cycle 
                       
                       n2 = rho%blkM(q1rho)%states(g1rho) 
                       n4 = rho%blkM(q1rho)%states(h1rho) 
                          
                       ! sum over elements of the block
                       do g2rho = 1,rho%map(q2rho) 
                          do h2rho = 1,rho%map(q2rho) 
                             den2 = rho%blkM(q2rho)%matrix(g2rho,h2rho)
                             if ( abs(den2) <1e-6) cycle       
                             
                             n22 = rho%blkM(q2rho)%states(g2rho) 
                             n44 = rho%blkM(q2rho)%states(h2rho) 
                             
                             ! calculated monopole matrix elements on the fly
                             x1 = N*N*(n2-1)+N*(n22-1)+n1
                             x2 = N*N*(n4-1)+N*(n44-1)+n3 
                             qmono = TB_MONO%hash(x1,1)
                             IImono = TB_MONO%hash(x1,2)
                             JJmono = TB_MONO%hash(x2,2)
                             
                             if ( IImono > JJmono ) then 
                                aux = bosonic_tp_index(JJmono,IImono,TB_MONO%dm(qmono))
                             else
                                aux = bosonic_tp_index(IImono,JJmono,TB_MONO%dm(qmono))
                             end if 
                             
                             IF (TB_MONO%mat(qmono)%RR(aux) < -99998.0) then
                                ! first step in here... (takes a long time) 
                                ! sum over allowed JJ values
                                sm_x = 0.d0 
                                
                                do JJ = abs(j1rho - j2rho),j1rho+j2rho,2                                
                                   do jtot = abs(JJ-jfoc) , JJ+jfoc , 2 
                                   
                                      sm_x = sm_x + (jtot + 1.d0) &
                                           * GetME_pn(JJ,JJ,jtot,n2,n22,n1,n4,n44,n3,THREEBOD,jbas)
                                      
                                   end do
                                end do
                                
                                TB_MONO%mat(qmono)%RR(aux) = sm_x
                             else
                                !additional steps out here. (fast)
                                sm_x = TB_MONO%mat(qmono)%RR(aux)
                             end if 
                             
                             sm = sm + 0.5*den1*den2*sm_x/(j1rho+1.d0)/(j2rho+1.d0) 
                             
                          end do
                       end do
                    end do
                 end do
              end do
           end do
           
           ! the matrix elements are multiplied by (2J+1)/(2j+1)/(2j'+1) 
           ! which amounts to averaging over J 
          
           gam%blkM(q)%matrix(i,j) = sm/(jfoc + 1.d0)
           gam%blkM(q)%matrix(j,i) = sm/(jfoc + 1.d0)

        end do 
     end do 

  end do 
!$OMP END PARALLEL DO  
end subroutine gamma_matrix_three_body
!===========================================================
!===========================================================
subroutine density_matrix(rho,D,DX,jbas)
  ! density matrix is scaled by (2j+1)
  implicit none 
  
  type(full_sp_block_mat) :: rho,D,DX 
  type(spd) :: jbas
  integer :: i,j,q,N
  
  
  do q = 1, rho%blocks
     N = rho%map(q)
     if (N == 0) cycle 
                     
    ! zero out eigenvectors which aren't holes.
    do i = 1, N
       DX%blkM(q)%matrix(:,i) = D%blkM(q)%matrix(:,i) *&
            jbas%con(rho%blkM(q)%states(i)) 
    end do 
    
    ! rho = Dx*Dx^T 
    call dgemm('N','T',N,N,N,al,Dx%blkM(q)%matrix,N,&
         Dx%blkM(q)%matrix,N,bet,rho%blkM(q)%matrix,N)
    rho%blkM(q)%matrix = rho%blkM(q)%matrix * (rho%blkM(q)%lmda(2)+1)
    
    
  end do 
 
end subroutine density_matrix
!===========================================================
!===========================================================
real(8) function e_HF(T,V,V3,jbas)
  ! calculate the hartree fock energy
  implicit none 
  
  real(8) :: sm
  integer :: q,i
  type(spd) :: jbas
  type(full_sp_block_mat) :: T,V,V3
  
 
  sm = 0.d0

  do q = 1, T%blocks
     ! sum over eigenvalues scaled by degeneracy
     do i = 1, T%map(q) 
        sm = sm + (T%blkM(q)%matrix(i,i) + &
        0.5d0 * V%blkM(q)%matrix(i,i)+ &
        (1.d0/3.d0) * V3%blkM(q)%matrix(i,i)) &
             * jbas%con(T%blkm(q)%states(i)) * &
             (T%blkM(q)%lmda(2) + 1.d0) 
                  
     end do   
 end do 

 
 e_HF = sm
      
end function e_HF
!===========================================================
!===========================================================
subroutine transform_1b_to_HF(D,Dx,F,H,jbas,T,V,V3) 
  ! typical transformation, remap to fancy array
  implicit none 
  
  type(sq_op) :: H
  type(spd) :: jbas
  type(full_sp_block_mat) :: D,F,Dx
  type(full_sp_block_mat),optional :: T,V,V3 ! for HF calculation
  integer :: q,dm,i,j,a,b,c1,c2,cx
  
  do q = 1, F%blocks
     
   
     dm = F%map(q)
     if (dm == 0)  cycle
     
     ! transform the Fock matrix
     call dgemm('N','N',dm,dm,dm,al,F%blkM(q)%matrix&
          ,dm,D%blkM(q)%matrix,dm,bet,Dx%blkM(q)%matrix,dm) 
     call dgemm('T','N',dm,dm,dm,al,D%blkM(q)%matrix&
          ,dm,Dx%blkM(q)%matrix,dm,bet,F%blkM(q)%matrix,dm) 
     
     if (present(T)) then 
     ! transform the KE matrix
     call dgemm('N','N',dm,dm,dm,al,T%blkM(q)%matrix&
          ,dm,D%blkM(q)%matrix,dm,bet,Dx%blkM(q)%matrix,dm) 
     call dgemm('T','N',dm,dm,dm,al,D%blkM(q)%matrix&
          ,dm,Dx%blkM(q)%matrix,dm,bet,T%blkM(q)%matrix,dm) 
     end if 

     if (present(V)) then 
     ! transform the KE matrix
     call dgemm('N','N',dm,dm,dm,al,V%blkM(q)%matrix&
          ,dm,D%blkM(q)%matrix,dm,bet,Dx%blkM(q)%matrix,dm) 
     call dgemm('T','N',dm,dm,dm,al,D%blkM(q)%matrix&
          ,dm,Dx%blkM(q)%matrix,dm,bet,V%blkM(q)%matrix,dm) 
     end if 

     if (present(V3)) then 
     ! transform the KE matrix
     call dgemm('N','N',dm,dm,dm,al,V3%blkM(q)%matrix&
          ,dm,D%blkM(q)%matrix,dm,bet,Dx%blkM(q)%matrix,dm) 
     call dgemm('T','N',dm,dm,dm,al,D%blkM(q)%matrix&
          ,dm,Dx%blkM(q)%matrix,dm,bet,V3%blkM(q)%matrix,dm) 
     end if 

     do a=1,dm
        do b=a,dm
           
           i = F%blkM(q)%states(a) 
           j = F%blkM(q)%states(b) 
           
           c1 = jbas%con(i) 
           c2 = jbas%con(j) 
           cx = c1 + c2 
           
           ! fancy array remap ( normal ordered now ) 
     select case (cx) 
        case(0) 
           H%fpp(i-hb4(i),j-hb4(j)) = &
                F%blkM(q)%matrix(a,b) 
           H%fpp(j-hb4(j),i-hb4(i)) = &
                F%blkM(q)%matrix(a,b) 
        case(1) 
           if (c2 > c1) then 
              H%fph(i-hb4(i),j-pb4(j)) = &
                F%blkM(q)%matrix(a,b)   
           else
              H%fph(j-hb4(j),i-pb4(i)) = &
                F%blkM(q)%matrix(b,a)
           end if 
        case(2) 
           H%fhh(i-pb4(i),j-pb4(j)) = &
                F%blkM(q)%matrix(a,b)
           H%fhh(j-pb4(j),i-pb4(i)) = &
                F%blkM(q)%matrix(a,b)
     end select
         
         end do 
     end do 
     
  end do 
    
end subroutine transform_1b_to_HF
!===========================================================
!===========================================================  
subroutine transform_2b_to_HF(D,H,jbas)
  ! map TBME into HF basis.  
  ! huge mess. It's way too much of a hastle if you don't
  ! put everything into large square matrices and work from there
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H
  type(full_sp_block_mat) :: D 
  integer :: q,i,j,a,b,nh,np,nb,nt,II,JJ,int1,int2
  integer :: j_i,j_j,j_a,j_b,li,lj,la,lb,ti,tj,ta,tb
  integer :: q1,q2,q3,q4
  integer,allocatable,dimension(:,:) :: qnbig
  real(8),allocatable,dimension(:,:) :: Vfull,Cfull,temp,Dsmall
  real(8),allocatable,dimension(:,:) :: V1,V2,V3,V4,crevfull
  
  print*, 'Transforming TBME to HF basis' 
  ! construct full D matrix (with zeros) 
  allocate(Dsmall(H%Nsp,H%Nsp)) 
  Dsmall = 0.d0 
  do q = 1, D%blocks
     if (D%map(q) == 0) cycle
     
     do i = 1, D%map(q) 
        do j = 1,D%map(q) 
           
           Dsmall(D%blkM(q)%states(i),D%blkM(q)%states(j)) = &
            D%blkM(q)%matrix(i,j) 
        end do 
     end do 
  end do 

  ! cycle over major blocks of the tp hamiltonian
  do q = 1, H%nblocks
     
     np = H%mat(q)%npp
     nh = H%mat(q)%nhh
     nb = H%mat(q)%nph
     nt = np+nh+nb
     if (nt == 0) cycle
     
     allocate(Vfull(nt,nt)) 
     allocate(V1(nt,nt),V2(nt,nt),V3(nt,nt),V4(nt,nt)) 
     allocate(Cfull(nt,nt)) 
     allocate(Crevfull(nt,nt))
     allocate(temp(nt,nt)) 
     allocate(qnbig(nt,2)) 
     
     ! mapping things out to a square matrix 
     
     Vfull(1:nh,1:nh) = H%mat(q)%gam(5)%X 
     
     Vfull(nh+1:nb+nh,1:nh) = H%mat(q)%gam(6)%X
     Vfull(1:nh,nh+1:nb+nh) = Transpose(H%mat(q)%gam(6)%X) 
     
     Vfull(nh+nb+1:nt,1:nh) = H%mat(q)%gam(3)%X 
     Vfull(1:nh,nh+nb+1:nt) = Transpose(H%mat(q)%gam(3)%X) 
     
     Vfull(nh+1:nb+nh,nh+1:nb+nh) = H%mat(q)%gam(4)%X
     
     Vfull(nh+nb+1:nt,nh+1:nh+nb) = H%mat(q)%gam(2)%X
     Vfull(nh+1:nh+nb,nh+nb+1:nt) = Transpose(H%mat(q)%gam(2)%X)
           
     Vfull(nh+nb+1:nt,nh+nb+1:nt) = H%mat(q)%gam(1)%X 
     
     qnbig(1:nh,:) = H%mat(q)%qn(3)%Y
     qnbig(nh+1:nh+nb,:) = H%mat(q)%qn(2)%Y
     qnbig(nh+nb+1:nt,:) = H%mat(q)%qn(1)%Y
     
     ! filling the C matrices
   
     do II = 1,nt
           
        i = qnbig(II,1)
        j = qnbig(II,2)
          
        do JJ = 1, nt
           
           a = qnbig(JJ,1)
           b = qnbig(JJ,2)
                              
           Cfull(JJ,II) = Dsmall(a,i)*Dsmall(b,j) * &
                sqrt(1.d0 + kron_del(a,b)) /sqrt(1.d0 + kron_del(i,j))  
           
           Crevfull(JJ,II) = Dsmall(b,i)*Dsmall(a,j) *&
                (1 - kron_del(a,b)) * &
           (-1)**( (jbas%jj(a) + jbas%jj(b)) /2 ) * &
           sqrt(1.d0 + kron_del(a,b)) /sqrt(1.d0 + kron_del(i,j)) 
         
           
        end do  
     end do 
        
     ! do the transformation
     call dgemm('N','N',nt,nt,nt,al,Vfull,nt,Cfull,nt,bet,temp,nt) 
     call dgemm('T','N',nt,nt,nt,al,Cfull,nt,temp,nt,bet,V1,nt) 
     
     call dgemm('N','N',nt,nt,nt,al,Vfull,nt,Crevfull,nt,bet,temp,nt) 
     call dgemm('T','N',nt,nt,nt,al,Cfull,nt,temp,nt,bet,V2,nt) 
     
     call dgemm('N','N',nt,nt,nt,al,Vfull,nt,Cfull,nt,bet,temp,nt) 
     call dgemm('T','N',nt,nt,nt,al,Crevfull,nt,temp,nt,bet,V3,nt) 
     
     call dgemm('N','N',nt,nt,nt,al,Vfull,nt,Crevfull,nt,bet,temp,nt) 
     call dgemm('T','N',nt,nt,nt,al,Crevfull,nt,temp,nt,bet,V4,nt)
     
     Vfull = V1 - (-1)**(H%mat(q)%lam(1)/2)*(V2+V3) + V4  
   
     ! put back into the fancy array
     H%mat(q)%gam(5)%X = Vfull(1:nh,1:nh) 
     H%mat(q)%gam(6)%X = Vfull(nh+1:nb+nh,1:nh) 
     H%mat(q)%gam(3)%X = Vfull(nh+nb+1:nt,1:nh) 
     H%mat(q)%gam(4)%X = Vfull(nh+1:nb+nh,nh+1:nb+nh)
     H%mat(q)%gam(2)%X = Vfull(nh+nb+1:nt,nh+1:nh+nb)      
     H%mat(q)%gam(1)%X = Vfull(nh+nb+1:nt,nh+nb+1:nt)
 
     ! deallocate everything
     deallocate(Vfull)
     deallocate(V1,V2,V3,V4)
     deallocate(Cfull)
     deallocate(Crevfull)
     deallocate(qnbig) 
     deallocate(temp)
     
  end do 
end subroutine     
!==============================================================
!==============================================================
subroutine transform_1b_to_HF_tensor(D,O1,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: O1
  type(full_sp_block_mat) :: D 
  integer :: q,i,j,a,b,nt,c1,c2,cx
  real(8),allocatable,dimension(:,:) :: Dfull,Ffull,Dx

  nt = O1%Nsp
  allocate(Dfull(nt,nt)) 
  Dfull = 0.d0 
  do q = 1, D%blocks
     if (D%map(q) == 0) cycle
     
     do i = 1, D%map(q) 
        do j = 1,D%map(q) 
           
           Dfull(D%blkM(q)%states(i),D%blkM(q)%states(j)) = &
            D%blkM(q)%matrix(i,j) 
        
        end do
     end do 
  end do 
  
  allocate(Ffull(nt,nt)) 
  allocate(Dx(nt,nt)) 
  do i = 1, nt
     do j = 1,nt      
        Ffull(i,j) = f_tensor_elem(i,j,O1,jbas) 
     end do 
  end do 
  
  ! transform the Fock matrix
  call dgemm('N','N',nt,nt,nt,al,Ffull,nt,Dfull,nt,bet,Dx,nt) 
  call dgemm('T','N',nt,nt,nt,al,Dfull,nt,Dx,nt,bet,Ffull,nt)
  
  do i=1,nt
     do j=1,nt
                
        c1 = jbas%con(i) 
        c2 = jbas%con(j) 
        cx = c1 + c2 
           
           ! fancy array remap ( normal ordered now ) 
        select case (cx) 
        case(0) 
           O1%fpp(i-hb4(i),j-hb4(j)) = &
                Ffull(i,j) 
        case(1) 
           if (c2 > c1) then 
              O1%fph(i-hb4(i),j-pb4(j)) = &
                   Ffull(i,j)  
           end if
        case(2) 
           O1%fhh(i-pb4(i),j-pb4(j)) = &
                Ffull(i,j)
        end select
        
     end do
  end do  


end subroutine transform_1b_to_HF_tensor
!=======================================================================
!=======================================================================
subroutine meanfield_2b(rho,H,TB,jbas) 
  ! normal ordering the three body force into the two body force. 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H
  type(full_sp_block_mat) :: rho
  type(three_body_force) :: TB
  integer :: q, II,JJ ,i,j,k,l,g,qrho,hrho,grho,n3,n6
  integer :: bigJ,jtot,jmin,jmax,jrho,c1,c2,jxstart
  real(8) :: sm,sm_x,den,pre1,pre2
  logical :: square

  print*, 'Normal ordering three-body force...'

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(rho,TB,jbas,H)
  do q = 1, H%nblocks
     bigJ = H%mat(q)%lam(1) 
     do g = 1, 6
       
        c1 = sea1(g) 
        c2 = sea2(g) 
        square = sqs(g) 
        jxstart = jst(g) 
          
        do II = 1,size(H%mat(q)%gam(g)%X(:,1))
           
           i = H%mat(q)%qn(c1)%Y(II,1)
           j = H%mat(q)%qn(c1)%Y(II,2)
           pre1 = 1.d0
           if (i == j ) pre1 = 1/sqrt(2.d0)
           
           do JJ = min(jxstart,II),size(H%mat(q)%gam(g)%X(1,:))  
              
              k = H%mat(q)%qn(c2)%Y(JJ,1)
              l = H%mat(q)%qn(c2)%Y(JJ,2)
              pre2 = 1.d0
              if (k == l ) pre2 = 1/sqrt(2.d0)
           
              sm = 0.d0 
                            
              do qrho =  1, rho%blocks
                 
                 jrho = rho%blkM(qrho)%lmda(2) 
                  
                 jmin = abs(bigJ - jrho) 
                 jmax = bigJ + jrho 
!                 if (jmin > jmax) cycle
                 ! sum over elements of the block
                 do grho = 1,rho%map(qrho) 
                    do hrho = 1,rho%map(qrho) 
                       
                       den = rho%blkM(qrho)%matrix(grho,hrho)
                       
                       if (abs(den) < 1e-6) cycle
                       
                       n3 = rho%blkM(qrho)%states(grho) 
                       n6 = rho%blkM(qrho)%states(hrho) 
                       
                       sm_x = 0.d0 
                       do jtot = jmin,jmax,2
                          sm_x = sm_x + (jtot+1.d0) * &
                               GetME_pn(bigJ,bigJ,jtot,i,j,n3,k,l,n6,TB,jbas) 
                       end do 
                       
                       sm = sm + den * sm_x /(jrho+1.d0) ! again we divide by jrho
                       ! because of how I've defined the density matrix.
                    end do
                 end do
              end do
              H%mat(q)%gam(g)%X(II,JJ)=H%mat(q)%gam(g)%X(II,JJ) + sm /(bigJ+1.d0)*pre1*pre2
              if (square) H%mat(q)%gam(g)%X(JJ,II) = H%mat(q)%gam(g)%X(II,JJ) 
           end do
        end do
     end do
  end do
!$OMP END PARALLEL DO
end subroutine meanfield_2b
!================================================================================================
!================================================================================================
end module HF_mod
                
              
                 
                 
           
            

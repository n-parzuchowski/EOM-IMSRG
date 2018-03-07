 module isospin_operators
  use basic_IMSRG
  
  TYPE :: iso_block
     integer :: lam(4) ! specifices J,Par,Tz1,Tz2 of the block
     integer :: Jpair(2) ! for tensor operators that have rank > 0  
     real(8),allocatable,dimension(:,:) :: Xpphh 
     integer :: npp1,nhh2,ntot! dimensions
     type(int_mat),dimension(2) :: qn
  END TYPE iso_block

  type iso_ladder
     type(iso_block),allocatable,dimension(:) :: tblck
     real(8),allocatable,dimension(:,:) :: fph
     integer :: nblocks,Aprot,Aneut,Nsp,herm,belowEF,neq
     integer :: Rank,dpar,dTz,eMax,lmax,xindx
     real(8) :: E0,hospace,lawson_beta,com_hw 
     logical :: pphh_ph
     character(2) :: trans_label
  end type iso_ladder

  TYPE :: iso_op_block
     integer :: lam(4) ! specifices J,Par,Tz1,Tz2 of the block
     integer :: Jpair(2) ! for tensor operators that have rank > 0  
     type(real_mat),dimension(9) :: tgam  
     integer :: npp1,npp2,nhh1,nph1,nph2,nhh2,ntot! dimensions
     type(int_mat),dimension(3,2) :: tensor_qn
  END TYPE iso_op_block 

  type iso_operator
     type(iso_op_block),allocatable,dimension(:) :: tblck
     real(8),allocatable,dimension(:,:) :: fock
     integer :: nblocks,Aprot,Aneut,Nsp,herm,belowEF,neq
     integer :: Rank,dpar,dTz,eMax,lmax,xindx
     real(8) :: E0,hospace,lawson_beta,com_hw 
     logical :: pphh_ph
     character(2) :: trans_label
  end type iso_operator
    
contains
!===============================================================
!===============================================================  
real(8) function iso_frob_norm(op) 
  implicit none 
  
  type(iso_operator) :: op
  integer :: q,g
  real(8) :: sm

  sm = sum(op%fock**2)
 
  do q = 1, op%nblocks
        do g = 1,9
           sm = sm + sum(op%tblck(q)%tgam(g)%X**2)
        end do
  end do 
  
  iso_frob_norm = sqrt(sm)
end function 
!=====================================================
!=====================================================
subroutine append_isospin_operator(BB,bx,AA) 
  ! make a copy of H onto op
  implicit none 
  
  type(iso_operator) :: AA,BB
  real(8) :: bx
  integer :: q,i,j,holes,parts,nh,np,nb
  
  AA%E0 = BB%E0*bx+AA%E0
  AA%fock = BB%fock*bx + AA%fock 
  
  do q = 1, AA%nblocks
     
     do i = 1,9
        AA%tblck(q)%tgam(i)%X = bx*BB%tblck(q)%tgam(i)%X + &
             AA%tblck(q)%tgam(i)%X
     end do
       
  end do
    
end subroutine append_isospin_operator
!=====================================================
!=====================================================
subroutine copy_isospin_operator(H,op) 
  ! make a copy of H onto op
  implicit none 
  
  type(iso_operator) :: H,op
  integer :: q,i,j,holes,parts,nh,np,nb
  
  op%herm = H%herm
  op%E0 = H%E0
  op%fock = H%fock
  
  do q = 1, op%nblocks
     
     do i = 1,9
        op%tblck(q)%tgam(i)%X = H%tblck(q)%tgam(i)%X
     end do
       
  end do
    
end subroutine copy_isospin_operator
!=====================================================
!=====================================================
subroutine clear_isospin_operator(op) 
  ! make a copy of H onto op
  implicit none 
  
  type(iso_operator) :: op
  integer :: q,i,j,holes,parts,nh,np,nb
  
  op%E0 = 0.d0
  op%fock = 0.d0
  
  do q = 1, op%nblocks
     
     do i = 1,9
        op%tblck(q)%tgam(i)%X = 0.d0
     end do
     
  end do
  
end subroutine clear_isospin_operator
!=========================================================
!=========================================================
subroutine duplicate_isospin_operator(A,B)
  ! I'm just making an empty copy of the operator, which carries all of it's parameters
  
  implicit none

  type(iso_operator) :: A,B
  integer :: q,i,n1,n2

  B%rank = A%rank
  B%dpar = A%dpar
  B%dTz = A%dTz
  B%neq = A%neq
  B%eMax = A%eMax
  B%lmax = A%lmax
  B%xindx = A%xindx
  B%herm = A%herm
  B%Aprot = A%Aprot
  B%Aneut = B%Aneut
  B%Nsp = A%Nsp
  B%belowEF = A%belowEF
  B%nblocks = A%nblocks
  B%trans_label = A%trans_label

  allocate(B%fock(B%Nsp, B%Nsp ))
  allocate(B%tblck(B%nblocks))
  B%fock = 0.d0
  
  do q = 1, B%nblocks
     B%tblck(q)%lam = A%tblck(q)%lam 
     B%tblck(q)%npp1 = A%tblck(q)%npp1
     B%tblck(q)%nph1 = A%tblck(q)%nph1
     B%tblck(q)%nhh1 = A%tblck(q)%nhh1
     B%tblck(q)%nhh2 = A%tblck(q)%nhh2
     B%tblck(q)%nph2 = A%tblck(q)%nph2
     B%tblck(q)%npp2 = A%tblck(q)%npp2
     B%tblck(q)%ntot = A%tblck(q)%ntot 
     B%tblck(q)%Jpair = A%tblck(q)%Jpair

     do i = 1, 9
        n1 = size( A%tblck(q)%tgam(i)%X(:,1))
        n2 = size( A%tblck(q)%tgam(i)%X(1,:))
        allocate(B%tblck(q)%tgam(i)%X(n1,n2))
        B%tblck(q)%tgam(i)%X = 0.d0
     end do
        
     do i = 1, 3
        allocate(B%tblck(q)%tensor_qn(i,1)%Y(size(A%tblck(q)%tensor_qn(i,1)%Y(:,1)),2))
        allocate(B%tblck(q)%tensor_qn(i,2)%Y(size(A%tblck(q)%tensor_qn(i,2)%Y(:,1)),2))        
        B%tblck(q)%tensor_qn(i,1)%Y = A%tblck(q)%tensor_qn(i,1)%Y
        B%tblck(q)%tensor_qn(i,2)%Y = A%tblck(q)%tensor_qn(i,2)%Y
     end do
  end do

end subroutine duplicate_isospin_operator
!=========================================================
!=========================================================
subroutine duplicate_isospin_ladder(A,B)
  implicit none

  type(iso_ladder) :: A,B
  integer :: q 

  B%rank = A%rank
  B%dpar = A%dpar
  B%dTz = A%dTz
  B%neq = A%neq
  B%eMax = A%eMax
  B%lmax = A%lmax
  B%xindx = A%xindx
  B%herm = A%herm
  B%Aprot = A%Aprot
  B%Aneut = B%Aneut
  B%Nsp = A%Nsp
  B%belowEF = A%belowEF
  B%nblocks = A%nblocks
  B%trans_label = A%trans_label

  allocate(B%fph(B%Nsp - B%belowEF , B%belowEF ))
  allocate(B%tblck(B%nblocks))
  B%fph = 0.d0
  
  do q = 1, B%nblocks
     B%tblck(q)%lam = A%tblck(q)%lam 
     B%tblck(q)%npp1 = A%tblck(q)%npp1
     B%tblck(q)%nhh2 = A%tblck(q)%nhh2
     B%tblck(q)%ntot = A%tblck(q)%ntot 
     B%tblck(q)%Jpair = A%tblck(q)%Jpair

     allocate(B%tblck(q)%qn(1)%Y(B%tblck(q)%npp1,2))
     allocate(B%tblck(q)%qn(2)%Y(B%tblck(q)%nhh2,2))
     allocate(B%tblck(q)%Xpphh(B%tblck(q)%npp1,B%tblck(q)%nhh2))
     B%tblck(q)%Xpphh = 0.d0 
     B%tblck(q)%qn(1)%Y = A%tblck(q)%qn(1)%Y
     B%tblck(q)%qn(2)%Y = A%tblck(q)%qn(2)%Y  
  end do

end subroutine duplicate_isospin_ladder
!=========================================================
!=========================================================
subroutine allocate_isospin_ladder(jbas,op,zerorank)
  implicit none
  
  type(spd) :: jbas
  type(iso_ladder) :: op
  type(sq_op) :: zerorank 
  integer :: N,AX,q,i,j,j1,j2,tz1,tz2,l1,l2,x,rank
  integer :: Jtot1,Jtot2,Par1,Par2,nph1,npp1,nhh1,nph2,npp2,nhh2
  integer :: CX,j_min,j_max,numJ,q1,q2,dTz , Tz
  real(8) :: d6ji,lwake,size,mem

  ! rank is multiplied by 2 as well
  rank = op%rank 
  op%hospace = zerorank%hospace
  if ( mod(rank,2) == 1 ) STOP 'RANK MUST BE MULTIPLIED BY 2'
  if ( mod(op%dpar,2) == 1) STOP 'dPAR MUST BE MULTIPLIED BY 2'
  dTz = op%dTz 
  mem = 0.d0 
  AX = sum(jbas%con) !number of shells under eF 
  op%belowEF = AX
  op%Nsp = jbas%total_orbits
  N = op%Nsp  !number of sp shells
  op%Aprot = zerorank%Aprot
  op%Aneut = zerorank%Aneut

  do i = 1,N
     do j = i,N
        
        j_min = abs(jbas%jj(i)-jbas%jj(j))
        j_max = jbas%jj(i) + jbas%jj(j) 
  
        numJ = (j_max - j_min)/2 + 2
  
        x = bosonic_tp_index(i,j,N) 

        if (.not. allocated(jbas%xmap_tensor(op%xindx,x)%Z)) then 
           allocate(jbas%xmap_tensor(op%xindx,x)%Z(numJ)) 
        end if 
        jbas%xmap_tensor(op%xindx,x)%Z = 0
        jbas%xmap_tensor(op%xindx,x)%Z(1) = j_min
        
     end do 
  end do 
  
  ! quantum numbers of the last block 
  Tz1 = 1 
  Par1 = 1
  Jtot1 = jbas%Jtotal_max*2 
  Jtot2 = Jtot1+ RANK 

  op%nblocks =  iso_ladder_block_index(Jtot1,Jtot2,RANK,Tz1,Par1)               
 
  allocate(op%tblck(op%nblocks)) 
 
  allocate(op%fph(N-AX,AX))
  mem = mem + sizeof(op%fph)
  op%fph = 0.d0
  
  q = 1
  do Jtot1 = 0,2*jbas%Jtotal_max,2 
     do Jtot2 = abs(Jtot1 - rank), Jtot1+rank, 2   
       
        do Tz1 = -1,1
           do Par1 = 0,1
                    
              if (mod(op%dpar/2,2) == 1) then ! this only works for EX transitions
                 Par2 = abs(Par1-1)  
              else
                 Par2 = Par1 
              end if 
              
              op%tblck(q)%Jpair(1) = Jtot1
              op%tblck(q)%Jpair(2) = Jtot2

              Tz2 = Tz1 - dTz 
              op%tblck(q)%lam(1) = 1 ! phase instead of J 
              op%tblck(q)%lam(2) = Par1 !just remember that they change if the operator has odd parity.
              op%tblck(q)%lam(3) = Tz1
              op%tblck(q)%lam(4) = Tz2

              ! we already figured this stuff out for 
              ! the rank 0 case, so lets re-use it

              
              if ( (abs(Tz2) > 1) .or. (Jtot2 > Jbas%Jtotal_max*2) ) then 
                 npp1=0
                 nhh2=0
              else
                 q1 = block_index(Jtot1,Tz1,Par1) 
                 q2 = block_index(Jtot2,Tz2,Par2) 
                 npp1 = zerorank%mat(q1)%npp
                 nhh2 = zerorank%mat(q2)%nhh
              end if 
              
              op%tblck(q)%npp1 = npp1
              op%tblck(q)%nhh2 = nhh2              
                            
              allocate(op%tblck(q)%qn(1)%Y(npp1,2)) !qnpp1
              allocate(op%tblck(q)%qn(2)%Y(nhh2,2)) !qnhh2
              
              if ( (abs(Tz2) .le. 1) .and. (Jtot2 .le. Jbas%Jtotal_max*2) ) then 
                 
              ! yeah this is a mess but it really facilitates the 
              ! commutators

              ! blocks such as pphh ALWAYS have pp first then hh
              ! looks scary but i'm just copying everything from the
                 ! rank zero operators we already have
                 op%tblck(q)%qn(1)%Y = zerorank%mat(q1)%qn(1)%Y
                 op%tblck(q)%qn(2)%Y = zerorank%mat(q2)%qn(3)%Y
              end if 

              allocate(op%tblck(q)%Xpphh(npp1,nhh2)) !Xpphh

              mem = mem + sizeof(op%tblck(q)%Xpphh) 
              
              op%tblck(q)%Xpphh = 0.d0
              
              q = q + 1
            
            end do
         end do
         
      end do
   end do

   
   print*, 'Isospin Ladder Operator Memory =', mem/1024./1024./1024. , 'GB' 

   ! i need to fill this last array for tensor_elem
   do q = 1, zerorank%nblocks
         Jtot1 = zerorank%mat(q)%lam(1)
         Tz = zerorank%mat(q)%lam(3)
         Par1 = zerorank%mat(q)%lam(2)

         npp1 = 0 ; nhh1 = 0
         do i = 1,N   !looping over sp states
            do j = i,N
               
               ! check if the two sp states can exist in this block 
               tz1 = jbas%itzp(i)
               tz2 = jbas%itzp(j) 
               if ( Tz .ne. (tz1 + tz2)/2 ) cycle
         
               l1 = jbas%ll(i) 
               l2 = jbas%ll(j)
               if ( Par1 .ne. (1 - (-1)**(l1 + l2))/2 ) cycle
         
               j1 = jbas%jj(i) 
               j2 = jbas%jj(j) 
               if (.not. triangle(j1,j2,Jtot1) ) cycle
                                                      
               cX = jbas%con(i) + jbas%con(j)
               
               x = bosonic_tp_index(i,j,N) 
               j_min = jbas%xmap_tensor(op%xindx,x)%Z(1) 
     
      
               select case (CX)
               case (0) 
                  npp1 = npp1 + 1
                  jbas%xmap_tensor(op%xindx,x)%Z((Jtot1-j_min)/2+2) = npp1 
               case (1)                      
                  nph1 = nph1 + 1
                  jbas%xmap_tensor(op%xindx,x)%Z((Jtot1-j_min)/2+2) = nph1
               case (2) 
                  nhh1 = nhh1 + 1
                  jbas%xmap_tensor(op%xindx,x)%Z((Jtot1-j_min)/2+2) = nhh1
               end select
        
            end do
         end do
      end do
      
   ! these are six j symbols that I don't already have,
   ! which the commutators need for this tensor.
   ! access with XXXsixj

      if (.not. allocated(half6j(op%xindx)%tp_mat)) then          
         call store_6j_3halfint(jbas,rank,op%xindx)
      end if
   !call divide_work_tensor(op) 

end subroutine allocate_isospin_ladder
!=========================================================
!=========================================================
subroutine allocate_isospin_operator(jbas,op,zerorank)
  implicit none
  
  type(spd) :: jbas
  type(iso_operator) :: op
  type(sq_op) :: zerorank 
  integer :: N,AX,q,i,j,j1,j2,tz1,tz2,l1,l2,x,rank
  integer :: Jtot1,Jtot2,Par1,Par2,nph1,npp1,nhh1,nph2,npp2,nhh2
  integer :: CX,j_min,j_max,numJ,q1,q2,dTz , Tz
  real(8) :: d6ji,lwake,size,mem

  ! rank is multiplied by 2 as well
  rank = op%rank 
  op%hospace = zerorank%hospace
  if ( mod(rank,2) == 1 ) STOP 'RANK MUST BE MULTIPLIED BY 2'
  if ( mod(op%dpar,2) == 1) STOP 'dPAR MUST BE MULTIPLIED BY 2'
  dTz = op%dTz 
  mem = 0.d0 
  AX = sum(jbas%con) !number of shells under eF 
  op%belowEF = AX
  op%Nsp = jbas%total_orbits
  N = op%Nsp  !number of sp shells
  op%Aprot = zerorank%Aprot
  op%Aneut = zerorank%Aneut

  do i = 1,N
     do j = i,N
        
        j_min = abs(jbas%jj(i)-jbas%jj(j))
        j_max = jbas%jj(i) + jbas%jj(j) 
  
        numJ = (j_max - j_min)/2 + 2
  
        x = bosonic_tp_index(i,j,N) 

        if (.not. allocated(jbas%xmap_tensor(op%xindx,x)%Z)) then 
           allocate(jbas%xmap_tensor(op%xindx,x)%Z(numJ)) 
        end if 
        jbas%xmap_tensor(op%xindx,x)%Z = 0
        jbas%xmap_tensor(op%xindx,x)%Z(1) = j_min
        
     end do 
  end do 
  
  ! quantum numbers of the last block 
  Tz1 = 1 
  Par1 = 1
  Jtot1 = jbas%Jtotal_max*2 
  Jtot2 = Jtot1+ RANK 

  op%nblocks =  iso_ladder_block_index(Jtot1,Jtot2,RANK,Tz1,Par1)               
 
  allocate(op%tblck(op%nblocks)) 
 
  allocate(op%fock(N,N))
  mem = mem + sizeof(op%fock)
  op%fock = 0.d0
  
  q = 1
  do Jtot1 = 0,2*jbas%Jtotal_max,2 
     do Jtot2 = abs(Jtot1 - rank), Jtot1+rank, 2   
       
        do Tz1 = -1,1
           do Par1 = 0,1
                    
              if (mod(op%dpar/2,2) == 1) then ! this only works for EX transitions
                 Par2 = abs(Par1-1)  
              else
                 Par2 = Par1 
              end if 
              
              op%tblck(q)%Jpair(1) = Jtot1
              op%tblck(q)%Jpair(2) = Jtot2

              Tz2 = Tz1 - dTz 
              op%tblck(q)%lam(1) = (-1)**((Jtot1-Jtot2)/2) ! phase instead of J 
              op%tblck(q)%lam(2) = Par1 !just remember that they change if the operator has odd parity.
              op%tblck(q)%lam(3) = Tz1
              op%tblck(q)%lam(4) = Tz2

              ! we already figured this stuff out for 
              ! the rank 0 case, so lets re-use it

              
              if ( (abs(Tz2) > 1) .or. (Jtot2 > Jbas%Jtotal_max*2) ) then 
                 npp1=0;nph1=0;nhh1=0
                 nhh2=0;nph2=0;npp2=0
              else
                 q1 = block_index(Jtot1,Tz1,Par1) 
                 q2 = block_index(Jtot2,Tz2,Par2) 
                 npp1 = zerorank%mat(q1)%npp
                 nph1 = zerorank%mat(q1)%nph
                 nhh1 = zerorank%mat(q1)%nhh

                 npp2 = zerorank%mat(q2)%npp
                 nph2 = zerorank%mat(q2)%nph
                 nhh2 = zerorank%mat(q2)%nhh
              end if 
              
              op%tblck(q)%npp1 = npp1
              op%tblck(q)%nph1 = nph1
              op%tblck(q)%nhh1 = nhh1
              op%tblck(q)%nhh2 = nhh2
              op%tblck(q)%nph2 = nph2
              op%tblck(q)%npp2 = npp2              
                            
                                          
              allocate(op%tblck(q)%tensor_qn(1,1)%Y(npp1,2)) !qnpp1
              allocate(op%tblck(q)%tensor_qn(2,1)%Y(nph1,2)) !qnph1
              allocate(op%tblck(q)%tensor_qn(3,1)%Y(nhh1,2)) !qnhh1
              
              allocate(op%tblck(q)%tensor_qn(1,2)%Y(npp2,2)) !qnpp2
              allocate(op%tblck(q)%tensor_qn(2,2)%Y(nph2,2)) !qnph2
              allocate(op%tblck(q)%tensor_qn(3,2)%Y(nhh2,2)) !qnhh2

              if ( (abs(Tz2) .le. 1) .and. (Jtot2 .le. Jbas%Jtotal_max*2) ) then 
                 
              ! yeah this is a mess but it really facilitates the 
              ! commutators

              ! blocks such as pphh ALWAYS have pp first then hh
              ! looks scary but i'm just copying everything from the
                 ! rank zero operators we already have
                 op%tblck(q)%tensor_qn(1,1)%Y = zerorank%mat(q1)%qn(1)%Y
                 op%tblck(q)%tensor_qn(2,1)%Y = zerorank%mat(q1)%qn(2)%Y
                 op%tblck(q)%tensor_qn(3,1)%Y = zerorank%mat(q1)%qn(3)%Y

                 op%tblck(q)%tensor_qn(1,2)%Y = zerorank%mat(q2)%qn(1)%Y
                 op%tblck(q)%tensor_qn(2,2)%Y = zerorank%mat(q2)%qn(2)%Y
                 op%tblck(q)%tensor_qn(3,2)%Y = zerorank%mat(q2)%qn(3)%Y
              end if 

              allocate(op%tblck(q)%tgam(3)%X(npp1,nhh2)) !Vpphh
              allocate(op%tblck(q)%tgam(7)%X(nhh1,npp2)) !Vhhpp              
              allocate(op%tblck(q)%tgam(1)%X(npp1,npp2)) !Vpppp
              allocate(op%tblck(q)%tgam(5)%X(nhh1,nhh2)) !Vhhhh              
              allocate(op%tblck(q)%tgam(4)%X(nph1,nph2)) !Vphph
              allocate(op%tblck(q)%tgam(2)%X(npp1,nph2)) !Vppph
              allocate(op%tblck(q)%tgam(6)%X(nph1,nhh2)) !Vphhh
              allocate(op%tblck(q)%tgam(8)%X(nph1,npp2)) !Vphpp
              allocate(op%tblck(q)%tgam(9)%X(nhh1,nph2)) !Vhhph

              mem = mem + sizeof(op%tblck(q)%tgam(3)%X) 
              mem = mem + sizeof(op%tblck(q)%tgam(7)%X)              
              mem = mem + sizeof(op%tblck(q)%tgam(1)%X)
              mem = mem + sizeof(op%tblck(q)%tgam(2)%X)
              mem = mem + sizeof(op%tblck(q)%tgam(4)%X)
              mem = mem + sizeof(op%tblck(q)%tgam(5)%X)
              mem = mem + sizeof(op%tblck(q)%tgam(6)%X)
              mem = mem + sizeof(op%tblck(q)%tgam(8)%X)
              mem = mem + sizeof(op%tblck(q)%tgam(9)%X)
                                        
              do i = 1,9
                 if (op%pphh_ph.and.(i.ne.3).and.(i.ne.7)) cycle 
                 op%tblck(q)%tgam(i)%X = 0.d0
              end do
              
              q = q + 1
            
            end do
         end do
         
      end do
   end do

   
   print*, 'Isospin Ladder Operator Memory =', mem/1024./1024./1024. , 'GB' 

   ! i need to fill this last array for tensor_elem
   do q = 1, zerorank%nblocks
         Jtot1 = zerorank%mat(q)%lam(1)
         Tz = zerorank%mat(q)%lam(3)
         Par1 = zerorank%mat(q)%lam(2)

         npp1 = 0 ; nhh1 = 0 ; nph1 = 0 
         do i = 1,N   !looping over sp states
            do j = i,N
               
               ! check if the two sp states can exist in this block 
               tz1 = jbas%itzp(i)
               tz2 = jbas%itzp(j) 
               if ( Tz .ne. (tz1 + tz2)/2 ) cycle
         
               l1 = jbas%ll(i) 
               l2 = jbas%ll(j)
               if ( Par1 .ne. (1 - (-1)**(l1 + l2))/2 ) cycle
         
               j1 = jbas%jj(i) 
               j2 = jbas%jj(j) 
               if (.not. triangle(j1,j2,Jtot1) ) cycle
                                                      
               cX = jbas%con(i) + jbas%con(j)
               
               x = bosonic_tp_index(i,j,N) 
               j_min = jbas%xmap_tensor(op%xindx,x)%Z(1) 
     
      
               select case (CX)
               case (0) 
                  npp1 = npp1 + 1
                  jbas%xmap_tensor(op%xindx,x)%Z((Jtot1-j_min)/2+2) = npp1 
               case (1)                      
                  nph1 = nph1 + 1
                  jbas%xmap_tensor(op%xindx,x)%Z((Jtot1-j_min)/2+2) = nph1
               case (2) 
                  nhh1 = nhh1 + 1
                  jbas%xmap_tensor(op%xindx,x)%Z((Jtot1-j_min)/2+2) = nhh1
               end select
        
            end do
         end do
      end do
      
   ! these are six j symbols that I don't already have,
   ! which the commutators need for this tensor.
   ! access with XXXsixj

      if (.not. allocated(half6j(op%xindx)%tp_mat)) then          
         call store_6j_3halfint(jbas,rank,op%xindx)
      end if
   !call divide_work_tensor(op) 

end subroutine allocate_isospin_operator
!==================================================================  
!==================================================================
real(8) function iso_ladder_elem(a,b,c,d,J1,J2,op,jbas) 
  ! grabs the matrix element you are looking for
  implicit none
  
  integer :: a,b,c,d,J1,J2,rank,T,P,q,qx,c1,c2,N,J1x,J2x
  integer :: int1,int2,i1,i2,j_min,x,k1,k2,ax,bx,cx,dx
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  ,phase
  logical :: fail_c,switch
  type(iso_ladder) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c


  if (jbas%con(a)+jbas%con(b) .ne. 0) then
     iso_ladder_elem = 0.d0
     return
  end if

  if (jbas%con(c)+jbas%con(d) .ne. 2) then
     iso_ladder_elem = 0.d0
     return

  end if
  
  !make sure the matrix element exists first
 
  rank = op%rank

  if ( .not. (triangle ( J1,J2,rank ))) then 
     iso_ladder_elem = 0.d0 
     return
  end if 
 
  
  fail_c = .true. 
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  
  if ( .not. ((triangle(ja,jb,J1)) .and. (triangle (jc,jd,J2))) ) then 
     iso_ladder_elem = 0.d0
     return
  end if
     
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)

  P = mod(la + lb,2) 
     
  if ( mod(lc + ld,2) .ne. abs(P - ((-1)**(op%dpar/2+1)+1)/2) ) then
    iso_ladder_elem = 0.d0 
    return
  end if 
        
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)
     
  T = (ta + tb)/2
     
  if ((tc+td) .ne. 2*(T-op%dTz)) then     
    iso_ladder_elem = 0.d0
    return
  end if 

  q = iso_ladder_block_index(J1,J2,rank,T,P) 
  
  ! see subroutine "allocate_blocks" for mapping from qx to each 
  ! of the 6 storage arrays
    
  pre = 1 
  
  N = op%Nsp
  ! get the indeces in the correct order
  if ( a > b )  then 
     x = bosonic_tp_index(b,a,N) 
     j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
     i1 = jbas%xmap_tensor(op%xindx,x)%Z( (J1-j_min)/2 + 2) 
     pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J1)/2 ) 
  else
     if (a == b) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(a,b,N)
     j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
     i1 = jbas%xmap_tensor(op%xindx,x)%Z( (J1-j_min)/2 + 2) 
  end if 

  if (c > d)  then     
     x = bosonic_tp_index(d,c,N) 
     j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
     i2 = jbas%xmap_tensor(op%xindx,x)%Z( (J2-j_min)/2 + 2) 
     pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J2)/2 ) 
  else 
     if (c == d) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(c,d,N) 
     j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
     i2 = jbas%xmap_tensor(op%xindx,x)%Z( (J2-j_min)/2 + 2)  
  end if
 
  ! grab the matrix element

  ! right now i1 and i2 still refer to where the pair is located
  ! in the rank zero qn storage

 ! if ((a==23).and.(b==28).and.(c==1).and.(d==5).and.(J1==8).and.(J2==4))then
!     print*, q,i1,i2,op%tblck(q)%Xpphh(i1,i2) * pre
 ! end if
 
  iso_ladder_elem = op%tblck(q)%Xpphh(i1,i2) * pre
    
 end function iso_ladder_elem
!==================================================================  
!==================================================================
real(8) function iso_op_elem(a,b,c,d,J1,J2,op,jbas) 
  ! grabs the matrix element you are looking for
  implicit none
  
  integer :: a,b,c,d,J1,J2,rank,T,P,q,qx,c1,c2,N,J1x,J2x
  integer :: int1,int2,i1,i2,j_min,x,k1,k2,ax,bx,cx,dx
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  ,phase
  logical :: fail_c,switch
  type(iso_operator) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c

  !make sure the matrix element exists first
 
  rank = op%rank

  if ( .not. (triangle ( J1,J2,rank ))) then 
     iso_op_elem = 0.d0 
     return
  end if 
   
  fail_c = .true. 
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  
  if ( .not. ((triangle(ja,jb,J1)) .and. (triangle (jc,jd,J2))) ) then 
     iso_op_elem = 0.d0
     return
  end if
     
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)

  P = mod(la + lb,2) 
     
  if ( mod(lc + ld,2) .ne. abs(P - ((-1)**(op%dpar/2+1)+1)/2) ) then
    iso_op_elem = 0.d0 
    return
  end if 
        
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)
     
  T = (ta + tb)/2
     
  if ((tc+td) .ne. 2*(T-op%dTz)) then     
    iso_op_elem = 0.d0
    return
  end if 

  q = iso_ladder_block_index(J1,J2,rank,T,P) 
  
  ! see subroutine "allocate_blocks" for mapping from qx to each 
  ! of the 6 storage arrays
      
  C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
  C2 = jbas%con(c)+jbas%con(d) + 1
    
  qx = C1*C2
  qx = qx + adjust_index(qx)   !Vpppp nature  

  pre = 1 
  
  N = op%Nsp
  ! get the indeces in the correct order
  if ( a > b )  then 
     x = bosonic_tp_index(b,a,N) 
     j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
     i1 = jbas%xmap_tensor(op%xindx,x)%Z( (J1-j_min)/2 + 2) 
     pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J1)/2 ) 
  else
     if (a == b) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(a,b,N)
     j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
     i1 = jbas%xmap_tensor(op%xindx,x)%Z( (J1-j_min)/2 + 2) 
  end if 

  if (c > d)  then     
     x = bosonic_tp_index(d,c,N) 
     j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
     i2 = jbas%xmap_tensor(op%xindx,x)%Z( (J2-j_min)/2 + 2) 
     pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J2)/2 ) 
  else 
     if (c == d) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(c,d,N) 
     j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
     i2 = jbas%xmap_tensor(op%xindx,x)%Z( (J2-j_min)/2 + 2)  
  end if
 
  ! grab the matrix element

   If (C1>C2) qx = qx + tensor_adjust(qx)       

   ! right now i1 and i2 still refer to where the pair is located
   iso_op_elem = op%tblck(q)%tgam(qx)%X(i1,i2) * pre
    
 end function iso_op_elem
!==================================================================  
!==================================================================
subroutine add_elem_to_ladder(V,a,b,c,d,J1,J2,op,jbas) 
  ! not safe. Don't use. 
  implicit none
  
  integer :: a,b,c,d,J1,J2,rank,T,P,q,qx,c1,c2,N,J1x,J2x
  integer :: int1,int2,i1,i2,j_min,x,k1,k2,ax,bx,cx,dx
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  ,phase
  logical :: fail_c,switch
  type(iso_ladder) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c,V


  !make sure the matrix element exists first

  rank = op%rank 
  
  fail_c = .true. 
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)

  P = mod(la + lb,2) 
        
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)
     
  T = (ta + tb)/2
  
  q = iso_ladder_block_index(J1,J2,rank,T,P) 
  
  ! see subroutine "allocate_blocks" for mapping from qx to each 
  ! of the 6 storage arrays
    
  pre = 1 
  
  N = op%Nsp
  ! get the indeces in the correct order

  if (a == b) pre = pre * sqrt( 2.d0 )
  x = bosonic_tp_index(a,b,N)
  j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
  i1 = jbas%xmap_tensor(op%xindx,x)%Z( (J1-j_min)/2 + 2) 
  
  if (c == d) pre = pre * sqrt( 2.d0 )
  x = bosonic_tp_index(c,d,N) 
  j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
  i2 = jbas%xmap_tensor(op%xindx,x)%Z( (J2-j_min)/2 + 2)  
  
 
  ! grab the matrix element

  ! right now i1 and i2 still refer to where the pair is located
  ! in the rank zero qn storage

  !$OMP ATOMIC
  op%tblck(q)%Xpphh(i1,i2) = op%tblck(q)%Xpphh(i1,i2) + V *pre 
    
end subroutine add_elem_to_ladder
!==================================================================  
!==================================================================
subroutine add_elem_to_iso_op(V,a,b,c,d,J1,J2,op,jbas) 
  ! not safe
  implicit none
  
  integer :: a,b,c,d,J1,J2,rank,T,P,q,qx,c1,c2,N,J1x,J2x
  integer :: int1,int2,i1,i2,j_min,x,k1,k2,ax,bx,cx,dx
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  ,phase
  logical :: fail_c,switch
  type(iso_operator) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c,V

  !make sure the matrix element exists first
 
  rank = op%rank
   
  fail_c = .true. 
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  
     
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)

  P = mod(la + lb,2) 
     
        
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)
     
  T = (ta + tb)/2
     

  q = iso_ladder_block_index(J1,J2,rank,T,P) 
  
  ! see subroutine "allocate_blocks" for mapping from qx to each 
  ! of the 6 storage arrays
      
  C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
  C2 = jbas%con(c)+jbas%con(d) + 1
    
  qx = C1*C2
  qx = qx + adjust_index(qx)   !Vpppp nature  

  pre = 1.d0 
  
  N = op%Nsp
  
  if (a == b) pre = pre * sqrt( 2.d0 )
  x = bosonic_tp_index(a,b,N)
  j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
  i1 = jbas%xmap_tensor(op%xindx,x)%Z( (J1-j_min)/2 + 2) 

  if (c == d) pre = pre * sqrt( 2.d0 )
  x = bosonic_tp_index(c,d,N) 
  j_min = jbas%xmap_tensor(op%xindx,x)%Z(1)  
  i2 = jbas%xmap_tensor(op%xindx,x)%Z( (J2-j_min)/2 + 2)  

 
  ! grab the matrix element

  If (C1>C2) qx = qx + tensor_adjust(qx)
  
   ! right now i1 and i2 still refer to where the pair is located
   !$OMP ATOMIC
   op%tblck(q)%tgam(qx)%X(i1,i2) = op%tblck(q)%tgam(qx)%X(i1,i2) + V * pre
    
 end subroutine add_elem_to_iso_op
!=================================================================     
!=================================================================
real(8) function f_iso_ladder_elem(a,b,op,jbas) 
  implicit none 
  
  integer :: a,b,x1,x2,c1,c2
  type(spd) :: jbas
  type(iso_ladder) :: op 
  
  ! are they holes or particles
  c1 = jbas%con(a)
  c2 = jbas%con(b) 
  
  select case(c1+c2) 
     case(0) 
        ! pp 
        f_iso_ladder_elem = 0.d0 
     case(1) 
        ! ph 
        if (c1 > c2) then 
           f_iso_ladder_elem = 0.d0 
        else 
           f_iso_ladder_elem = op%fph(a-hb4(a),b-pb4(b)) 
        end if
     case(2) 
        ! hh 
        f_iso_ladder_elem = 0.d0          
  end select

end function f_iso_ladder_elem
!=================================================================     
!=================================================================
real(8) function f_iso_op_elem(a,b,op,jbas) 
  implicit none 
  
  integer :: a,b,x1,x2,c1,c2
  type(spd) :: jbas
  type(iso_operator) :: op 
  
  f_iso_op_elem = op%fock(a,b)
  
end function f_iso_op_elem
!=================================================================     
!=================================================================
integer function iso_ladder_block_index(J1,J2,RANK,T,P) 
  ! input 2*J1 2*J2,rank,Tz,and Parity to get block index
  integer :: J1,J2,T,P,RANK,MORE,JX
  
  IF ( J1 < RANK ) THEN
     iso_ladder_block_index = 6*(J1*J1/4+(J1+J2-rank)/2) &
          + 2*(T+1) + P + 1      
  ELSE
     iso_ladder_block_index = 6*((J2+rank*J1)/2 -rank*rank/4) &
          + 2*(T+1) + P + 1 
  END IF
     
end function iso_ladder_block_index
!=================================================================     
!=================================================================
real(8) function count_dTz_configs(J1,Tz,PAR,jbas,ph,qn)
   implicit none
   
   type(spd) :: jbas
   integer :: J1,Tz,PAR,i,j,ji,jj,NX,r1,A,ix,jx 
   integer,optional,dimension(:,:) :: qn 
   logical,intent(in) :: ph 
   
   NX = jbas%total_orbits
   A = sum(jbas%con)
   r1 = 0
   do ix = 1, NX-A
      i = jbas%parts(ix) 
      do jx = 1,A 
         j = jbas%holes(jx)
         
         ji = jbas%jj(i) 
         jj = jbas%jj(j) 

         if ( abs(jbas%itzp(i) - jbas%itzp(j)) .ne. Tz ) cycle 

         if ( mod(jbas%ll(i) + jbas%ll(j),2) == PAR ) then
            if (triangle(ji,jj,J1)) then 

               r1 = r1 + 1 
               if (present(qn)) then
                  qn(r1,1) = i
                  qn(r1,2) = j
               end if                              
             
            end if
         end if
      end do
   end do

   count_dTz_configs = r1 
 end function count_dTz_configs
!===================================================================================
!===================================================================================
 real(8) function count_iso_op_configs(J1,Tz,PAR,jbas,ph,qn)
  implicit none
   
   type(spd) :: jbas
   integer :: J1,Tz,PAR,i,j,ji,jj,NX,r1
   integer,optional,dimension(:,:) :: qn 
   logical,intent(in) :: ph 
   
   NX = jbas%total_orbits
   r1 = 0
   do i = 1, NX
      do j = 1,NX 
         
         ji = jbas%jj(i) 
         jj = jbas%jj(j) 
         if (abs(jbas%itzp(i) - jbas%itzp(j)) .ne. Tz ) cycle 
         if ( mod(jbas%ll(i) + jbas%ll(j),2) == PAR ) then
            if (triangle(ji,jj,J1)) then 
               
               if (ph) then 
                  if ( (jbas%con(i) - jbas%con(j) .ne. 1)) cycle 
               end if
                              
               r1 = r1+1                       
               if (present(qn)) then
                  qn(r1,1) = i
                  qn(r1,2) = j
               end if
            end if
         end if
      end do
   end do

   count_iso_op_configs = r1 
 end function count_iso_op_configs
!=======================================================================
!=======================================================================
 subroutine fill_generalized_oppandya_matrix(J1,J2,MAT,qn1,qn2,OP,jbas,hp)
   ! CALCULATES THE CROSS GENERALIZED PANDYA MATRIX ELEMENTS OF
   ! OP FOR A GIVEN CHANNEL
   implicit none

   type(spd) :: jbas
   type(iso_operator) :: OP 
   integer,dimension(:,:) :: qn1,qn2
   real(8),dimension(:,:) :: MAT
   integer :: a,b,c,d,J1,J2,N1,N2,II,JJ,p,h
   logical :: hp

   N1  = size(MAT(:,1))
   N2  = size(MAT(1,:))   

   do JJ = 1,N2
      do II = 1,N1

         a=qn1(II,1)
         b=qn1(II,2)

         if ( jbas%con(a) == 1) then
            h = a
            p = b
         else
            h = b
            p = a
         end if
         
         c=qn2(JJ,1)
         d=qn2(JJ,2)

         if (hp) then 
            MAT(II,JJ) = Voppandya(h,p,c,d,J1,J2,Op,jbas)
         else
            MAT(II,JJ) = Voppandya(p,h,c,d,J1,J2,Op,jbas)
         end if
      end do
   end do
   
 end subroutine fill_generalized_oppandya_matrix
!=======================================================================
!=======================================================================
 subroutine fill_generalized_isopandya_matrix(J1,J2,MAT,qn1,qn2,OP,jbas)
   ! CALCULATES THE CROSS GENERALIZED PANDYA MATRIX ELEMENTS OF
   ! OP FOR A GIVEN CHANNEL
   implicit none

   type(spd) :: jbas
   type(iso_ladder) :: OP 
   integer,dimension(:,:) :: qn1,qn2
   real(8),dimension(:,:) :: MAT
   integer :: a,b,c,d,J1,J2,N1,N2,II,JJ

   N1  = size(MAT(:,1))
   N2  = size(MAT(1,:))   

   do JJ = 1,N2
      do II = 1,N1

         a=qn1(II,1)
         b=qn1(II,2)
         c=qn2(JJ,2)
         d=qn2(JJ,1)

         MAT(II,JJ) = Visopandya(a,b,c,d,J1,J2,Op,jbas)

      end do
   end do
   
 end subroutine fill_generalized_isopandya_matrix
!=======================================================================
!=======================================================================
 subroutine fill_cc_matrix(J1,MAT,qn1,OP,jbas)
   ! CALCULATES THE CROSS GENERALIZED PANDYA MATRIX ELEMENTS OF
   ! OP FOR A GIVEN CHANNEL
   implicit none

   type(spd) :: jbas
   type(sq_op) :: OP 
   integer,dimension(:,:) :: qn1
   real(8),dimension(:,:) :: MAT
   integer :: a,b,c,d,J1,N1,N2,II,JJ

   N1  = size(MAT(:,1))
   N2  = size(MAT(1,:))   

   do JJ = 1,N2
      do II = 1,N1

         a=qn1(II,1)
         b=qn1(II,2)
         c=qn1(JJ,2)
         d=qn1(JJ,1)

         MAT(II,JJ) = VCC(a,b,c,d,J1,Op,jbas)

      end do
   end do
   
 end subroutine fill_cc_matrix
!=========================================================
!=========================================================
 subroutine fill_rectangle_cc_matrix(J1,MAT,qn1,qn2,OP,jbas,ph)
   ! CALCULATES THE CROSS GENERALIZED PANDYA MATRIX ELEMENTS OF
   ! OP FOR A GIVEN CHANNEL
   ! ph is true if the summing channel is arranged "ph", else
   ! it should be "hp" 
   implicit none

   type(spd) :: jbas
   type(sq_op) :: OP 
   integer,dimension(:,:) :: qn1,qn2
   real(8),dimension(:,:) :: MAT
   integer :: a,b,c,d,J1,N1,N2,II,JJ,p,h
   logical :: ph 

   N1  = size(MAT(:,1))
   N2  = size(MAT(1,:))   

   do JJ = 1,N2
      do II = 1,N1

         a=qn1(II,1)
         b=qn1(II,2)
         c=qn2(JJ,2)
         d=qn2(JJ,1)

         if ( jbas%con(c) == 1) then
            h = c
            p = d
         else
            h = d
            p = c
         end if


         if (ph) then
            MAT(II,JJ) = VCC(a,b,p,h,J1,Op,jbas)
         else
            MAT(II,JJ) = VCC(a,b,h,p,J1,Op,jbas)
         end if
         
      end do
   end do
   
 end subroutine fill_rectangle_cc_matrix
!=========================================================
!=========================================================
 real(8) function Visopandya(a,d,c,b,J1,J2,Op,jbas)
  ! \overbar{V}^J_{ a \bar{d} c \bar{b} } 
  implicit none 
  
  integer :: a,b,c,d,J1,J2,rank,J4min,J4max
  integer :: ja,jb,jc,jd,J3,J4,j3min,j3max
  type(spd) :: jbas
  type(iso_ladder) :: Op
  real(8) :: sm ,coef9
  
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  rank = Op%rank 

  j3min = abs(ja-jb)
  j4min = abs(jc-jd) 
  j3max = ja+jb
  j4max = jc+jd  
  
  sm = 0.d0 
  
  do J3 = j3min,j3max,2
     do J4 = j4min,j4max,2 
     sm = sm - sqrt((J1+1.d0)*(J2+1.d0) &
          *(J3+1.d0)*(J4+1.d0)) * &
          ninej(OP%xindx,ja,jd,J1,jb,jc,J2,J3,J4,rank) * &
          iso_ladder_elem(a,b,c,d,J3,J4,Op,jbas) * &
          (-1)**((jb+jd+J2+J4)/2) 
     end do 
  end do 
  
  Visopandya = sm 
end function Visopandya
!=========================================================
!=========================================================
real(8) function Voppandya(a,d,c,b,J1,J2,Op,jbas)
  ! \overbar{V}^J_{ a \bar{d} c \bar{b} } 
  implicit none 
  
  integer :: a,b,c,d,J1,J2,rank,J4min,J4max
  integer :: ja,jb,jc,jd,J3,J4,j3min,j3max
  type(spd) :: jbas
  type(iso_operator) :: Op
  real(8) :: sm ,coef9
  
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  rank = Op%rank 

  j3min = abs(ja-jb)
  j4min = abs(jc-jd) 
  j3max = ja+jb
  j4max = jc+jd  
  
  sm = 0.d0 

  do J3 = j3min,j3max,2
     do J4 = j4min,j4max,2 
     sm = sm - sqrt((J1+1.d0)*(J2+1.d0) &
          *(J3+1.d0)*(J4+1.d0)) * &
          ninej(OP%xindx,ja,jd,J1,jb,jc,J2,J3,J4,rank) * &
          iso_op_elem(a,b,c,d,J3,J4,Op,jbas) * &
          (-1)**((jb+jd+J2+J4)/2) 
     end do 
  end do 
  
  Voppandya = sm 
end function Voppandya
!==================================================================
!==================================================================
real(8) function iso_ladder_mscheme(a,ma,b,mb,c,mc,d,md,MU,Op,jbas) 
  implicit none
  
  integer :: a,ma,b,mb,c,mc,d,md,J1min,J1max,mu,rank
  integer :: ja,jb,jc,jd,J1,M1,J2,M2,J2min,J2max
  type(spd) :: jbas
  type(iso_ladder) :: Op
  real(8) :: dcgi,sm
  
  M1 = ma+mb 
  M2 = mc+md 
  rank = Op%rank
  if (MU .ne. M1-M2) then
     iso_ladder_mscheme = 0.d0
     return
  end if
  
  sm = 0.d0 
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c) 
  jd = jbas%jj(d) 

  J1min = abs(ja-jb)
  J1max = ja+jb
  J2min = abs(jc-jd)
  J2max = jc+jd
  
  do J1 = J1min,J1max,2
     do J2 = J2min,J2max,2 
        
        sm = sm + iso_ladder_elem(a,b,c,d,J1,J2,op,jbas) *&
             dcgi(ja,ma,jb,mb,J1,M1)*&
             dcgi(jc,mc,jd,md,J2,M2)*&
             dcgi(J2,M2,RANK,MU,J1,M1)/&
             sqrt(J1 + 1.d0)
     end do
  end do
  
  iso_ladder_mscheme = sm 
end function iso_ladder_mscheme
!==================================================================
!==================================================================
real(8) function f_iso_ladder_mscheme(a,ma,b,mb,mu,Op,jbas)
  implicit none

  integer :: a,ma,b,mb,ja,jb,rank,mu
  type(spd) :: jbas
  type(iso_ladder) :: Op
  real(8) :: dcgi,sm 

  ja = jbas%jj(a)
  jb = jbas%jj(b)
  rank = op%rank
  if (mu+mb.ne.ma) then
     f_iso_ladder_mscheme= 0.d0 
  else
     f_iso_ladder_mscheme = dcgi(jb,mb,rank,mu,ja,ma)&
          /sqrt(ja+1.d0)*f_iso_ladder_elem(a,b,Op,jbas)
  end if
end function f_iso_ladder_mscheme
!===================================================================
!===================================================================
end module isospin_operators

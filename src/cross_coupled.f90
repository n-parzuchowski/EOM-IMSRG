module cross_coupled
  use basic_IMSRG
  
   type :: ph_mat
      type(int_vec),allocatable,dimension(:) :: qmap,nbmap
      type(real_mat),allocatable,dimension(:) :: CCX 
      integer :: nblocks,Nsp,herm
      integer,allocatable,dimension(:) :: Jval
   end type ph_mat 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  type,extends(ph_mat) :: cc_mat
     type(int_vec),allocatable,dimension(:) :: rmap
     integer,allocatable,dimension(:) :: nph,rlen
!     type(int_mat), allocatable,dimension(:) :: qn1,qn2
  end type cc_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type,extends(ph_mat) :: pandya_mat
     type(int_vec),allocatable,dimension(:) :: rmap
     type(real_mat),allocatable,dimension(:) :: CCR
     type(int_mat), allocatable,dimension(:) :: qn1,qn2
     integer,allocatable,dimension(:) :: Jval2,nb1,nb2
     integer :: rank,dpar
  end type pandya_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
   type,extends(ph_mat) :: ex_cc_mat ! excitation operator
      integer,allocatable,dimension(:) :: nph
   end type ex_cc_mat
! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   type,extends(ph_mat) ::  ex_pandya_mat ! excitation operator
      type(real_mat),allocatable,dimension(:) :: CCR
      integer,allocatable,dimension(:) :: Jval2
      integer :: rank,dpar
   end type ex_pandya_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  interface init_ph_mat 
     module procedure allocate_CCMAT,allocate_tensor_CCMAT,&
          allocate_ex_CCMAT,allocate_ex_tensor_CCMAT
  end interface

  interface duplicate_ph_mat
     module procedure dup_CC,dup_pandya,dup_ex_cc,dup_ex_pandya
  end interface
  
  interface init_ph_wkspc
     module procedure allocate_pandya_wkspc,allocate_CC_wkspc,&
          allocate_ex_pandya_wkspc,allocate_ex_CC_wkspc
  end interface
  
  interface fetch_rval
     module procedure pandya_rval,specific_rval ,ex_specific_rval,ex_pandya_rval
  end interface
  
  interface ph_rval
     module procedure ex_ph_rval,reg_ph_rval
  end interface
  
contains
!=======================================================  
!=======================================================
subroutine allocate_CCMAT(OP,CCME,jbas) 
  ! allocates a cross-coupled ME storage structure
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  
  type(spd) :: jbas
  type(sq_op) :: OP
  type(cc_mat) :: CCME
  integer :: JT,ji,jp,jj,jh,JC,q1,q2,g,li,lj,ti,tj
  integer :: a,b,p,h,i,j,r,Jmin,Jmax,NX,TZ,PAR,x,JTM
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,numJ,mes 
  real(8) :: sm,sm2

  mes = 0
  NX = OP%Nsp
  CCME%Nsp = NX
  CCME%herm = OP%herm
  JTM = jbas%Jtotal_max
  CCME%nblocks = (JTM + 1) * 2 * 2
  ! 2 dof for parity, 2 for Tz (i'm only worried about abs(Tz) ) 
  allocate(CCME%CCX(CCME%nblocks)) ! < h a |v | p b> 
  allocate(CCME%nph(CCME%nblocks)) ! number of ph pairs in block 
  allocate(CCME%rlen(CCME%nblocks)) ! number of total pairs in block
  allocate(CCME%Jval(CCME%nblocks)) ! J value for block
  allocate(CCME%rmap(NX*NX))  ! map of pair index to r indeces
  allocate(CCME%qmap(NX*NX))  ! map of pair index to q indeces
  allocate(CCME%nbmap(NX*NX)) 
  
  do i = 1, NX
     do j = 1,NX
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j) 
        
        numJ = (ji + jj - abs(ji-jj))/2 + 1
        
        x = CCindex(i,j,NX) 
        allocate(CCME%rmap(x)%Z(numJ)) 
        allocate(CCME%qmap(x)%Z(numJ))
        allocate(CCME%nbmap(x)%Z(numJ))
        CCME%qmap(x)%Z = 0
        CCME%nbmap(x)%Z = 0
     end do 
  end do 
  
  
  do q1 = 1, CCME%nblocks
     
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
     ! fastest changing quantity : JC
     ! slowest: PAR 
     
     nb = 0 
     r = 0 
     do i = 1, NX
        do j = 1,NX 
           
           ji = jbas%jj(i) 
           jj = jbas%jj(j) 
           
           if (.not. (triangle(ji,jj,JC))) cycle 
           if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
           if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 
           
           if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
              nb = nb + 1 
           end if 
           
           r = r+1
        end do 
     end do 
     
     allocate( CCME%CCX(q1)%X(r,nb) ) 
     mes = mes + r*nb
     
     nb = 0 
     r = 0 
     do i = 1, NX
        do j = 1,NX 
           
           ji = jbas%jj(i) 
           jj = jbas%jj(j) 
           if (.not. (triangle(ji,jj,JC))) cycle 
           if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
           if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 
          
           x = CCindex(i,j,NX) 
           
           g = 1
           do while (CCME%qmap(x)%Z(g) .ne. 0) 
              g = g + 1
           end do 
         
           if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then
              nb = nb + 1
              CCME%nbmap(x)%Z(g) = nb 
           end if
           
           r = r+1
           
           CCME%qmap(x)%Z(g) = q1
           CCME%rmap(x)%Z(g) = r

        end do 
     end do 
               
     CCME%nph(q1) = nb
     CCME%rlen(q1) = r
     CCME%Jval(q1) = JC 
  end do 
!  print*, "CC mat storage: ", mes * 8.d0 / 1024.d0/1024.d0, "MB" 
end subroutine allocate_CCMAT
!=======================================================  
!=======================================================
subroutine allocate_ex_CCMAT(OP,CCME,jbas) 
  ! allocates a cross-coupled ME storage structure
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  
  type(spd) :: jbas
  type(sq_op) :: OP
  type(ex_cc_mat) :: CCME
  integer :: JT,ji,jp,jj,jh,JC,q1,q2,g,li,lj,ti,tj,holes,parts
  integer :: a,b,p,h,i,j,r,Jmin,Jmax,NX,TZ,PAR,x,JTM
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,numJ
  real(8) :: sm,sm2
  
  NX = OP%Nsp
  CCME%Nsp = NX
  CCME%herm = OP%herm
  JTM = jbas%Jtotal_max
  CCME%nblocks = (JTM + 1) * 2 * 2
  holes = sum(jbas%con)
  parts = size(jbas%con)-holes
  ! 2 dof for parity, 2 for Tz (i'm only worried about abs(Tz) ) 
  allocate(CCME%CCX(CCME%nblocks)) ! < h a |v | p b> 
  allocate(CCME%nph(CCME%nblocks)) ! number of ph pairs in block 
  allocate(CCME%Jval(CCME%nblocks)) ! J value for block
  allocate(CCME%qmap(holes*parts))  ! map of pair index to q indeces
  allocate(CCME%nbmap(holes*parts)) 
  
  do i = 1, NX
     do j = 1,NX
        
        if ((jbas%con(j)-jbas%con(i)).ne.1) cycle
        ji = jbas%jj(i) 
        jj = jbas%jj(j) 
        
        numJ = (ji + jj - abs(ji-jj))/2 + 1
        
        x = PHindex(i,j,holes) 
        allocate(CCME%qmap(x)%Z(numJ))
        allocate(CCME%nbmap(x)%Z(numJ))
        CCME%qmap(x)%Z = 0
        CCME%nbmap(x)%Z = 0
     end do 
  end do 
  
  
  do q1 = 1, CCME%nblocks
     
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
     ! fastest changing quantity : JC
     ! slowest: PAR 
     
     nb = 0 
     r = 0 
     do i = 1, NX
        do j = 1,NX 
           
           ji = jbas%jj(i) 
           jj = jbas%jj(j) 
           
           if (.not. (triangle(ji,jj,JC))) cycle 
           if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
           if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 
           
           if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
              nb = nb + 1 
           end if 
           
           r = r+1
        end do 
     end do 
     
     allocate( CCME%CCX(q1)%X(nb,nb) ) 
    
     nb = 0 
     r = 0 
     do i = 1, NX
        do j = 1,NX 
           
           if ((jbas%con(j)-jbas%con(i)).ne.1)cycle
           ji = jbas%jj(i) 
           jj = jbas%jj(j) 
           if (.not. (triangle(ji,jj,JC))) cycle 
           if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
           if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 
          
           x = PHindex(i,j,holes) 
           
           g = 1
           do while (CCME%qmap(x)%Z(g) .ne. 0) 
              g = g + 1
           end do 
         
           if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then
              nb = nb + 1
              CCME%nbmap(x)%Z(g) = nb 
           end if
           
           r = r+1
           
           CCME%qmap(x)%Z(g) = q1

        end do 
     end do 
               
     CCME%nph(q1) = nb
     CCME%Jval(q1) = JC 
  end do 
end subroutine allocate_ex_CCMAT
!=======================================================  
!=======================================================
subroutine allocate_tensor_CCMAT(OP,CCME,jbas) 
  ! allocates a cross-coupled ME storage structure
  ! currently the only CCME of interest are phab terms    |--->--| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(bh)J>
  !                                                      |---<----|
  implicit none 
  
  
  type(spd) :: jbas
  type(sq_op) :: OP
  type(pandya_mat) :: CCME
  integer :: Jtot1,Jtot2,ji,jp,jj,jh,JC,q1,q2,g,li,lj,ti,tj,q
  integer :: a,b,p,h,i,j,Jmin,Jmax,NX,TZ,PAR,x,JTM,RANK,Jold
  integer :: int1,int2,IX,JX,i1,i2,nb1,nb2,r1,r2,nh,np,numJ
  real(8) :: sm,sm2
  
  NX = OP%Nsp
  RANK = OP%rank
  CCME%rank = OP%rank
  CCME%dpar = OP%dpar
  CCME%Nsp = NX
  CCME%herm = OP%herm  
  JTM = jbas%Jtotal_max
  Jold = 1
! quantum numbers of the last block 

  Tz = 1 
  Par = 1
  Jtot1 = jbas%Jtotal_max*2 
  Jtot2 = Jtot1+ RANK 
  CCME%nblocks  =  tensor_block_index(Jtot1,Jtot2,RANK,Tz,Par)/3*2 
  ! 2 dof for parity, 2 for Tz (i'm only worried about abs(Tz) ) 

  allocate(CCME%CCX(CCME%nblocks)) ! h(p)b(a)
  allocate(CCME%CCR(CCME%nblocks)) ! p(h)a(b)
  allocate(CCME%Jval(CCME%nblocks)) ! J value for block
  allocate(CCME%Jval2(CCME%nblocks)) ! J2 value for the block
  allocate(CCME%rmap(NX*NX))  ! map of pair index to r indeces
  allocate(CCME%qmap(NX*NX))  ! map of pair index to q indeces
  allocate(CCME%nbmap(NX*NX)) 
  allocate(CCME%qn1(NX*NX),CCME%qn2(NX*NX))

  do i = 1, NX
     do j = 1,NX
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j) 
        
        numJ = (ji + jj - abs(ji-jj))/2 + 1
        
        x = CCindex(i,j,NX) 
        allocate(CCME%rmap(x)%Z(numJ)) 
        allocate(CCME%qmap(x)%Z(numJ))
        allocate(CCME%nbmap(x)%Z(numJ))
        CCME%qmap(x)%Z = 0
        CCME%nbmap(x)%Z = 0
     end do 
  end do 
  
  q1 = 0
  do Jtot1 = 0,2*jbas%jtotal_max,2 
     do Jtot2 = max(abs(Jtot1 - rank),Jtot1),Jtot1+rank,2
        do Tz = 0, 1    
           do PAR = 0,1 
              
              q1 = q1 + 1
              CCME%Jval(q1) = Jtot1
              CCME%Jval2(q1) = Jtot2
              
              if (jtot2 > 2*jbas%jtotal_max) cycle 

              ! fastest changing quantity : JC
              ! slowest: PAR 

              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              do i = 1, NX
                 do j = 1,NX 

                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 

                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == PAR ) then
                       if (triangle(ji,jj,Jtot1)) then 

                          if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
                             nb1 = nb1 + 1 
                          end if
                          r1 = r1+1                       
                          
                       end if
                    end if
                    
                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == mod(PAR+op%dpar/2,2) ) then
                       if (triangle(ji,jj,Jtot2)) then 
                          if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
                             nb2 = nb2 + 1 
                          end if
                          
                          r2 = r2+1
                       end if
                    end if
                 end do
              end do

              allocate( CCME%CCX(q1)%X(r1,nb2) ) 
              allocate( CCME%CCR(q1)%X(nb1,r2) ) 

              allocate( CCME%qn1(q1)%Y(r1,2),CCME%qn2(q1)%Y(r2,2)) 
              
              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              ! I only need one of these arrays per J, so I use the same shape as before. 
              do i = 1, NX
                 do j = 1,NX 
                   
                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 
                                        
                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == PAR ) then
                       if (triangle(ji,jj,Jtot1)) then 

                          r1 = r1+1                       
                          
                          CCME%qn1(q1)%Y(r1,1) = i 
                          CCME%qn1(q1)%Y(r1,2) = j

                       end if
                    end if
                    

                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == mod(PAR+op%dpar/2,2) ) then
                       if (triangle(ji,jj,Jtot2)) then 
                          
                          r2 = r2+1
                          CCME%qn2(q1)%Y(r2,1) = i 
                          CCME%qn2(q1)%Y(r2,2) = j 
 
                       end if
                    end if
                   
                    
                 end do
              end do
              
              if (max(abs(Jtot1 - rank),Jtot1) .ne. Jtot2) cycle

              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              
              ! I only need one of these arrays per J, so I use the same shape as before. 
              do i = 1, NX
                 do j = 1,NX 
                   

                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 

                    if (.not. (triangle(ji,jj,Jtot1))) cycle 
                    if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    x = CCindex(i,j,NX) 
                    
                    g = 1
                    do while (CCME%qmap(x)%Z(g) .ne. 0) 
                       g = g + 1
                    end do

                    if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then
                       nb1 = nb1 + 1
                       CCME%nbmap(x)%Z(g) = nb1 
                    end if

                    r1 = r1+1
                    q = block_index(Jtot1,Tz,Par)          
                    CCME%qmap(x)%Z(g) = q
                    CCME%rmap(x)%Z(g) = r1
                                        
                 end do
              end do
              
           end do
        end do
     end do
  end do
  
end subroutine allocate_tensor_CCMAT        
!=======================================================  
!=======================================================
subroutine allocate_small_tensor_CCMAT(OP,CCME,jbas) 
  ! allocates a cross-coupled ME storage structure
  ! currently the only CCME of interest are phab terms    |--->--| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(bh)J>
  !                                                      |---<----|
  implicit none 
  
  
  type(spd) :: jbas
  type(sq_op) :: OP
  type(pandya_mat) :: CCME
  integer :: Jtot1,Jtot2,ji,jp,jj,jh,JC,q1,q2,g,li,lj,ti,tj,q
  integer :: a,b,p,h,i,j,Jmin,Jmax,NX,TZ,PAR,x,JTM,RANK,Jold
  integer :: int1,int2,IX,JX,i1,i2,nb1,nb2,r1,r2,nh,np,numJ
  real(8) :: sm,sm2
  
  NX = OP%Nsp
  RANK = OP%rank
  CCME%rank = OP%rank
  CCME%dpar = OP%dpar
  CCME%Nsp = NX
  CCME%herm = OP%herm  
  JTM = jbas%Jtotal_max
  Jold = 1
! quantum numbers of the last block 

  Tz = 1 
  Par = 1
  Jtot1 = jbas%Jtotal_max*2 
  Jtot2 = Jtot1+ RANK 
  CCME%nblocks  =  tensor_block_index(Jtot1,Jtot2,RANK,Tz,Par)/3*2 
  ! 2 dof for parity, 2 for Tz (i'm only worried about abs(Tz) ) 

  allocate(CCME%CCX(CCME%nblocks)) ! h(p)b(a)
  allocate(CCME%CCR(CCME%nblocks)) ! p(h)a(b)
  allocate(CCME%nb1(CCME%nblocks),CCME%nb2(CCME%nblocks))
  allocate(CCME%Jval(CCME%nblocks)) ! J value for block
  allocate(CCME%Jval2(CCME%nblocks)) ! J2 value for the block
  allocate(CCME%rmap(NX*NX))  ! map of pair index to r indeces
  allocate(CCME%qmap(NX*NX))  ! map of pair index to q indeces
  allocate(CCME%nbmap(NX*NX)) 
  allocate(CCME%qn1(NX*NX),CCME%qn2(NX*NX))

  do i = 1, NX
     do j = 1,NX
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j) 
        
        numJ = (ji + jj - abs(ji-jj))/2 + 1
        
        x = CCindex(i,j,NX) 
        allocate(CCME%rmap(x)%Z(numJ)) 
        allocate(CCME%qmap(x)%Z(numJ))
        allocate(CCME%nbmap(x)%Z(numJ))
        CCME%qmap(x)%Z = 0
        CCME%nbmap(x)%Z = 0
     end do 
  end do 
  
  q1 = 0
  do Jtot1 = 0,2*jbas%jtotal_max,2 
     do Jtot2 = max(abs(Jtot1 - rank),Jtot1),Jtot1+rank,2
        do Tz = 0, 1    
           do PAR = 0,1 
              
              q1 = q1 + 1
              CCME%Jval(q1) = Jtot1
              CCME%Jval2(q1) = Jtot2
              
              if (jtot2 > 2*jbas%jtotal_max) cycle 

              ! fastest changing quantity : JC
              ! slowest: PAR 

              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              do i = 1, NX
                 do j = 1,NX 

                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 

                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == PAR ) then
                       if (triangle(ji,jj,Jtot1)) then 

                          if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
                             nb1 = nb1 + 1 
                          end if
                          r1 = r1+1                       
                          
                       end if
                    end if
                    
                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == mod(PAR+op%dpar/2,2) ) then
                       if (triangle(ji,jj,Jtot2)) then 
                          if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
                             nb2 = nb2 + 1 
                          end if
                          
                          r2 = r2+1
                       end if
                    end if
                 end do
              end do

!              allocate( CCME%CCX(q1)%X(r1,nb2) ) 
 !             allocate( CCME%CCR(q1)%X(nb1,r2) ) 
              CCME%nb1(q1) = nb1
              CCME%nb2(q1) = nb2 
              allocate( CCME%qn1(q1)%Y(r1,2),CCME%qn2(q1)%Y(r2,2)) 
              
              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              ! I only need one of these arrays per J, so I use the same shape as before. 
              do i = 1, NX
                 do j = 1,NX 
                   
                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 
                                        
                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == PAR ) then
                       if (triangle(ji,jj,Jtot1)) then 

                          r1 = r1+1                       
                          
                          CCME%qn1(q1)%Y(r1,1) = i 
                          CCME%qn1(q1)%Y(r1,2) = j

                       end if
                    end if
                    

                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == mod(PAR+op%dpar/2,2) ) then
                       if (triangle(ji,jj,Jtot2)) then 
                          
                          r2 = r2+1
                          CCME%qn2(q1)%Y(r2,1) = i 
                          CCME%qn2(q1)%Y(r2,2) = j 
 
                       end if
                    end if
                   
                    
                 end do
              end do
              
              if (max(abs(Jtot1 - rank),Jtot1) .ne. Jtot2) cycle

              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              
              ! I only need one of these arrays per J, so I use the same shape as before. 
              do i = 1, NX
                 do j = 1,NX 
                   

                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 

                    if (.not. (triangle(ji,jj,Jtot1))) cycle 
                    if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    x = CCindex(i,j,NX) 
                    
                    g = 1
                    do while (CCME%qmap(x)%Z(g) .ne. 0) 
                       g = g + 1
                    end do

                    if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then
                       nb1 = nb1 + 1
                       CCME%nbmap(x)%Z(g) = nb1 
                    end if

                    r1 = r1+1
                    q = block_index(Jtot1,Tz,Par)          
                    CCME%qmap(x)%Z(g) = q
                    CCME%rmap(x)%Z(g) = r1
                                        
                 end do
              end do
              
           end do
        end do
     end do
  end do
  
end subroutine allocate_small_tensor_CCMAT
!=======================================================  
!=======================================================
subroutine allocate_ex_tensor_CCMAT(OP,CCME,jbas) 
  ! allocates a cross-coupled ME storage structure
  ! currently the only CCME of interest are phab terms    |--->--| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(bh)J>
  !                                                      |---<----|
  implicit none 
  
  
  type(spd) :: jbas
  type(sq_op) :: OP
  type(ex_pandya_mat) :: CCME
  integer :: Jtot1,Jtot2,ji,jp,jj,jh,JC,q1,q2,g,li,lj,ti,tj,q
  integer :: a,b,p,h,i,j,Jmin,Jmax,NX,TZ,PAR,x,JTM,RANK,Jold
  integer :: int1,int2,IX,JX,i1,i2,nb1,nb2,r1,r2,nh,np,numJ,holes,parts
  real(8) :: sm,sm2
  
  NX = OP%Nsp
  RANK = OP%rank
  CCME%rank = OP%rank
  CCME%dpar = OP%dpar
  CCME%Nsp = NX
  CCME%herm = OP%herm  
  JTM = jbas%Jtotal_max
  Jold = 1
! quantum numbers of the last block 
  
  holes = sum(jbas%con) 
  parts = size(jbas%con)-holes
  
  Tz = 1 
  Par = 1
  Jtot1 = jbas%Jtotal_max*2 
  Jtot2 = Jtot1+ RANK 
  CCME%nblocks  =  tensor_block_index(Jtot1,Jtot2,RANK,Tz,Par)/3*2 
  ! 2 dof for parity, 2 for Tz (i'm only worried about abs(Tz) ) 

  allocate(CCME%CCX(CCME%nblocks)) ! h(p)b(a)
  allocate(CCME%CCR(CCME%nblocks)) ! p(h)a(b)
  allocate(CCME%Jval(CCME%nblocks)) ! J value for block
  allocate(CCME%Jval2(CCME%nblocks)) ! J2 value for the block
  allocate(CCME%qmap(parts*holes))  ! map of pair index to q indeces
  allocate(CCME%nbmap(parts*holes)) 
  
  do i = 1, NX
     do j = 1,NX
        
        if ((jbas%con(j)-jbas%con(i)).ne.1) cycle
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j) 
        numJ = (ji + jj - abs(ji-jj))/2 + 1        
        x = PHindex(i,j,holes) 
        allocate(CCME%qmap(x)%Z(numJ))
        allocate(CCME%nbmap(x)%Z(numJ))
        CCME%qmap(x)%Z = 0
        CCME%nbmap(x)%Z = 0
     end do 
  end do 

  q1 = 0
  do Jtot1 = 0,2*jbas%jtotal_max,2 
     do Jtot2 = max(abs(Jtot1 - rank),Jtot1),Jtot1+rank,2
        do Tz = 0, 1    
           do PAR = 0,1 
              
              q1 = q1 + 1
              CCME%Jval(q1) = Jtot1
              CCME%Jval2(q1) = Jtot2
              
              if (jtot2 > 2*jbas%jtotal_max) then 
                 allocate( CCME%CCX(q1)%X(0,0) ) 
                 allocate( CCME%CCR(q1)%X(0,0) ) 
                 cycle 
              end if 
              ! fastest changing quantity : JC
              ! slowest: PAR 

              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              do i = 1, NX
                 do j = 1,NX 

                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 


                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == PAR ) then
                       if (triangle(ji,jj,Jtot1)) then 
                          if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
                             nb1 = nb1 + 1 
                          end if
                          r1 = r1 + 1
                       end if
                    end if
                    
                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == mod(PAR+op%dpar/2,2) ) then
                       if (triangle(ji,jj,Jtot2)) then 
                          if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
                             nb2 = nb2 + 1 
                          end if                         
                          r2 = r2 + 1
                       end if
                    end if
                 end do
              end do

              allocate( CCME%CCX(q1)%X(nb1,nb2) ) 
              allocate( CCME%CCR(q1)%X(nb1,nb2) ) 

              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              if (max(abs(Jtot1 - rank),Jtot1) .ne. Jtot2) cycle
              ! I only need one of these arrays per J, so I use the same shape as before. 
              do i = 1, NX
                 do j = 1,NX 
                   
                    if ((jbas%con(j)-jbas%con(i)).ne.1) cycle
                    
                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 

                    if (.not. (triangle(ji,jj,Jtot1))) cycle 
                    if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    x = PHindex(i,j,holes) 

                    g = 1
                    do while (CCME%qmap(x)%Z(g) .ne. 0) 
                       g = g + 1
                    end do

                    if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then
                       nb1 = nb1 + 1
                       CCME%nbmap(x)%Z(g) = nb1 
                    end if

                    r1 = r1+1
                    q = block_index(Jtot1,Tz,Par)          
                    CCME%qmap(x)%Z(g) = q
                    
                 end do
              end do
              
           end do
        end do
     end do
  end do
  
end subroutine allocate_ex_tensor_CCMAT
!=======================================================  
!=======================================================          
subroutine dup_CC(C1,CCME) 
  ! makes a copy of C1 onto CCME
  implicit none 
  
  type(cc_mat) :: C1,CCME
  integer :: JT,ji,jp,jj,jh,JC,q1,q2
  integer :: a,b,p,h,i,j,r,Jmin,Jmax,NX,TZ,PAR
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,numJ
  real(8) :: sm,sm2
  
  NX = C1%Nsp
  CCME%Nsp = NX
  CCME%nblocks = C1%nblocks
  CCME%herm = C1%herm
 
  allocate(CCME%CCX(C1%nblocks))
  allocate(CCME%nph(C1%nblocks))
  allocate(CCME%rlen(C1%nblocks)) 
  allocate(CCME%Jval(C1%nblocks)) 
  allocate(CCME%rmap(NX*NX))  ! map of pair index to r indeces
  allocate(CCME%qmap(NX*NX))  ! map of pair index to q indeces
  allocate(CCME%nbmap(NX*NX)) ! map of ph pair index to nb indeces
  
  do i = 1, NX
     do j = 1,NX
        
        r = CCindex(i,j,NX)
        
        numJ = size(C1%rmap(r)%Z) 
        
        allocate(CCME%rmap(r)%Z(numJ)) 
        allocate(CCME%qmap(r)%Z(numJ))
        allocate(CCME%nbmap(r)%Z(numJ)) 
        CCME%rmap(r)%Z = c1%rmap(r)%Z
        CCME%qmap(r)%Z = c1%qmap(r)%Z
        CCME%nbmap(r)%Z = c1%nbmap(r)%Z
       
     end do 
  end do 
   
  CCME%nph = C1%nph
  CCME%rlen = C1%rlen
  CCME%Jval = C1%Jval
  
  do q1 = 1, CCME%nblocks
     
     r = CCME%rlen(q1)
     nb = CCME%nph(q1) 
     allocate( CCME%CCX(q1)%X(r,nb) ) 
     
  end do 
end subroutine dup_CC
!=======================================================  
!=======================================================
subroutine dup_ex_CC(C1,CCME) 
  ! makes a copy of C1 onto CCME
  implicit none 
  
  type(ex_cc_mat) :: C1,CCME
  integer :: JT,ji,jp,jj,jh,JC,q1,q2,hp
  integer :: a,b,p,h,i,j,r,Jmin,Jmax,NX,TZ,PAR
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,numJ
  real(8) :: sm,sm2
  
  NX = C1%Nsp
  CCME%Nsp = NX
  CCME%nblocks = C1%nblocks
  CCME%herm = C1%herm
 
  hp = size(C1%qmap)
  allocate(CCME%CCX(C1%nblocks))
  allocate(CCME%nph(C1%nblocks))
  allocate(CCME%Jval(C1%nblocks)) 
  allocate(CCME%qmap(hp))  ! map of pair index to q indeces
  allocate(CCME%nbmap(hp)) ! map of ph pair index to nb indeces
  
  do r = 1, hp 

        numJ = size(C1%qmap(r)%Z) 
        allocate(CCME%qmap(r)%Z(numJ))
        allocate(CCME%nbmap(r)%Z(numJ)) 
        CCME%qmap(r)%Z = c1%qmap(r)%Z
        CCME%nbmap(r)%Z = c1%nbmap(r)%Z
           
  end do 
   
  CCME%nph = C1%nph
  CCME%Jval = C1%Jval
  
  do q1 = 1, CCME%nblocks
     
     nb = CCME%nph(q1) 
     allocate( CCME%CCX(q1)%X(nb,nb) ) 
     
  end do 
end subroutine dup_ex_CC
!=======================================================  
!=======================================================
subroutine dup_pandya(C1,CCME) 
  ! makes a copy of C1 onto CCME
  implicit none 
  
  type(pandya_mat) :: C1,CCME
  integer :: JT,ji,jp,jj,jh,JC,q1,q2
  integer :: a,b,p,h,i,j,r,Jmin,Jmax,NX,TZ,PAR
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,numJ
  real(8) :: sm,sm2
    
  NX = C1%Nsp
  CCME%Nsp = NX
  CCME%rank = C1%rank
  CCME%dpar = C1%dpar
  CCME%nblocks = C1%nblocks
  CCME%herm = C1%herm
 
  allocate(CCME%CCX(C1%nblocks))
  allocate(CCME%CCR(C1%nblocks))
  allocate(CCME%Jval(C1%nblocks)) 
  allocate(CCME%rmap(NX*NX))  ! map of pair index to r indeces
  allocate(CCME%qmap(NX*NX))  ! map of pair index to q indeces
  allocate(CCME%nbmap(NX*NX)) ! map of ph pair index to nb indeces
  
  do i = 1, NX
     do j = 1,NX
        
        r = CCindex(i,j,NX)
        
        numJ = size(C1%rmap(r)%Z) 
        
        allocate(CCME%rmap(r)%Z(numJ)) 
        allocate(CCME%qmap(r)%Z(numJ))
        allocate(CCME%nbmap(r)%Z(numJ)) 
        CCME%rmap(r)%Z = c1%rmap(r)%Z
        CCME%qmap(r)%Z = c1%qmap(r)%Z
        CCME%nbmap(r)%Z = c1%nbmap(r)%Z
       
     end do 
  end do 
   
  CCME%Jval = C1%Jval
  allocate(CCME%Jval2(C1%nblocks)) 
  CCME%Jval2 = C1%Jval2
  
  
  do q1 = 1, CCME%nblocks
     
     r = size(C1%CCX(q1)%X(:,1))    
     nb = size(C1%CCX(q1)%X(1,:))
     allocate( CCME%CCX(q1)%X(r,nb) ) 

     nb = size(C1%CCR(q1)%X(:,1))    
     r = size(C1%CCR(q1)%X(1,:))
     allocate( CCME%CCR(q1)%X(nb,r) ) 
     
  end do 

end subroutine dup_pandya
!=======================================================  
!=======================================================
subroutine dup_ex_pandya(C1,CCME) 
  ! makes a copy of C1 onto CCME
  implicit none 
  
  type(ex_pandya_mat) :: C1,CCME
  integer :: JT,ji,jp,jj,jh,JC,q1,q2,hp
  integer :: a,b,p,h,i,j,r,Jmin,Jmax,NX,TZ,PAR
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,numJ
  real(8) :: sm,sm2
    
  NX = C1%Nsp
  CCME%Nsp = NX
  CCME%rank = C1%rank
  CCME%dpar = C1%dpar
  CCME%nblocks = C1%nblocks
  CCME%herm = C1%herm
  hp = size(C1%qmap)
  
  allocate(CCME%CCX(C1%nblocks))
  allocate(CCME%CCR(C1%nblocks))
  allocate(CCME%Jval(C1%nblocks)) 
  allocate(CCME%qmap(hp))  ! map of pair index to q indeces
  allocate(CCME%nbmap(hp)) ! map of ph pair index to nb indeces
  
  do r = 1, hp
        
     numJ = size(C1%qmap(r)%Z) 
     allocate(CCME%qmap(r)%Z(numJ))
     allocate(CCME%nbmap(r)%Z(numJ)) 
     CCME%qmap(r)%Z = c1%qmap(r)%Z
     CCME%nbmap(r)%Z = c1%nbmap(r)%Z
       
  end do 
   
  CCME%Jval = C1%Jval
  allocate(CCME%Jval2(C1%nblocks)) 
  CCME%Jval2 = C1%Jval2
  
  do q1 = 1, CCME%nblocks
     
     r = size(C1%CCX(q1)%X(:,1))    
     nb = size(C1%CCX(q1)%X(1,:))
     allocate( CCME%CCX(q1)%X(r,nb) ) 

     nb = size(C1%CCR(q1)%X(:,1))    
     r = size(C1%CCR(q1)%X(1,:))
     allocate( CCME%CCR(q1)%X(nb,r) ) 
     
  end do 

end subroutine dup_ex_pandya
!=======================================================
!=======================================================
subroutine allocate_CC_wkspc(CCOP,WCC)
  implicit none 
  
  type(cc_mat) :: CCOP,WCC 
  integer :: q,r,mes

  mes= 0 
  allocate(WCC%CCX(CCOP%nblocks))

  do q = 1,CCOP%nblocks
     
     r = CCOP%rlen(q)      
     allocate(WCC%CCX(q)%X(r,r)) 
     mes = mes +r*r
     WCC%CCX(q)%X = 0.d0
     
  end do

!  print*, "CC workspace:", mes*8.d0/1024.d0/1024.d0, "MB"
  

end subroutine allocate_CC_wkspc
!===========================================================
!===========================================================  
subroutine allocate_pandya_wkspc(CCOP,WCC)
  implicit none 
  
  type(pandya_mat) :: CCOP,WCC 
  integer :: q,r1,r2
  
  allocate(WCC%CCX(CCOP%nblocks))
  allocate(WCC%CCR(CCOP%nblocks))
  
  do q = 1,CCOP%nblocks
     
     r1 = size(CCOP%CCX(q)%X(:,1)) 
     r2 = size(CCOP%CCR(q)%X(1,:))
     
     allocate(WCC%CCX(q)%X(r1,r2)) 
     allocate(WCC%CCR(q)%X(r1,r2)) 
     WCC%CCX(q)%X = 0.d0
     WCC%CCR(q)%X = 0.d0
     
  end do

end subroutine allocate_pandya_wkspc
!=======================================================
!=======================================================
subroutine allocate_ex_CC_wkspc(CCOP,WCC)
  implicit none 
  
  type(ex_cc_mat) :: CCOP,WCC 
  integer :: q,r
  
  allocate(WCC%CCX(CCOP%nblocks))

  do q = 1,CCOP%nblocks
     
     r = CCOP%nph(q)      
     allocate(WCC%CCX(q)%X(r,r)) 
     WCC%CCX(q)%X = 0.d0
     
  end do
  

end subroutine allocate_ex_CC_wkspc
!===========================================================
!===========================================================  
subroutine allocate_ex_pandya_wkspc(CCOP,WCC)
  implicit none 
  
  type(ex_pandya_mat) :: CCOP,WCC 
  integer :: q,r1,r2
  
  allocate(WCC%CCX(CCOP%nblocks))
  allocate(WCC%CCR(CCOP%nblocks))
  
  do q = 1,CCOP%nblocks
     
     r1 = size(CCOP%CCX(q)%X(:,1)) 
     r2 = size(CCOP%CCR(q)%X(1,:))
     
     allocate(WCC%CCX(q)%X(r1,r2)) 
     allocate(WCC%CCR(q)%X(r1,r2)) 
     WCC%CCX(q)%X = 0.d0
     WCC%CCR(q)%X = 0.d0
     
  end do

end subroutine allocate_ex_pandya_wkspc
!===========================================================
!===========================================================     
integer function CCindex(a,b,N)
  implicit none 
  
  integer :: a,b,N 
  
  CCindex = N*(a-1) + b
end function CCindex
!===========================================================
!===========================================================     
integer function PHindex(a,b,holes)
  ! here N are the holes
  implicit none 
  
  integer :: a,b,holes,p,h  
  
  p = pb4(a)+1
  h = hb4(b)+1
  
  PHindex = holes*(p-1) + h

end function PHindex
!===========================================================
!===========================================================     
subroutine calculate_cross_coupled(HS,CCME,jbas) 
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: HS
  type(cc_mat) :: CCME
  integer :: JT,ja,jp,jb,jh,JC,q1,q2,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin,Jmax,Rindx,Gindx,g,ta,tb,Atot,hg,pg
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx,x,JTM
  real(8) :: sm,sm2,pre

  Atot = HS%belowEF
  Ntot = HS%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  CCME%herm = HS%Herm
  
!$omp parallel do default(firstprivate),shared(CCME,HS,jbas) 
  do q1 = 1, CCME%nblocks
      
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
     
     CCME%CCX(q1)%X = 0.d0
         
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. triangle(jp,jh,JC) )  cycle
           if ( mod(lp + lh,2) .ne. PAR ) cycle
           if (abs(tp - th)/2 .ne. Tz ) cycle 
           
           x = CCindex(p,h,HS%Nsp)
           gnb = 1
           do while (CCME%qmap(x)%Z(gnb) .ne. q1 )
              gnb = gnb + 1
           end do
              
           NBindx = CCME%nbmap(x)%Z(gnb) 

           ! for the ph  channel 2body derivative 
        
           do a = 1, HS%nsp
              do b = 1, HS%nsp
  
                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 
                 if (.not. triangle(ja,jb,JC) )  cycle
                 if ( mod(la + lb,2) .ne. PAR ) cycle
                 if (abs(ta - tb)/2 .ne. Tz ) cycle 
                 
                 x = CCindex(a,b,HS%Nsp) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Rindx = CCME%rmap(x)%Z(g)
                    
                 sm = 0.d0 
               
!               
                 if ( (mod(la + lh,2) == mod(lb + lp,2)) .and. &
                      ( (ta + th) == (tb + tp) ) ) then  
               
                    ! hapb 
                    Jmin = max(abs(jp - jb),abs(ja - jh)) 
                    Jmax = min(jp+jb,ja+jh) 
                    
                    sm = 0.d0 
                    do JT = Jmin,Jmax,2
                       sm = sm + (-1)**(JT/2) * (JT + 1) * &
                            sixj(jp,jh,JC,ja,jb,JT)  * &
                            v_elem(h,a,p,b,JT,HS,jbas) 
                    end do
                 
                    ! store  < h a | v | p b>   ( V )_a(b)h(p)
                    CCME%CCX(q1)%X(Rindx,NBindx) = sm * &
                         (-1) **( (jh + jb + JC) / 2) * pre * sqrt(JC + 1.d0)
                    ! scaled by sqrt(JC + 1) for convience in ph derivative
                                     
                 end if
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine calculate_cross_coupled
!=======================================================  
!=======================================================          
subroutine calculate_cross_coupled_pphh(HS,CCME,jbas) 
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: HS
  type(ex_cc_mat) :: CCME
  integer :: JT,ja,jp,jb,jh,JC,q1,q2,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin,Jmax,Rindx,Gindx,g,ta,tb,Atot,hg,pg
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx,x,JTM,ax,bx
  real(8) :: sm,sm2,pre,horse

  Atot = HS%belowEF
  Ntot = HS%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  CCME%herm = HS%Herm
  
!$omp parallel do default(firstprivate),shared(CCME,HS,jbas) 
  do q1 = 1, CCME%nblocks
      
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
     
     CCME%CCX(q1)%X = 0.d0
         
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. triangle(jp,jh,JC) )  cycle
           if ( mod(lp + lh,2) .ne. PAR ) cycle
           if (abs(tp - th)/2 .ne. Tz ) cycle 
           
           x = PHindex(p,h,Atot)
           gnb = 1
           do while (CCME%qmap(x)%Z(gnb) .ne. q1 )
              gnb = gnb + 1
           end do
              
           NBindx = CCME%nbmap(x)%Z(gnb) 
        
           do ax = 1, HS%belowEF
              a = jbas%holes(ax) 
              do bx = 1, HS%nsp-HS%belowEF
                 b = jbas%parts(bx) 
                 
                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 
                 if (.not. triangle(ja,jb,JC) )  cycle
                 if ( mod(la + lb,2) .ne. PAR ) cycle
                 if (abs(ta - tb)/2 .ne. Tz ) cycle 
                 
                 x = PHindex(b,a,Atot) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Gindx = CCME%nbmap(x)%Z(g)
                 
                 sm = 0.d0 
               
!                 horse = 0.d0 
                 if ( (mod(lb + lh,2) == mod(la + lp,2)) .and. &
                      ( (tb + th) == (ta + tp) ) ) then  
               
                    ! hapb 
                    Jmin = max(abs(jp - ja),abs(jb - jh)) 
                    Jmax = min(jp+ja,jb+jh) 
                    
                    sm = 0.d0 
                    do JT = Jmin,Jmax,2
                       sm = sm + (-1)**(JT/2) * (JT + 1) * &
                            sixj(jp,jh,JC,jb,ja,JT)  * &
                            v_elem(h,b,p,a,JT,HS,jbas) 
                    end do
                 
                    ! store  < h a | v | p b>    Pandya ( V )_h(p)b(a)
                    CCME%CCX(q1)%X(Gindx,NBindx) = sm * &
                         (-1) **( (jh + ja + JC) / 2) * pre * sqrt(JC + 1.d0)
                    ! scaled by sqrt(JC + 1) for convience in ph derivative

                 end if
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine calculate_cross_coupled_pphh
!=======================================================  
!=======================================================          
subroutine calculate_single_pandya(OP,CCME,jbas,q1,term) 
  ! currently the only CCME of interest are phab terms    |---<--|  J1 
  ! coupling in the 3-1 channel                        <(pa)|V|(hb)> rank
  !                                                      |---<--| J2 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: OP 
  type(pandya_mat) :: CCME
  integer :: Jtot1,Jtot2,ja,jp,jb,jh,JC,q2,q,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin1,Jmax1,Rindx,Gindx,g,ta,tb,Atot,hg,pg,J3,J4,NBindx2,qONE,qTWO
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx1,x,JTM,rank,Jmin2,Jmax2
  integer,intent(in) :: q1,term
  real(8) :: sm,sm2,pre,horse
  logical :: parflip

  Atot = OP%belowEF
  Ntot = OP%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  rank = OP%rank

  parflip = .false. 
  if ( mod(op%dpar/2,2) == 1) parflip = .true. 

  CCME%herm = OP%Herm
!!$omp parallel do default(firstprivate),shared(CCME,OP,jbas) 
      
  Jtot1 = CCME%Jval(q1)
  Jtot2 = CCME%Jval2(q1)
!  if (Jtot2 > 2*JTM) cycle

  PAR = mod(q1-1,2)
  Tz = mod((q1-1)/2,2) 

  if ( term == 1) then 
     CCME%CCX(q1)%X = 0.d0
  else
     CCME%CCR(q1)%X = 0.d0
  end if 

  qONE = block_index(Jtot1,Tz,Par)
  qTWO = block_index(Jtot2,Tz,mod(Par+op%dpar/2,2))

  ! ab = ph 
  do hg = 1, Atot
     do pg = 1, Ntot - Atot 

        h = jbas%holes(hg) 
        p = jbas%parts(pg) 

        jp = jbas%jj(p) 
        jh = jbas%jj(h)
        lp = jbas%ll(p) 
        lh = jbas%ll(h)
        tp = jbas%itzp(p) 
        th = jbas%itzp(h)

        if (.not. (parflip)) then 
           if ( mod(lp + lh,2) .ne. PAR ) cycle
        end if

        if (abs(tp - th)/2 .ne. Tz ) cycle 

        NBindx1 = 0
        NBindx2 = 0 
        if ( triangle(jp,jh,Jtot1) )  then 


           if ( mod(lp + lh,2) ==  PAR ) then   
              x = CCindex(p,h,OP%Nsp)
              gnb = 1

              do while (CCME%qmap(x)%Z(gnb) .ne. qONE )
                 gnb = gnb + 1 
              end do

              NBindx1 = CCME%nbmap(x)%Z(gnb) 

           end if

           if  ( triangle(jp,jh,Jtot2) )  then 

              if ( mod(lp+lh+op%dpar/2,2) == PAR) then 

                 x = CCindex(p,h,OP%Nsp)
                 gnb = 1

                 do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                    gnb = gnb + 1 
                 end do

                 NBindx2 = CCME%nbmap(x)%Z(gnb)
              end if
           end if

        else if  ( triangle(jp,jh,Jtot2) )  then 

           if ( mod(lp+lh+op%dpar/2,2) == PAR) then 

              x = CCindex(p,h,OP%Nsp)
              gnb = 1

              do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                 gnb = gnb + 1 
              end do

              NBindx2 = CCME%nbmap(x)%Z(gnb)
           end if
        else 
           cycle
        end if


        ! for the ph  channel 2body derivative 

        do a = 1, OP%nsp
           do b = 1, OP%nsp

              ja = jbas%jj(a) 
              jb = jbas%jj(b)
              la = jbas%ll(a) 
              lb = jbas%ll(b)
              ta = jbas%itzp(a) 
              tb = jbas%itzp(b)


              if (.not. (parflip)) then 
                 if ( mod(la + lb,2) .ne. PAR ) cycle
              end if

              if (abs(ta - tb)/2 .ne. Tz ) cycle 
              
              IF (term ==1 ) then 
                 if ( (triangle(ja,jb,Jtot1)) .and. (NBindx2 .ne. 0) ) then 
                    
                    if ( mod(la+lb,2) == PAR ) then 
                       x = CCindex(a,b,OP%Nsp) 

                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qONE )
                          g = g + 1
                       end do

                       Rindx = CCME%rmap(x)%Z(g)

                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  

                          ! hapb 
                          Jmin1 = abs(ja - jh) 
                          Jmax1 = ja+jh 
                          Jmin2 = abs(jp - jb)
                          Jmax2 = jp+jb 

                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(OP%xindx,ja,jb,Jtot1,jh,jp,Jtot2,J3,J4,rank) * &
                                     tensor_elem(a,h,p,b,J3,J4,OP,jbas) 

                             end do
                          end do

                          CCME%CCX(q1)%X(Rindx,NBindx2) = sm * & 
                               (-1) **( (jh+jb+Jtot2) / 2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))
                          ! NOTE THAT the ph belongs to the LARGER J
                       end if
                    end if
                 end if
              else
                 
                 if ((triangle(ja,jb,Jtot2)) .and. (NBindx1 .ne. 0) ) then 

                    if ( mod(la+lb+op%dpar/2,2) == PAR ) then 
                       x = CCindex(b,a,OP%Nsp) 
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qTWO )
                          g = g + 1
                       end do

                       Gindx = CCME%rmap(x)%Z(g)


                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  

                          ! hapb 
                          Jmin1 = abs(ja - jh) 
                          Jmax1 = ja+jh 
                          Jmin2 = abs(jp - jb)
                          Jmax2 = jp+jb 

                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(OP%xindx,jh,jp,Jtot1,ja,jb,Jtot2,J3,J4,rank) * &
                                     tensor_elem(h,a,b,p,J3,J4,OP,jbas) 

                             end do
                          end do

                          ! store  ( V )_h(p)b(a)
                          ! NOTE THAT the ph belongs to the SMALLER J
                          CCME%CCR(q1)%X(NBindx1,Gindx) = sm * &
                               (-1) **( (jp+ja+Jtot2)/2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))&
                               * (-1)**((jp+jh)/2) 
                       end if
                    end if
                 end if
              end if
           end do
        end do
     end do
  end do

end subroutine calculate_single_pandya
!=======================================================  
!=======================================================          
subroutine calculate_generalized_pandya(OP,CCME,jbas) 
  ! currently the only CCME of interest are phab terms    |---<--|  J1 
  ! coupling in the 3-1 channel                        <(pa)|V|(hb)> rank
  !                                                      |---<--| J2 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: OP 
  type(pandya_mat) :: CCME
  integer :: Jtot1,Jtot2,ja,jp,jb,jh,JC,q1,q2,q,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin1,Jmax1,Rindx,Gindx,g,ta,tb,Atot,hg,pg,J3,J4,NBindx2,qONE,qTWO
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx1,x,JTM,rank,Jmin2,Jmax2
  real(8) :: sm,sm2,pre,horse
  logical :: parflip

  Atot = OP%belowEF
  Ntot = OP%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  rank = OP%rank

  parflip = .false. 
  if ( mod(op%dpar/2,2) == 1) parflip = .true. 

  CCME%herm = OP%Herm
!!$omp parallel do default(firstprivate),shared(CCME,OP,jbas) 
  do q1 = 1, CCME%nblocks
      
     Jtot1 = CCME%Jval(q1)
     Jtot2 = CCME%Jval2(q1)
     if (Jtot2 > 2*JTM) cycle
     
     PAR = mod(q1-1,2)
     Tz = mod((q1-1)/2,2) 
     
     CCME%CCR(q1)%X = 0.d0
     CCME%CCX(q1)%X = 0.d0
         
     qONE = block_index(Jtot1,Tz,Par)
     qTWO = block_index(Jtot2,Tz,mod(Par+op%dpar/2,2))
     
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. (parflip)) then 
              if ( mod(lp + lh,2) .ne. PAR ) cycle
           end if 
           
           if (abs(tp - th)/2 .ne. Tz ) cycle 
        
           NBindx1 = 0
           NBindx2 = 0 
           if ( triangle(jp,jh,Jtot1) )  then 
              
              
              if ( mod(lp + lh,2) ==  PAR ) then   
                 x = CCindex(p,h,OP%Nsp)
                 gnb = 1
              
                 do while (CCME%qmap(x)%Z(gnb) .ne. qONE )
                    gnb = gnb + 1 
                 end do
              
                 NBindx1 = CCME%nbmap(x)%Z(gnb) 
             
              end if 
              
              if  ( triangle(jp,jh,Jtot2) )  then 
              
                 if ( mod(lp+lh+op%dpar/2,2) == PAR) then 
        
                    x = CCindex(p,h,OP%Nsp)
                    gnb = 1
                    
                    do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                       gnb = gnb + 1 
                    end do
                 
                    NBindx2 = CCME%nbmap(x)%Z(gnb)
                 end if 
              end if
 
           else if  ( triangle(jp,jh,Jtot2) )  then 
              
              if ( mod(lp+lh+op%dpar/2,2) == PAR) then 
        
                 x = CCindex(p,h,OP%Nsp)
                 gnb = 1
              
                 do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                    gnb = gnb + 1 
                 end do
              
                 NBindx2 = CCME%nbmap(x)%Z(gnb)
              end if
           else 
              cycle
           end if 
          
           
           ! for the ph  channel 2body derivative 
        
           do a = 1, OP%nsp
              do b = 1, OP%nsp
           
                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 

                 if (.not. (parflip)) then 
                    if ( mod(la + lb,2) .ne. PAR ) cycle
                 end if 

                 if (abs(ta - tb)/2 .ne. Tz ) cycle 

       
                 if ( (triangle(ja,jb,Jtot1)) .and. (NBindx2 .ne. 0) ) then 
       
                    if ( mod(la+lb,2) == PAR ) then 
                       x = CCindex(a,b,OP%Nsp) 
                       
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qONE )
                          g = g + 1
                       end do
              
                       Rindx = CCME%rmap(x)%Z(g)
                 
                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  
                         
                          ! hapb 
                          Jmin1 = abs(ja - jh) 
                          Jmax1 = ja+jh 
                          Jmin2 = abs(jp - jb)
                          Jmax2 = jp+jb 
                          
                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(OP%xindx,ja,jb,Jtot1,jh,jp,Jtot2,J3,J4,rank) * &
                                     tensor_elem(a,h,p,b,J3,J4,OP,jbas) 
                                
                             end do
                          end do
                                               
                          CCME%CCX(q1)%X(Rindx,NBindx2) = sm * & 
                               (-1) **( (jh+jb+Jtot2) / 2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))
                          ! NOTE THAT the ph belongs to the LARGER J
                       end if
                    end if
                 end if 
                 
                 if ((triangle(ja,jb,Jtot2)) .and. (NBindx1 .ne. 0) ) then 
                    
                    if ( mod(la+lb+op%dpar/2,2) == PAR ) then 
                       x = CCindex(b,a,OP%Nsp) 
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qTWO )
                          g = g + 1
                       end do
                    
                       Gindx = CCME%rmap(x)%Z(g)
                 

                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  
               
                          ! hapb 
                          Jmin1 = abs(ja - jh) 
                          Jmax1 = ja+jh 
                          Jmin2 = abs(jp - jb)
                          Jmax2 = jp+jb 
                    
                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(OP%xindx,jh,jp,Jtot1,ja,jb,Jtot2,J3,J4,rank) * &
                                     tensor_elem(h,a,b,p,J3,J4,OP,jbas) 
                       
                             end do
                          end do

                          ! store  ( V )_h(p)b(a)
                          ! NOTE THAT the ph belongs to the SMALLER J
                          CCME%CCR(q1)%X(NBindx1,Gindx) = sm * &
                               (-1) **( (jp+ja+Jtot2)/2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))&
                               * (-1)**((jp+jh)/2) 
                       end if
                    end if
                 end if 
                 
              end do
           end do
        end do
     end do

  end do 
!!$omp end parallel do

end subroutine calculate_generalized_pandya
!=======================================================  
!=======================================================          
subroutine EOM_generalized_pandya(OP,CCME,jbas) 
  ! currently the only CCME of interest are phab terms    |---<--|  J1 
  ! coupling in the 3-1 channel                        <(pa)|V|(hb)> rank
  !                                                      |---<--| J2 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: OP 
  type(ex_pandya_mat) :: CCME
  integer :: Jtot1,Jtot2,ja,jp,jb,jh,JC,q1,q2,q,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin1,Jmax1,Rindx,Gindx,g,ta,tb,Atot,hg,pg,J3,J4,NBindx2,qONE,qTWO
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx1,x,JTM,rank,Jmin2,Jmax2,bx,ax
  real(8) :: sm,sm2,pre,horse
  logical :: phase,parflip

  Atot = OP%belowEF
  Ntot = OP%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  rank = OP%rank

  parflip = .false. 
  if ( mod(op%dpar/2,2) == 1) parflip = .true. 

  CCME%herm = OP%Herm
!$omp parallel do default(firstprivate),shared(CCME,OP,jbas) 
  do q1 = 1, CCME%nblocks
      
     Jtot1 = CCME%Jval(q1)
     Jtot2 = CCME%Jval2(q1)
     if (Jtot2 > 2*JTM) cycle
     
     PAR = mod(q1-1,2)
     Tz = mod((q1-1)/2,2) 
     
     CCME%CCR(q1)%X = 0.d0
     CCME%CCX(q1)%X = 0.d0
         
     qONE = block_index(Jtot1,Tz,Par)
     qTWO = block_index(Jtot2,Tz,mod(Par+op%dpar/2,2))
     
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. (parflip)) then 
              if ( mod(lp + lh,2) .ne. PAR ) cycle
           end if 
           
           if (abs(tp - th)/2 .ne. Tz ) cycle 
        
           NBindx1 = 0
           NBindx2 = 0 
           if ( triangle(jp,jh,Jtot1) )  then 
              
              
              if ( mod(lp + lh,2) ==  PAR ) then   
                 x = PHindex(p,h,Atot)
                 gnb = 1
              
                 do while (CCME%qmap(x)%Z(gnb) .ne. qONE )
                    gnb = gnb + 1 
                 end do
              
                 NBindx1 = CCME%nbmap(x)%Z(gnb) 
             
              end if 
              
              if  ( triangle(jp,jh,Jtot2) )  then 
              
                 if ( mod(lp+lh+op%dpar/2,2) == PAR) then 
        
                    x = PHindex(p,h,Atot)
                    gnb = 1
                    
                    do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                       gnb = gnb + 1 
                    end do
                 
                    NBindx2 = CCME%nbmap(x)%Z(gnb)
                 end if 
              end if
 
           else if  ( triangle(jp,jh,Jtot2) )  then 
              
              if ( mod(lp+lh+op%dpar/2,2) == PAR) then 
        
                 x = PHindex(p,h,Atot)
                 gnb = 1
              
                 do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                    gnb = gnb + 1 
                 end do
              
                 NBindx2 = CCME%nbmap(x)%Z(gnb)
              end if
           else 
              cycle
           end if 
          
           
           
           do ax = 1, OP%belowEF
              a = jbas%holes(ax)
              do bx = 1, OP%nsp - OP%belowEF
                 b = jbas%parts(bx) 
                 
                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 

                 if (.not. (parflip)) then 
                    if ( mod(la + lb,2) .ne. PAR ) cycle
                 end if 

                 if (abs(ta - tb)/2 .ne. Tz ) cycle 

       
                 if ( (triangle(ja,jb,Jtot1)) .and. (NBindx2 .ne. 0) ) then 
       
                    if ( mod(la+lb,2) == PAR ) then 
                       x = PHindex(b,a,Atot) 
                       
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qONE )
                          g = g + 1
                       end do
              
                       Rindx = CCME%nbmap(x)%Z(g)
                 
                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  
                         
                          ! hapb 
                          Jmin2 = abs(ja - jh) 
                          Jmax2 = ja+jh 
                          Jmin1 = abs(jp - jb)
                          Jmax1 = jp+jb 
                          
                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(OP%xindx,jp,jh,Jtot2,jb,ja,Jtot1,J3,J4,rank) * &
                                     tensor_elem(p,b,a,h,J3,J4,OP,jbas) 
                                
                             end do
                          end do
                          
                          ! STORED SUCH THAT THE PH JTOT is GREATER. 
                          CCME%CCX(q1)%X(Rindx,NBindx2) = sm * & 
                               (-1) **( (jh+jb+Jtot1) / 2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))
                          ! stored backwards because. 

                       end if
                    end if
                 end if 
                 
                 if ((triangle(ja,jb,Jtot2)) .and. (NBindx1 .ne. 0) ) then 
                    
                    if ( mod(la+lb+op%dpar/2,2) == PAR ) then 
                       x = PHindex(b,a,Atot) 
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qTWO )
                          g = g + 1
                       end do
                    
                       Gindx = CCME%nbmap(x)%Z(g)
                 

                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  
               
                          ! p(h)a(b)  
                          Jmin2 = abs(ja - jh) 
                          Jmax2 = ja+jh 
                          Jmin1 = abs(jp - jb)
                          Jmax1 = jp+jb 
                    
                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(OP%xindx,jp,jh,Jtot1,jb,ja,Jtot2,J3,J4,rank) * &
                                     tensor_elem(p,b,a,h,J3,J4,OP,jbas) 
                       
                             end do
                          end do

                          ! store  ( V )_p(h)a(b)

                          ! PH JTOT IS LESS THAN OR EQUAL TO.
                          CCME%CCR(q1)%X(NBindx1,Gindx) = sm * &
                               (-1) **( (jb+jh+Jtot2)/2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))
                       end if
                    end if
                 end if 
                 
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine EOM_generalized_pandya
!=======================================================  
!=======================================================          
subroutine EOM_scalar_cross_coupled(HS,CCME,jbas) 
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: HS
  type(cc_mat) :: CCME
  integer :: JT,ja,jp,jb,jh,JC,q1,q2,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin,Jmax,Rindx,Gindx,g,ta,tb,Atot,hg,pg
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx,x,JTM
  real(8) :: sm,sm2,pre,horse

  Atot = HS%belowEF
  Ntot = HS%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 

!$omp parallel do default(firstprivate),shared(CCME,HS,jbas) 
  do q1 = 1, CCME%nblocks
      
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
    
     CCME%CCX(q1)%X = 0.d0
         
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. triangle(jp,jh,JC) )  cycle
           if ( mod(lp + lh,2) .ne. PAR ) cycle
           if (abs(tp - th)/2 .ne. Tz ) cycle 
           
           x = CCindex(p,h,HS%Nsp)
           gnb = 1
           do while (CCME%qmap(x)%Z(gnb) .ne. q1 )
              gnb = gnb + 1
           end do
              
           NBindx = CCME%nbmap(x)%Z(gnb) 

           ! for the ph  channel 2body derivative 
        
           do a = 1, HS%nsp
              if ( jbas%con(a) .ne. 1 )  cycle

              do b = 1, HS%nsp
                 if (jbas%con(b) .ne. 0) cycle

                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 
                 if (.not. triangle(ja,jb,JC) )  cycle
                 if ( mod(la + lb,2) .ne. PAR ) cycle
                 if (abs(ta - tb)/2 .ne. Tz ) cycle 
                 
                 x = CCindex(a,b,HS%Nsp) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Rindx = CCME%rmap(x)%Z(g)
                
                 sm = 0.d0 
               
!                 horse = 0.d0 
                 if ( (mod(la + lh,2) == mod(lb + lp,2)) .and. &
                      ( (ta + th) == (tb + tp) ) ) then  
               
                    ! hapb 
                    Jmin = max(abs(jp - jb),abs(ja - jh)) 
                    Jmax = min(jp+jb,ja+jh) 
                    
                    sm = 0.d0 
                    do JT = Jmin,Jmax,2
                       sm = sm + (-1)**(JT/2) * (JT + 1) * &
                            sixj(jp,jh,JC,ja,jb,JT)  * &
                            v_elem(h,a,p,b,JT,HS,jbas) 
                    end do
                 
                    ! store  < p b | v | h a> 
                    ! scaled by sqrt(JC + 1) for convience in ph derivative
                             
                    CCME%CCX(q1)%X(Rindx,NBindx) = sm * HS%herm * &
                         (-1) **( (jh - ja + JC) / 2) * sqrt(JC + 1.d0)
                 
                 end if
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine EOM_scalar_cross_coupled
!=====================================================
!=====================================================      
integer function specific_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(cc_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1
  do while (LCC%qmap(x)%Z(g) .ne. q )
  
     g = g + 1
  end do
  
  specific_rval = LCC%rmap(x)%Z(g)
end function specific_rval
!=====================================================
!=====================================================      
integer function ex_specific_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(ex_cc_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = PHindex(i,l,Ntot)
  g = 1
  do while (LCC%qmap(x)%Z(g) .ne. q )
  
     g = g + 1
  end do
  
  ex_specific_rval = LCC%nbmap(x)%Z(g)
end function ex_specific_rval
!============================================
!============================================
integer function reg_ph_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(cc_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1

  do while (LCC%qmap(x)%Z(g) .ne. q )
     g = g + 1
  end do
  
  reg_ph_rval = LCC%nbmap(x)%Z(g)
end function reg_ph_rval
!=====================================================
!=====================================================   
integer function pandya_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(pandya_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1
  do while (LCC%qmap(x)%Z(g) .ne. q )  
     g = g + 1
  end do
  
  pandya_rval = LCC%rmap(x)%Z(g)
end function pandya_rval
!=====================================================
!=====================================================   
integer function ex_ph_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(ex_cc_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1

  do while (LCC%qmap(x)%Z(g) .ne. q )
     g = g + 1
  end do
  
  ex_ph_rval = LCC%nbmap(x)%Z(g)
end function ex_ph_rval
!=====================================================
!=====================================================   
integer function ex_pandya_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(ex_pandya_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = PHindex(i,l,Ntot)
  g = 1
  do while (LCC%qmap(x)%Z(g) .ne. q )  
     g = g + 1
  end do
  
  ex_pandya_rval = LCC%nbmap(x)%Z(g)
end function ex_pandya_rval
!=================================================================     
!=================================================================
real(8) function WCCX(a,d,b,c,J1,J2,WX,RCC,jbas)
  implicit none
  
  integer :: a,b,c,d,rad,rbc,Ntot,q1,q2
  integer :: J1,J2,PAR,PAR2,Tz,qx,rank  
  type(spd) :: jbas
  type(pandya_mat) :: RCC 
  real(8),dimension(:,:) :: WX 
  
  Ntot = jbas%total_orbits
  rank = RCC%rank
  
  PAR = mod(jbas%ll(a) + jbas%ll(d) ,2) 
  Tz = abs(jbas%itzp(a) - jbas%itzp(d))/2 
  
  PAR2 = mod(PAR + RCC%dpar/2,2)
  if ( abs(jbas%itzp(b) - jbas%itzp(c)) .ne. 2*Tz) then 
     WCCX = 0.d0 
     return
  end if 
  
  if (mod(jbas%ll(b) + jbas%ll(c) ,2)  .ne. PAR2) then 
     WCCX = 0.d0 
     return
  end if
  
  q1 = block_index(J1,Tz,PAR)
  q2 = block_index(J2,Tz,PAR2)
  
  rad = fetch_rval(a,d,Ntot,q1,RCC)
  rbc = fetch_rval(b,c,Ntot,q2,RCC)
  
  WCCX = WX(rad,rbc)

end function WCCX
!============================================================
!============================================================
real(8) function WCCR(a,d,b,c,J1,J2,WCC,RCC,jbas)
  implicit none
  
  integer :: a,b,c,d,rad,rbc,Ntot,q1,q2
  integer :: J1,J2,PAR,PAR2,Tz,qx,rank  
  type(spd) :: jbas
  type(pandya_mat) :: WCC,RCC 
  
  Ntot = jbas%total_orbits
  rank = RCC%rank
  
  PAR = mod(jbas%ll(a) + jbas%ll(d) ,2) 
  Tz = abs(jbas%itzp(a) - jbas%itzp(d))/2 
  
  PAR2 = mod(PAR + RCC%dpar/2,2)
  if ( abs(jbas%itzp(b) - jbas%itzp(c)) .ne. 2*Tz) then 
     WCCR = 0.d0 
     return
  end if 
  
  if (mod(jbas%ll(b) + jbas%ll(c) ,2)  .ne. PAR2) then 
     WCCR = 0.d0 
     return
  end if
  
  q1 = block_index(J1,Tz,PAR)
  q2 = block_index(J2,Tz,PAR2)
  
  rad = fetch_rval(a,d,Ntot,q1,RCC)
  rbc = fetch_rval(b,c,Ntot,q2,RCC)
  
  qx = CCtensor_block_index(J1,J2,rank,Tz,PAR)
  
  WCCR = WCC%CCR(qx)%X(rad,rbc)

end function WCCR
!============================================================
!============================================================
end module

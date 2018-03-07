module brute_force_testing
  use cross_coupled
  use isospin_operators
  use commutators
  use TS_commutators
  use EOM_TS_commutators
  use EOM_dTZ_commutators
  use operator_commutators
  use EOM_scalar_commutators
  
  
contains
!============================================================
!============================================================
subroutine construct_random_rank0(OP,HERM,jbas) 
  implicit none 
  
  integer,intent(in):: HERM
  type(sq_op) :: OP 
  type(spd) :: jbas
  integer :: i,j,k,l,ji,jj,jk,jl
  integer :: J1,J2,IX,JX,q,ig,jg 
  real(8) :: x
  
  OP%rank = 0
  OP%herm = HERM 
  
  do ig = 1, OP%belowEF
     do jg = ig, OP%belowEF
        
        i = jbas%holes(ig)
        j = jbas%holes(jg)
 
        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i) .ne. jbas%ll(j) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fhh(ig,jg) = x
        OP%fhh(jg,ig) = OP%herm * x
        
        
     end do
     
     OP%fhh(ig,ig) = (OP%herm+1)*OP%fhh(ig,ig)/2.d0
     
  end do

  do ig = 1, OP%nsp-OP%belowEF
     do jg = 1, OP%belowEF

        i = jbas%parts(ig)
        j = jbas%holes(jg)
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i) .ne. jbas%ll(j) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fph(ig,jg) = x

     end do
  end do

  do ig = 1,OP%nsp- OP%belowEF
     do jg = ig,OP%nsp- OP%belowEF
        
        i = jbas%parts(ig)
        j = jbas%parts(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i) .ne. jbas%ll(j) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fpp(ig,jg) = x 
        OP%fpp(jg,ig) = OP%herm * x
     end do
     OP%fpp(ig,ig) = (OP%herm+1)*OP%fpp(ig,ig)/2.d0
  end do
  
  do q=1,OP%nblocks
     do i = 1, 6
        do IX = 1,size(OP%mat(q)%gam(i)%X(:,1))
           do JX = 1,size(OP%mat(q)%gam(i)%X(1,:))
              
              call random_number(x) 
              x = 10.*x-5.

              ! PAULI PRINCIPLE 
              if ( OP%mat(q)%qn(sea1(i))%Y(IX,1) == OP%mat(q)%qn(sea1(i))%Y(IX,2) ) then 
                 if ( mod( OP%mat(q)%lam(1)/2, 2) == 1) x = 0.d0 
              end if 

              if ( OP%mat(q)%qn(sea2(i))%Y(JX,1) == OP%mat(q)%qn(sea2(i))%Y(JX,2) ) then 
                 if ( mod( OP%mat(q)%lam(1)/2, 2) == 1) x = 0.d0 
              end if 
              
              OP%mat(q)%gam(i)%X(IX,JX) =  x * (1 + OP%herm *kron_del(IX,JX) ) 
              if (sqs(i)) OP%mat(q)%gam(i)%X(JX,IX) = OP%herm*x*  (1 + OP%herm *kron_del(IX,JX) )  
              
              
           end do 
        end do 
     end do
  end do 

end subroutine 
!============================================================
!============================================================
subroutine construct_random_rankX(OP,HERM,jbas) 
  implicit none 
  
  integer,intent(in) :: HERM
  type(sq_op) :: OP 
  type(spd) :: jbas
  integer :: i,j,k,l,ji,jj,jk,jl,ig,jg
  integer :: J1,J2,IX,JX,q,qx,rank
  real(8) :: x

  rank = OP%rank
  OP%herm = HERM 

  do ig = 1, OP%belowEF
     do jg = ig, OP%belowEF
        
        i = jbas%holes(ig)
        j = jbas%holes(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,OP%rank)) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j)+op%dpar/2,2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fhh(ig,jg) = x
        OP%fhh(jg,ig) = (-1)**((ji-jj)/2) * x* OP%herm
         
     end do
     OP%fhh(ig,ig) = (OP%herm+1)*OP%fhh(ig,ig)/2.d0      
  end do

  do ig = 1, OP%nsp- OP%belowEF
     do jg = ig, OP%nsp-OP%belowEF
        
        i = jbas%parts(ig)
        j = jbas%parts(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,OP%rank)) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j)+op%dpar/2,2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fpp(ig,jg) = x
        OP%fpp(jg,ig) = (-1)**((ji-jj)/2) * x* OP%herm 
     end do 
     OP%fpp(ig,ig) = (OP%herm+1)*OP%fpp(ig,ig)/2.d0 
  end do 
  
  do ig = 1, OP%nsp- OP%belowEF
     do jg = 1, OP%belowEF
        
        i = jbas%parts(ig)
        j = jbas%holes(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,OP%rank)) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j)+op%dpar/2,2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fph(ig,jg) = x
     end do 
  end do 

  do q = 1, OP%nblocks
     do i = 1, 9
         
        do IX = 1, size(OP%tblck(q)%tgam(i)%X(:,1) )
           do JX = 1, size(OP%tblck(q)%tgam(i)%X(1,:) )
              
              call random_number(x)

              x = 10.*x-5.              
              ! PAULI PRINCIPLE 
              if ( OP%tblck(q)%tensor_qn(sea1(i),1)%Y(IX,1) == OP%tblck(q)%tensor_qn(sea1(i),1)%Y(IX,2) ) then 
                 if ( mod( OP%tblck(q)%Jpair(1)/2, 2) == 1) x = 0.d0 
              end if
              
              if ( OP%tblck(q)%tensor_qn(sea2(i),2)%Y(JX,1) == OP%tblck(q)%tensor_qn(sea2(i),2)%Y(JX,2) ) then 
                 if ( mod( OP%tblck(q)%Jpair(2)/2, 2) == 1) x = 0.d0 
              end if
              
              OP%tblck(q)%tgam(i)%X(IX,JX) = x
              
           end do
        end do
        
        

        if ( (OP%tblck(q)%Jpair(1) == OP%tblck(q)%Jpair(2)).and.(sqs(i)) ) then
           if (mod(OP%dpar/2,2) == 0 ) then
              OP%tblck(q)%tgam(i)%X = &
                   0.5*(OP%tblck(q)%tgam(i)%X + OP%herm * Transpose(OP%tblck(q)%tgam(i)%X))
           
           else
              qx = tensor_block_index(OP%tblck(q)%Jpair(1),OP%tblck(q)%Jpair(2)&
                   ,rank,OP%tblck(q)%lam(3),mod(OP%tblck(q)%lam(2)+1,2))
              
              OP%tblck(qx)%tgam(i)%X = Transpose( OP%tblck(q)%tgam(i)%X ) * OP%herm 
              
           end if 
       end if 
        
     end do
     
     
     if ( OP%tblck(q)%Jpair(1) == OP%tblck(q)%Jpair(2) ) then
        if (mod(OP%dpar/2,2) == 0 ) then 

           OP%tblck(q)%tgam(7)%X = Transpose(OP%tblck(q)%tgam(3)%X) * OP%herm
           OP%tblck(q)%tgam(8)%X = Transpose(OP%tblck(q)%tgam(2)%X) * OP%herm
           OP%tblck(q)%tgam(9)%X = Transpose(OP%tblck(q)%tgam(6)%X) * OP%herm
        else
           

           qx = tensor_block_index(OP%tblck(q)%Jpair(1),OP%tblck(q)%Jpair(2)&
                   ,rank,OP%tblck(q)%lam(3),mod(OP%tblck(q)%lam(2)+1,2))
           
           OP%tblck(qx)%tgam(7)%X = Transpose(OP%tblck(q)%tgam(3)%X)*OP%herm
           OP%tblck(qx)%tgam(8)%X = Transpose(OP%tblck(q)%tgam(2)%X)*OP%herm 
           OP%tblck(qx)%tgam(9)%X = Transpose(OP%tblck(q)%tgam(6)%X)*OP%herm
           OP%tblck(qx)%tgam(3)%X = Transpose(OP%tblck(q)%tgam(7)%X)*OP%herm
           OP%tblck(qx)%tgam(2)%X = Transpose(OP%tblck(q)%tgam(8)%X)*OP%herm 
           OP%tblck(qx)%tgam(6)%X = Transpose(OP%tblck(q)%tgam(9)%X)*OP%herm
        
        end if
      
     end if 
     
  end do

end subroutine construct_random_rankX
!============================================================
!============================================================
subroutine construct_random_opX(OP,HERM,jbas) 
  implicit none 
  
  integer,intent(in) :: HERM
  type(iso_operator) :: OP 
  type(spd) :: jbas
  integer :: i,j,k,l,ji,jj,jk,jl,ig,jg
  integer :: J1,J2,IX,JX,q,qx,rank
  real(8) :: x

  rank = OP%rank
  OP%herm = HERM 
  
  do i = 1, OP%nsp
     do j = 1, OP%nsp
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,OP%rank)) cycle        
        if (jbas%itzp(i)-OP%dTz*2 .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j)+op%dpar/2,2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fock(i,j) = x
        if (Op%dTZ == 0 ) OP%fock(j,i) = x*(-1)**((ji-jj)/2)*herm 
     end do 
  end do 

  do q = 1, OP%nblocks
     if (OP%tblck(q)%Jpair(2) > jbas%jtotal_max*2) cycle
     if (abs(OP%tblck(q)%lam(4)) > 1) cycle
     do i = 1, 9
        
        do IX = 1, size(OP%tblck(q)%tgam(i)%X(:,1) )
           do JX = 1, size(OP%tblck(q)%tgam(i)%X(1,:) )
              
              call random_number(x)

              x = 10.*x-5.              
              ! PAULI PRINCIPLE 
              if ( OP%tblck(q)%tensor_qn(sea1(i),1)%Y(IX,1) == OP%tblck(q)%tensor_qn(sea1(i),1)%Y(IX,2) ) then 
                 if ( mod( OP%tblck(q)%Jpair(1)/2, 2) == 1) x = 0.d0 
              end if
              
              if ( OP%tblck(q)%tensor_qn(sea2(i),2)%Y(JX,1) == OP%tblck(q)%tensor_qn(sea2(i),2)%Y(JX,2) ) then 
                 if ( mod( OP%tblck(q)%Jpair(2)/2, 2) == 1) x = 0.d0 
              end if
              
              OP%tblck(q)%tgam(i)%X(IX,JX) = x
                                             
           end do
        end do
     end do

     if (OP%dTZ == 0) then
        ! need to make sure tranposes look okay
         qx = iso_ladder_block_index(OP%tblck(q)%Jpair(2),OP%tblck(q)%Jpair(1)&
                   ,rank,OP%tblck(q)%lam(3),mod(OP%tblck(q)%lam(2)+OP%dpar/2,2))
         if (qx > q) cycle
         if (qx == q) then
            OP%tblck(q)%tgam(4)%X = 0.5*OP%tblck(q)%tgam(4)%X + 0.5*transpose(OP%tblck(q)%tgam(4)%X)
            OP%tblck(q)%tgam(1)%X = 0.5*OP%tblck(q)%tgam(1)%X + 0.5*transpose(OP%tblck(q)%tgam(1)%X)
            OP%tblck(q)%tgam(5)%X = 0.5*OP%tblck(q)%tgam(5)%X + 0.5*transpose(OP%tblck(q)%tgam(5)%X)
         end if
         OP%tblck(qx)%tgam(7)%X = Transpose(OP%tblck(q)%tgam(3)%X)*OP%herm*OP%tblck(q)%lam(1)
         OP%tblck(qx)%tgam(8)%X = Transpose(OP%tblck(q)%tgam(2)%X)*OP%herm*OP%tblck(q)%lam(1) 
         OP%tblck(qx)%tgam(9)%X = Transpose(OP%tblck(q)%tgam(6)%X)*OP%herm*OP%tblck(q)%lam(1)
         OP%tblck(qx)%tgam(3)%X = Transpose(OP%tblck(q)%tgam(7)%X)*OP%herm*OP%tblck(q)%lam(1)
         OP%tblck(qx)%tgam(2)%X = Transpose(OP%tblck(q)%tgam(8)%X)*OP%herm*OP%tblck(q)%lam(1) 
         OP%tblck(qx)%tgam(6)%X = Transpose(OP%tblck(q)%tgam(9)%X)*OP%herm*OP%tblck(q)%lam(1)
         OP%tblck(qx)%tgam(1)%X = Transpose(OP%tblck(q)%tgam(1)%X)*OP%herm*OP%tblck(q)%lam(1)
         OP%tblck(qx)%tgam(4)%X = Transpose(OP%tblck(q)%tgam(4)%X)*OP%herm*OP%tblck(q)%lam(1)
         OP%tblck(qx)%tgam(5)%X = Transpose(OP%tblck(q)%tgam(5)%X)*OP%herm*OP%tblck(q)%lam(1)
      end if
  end do

end subroutine construct_random_opX
!============================================================
!============================================================
subroutine construct_random_isoX(OP,HERM,jbas) 
  implicit none 
  
  integer,intent(in) :: HERM
  type(iso_ladder) :: OP 
  type(spd) :: jbas
  integer :: i,j,k,l,ji,jj,jk,jl,ig,jg
  integer :: J1,J2,IX,JX,q,qx,rank
  real(8) :: x

  rank = OP%rank
  OP%herm = HERM 
  
  do ig = 1, OP%nsp- OP%belowEF
     do jg = 1, OP%belowEF
        
        i = jbas%parts(ig)
        j = jbas%holes(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,OP%rank)) cycle        
        if (jbas%itzp(i)-OP%dTz*2 .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j)+op%dpar/2,2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fph(ig,jg) = x
     end do 
  end do 

  do q = 1, OP%nblocks
         
        do IX = 1, OP%tblck(q)%npp1 
           do JX = 1, OP%tblck(q)%nhh2
              
              call random_number(x)

              x = 10.*x-5.              
              ! PAULI PRINCIPLE 
              if ( OP%tblck(q)%qn(1)%Y(IX,1) == OP%tblck(q)%qn(1)%Y(IX,2) ) then 
                 if ( mod( OP%tblck(q)%Jpair(1)/2, 2) == 1) x = 0.d0 
              end if
              
              if ( OP%tblck(q)%qn(2)%Y(JX,1) == OP%tblck(q)%qn(2)%Y(JX,2) ) then 
                 if ( mod( OP%tblck(q)%Jpair(2)/2, 2) == 1) x = 0.d0 
              end if

              OP%tblck(q)%Xpphh(IX,JX) = x

        end do
        
                
     end do
     
     
  end do

end subroutine construct_random_isoX
!============================================================
!============================================================
subroutine construct_delta_rankX(OP,HERM,jbas) 
  implicit none 
  
  integer,intent(in) :: HERM
  type(sq_op) :: OP 
  type(spd) :: jbas
  integer :: i,j,k,l,ji,jj,jk,jl,ig,jg
  integer :: J1,J2,IX,JX,q,qx,rank,m
  real(8) :: x,sm,dcgi

  rank = OP%rank
  OP%herm = HERM 

  do ig = 1, OP%belowEF
     i = jbas%holes(ig)
     ji = jbas%jj(i)
     sm = 0.d0
     do m = -1*ji,ji,2
        sm = sm + dcgi(ji,m,rank,0,ji,m)
     end do
     IF (ig == 1) then
        print*,i, sm , sqrt(ji+1.d0)
     end if 
     OP%fhh(ig,ig) = sm/sqrt(ji+1.d0)
!     OP%fhh(ig,ig) = (OP%herm+1)*OP%fhh(ig,ig)/2.d0      
  end do

  do ig = 1, OP%nsp- OP%belowEF
     i = jbas%parts(ig)
     ji = jbas%jj(i)
     sm = 0.d0
     do m = -1*ji,ji,2
        sm = sm + dcgi(ji,m,rank,0,ji,m)
     end do 
     OP%fpp(ig,ig) = sm/sqrt(ji+1.d0)
 !    OP%fpp(ig,ig) = (OP%herm+1)*OP%fpp(ig,ig)/2.d0 
  end do 

end subroutine
!============================================================
!============================================================
subroutine seed_random_number
  implicit none 
  
  integer :: i,a
  real(8) :: x
  logical :: ext
  
  inquire(file='seed',exist=ext)
  
  if (ext) then 
     open(unit=41,file='seed')
     read(41,*) a
     close(41)
  else 
     a = 2934
  end if 
  
  
  do i = 1, a 
     call random_number(x)
  end do 
  
  a = nint(12349.d0*x)
  
  open(unit=41,file='seed')
  write(41,*)  a
  close(41)

end subroutine
!============================================================
!============================================================
subroutine test_scalar_scalar_commutator(jbas,h1,h2) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2
  type(cc_mat) :: AACC,BBCC,WCC
  integer :: a,b,c,d,ja,jb,jc,jd,jmin,jmax,PAR,TZ,Jtot
  integer :: hole,part,iii
  integer,intent(in) :: h1,h2
  real(8) :: val,sm
  real(8) :: vv,xx,yy,zz
    
  call seed_random_number
  
  call allocate_blocks(jbas,AA)
  call duplicate_sq_op(AA,BB)
  call duplicate_sq_op(AA,OUT)
  call duplicate_sq_op(AA,w1) !workspace
  call duplicate_sq_op(AA,w2) !workspace
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
  call duplicate_ph_mat(AACC,BBCC) !cross coupled ME
  call init_ph_wkspc(BBCC,WCC) !
  
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_rank0(BB,h2,jbas) 

  OUT%herm = -1* AA%herm * BB%herm 
  
  print*, 'TESTING SCALAR-SCALAR COMMUTATORS' 
  
  OUT%E0 = commutator_110(AA,BB,jbas) + commutator_220(AA,BB,jbas)
  
  val = scalar_scalar_0body_comm(AA,BB,jbas)
 
  if ( abs( val - OUT%E0) > 1e-10 ) then 
     print*, val, OUT%E0
     STOP 'ZERO BODY FAILURE' 
  end if 
  
  call calculate_cross_coupled(BB,BBCC,jbas)
  call calculate_cross_coupled(AA,AACC,jbas)
  
  call commutator_111(AA,BB,OUT,jbas) 
  call commutator_121(AA,BB,OUT,jbas)
  call commutator_122(AA,BB,OUT,jbas)
 
  call commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)
 
  call commutator_221(AA,BB,OUT,w1,w2,jbas)
  call commutator_222_ph(AACC,BBCC,OUT,WCC,jbas)
 
 
!  do a = 1, jbas%total_orbits
 !    do b = 1, jbas%total_orbits
  do iii = 1, 50   
     call random_number(vv)
     call random_number(yy)
   
     a = ceiling(vv*(AA%Nsp))
     b = ceiling(yy*(AA%Nsp))
        
        val = scalar_scalar_1body_comm(AA,BB,a,b,jbas) 
        
        if (abs(val-f_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if 
        
        print*, 'success:', a,b
  !   end do 
  end do 
 
  !do a = 1, jbas%total_orbits
  iii = 0 
  do while (iii < 15)  
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)
   
     a = ceiling(vv*AA%Nsp)
     b = ceiling(xx*AA%Nsp)
     c = ceiling(yy*AA%Nsp)
     d = ceiling(zz*AA%Nsp)
     
     ja = jbas%jj(a) 
   !  do b = 1, jbas%total_orbits
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
    !    do c = 1, jbas%total_orbits
           jc = jbas%jj(c)
     !      do d = 1, jbas%total_orbits
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d),2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              iii = iii + 1
              jmin = max( abs(ja-jb) , abs(jc-jd) )
              jmax = min( ja+jb , jc+jd) 
              
              do Jtot = jmin,jmax,2
                 
                 val = scalar_scalar_2body_comm(AA,BB,a,b,c,d,Jtot,jbas)
                 
                 if (abs(val-v_elem(a,b,c,d,Jtot,OUT,jbas)) > 1e-6) then
                    print*, 'at:',a,b,c,d, 'J:', Jtot ,val,v_elem(a,b,c,d,Jtot,OUT,jbas)                
                    print*, val,v_elem(a,b,c,d,Jtot,OUT,jbas)
                    STOP 'TWO BODY FAILURE'  
                 end if
                 
              end do 
              
              print*, 'success:', a,b,c,d
!           end do
 !       end do
  !   end do
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine
!============================================================
!============================================================
subroutine test_EOM_scalar_scalar_commutator(jbas,h1,h2) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2
  type(cc_mat) :: AACC,BBCC,WCC
  integer :: a,b,c,d,q,g,ja,jb,jc,jd,jmin
  integer :: ax,bx,cx,dx,jmax,PAR,TZ,Jtot
  integer,intent(in) :: h1,h2
  real(8) :: val
  
  
  call seed_random_number
  
  call allocate_blocks(jbas,AA)
  call duplicate_sq_op(AA,BB)
  call duplicate_sq_op(AA,OUT)
  call duplicate_sq_op(AA,w1) !workspace
  call duplicate_sq_op(AA,w2) !workspace
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
  call duplicate_ph_mat(AACC,BBCC) !cross coupled ME
  call init_ph_wkspc(BBCC,WCC) !
  
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_rank0(BB,h2,jbas) 
  
  ! make BB look like an excitation operator. 
  BB%fhh = 0.d0 
  BB%fpp = 0.d0 
  
  do q = 1, BB%nblocks
     do g = 1, 6
        if (g==3) cycle
        BB%mat(q)%gam(g)%X = 0.d0 
     end do 
  end do 

  OUT%herm = -1* AA%herm * BB%herm 
  
  print*, 'TESTING EOM SCALAR-SCALAR COMMUTATORS' 
  
  OUT%E0 = EOM_scalar_commutator_110(AA,BB,jbas) + &
       EOM_scalar_commutator_220(AA,BB,jbas)
  
  val = EOM_scalar_scalar_0body_comm(AA,BB,jbas)
 
  if ( abs( val - OUT%E0) > 1e-10 ) then 
     print*, val, OUT%E0
     STOP 'ZERO BODY FAILURE' 
  end if 
  
  call EOM_scalar_cross_coupled(BB,BBCC,jbas)
  call calculate_cross_coupled(AA,AACC,jbas) 
  
  call EOM_scalar_commutator_111(AA,BB,OUT,jbas) 
  call EOM_scalar_commutator_121(AA,BB,OUT,jbas)
  call EOM_scalar_commutator_122(AA,BB,OUT,jbas)
  
  call EOM_scalar_commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)
  
  call EOM_scalar_commutator_221(AA,BB,OUT,w1,w2,jbas)
  call EOM_scalar_commutator_222_ph(AACC,BBCC,OUT,WCC,jbas)
  
  do ax = 1, AA%Nsp-AA%belowEF
     a = jbas%parts(ax)
     do bx = 1, AA%belowEF
        b= jbas%holes(bx) 
        
        val = EOM_scalar_scalar_1body_comm(AA,BB,a,b,jbas) 
        
        if (abs(val-f_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if 
        
        print*, 'success:', a,b
     end do 
  end do 
 
  do ax = 1, AA%Nsp-AA%belowEF
     a = jbas%parts(ax)
     ja = jbas%jj(a) 
     do bx = 1, AA%Nsp-AA%belowEF
        b = jbas%parts(bx)
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
        do cx = 1, AA%belowEF
           c = jbas%holes(cx) 
           jc = jbas%jj(c)
           do dx = 1, AA%belowEF
              d = jbas%holes(dx) 
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d),2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              
              jmin = max( abs(ja-jb) , abs(jc-jd) )
              jmax = min( ja+jb , jc+jd) 
              
              do Jtot = jmin,jmax,2
                 
                 val = EOM_scalar_scalar_2body_comm(AA,BB,a,b,c,d,Jtot,jbas)
                 
                 if (abs(val-v_elem(a,b,c,d,Jtot,OUT,jbas)) > 1e-6) then
                    print*, 'at:',a,b,c,d, 'J:', Jtot ,val,v_elem(a,b,c,d,Jtot,OUT,jbas)                
                    print*, val,v_elem(a,b,c,d,Jtot,OUT,jbas)
                    STOP 'TWO BODY FAILURE'  
                 end if
              end do 
              
              print*, 'success:', a,b,c,d
           end do
        end do
     end do
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine
!============================================================
subroutine compare_tensor_scalar_commutator(jbas,h1,h2) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,OUTs,w1,w2,BBY,w1s,w2s
  type(pandya_mat) :: BBCC,WCC
  type(cc_mat) :: AACC,BBYC,WCCs
  integer :: a,b,c,d,ja,jb,jc,jd,j1min,j1max,rank,q,Jab,Jdi
  integer :: j2min,j2max,PAR,TZ,J1,J2,dpar,iii,i,j,jtot
  integer,intent(in) :: h1,h2
  real(8) :: val,t1,t2,t3,t4,omp_get_wtime
  real(8) :: vv,xx,yy,zz,x
  
  call seed_random_number
  
  BB%rank = 0
  BB%dpar = 0
  AA%rank = 0
  BB%pphh_ph = .false.
  
  AA%herm = h1 
  BB%herm = h2 
  BBY%herm = h2 
 
  call allocate_blocks(jbas,AA)
  call allocate_tensor(jbas,BB,AA)
  call duplicate_sq_op(AA,BBY) 
  call duplicate_sq_op(AA,outS) 


  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_rank0(BBY,h2,jbas) 
  
  call copy_rank0_to_tensor_format(BBY,BB,jbas) 
  
  call duplicate_sq_op(BB,OUT)
  call duplicate_sq_op(BB,w1) !workspace
  call duplicate_sq_op(BB,w2) !workspace
  call duplicate_sq_op(BBY,w1s) !workspace
  call duplicate_sq_op(BBY,w2s) !workspace
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
  CALL allocate_small_tensor_CCMAT(BB,BBCC,jbas) 
  ! call init_ph_mat(BB,BBCC,jbas) !cross coupled ME
  call init_ph_mat(BBy,BBYC,jbas)
  call init_ph_wkspc(BBYC,WCCs)  
  
  OUT%herm = -1* AA%herm * BB%herm 
  OUTs%herm = -1* AA%herm * BBY%herm 
 
  print*, 'TESTING SCALAR-TENSOR COMMUTATORS' 

!  call calculate_generalized_pandya(BB,BBCC,jbas) 
  call calculate_cross_coupled(AA,AACC,jbas) 
  
  call TS_commutator_111(AA,BB,OUT,jbas) 
  call TS_commutator_121(AA,BB,OUT,jbas)
  call TS_commutator_211(AACC,BB,OUT,jbas) 
  call TS_commutator_122(AA,BB,OUT,jbas)
  call TS_commutator_212(AA,BB,OUT,jbas)
  
  call TS_commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)  
  call TS_commutator_221(w1,w2,AA%herm*BB%herm,OUT,jbas)
  call TS_commutator_222_ph(AACC,BBCC,BB,OUT,jbas)
q=1 


  
  !do a = 7,30
   !  do b = 1,30


do q = 1,1000 
   call random_number(x) 
   a = ceiling(x*jbas%total_orbits)
   call random_number(x) 
   b = ceiling(x*jbas%total_orbits)
   call random_number(x) 
   c = ceiling(x*jbas%total_orbits)
   call random_number(x) 
   d = ceiling(x*jbas%total_orbits)
   call random_number(x) 
   i = ceiling(x*jbas%total_orbits)
   call random_number(x) 
   j = ceiling(x*jbas%total_orbits)
   
   do jtot = 1,19,2 
      do Jab = abs(jbas%jj(c)-jtot), jbas%jj(c) +jtot, 2 
         do Jdi = abs(jbas%jj(j)-jtot), jbas%jj(j) +jtot, 2 

                    t1 = commutator_223_single(AA,BBY,a,b,c,d,i,j,jtot,Jab,Jdi,jbas)*sqrt(jtot+1.d0)

                    t2 = TS_commutator_223_single(AA,BB,a,b,c,d,i,j,jtot,jtot,Jab,Jdi,jbas)         


                    if ( abs(t1 - t2) > 1e-6) then 
                       print*, t1,t2,'abcdef',a,b,c,d,i,j,'J',Jab,Jdi,jtot 
                      stop
                    else
                       if (abs(t1) > 1e-6 ) print*, t1,t2 
                    end if
           end do
        end do
     end do
  end do 
  
  stop
!  print*, 'time:', t3-t1,t2-t1,t3-t4
 
  call calculate_cross_coupled(BBY,BBYC,jbas)
  call calculate_cross_coupled(AA,AACC,jbas)
  
  call commutator_111(AA,BBY,OUTs,jbas) 
  call commutator_121(AA,BBY,OUTs,jbas)
  call commutator_122(AA,BBY,OUTs,jbas)
 
  call commutator_222_pp_hh(AA,BBY,OUTs,w1s,w2s,jbas)
 
  call commutator_221(AA,BBY,OUTs,w1s,w2s,jbas)
  call commutator_222_ph(AACC,BBYC,OUTs,WCCs,jbas)

  do a = 1, jbas%total_orbits
     do b = 1, jbas%total_orbits
        
        val = f_elem(a,b,OUTs,jbas)*sqrt(jbas%jj(a)+1.d0) 
        
        if (abs(val-f_tensor_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_tensor_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if 
        
        print*, 'success:', a,b
     end do 
  end do 

 do a = 1, jbas%total_orbits
     
     ja = jbas%jj(a) 
     do b = 1, jbas%total_orbits
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
        do c = 1, jbas%total_orbits
           jc = jbas%jj(c)
           do d = 1, jbas%total_orbits
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+BB%dpar/2,2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
  
              j1min = max(abs(ja-jb),abs(jc-jd))  
              j1max = min(ja+jb,jc+jd) 
              
              do J1 = j1min,j1max,2
                    
                    if (.not. (triangle(J1,J2,rank))) cycle
                    
                    val = v_elem(a,b,c,d,J1,OUTs,jbas)*Sqrt(J1+1.d0)
                  
                    if (abs(val-tensor_elem(a,b,c,d,J1,J1,OUT,jbas)) > 1e-8) then
                       print*, 'at:',a,b,c,d, 'J:', J1 ,val,tensor_elem(a,b,c,d,J1,J1,OUT,jbas)                
            
                       STOP 'TWO BODY FAILURE'  
                    end if 
              end do
              
              print*, 'success:', a,b,c,d
           end do
        end do
     end do
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine
!============================================================
subroutine test_scalar_tensor_commutator(jbas,h1,h2,rank,dpar,AAX,BBX) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2
  type(sq_op),optional :: AAX,BBX 
  type(pandya_mat) :: BBCC,WCC
  type(cc_mat) :: AACC
  integer :: a,b,c,d,ja,jb,jc,jd,j1min,j1max
  integer :: j2min,j2max,PAR,TZ,J1,J2,dpar,iii
  integer,intent(in) :: h1,h2,rank
  real(8) :: val,t1,t2,t3,t4,t5,t6,t7,omp_get_wtime
  real(8) :: vv,xx,yy,zz
  
  
  if (present(AAX)) then 
     call duplicate_sq_op(AAX,AA)
     call duplicate_sq_op(BBX,BB)
     call copy_sq_op(AAX,AA)
     call copy_sq_op(BBX,BB) 
  else
     call seed_random_number
  
     BB%rank = rank
     BB%dpar = dpar
     AA%rank = 0
     BB%pphh_ph = .false.
 
     call allocate_blocks(jbas,AA)
     call allocate_tensor(jbas,BB,AA)
     call construct_random_rank0(AA,h1,jbas) 
     call construct_random_rankX(BB,h2,jbas) 
  end if 
  
  call duplicate_sq_op(BB,OUT)
  call duplicate_sq_op(BB,w1) !workspace
  call duplicate_sq_op(BB,w2) !workspace
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
!  call init_ph_mat(BB,BBCC,jbas) !cross coupled ME
  CALL allocate_small_tensor_CCMAT(BB,BBCC,jbas) 

  OUT%herm = -1* AA%herm * BB%herm 
  
  print*, 'TESTING SCALAR-TENSOR COMMUTATORS rank:' ,BB%rank,'parity:',BB%dpar
!  t1 = OMP_get_wtime()
 ! call calculate_generalized_pandya(BB,BBCC,jbas)
!  t2 = OMP_get_wtime()
  call calculate_cross_coupled(AA,AACC,jbas) 
!  t3 = OMP_get_wtime() 
  call TS_commutator_111(AA,BB,OUT,jbas) 
  call TS_commutator_121(AA,BB,OUT,jbas)
  call TS_commutator_211(AACC,BB,OUT,jbas) 
!  t4 = OMP_get_wtime()
  call TS_commutator_122(AA,BB,OUT,jbas)
  call TS_commutator_212(AA,BB,OUT,jbas)
!  t5 = OMP_get_wtime()
  call TS_commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)   
  call TS_commutator_221(w1,w2,AA%herm*BB%herm,OUT,jbas)
!  t6= OMP_get_wtime() 

  call TS_commutator_222_ph(AACC,BBCC,BB,OUT,jbas)
!  t7 = OMP_get_wtime()
  
  print*, 'TIMES:'
  write(*,'(6(f14.7))') t2-t1,t3-t1,t4-t1,t5-t1,t6-t1,t7-t1
!  'porkchop'
!goto 12
do a =  1, jbas%total_orbits
    do b = 1, jbas%total_orbits
    !   do iii = 1, 50   
     !     call random_number(vv)
      !    call random_number(yy)
          
!          a = ceiling(vv*(AA%Nsp))
 !         b = ceiling(yy*(AA%Nsp))
          
          val = scalar_tensor_1body_comm(AA,BB,a,b,jbas) 
        
          if (abs(val-f_tensor_elem(a,b,OUT,jbas)) > 1e-10) then
             print*, 'at: ',a,b
             print*, val, f_tensor_elem(a,b,OUT,jbas)
             STOP 'ONE BODY FAILURE'  
          end if
          
          print*, 'success:', a,b
       end do
    end do
       
 !do a = 12, jbas%total_orbits
     
  iii = 0 
  do while (iii < 55)  
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)
        
     a = ceiling(vv*AA%Nsp)
     b = ceiling(xx*AA%Nsp)
     c = ceiling(yy*AA%Nsp)
     d = ceiling(zz*AA%Nsp)
     
     ja = jbas%jj(a) 
!     do b = 7, jbas%total_orbits
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
 !       do c = 1, jbas%total_orbits
           jc = jbas%jj(c)
  !         do d = 1, jbas%total_orbits
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+BB%dpar/2,2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              iii = iii+1 
              j1min = abs(ja-jb) 
              j1max = ja+jb 
              j2min = abs(jc-jd) 
              j2max = jc+jd
              
              do J1 = j1min,j1max,2
                 do J2 = j2min,j2max,2
                    
                    if (.not. (triangle(J1,J2,rank))) cycle
                    
                    val = scalar_tensor_2body_comm(AA,BB,a,b,c,d,J1,J2,jbas)
                                  
                    if (abs(val-tensor_elem(a,b,c,d,J1,J2,OUT,jbas)) > 1e-8) then
                       print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,tensor_elem(a,b,c,d,J1,J2,OUT,jbas)                       
                       STOP 'TWO BODY FAILURE'  
                    end if
                 end do 
              end do
              
              print*, 'success:', a,b,c,d
   !        end do
    !    end do
     !end do
  end do

  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine test_scalar_tensor_commutator
!=========================================================================
!=========================================================================
subroutine test_scalar_iso_commutator(jbas,h1,h2,rank,dpar,dTz) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA
  type(iso_operator) :: BB,OUT
  type(cc_mat) :: AACC
  integer :: a,b,c,d,ja,jb,jc,jd,j1min,j1max,N
  integer :: j2min,j2max,PAR,TZ,J1,J2,dpar,iii,dTz
  integer,intent(in) :: h1,h2,rank
  real(8) :: val,t1,t2,t3,t4,t5,t6,t7,omp_get_wtime
  real(8) :: vv,xx,yy,zz
    
  call seed_random_number
  
  BB%rank = rank
  BB%dpar = dpar
  BB%dTz = dTz
  BB%xindx = 1 
  AA%rank = 0
  N = jbas%total_orbits
  call allocate_blocks(jbas,AA) 
  allocate(jbas%xmap_tensor(1,N*(N+1)/2)) 
  allocate(half6j(1))
  call allocate_isospin_operator(jbas,BB,AA)

  call duplicate_isospin_operator(BB,OUT)
  BB%herm = 1 
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_opX(BB,h2,jbas) 
    
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME

  OUT%herm = -1* AA%herm * BB%herm 
  
  print*, 'TESTING SCALAR-TENSOR COMMUTATORS rank:' ,BB%rank,'parity:',BB%dpar
!  call calculate_cross_coupled(AA,AACC,jbas) 

  ! siete 
  call operator_commutator_111(AA,BB,OUT,jbas) 
  call operator_commutator_121(AA,BB,OUT,jbas)
  call operator_commutator_211(AA,BB,OUT,jbas) 
  
  call operator_commutator_122(AA,BB,OUT,jbas)
  call operator_commutator_212(AA,BB,OUT,jbas)
  
  call operator_commutator_222_pp_hh(AA,BB,OUT,jbas)   

  call operator_commutator_222_ph(AA,BB,OUT,jbas)

!goto 12
  do a =  1, jbas%total_orbits
     do b = 1, jbas%total_orbits
        !   do iii = 1, 50   
     !     call random_number(vv)
      !    call random_number(yy)
          
!          a = ceiling(vv*(AA%Nsp))
 !         b = ceiling(yy*(AA%Nsp))
          
        val = scalar_tensor_iso1body_comm(AA,BB,a,b,jbas) 
        
        if (abs(val-f_iso_op_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_iso_op_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if
        
        print*, 'success:', a,b,val
     end do
  end do
  
! ! 
!  do a = 12, jbas%total_orbits
     
  iii = 0 
  do while (iii < 55)  
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)
        
     a = ceiling(vv*AA%Nsp)
     b = ceiling(xx*AA%Nsp)
     c = ceiling(yy*AA%Nsp)
     d = ceiling(zz*AA%Nsp)
     ! a = 4!ceiling(vv*AA%Nsp)
     ! b = 7!ceiling(xx*AA%Nsp)
     ! c = 4!ceiling(yy*AA%Nsp)
     ! d = 7!ceiling(zz*AA%Nsp)
     
     ja = jbas%jj(a) 
!     do b = 7, jbas%total_orbits
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
 !       do c = 1, jbas%total_orbits
           jc = jbas%jj(c)
  !         do d = 1, jbas%total_orbits
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+BB%dpar/2,2)) cycle 
              if ( TZ - BB%dTz*2 .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              iii = iii+1 
              j1min = abs(ja-jb) 
              j1max = ja+jb 
              j2min = abs(jc-jd) 
              j2max = jc+jd
              
              do J1 = j1min,j1max,2
                 do J2 = j2min,j2max,2
                    
                    if (.not. (triangle(J1,J2,rank))) cycle
                    
                    val = scalar_tensor_iso2body_comm(AA,BB,a,b,c,d,J1,J2,jbas)
                    print*, a,b,c,d, 'J:', J1,J2 ,val,iso_op_elem(a,b,c,d,J1,J2,OUT,jbas)                       
                    if (abs(val-iso_op_elem(a,b,c,d,J1,J2,OUT,jbas)) > 1e-8) then
                       print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,iso_op_elem(a,b,c,d,J1,J2,OUT,jbas)                       
!                       print*, 'FAIL'
                       STOP 'TWO BODY FAILURE'  
                    end if
                 end do 
              end do
              
              print*, 'success:', a,b,c,d
   !        end do
    !    end do
     !end do
  end do

  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine test_scalar_iso_commutator
!============================================================
subroutine test_tensor_product(jbas,h1,h2,rank_a,rank_b,rank_c,dpar_a,dpar_b,dpar_c) 
  use operators
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2,HS
  type(pandya_mat) :: BBCC,WCC
  type(cc_mat) :: AACC
  integer :: a,b,c,d,ja,jb,jc,jd,j1min,j1max,g
  integer :: j2min,j2max,PAR,TZ,J1,J2,iii,q,N
  integer,intent(in) :: h1,h2,rank_a,rank_b,rank_c,dpar_A,dpar_B,dpar_c
  real(8) :: val,t1,t2,t3,t4,t5,t6,t7,omp_get_wtime
  real(8) :: vv,xx,yy,zz,dcgi00
   
!  call seed_random_number
  vv = dcgi00()
  BB%rank = rank_b 
  BB%dpar = dpar_b
  AA%rank = rank_a
  AA%dpar = dpar_a
  OUT%rank = rank_c
  OUT%dpar = dpar_c 
  BB%pphh_ph = .false.
  AA%pphh_ph = .false.
  OUT%pphh_ph = .false.

  HS%rank = 0
  HS%herm = 1

  N = jbas%total_orbits
  allocate(jbas%xmap_tensor(3,N*(N+1)/2)) 
  allocate(half6j(3))
  call allocate_blocks(jbas,HS)
 
  AA%xindx = 1
  call allocate_tensor(jbas,AA,HS)
  BB%xindx = 2
  deallocate(phase_hh,phase_pp)
  call allocate_tensor(jbas,BB,HS)
  OUT%xindx = 3
  deallocate(phase_hh,phase_pp)
  call allocate_tensor(jbas,OUT,HS)
  
  call construct_random_rankX(AA,h1,jbas) 
  call construct_random_rankX(BB,h2,jbas)   
  
  BB%fpp = 0.d0
  BB%fhh = 0.d0

  
  do q = 1, BB%nblocks
     BB%tblck(q)%lam(1) = 1
     do g = 1, 9
        if ( g == 3) cycle
        if ( g == 7) cycle
        BB%tblck(q)%tgam(g)%X = 0.d0 
     end do 
  end do
       
  do q = 1, OUT%nblocks
     OUT%tblck(q)%lam(1) = 1
  end do
  OUT%herm = 1 
  
  print*, 'TESTING TENSOR PRODUCT' 
  t1 = omp_get_wtime()
  call tensor_product(AA,BB,OUT,jbas) 
  t2 = omp_get_wtime()  
  print*, 'TIME:',t2-t1

  
  ! do a =  1, jbas%total_orbits
  !    if (jbas%con(a) == 1) cycle
  !    do b = 1, jbas%total_orbits
  !       if (jbas%con(b) == 0) cycle

  !       val = EOM_tensor_prod_1body(AA,BB,a,b,rank_c,jbas)
  !       if (abs(val-f_tensor_elem(a,b,OUT,jbas)) > 1e-10) then
  !          print*, 'at: ',a,b
  !          print*, val, f_tensor_elem(a,b,OUT,jbas),f_tensor_elem(a,b,OUT,jbas)/val  
  !          STOP 'ONE BODY FAILURE'  
  !       end if

  !       print*, 'success:', a,b , val 
  !    end do
  ! end do
       
!  !do a = 12, jbas%total_orbits
     
  iii = 0 
  do while (iii < 55)  
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)

     a = jbas%parts(ceiling(vv*(AA%Nsp-AA%belowEF)))
     b = jbas%parts(ceiling(xx*(AA%Nsp-AA%belowEF)))
     c = jbas%holes(ceiling(yy*AA%belowEF))
     d = jbas%holes(ceiling(zz*AA%belowEF))
     
     ja = jbas%jj(a) 
!     do b = 7, jbas%total_orbits
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
 !       do c = 1, jbas%total_orbits
           jc = jbas%jj(c)
  !         do d = 1, jbas%total_orbits
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+OUT%dpar/2,2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              iii = iii+1 
              j1min = abs(ja-jb) 
              j1max = ja+jb 
              j2min = abs(jc-jd) 
              j2max = jc+jd
              
              do J1 = j1min,j1max,2
                 do J2 = j2min,j2max,2
                    
                    if (.not. (triangle(J1,J2,rank_c))) cycle
                    print*, a,b,c,d,J1,J2
                    val = EOM_tensor_prod_2body(AA,BB,a,b,c,d,J1,J2,rank_c,jbas)
                    print*, 'at:', a, b, c, d, 'J:', J1, J2, val, tensor_elem(a,b,c,d,J1,J2,OUT,jbas)                       
                    if (abs(val-tensor_elem(a,b,c,d,J1,J2,OUT,jbas)) > 1e-8) then
                       !print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,tensor_elem(a,b,c,d,J1,J2,OUT,jbas)                       
                       STOP 'TWO BODY FAILURE'  
                    end if
                 end do 
              end do
              
              print*, 'success:', a,b,c,d
   !        end do
    !    end do
     !end do
  end do

  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine test_tensor_product
!============================================================
!============================================================
subroutine test_tensor_dTZ_product(jbas,h1,h2,rank_a,rank_b,rank_c,dpar_a,dpar_b,dpar_c,DtZ) 
  use operators
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,HS
  type(iso_ladder) :: BB,OUT
!  type(pandya_mat) :: BBCC,WCC
 ! type(cc_mat) :: AACC
  integer :: a,b,c,d,ja,jb,jc,jd,j1min,j1max,g
  integer :: j2min,j2max,PAR,TZ,J1,J2,iii,q,N,dtz
  integer,intent(in) :: h1,h2,rank_a,rank_b,rank_c,dpar_A,dpar_B,dpar_c
  real(8) :: val,t1,t2,t3,t4,t5,t6,t7,omp_get_wtime
  real(8) :: vv,xx,yy,zz,dcgi00
   
  call seed_random_number
  vv = dcgi00()
  BB%rank = rank_b 
  BB%dpar = dpar_b
  AA%rank = rank_a
  AA%dpar = dpar_a
  OUT%rank = rank_c
  OUT%dpar = dpar_c 
  BB%pphh_ph = .false.
  AA%pphh_ph = .false.
  OUT%pphh_ph = .false.

  BB%Dtz = dTz
  OUT%dTz = dTz 
  HS%rank = 0
  HS%herm = 1

  N = jbas%total_orbits
  allocate(jbas%xmap_tensor(3,N*(N+1)/2)) 
  allocate(half6j(3))
  call allocate_blocks(jbas,HS)
 
  AA%xindx = 1

  call allocate_tensor(jbas,AA,HS)
  BB%xindx = 2
  call allocate_isospin_ladder(jbas,BB,HS)
  OUT%xindx = 3
  call allocate_isospin_ladder(jbas,OUT,HS)
  
  call construct_random_rankX(AA,h1,jbas) 
  call construct_random_isoX(BB,h2,jbas)   
         
  OUT%herm = 1 

  print*, 'TESTING TENSOR PRODUCT' 
  t1 = omp_get_wtime()
  call tensor_dTz_tensor_product(AA,BB,OUT,jbas) 
  t2 = omp_get_wtime()  
  print*, 'TIME:',t2-t1

  
  ! do a =  1, jbas%total_orbits
  !    if (jbas%con(a) == 1) cycle
  !    do b = 1, jbas%total_orbits
  !       if (jbas%con(b) == 0) cycle

  !       val = EOM_dTz_tensor_prod_1body(AA,BB,a,b,rank_c,jbas)
  !       if (abs(val-f_iso_ladder_elem(a,b,OUT,jbas)) > 1e-10) then
  !          print*, 'at: ',a,b
  !          print*, val, f_iso_ladder_elem(a,b,OUT,jbas),f_iso_ladder_elem(a,b,OUT,jbas)/val  
  !          STOP 'ONE BODY FAILURE'  
  !       end if

  !       print*, 'success:', a,b , val 
  !    end do
  ! end do

! TIT
  iii = 0 
  do while (iii < 55)  
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)

     a = jbas%parts(ceiling(vv*(AA%Nsp-AA%belowEF)))
     b = jbas%parts(ceiling(xx*(AA%Nsp-AA%belowEF)))
     c = jbas%holes(ceiling(yy*AA%belowEF))
     d = jbas%holes(ceiling(zz*AA%belowEF))

     ja = jbas%jj(a) 
     jb = jbas%jj(b)

     PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
     TZ = jbas%itzp(a) + jbas%itzp(b) 

     jc = jbas%jj(c)
     jd = jbas%jj(d) 

     if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+OUT%dpar/2,2)) cycle 
     if ( TZ - 2*DTz  .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
     iii = iii+1 
     j1min = abs(ja-jb) 
     j1max = ja+jb 
     j2min = abs(jc-jd) 
     j2max = jc+jd

     do J1 = j1min,j1max,2
        do J2 = j2min,j2max,2
           print*, J1,J2,rank_c
           if (.not. (triangle(J1,J2,rank_c))) cycle
           val = EOM_dTz_tensor_prod_2body(AA,BB,a,b,c,d,J1,J2,rank_c,jbas)

           print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,iso_ladder_elem(a,b,c,d,J1,J2,OUT,jbas)                       
           if (abs(val-iso_ladder_elem(a,b,c,d,J1,J2,OUT,jbas)) > 1e-8) then
              print*, iso_ladder_elem(a,b,c,d,J1,J2,OUT,jbas)/val
              !print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,iso_ladder_elem(a,b,c,d,J1,J2,OUT,jbas)                       
              STOP 'TWO BODY FAILURE'  
           end if
        end do
     end do

     print*, 'success:', a,b,c,d
  end do

  print*, ' TENSOR PRODUCT EXPRESSIONS CONFIRMED '
  
end subroutine test_tensor_dTZ_product
!============================================================
!============================================================
subroutine test_EOM_scalar_tensor_commutator(jbas,h1,h2,rank,dpar) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2
  type(ex_pandya_mat) :: BBCC,WCC
  type(ex_cc_mat) :: AACC 
  integer :: a,b,c,d,g,q,ja,jb,jc,jd,j1min,j1max,dpar
  integer :: j2min,j2max,PAR,TZ,J1,J2,ax,bx,cx,dx,iii
  integer,intent(in) :: h1,h2,rank
  real(8) :: val,t1,t2,t3,t4,omp_get_wtime
  real(8) :: vv,xx,yy,zz
  
!  call seed_random_number
  
  BB%rank = rank
  BB%dpar = dpar
  BB%pphh_ph = .false.
  AA%rank = 0
  call allocate_blocks(jbas,AA)
  call allocate_tensor(jbas,BB,AA)

  BB%herm = 1 
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_rankX(BB,h2,jbas) 
 
  call duplicate_sq_op(BB,OUT)
  call duplicate_sq_op(BB,w1) !workspace
  call duplicate_sq_op(BB,w2) !workspace
 
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
  call init_ph_mat(BB,BBCC,jbas) !cross coupled ME
  call init_ph_wkspc(BBCC,WCC) !
  

  do q = 1, BB%nblocks
     BB%tblck(q)%lam(1) = 1
     OUT%tblck(q)%lam(1) = 1
     do g = 1, 9
        if ( g == 3) cycle
        if ( g == 7) cycle
        BB%tblck(q)%tgam(g)%X = 0.d0 
     end do 
  end do 
  
  OUT%herm = 1
  
  print*, 'TESTING EOM SCALAR-TENSOR COMMUTATORS' 
 ! t1 = OMP_get_wtime()
  call EOM_generalized_pandya(BB,BBCC,jbas)
 ! t2 = OMP_get_wtime()
  call calculate_cross_coupled_pphh(AA,AACC,jbas) 
  
  call EOM_TS_commutator_111(AA,BB,OUT,jbas) 
  call EOM_TS_commutator_121(AA,BB,OUT,jbas)
  call EOM_TS_commutator_211(AACC,BB,OUT,jbas) 
  call EOM_TS_commutator_122(AA,BB,OUT,jbas)
  call EOM_TS_commutator_212(AA,BB,OUT,jbas)
  
  call EOM_TS_commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)
  
  call EOM_TS_commutator_221(w1,w2,AA%herm*BB%herm,OUT,jbas)
  ! t4 = OMP_get_wtime()
  call EOM_TS_commutator_222_ph(AACC,BBCC,OUT,WCC,jbas)
  ! t3 = OMP_get_wtime()
  
  print*, 'time:', t3-t1,t2-t1,t3-t4
!goto 12
 ! do ax = 1, AA%Nsp-AA%belowEF
   
  do iii = 1, 50   
     call random_number(vv)
     call random_number(yy)
   
     ax = ceiling(vv*(AA%Nsp-AA%belowEF))
     bx = ceiling(yy*(AA%belowEF))
     
     a = jbas%parts(ax)
  !   do bx = 1, AA%belowEF
        b = jbas%holes(bx)
        
        val = EOM_scalar_tensor_1body_comm(AA,BB,a,b,jbas) 
        

        if (abs(val-f_tensor_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_tensor_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if 
        
        print*, 'success:', a,b
!     end do 
  end do 

  
  iii = 0 
  do while (iii < 15) 
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)
   
     ax = ceiling(vv*(AA%Nsp-AA%belowEF))
     bx = ceiling(xx*(AA%Nsp-AA%belowEF))
     cx = ceiling(yy*(AA%belowEF))
     dx = ceiling(zz*(AA%belowEF))
     
     ! do ax = 1, AA%Nsp-AA%belowEF
     a = jbas%parts(ax)
     ja = jbas%jj(a) 
  !   do bx = 1, AA%Nsp-AA%belowEF
        b = jbas%parts(bx)
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
   !     do cx = 1,AA%belowEF
           c = jbas%holes(cx)
           jc = jbas%jj(c)
    !       do dx = 1,AA%belowEF
              d = jbas%holes(dx)
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+BB%dpar/2,2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              iii = iii + 1
              
              j1min = abs(ja-jb) 
              j1max = ja+jb 
              j2min = abs(jc-jd) 
              j2max = jc+jd
              
              do J1 = j1min,j1max,2
                 do J2 = j2min,j2max,2
                    
                    if (.not. (triangle(J1,J2,rank))) cycle
                    
                    val = EOM_scalar_tensor_2body_comm(AA,BB,a,b,c,d,J1,J2,jbas)
!                    print*, a,b,c,d, 'J:', J1,J2 ,val,tensor_elem(a,b,c,d,J1,J2,OUT,jbas)  
                    if (abs(val-tensor_elem(a,b,c,d,J1,J2,OUT,jbas)) > 1e-8) then
                       print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,tensor_elem(a,b,c,d,J1,J2,OUT,jbas)  
                       print*, tensor_elem(c,d,a,b,J2,J1,OUT,jbas)  
            
                       STOP 'TWO BODY FAILURE'  
                    end if 
                 end do 
              end do
              
              print*, 'success:', a,b,c,d
        !   end do
       ! end do
     !end do
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine test_EOM_scalar_tensor_commutator
!============================================================
!============================================================
subroutine test_EOM_iso_commutator(jbas,h1,h2,rank,dpar,dTz) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA
  type(iso_ladder) :: BB,OUT
  integer :: a,b,c,d,g,q,ja,jb,jc,jd,j1min,j1max,dpar,dTz
  integer :: j2min,j2max,PAR,TZ,J1,J2,ax,bx,cx,dx,iii,N
  integer,intent(in) :: h1,h2,rank
  real(8) :: val,t1,t2,t3,t4,omp_get_wtime
  real(8) :: vv,xx,yy,zz
  
  call seed_random_number

  N = jbas%total_orbits
  BB%rank = rank
  BB%dpar = dpar
  BB%dTz = dTz
  BB%xindx = 1
  BB%pphh_ph = .false.
  AA%rank = 0
  call allocate_blocks(jbas,AA)
  allocate(jbas%xmap_tensor(1,N*(N+1)/2)) 
  allocate(half6j(1))
  call allocate_isospin_ladder(jbas,BB,AA)
  call duplicate_isospin_ladder(BB,OUT)
  BB%herm = 1 
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_isoX(BB,h2,jbas) 
 
 
  do q = 1, BB%nblocks
     BB%tblck(q)%lam(1) = 1
     OUT%tblck(q)%lam(1) = 1
  end do 
  
  OUT%herm = 1
  
  print*, 'TESTING EOM ISOSPIN-CHANGING-TENSOR COMMUTATORS' 
  
  call EOM_dTZ_commutator_111(AA,BB,OUT,jbas) 
  call EOM_dTZ_commutator_121(AA,BB,OUT,jbas)
  call EOM_dTZ_commutator_211(AA,BB,OUT,jbas) 
  call EOM_dTZ_commutator_122(AA,BB,OUT,jbas)
  call EOM_dTZ_commutator_212(AA,BB,OUT,jbas)
  
  call EOM_dTZ_commutator_222_pp_hh(AA,BB,OUT,jbas)
  
  call EOM_dTZ_commutator_222_ph(AA,BB,OUT,jbas)

  

   
  do iii = 1, 50   
     call random_number(vv)
     call random_number(yy)
   
     ax = ceiling(vv*(AA%Nsp-AA%belowEF))
     bx = ceiling(yy*(AA%belowEF))
     
     a = jbas%parts(ax)
     
     
     b = jbas%holes(bx)

     val = EOM_scalar_tensor_iso1body_comm(AA,BB,a,b,jbas) 


     if (abs(val-f_iso_ladder_elem(a,b,OUT,jbas)) > 1e-10) then
        print*, 'at: ',a,b
        print*, val, f_iso_ladder_elem(a,b,OUT,jbas)
        print*, 'fail', f_iso_ladder_elem(a,b,OUT,jbas)/val
        !STOP 'ONE BODY FAILURE'  
     end if

     print*, 'success:', a,b,val     
  end do 
  
  
  iii = 0 
  do while (iii < 50) 
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)

     ax = ceiling(vv*(AA%Nsp-AA%belowEF))
     bx = ceiling(xx*(AA%Nsp-AA%belowEF))
     cx = ceiling(yy*(AA%belowEF))
     dx = ceiling(zz*(AA%belowEF))


     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     b = jbas%parts(bx)
     jb = jbas%jj(b)

     PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
     TZ = jbas%itzp(a) + jbas%itzp(b) 


     c = jbas%holes(cx)
     jc = jbas%jj(c)

     d = jbas%holes(dx)
     jd = jbas%jj(d) 

     if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+BB%dpar/2,2)) cycle 
     if ( TZ - BB%dTz*2 .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
     iii = iii + 1

     j1min = abs(ja-jb) 
     j1max = ja+jb 
     j2min = abs(jc-jd) 
     j2max = jc+jd

     do J1 = j1min,j1max,2
        do J2 = j2min,j2max,2

           if (.not. (triangle(J1,J2,rank))) cycle

           val = EOM_scalar_tensor_iso2body_comm(AA,BB,a,b,c,d,J1,J2,jbas)
           if (abs(val-iso_ladder_elem(a,b,c,d,J1,J2,OUT,jbas)) > 1e-8) then
              print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,iso_ladder_elem(a,b,c,d,J1,J2,OUT,jbas)                            
              STOP 'TWO BODY FAILURE'  
           end if
        end do
     end do

     print*, 'success:', a,b,c,d
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine test_EOM_iso_commutator
!============================================================
!============================================================
real(8) function scalar_scalar_1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^0]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k
  integer :: ja,jb,jj,ji,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm 
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * f_elem(i,b,BB,jbas) &
          - f_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  
  
  If (ja == jb) then 
     ! this if statement is built into the expression.
     do i = 1, totorb
        do j = 1, totorb
        
           do Jtot = 0,JTM,2
              
              sm = sm + (jbas%con(i) -jbas%con(j) ) * (Jtot+1.d0) /(ja +1.d0) * &
                   ( f_elem(i,j,AA,jbas) * v_elem(j,a,i,b,Jtot,BB,jbas) - &
                   f_elem(i,j,BB,jbas) * v_elem(j,a,i,b,Jtot,AA,jbas) )
           end do
        end do
     end do
  
  
     do i =  1, totorb
        do j =  1, totorb
           do k =  1, totorb
              do Jtot = 0,JTM,2 

                 sm = sm + (jbas%con(i)*jbas%con(j)*(1-jbas%con(k)) + &
                      (1-jbas%con(i))*(1-jbas%con(j))*jbas%con(k)) * (Jtot + 1.d0) &
                      / (ja + 1.d0) * ( v_elem(a,k,i,j,Jtot,AA,jbas)*v_elem(i,j,b,k,Jtot,BB,jbas) &
                      - v_elem(a,k,i,j,Jtot,BB,jbas)*v_elem(i,j,b,k,Jtot,AA,jbas))/2.d0
           
              end do
           end do
        end do
     end do
  end if 
  scalar_scalar_1body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function EOM_scalar_scalar_1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^0]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k
  integer :: ja,jb,jj,ji,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm 
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * ph_elem(i,b,BB,jbas) &
          - ph_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  
  
  If (ja == jb) then 
     ! this if statement is built into the expression.
     do i = 1, totorb
        do j = 1, totorb
        
           do Jtot = 0,JTM,2
              
              sm = sm + (jbas%con(i) -jbas%con(j) ) * (Jtot+1.d0) /(ja +1.d0) * &
                   ( f_elem(i,j,AA,jbas) * pphh_elem(j,a,i,b,Jtot,BB,jbas) - &
                   ph_elem(i,j,BB,jbas) * v_elem(j,a,i,b,Jtot,AA,jbas) )
           end do
        end do
     end do
  
  
     do i =  1, totorb
        do j =  1, totorb
           do k =  1, totorb
              do Jtot = 0,JTM,2 

                 sm = sm + (jbas%con(i)*jbas%con(j)*(1-jbas%con(k)) + &
                      (1-jbas%con(i))*(1-jbas%con(j))*jbas%con(k)) * (Jtot + 1.d0) &
                      / (ja + 1.d0) * ( v_elem(a,k,i,j,Jtot,AA,jbas)*pphh_elem(i,j,b,k,Jtot,BB,jbas) &
                      - pphh_elem(a,k,i,j,Jtot,BB,jbas)*v_elem(i,j,b,k,Jtot,AA,jbas))/2.d0
           
              end do
           end do
        end do
     end do
  end if 
  EOM_scalar_scalar_1body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function scalar_scalar_0body_comm(AA,BB,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,l
  integer :: ja,jb,jj,ji,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm 
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
    
  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        
        sm = sm + (jbas%con(i) - jbas%con(j)) * (ji + 1.d0) * &
             f_elem(i,j,AA,jbas) * f_elem(j,i,BB,jbas) 
     end do 
  end do 

  do i = 1, totorb
     do j = 1, totorb
        do k = 1, totorb
           do l = 1, totorb
              do Jtot = 0,JTM,2
                 
                 sm = sm + (jbas%con(i) * jbas%con(j) * (1-jbas%con(k)) * (1-jbas%con(l)) &
                      - jbas%con(k) * jbas%con(l) * (1-jbas%con(i)) * (1-jbas%con(j)) ) * &
                      (Jtot+1.d0) * v_elem(i,j,k,l,Jtot,AA,jbas)*v_elem(k,l,i,j,Jtot,BB,jbas)/4.d0
              
              end do
           end do 
        end do
     end do
  end do
  scalar_scalar_0body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function EOM_scalar_scalar_0body_comm(AA,BB,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,l
  integer :: ja,jb,jj,ji,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm 
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
    
  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        
        sm = sm + (jbas%con(i) - jbas%con(j)) * (ji + 1.d0) * &
             f_elem(i,j,AA,jbas) * ph_elem(j,i,BB,jbas) 
     end do 
  end do 

  do i = 1, totorb
     do j = 1, totorb
        do k = 1, totorb
           do l = 1, totorb
              do Jtot = 0,JTM,2
                 
                 sm = sm + (jbas%con(i) * jbas%con(j) * (1-jbas%con(k)) * (1-jbas%con(l)) &
                      - jbas%con(k) * jbas%con(l) * (1-jbas%con(i)) * (1-jbas%con(j)) ) * &
                      (Jtot+1.d0) * v_elem(i,j,k,l,Jtot,AA,jbas)*pphh_elem(k,l,i,j,Jtot,BB,jbas)/4.d0
              
              end do
           end do 
        end do
     end do
  end do
  EOM_scalar_scalar_0body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function scalar_scalar_2body_comm(AA,BB,a,b,c,d,Jtot,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm,coef9
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

  do i = 1, totorb
     sm = sm + f_elem(a,i,AA,jbas) * v_elem( i,b,c,d,Jtot,BB,jbas) &
          + f_elem(b,i,AA,jbas) * v_elem( a,i,c,d,Jtot,BB,jbas) &
          - f_elem(i,c,AA,jbas) * v_elem( a,b,i,d,Jtot,BB,jbas) &
          - f_elem(i,d,AA,jbas) * v_elem( a,b,c,i,Jtot,BB,jbas) 
     
     sm = sm - f_elem(a,i,BB,jbas) * v_elem( i,b,c,d,Jtot,AA,jbas) &
          - f_elem(b,i,BB,jbas) * v_elem( a,i,c,d,Jtot,AA,jbas) &
          + f_elem(i,c,BB,jbas) * v_elem( a,b,i,d,Jtot,AA,jbas) &
          + f_elem(i,d,BB,jbas) * v_elem( a,b,c,i,Jtot,AA,jbas) 
  
  end do
  

  do i = 1, totorb
     do j = 1, totorb
        
        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,Jtot,AA,jbas)*v_elem(i,j,c,d,Jtot,BB,jbas)   &
             - v_elem(a,b,i,j,Jtot,BB,jbas)*v_elem(i,j,c,d,Jtot,AA,jbas)) 
     end do
  end do

  do i = 1, totorb
     ji =jbas%jj(i)
     do j = 1,totorb
        jj = jbas%jj(j) 
        
        if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
        
        do J1 = 0, JTM,2
           do J2 = 0, JTM,2 
              
              sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                   (-1)** ((J1+J2 + ja-jc)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,ja,J2,ji,jb,jc,jd,Jtot) * v_elem(j,a,i,d,J1,AA,jbas) &
                   * v_elem(i,b,j,c,J2,BB,jbas) &
                   
                   - (-1)** ((J1+J2 + jb-jc)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,jb,J2,ji,ja,jc,jd,Jtot) * v_elem(j,b,i,d,J1,AA,jbas) &
                   * v_elem(i,a,j,c,J2,BB,jbas) *(-1)**((ja+jb-Jtot)/2) &
                   
                   - (-1)** ((J1+J2 + ja-jd)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,ja,J2,ji,jb,jd,jc,Jtot) * v_elem(j,a,i,c,J1,AA,jbas) &
                   * v_elem(i,b,j,d,J2,BB,jbas) * (-1)**((jc+jd-Jtot)/2) &
                   
                    + (-1)** ((J1+J2 + jb-jd)/2) * (J1+1.d0) * (J2+1.d0) &
                    * coef9(jj,J1,jb,J2,ji,ja,jd,jc,Jtot) * v_elem(j,b,i,c,J1,AA,jbas) &
                   * v_elem(i,a,j,d,J2,BB,jbas) *(-1)**((ja+jb+jc+jd)/2)  )
           end do
        end do
        
       
        
     end do
  end do
  
  
  ! do J2= 0,4,2
  !    do i = 1, totorb
  !       ji =jbas%jj(i)
  !       do j = 1,totorb
  !          jj = jbas%jj(j) 
           
  !          if ((jbas%con(i)-jbas%con(j)) == 0) cycle 

  !          sm = sm + (jbas%con(i)-jbas%con(j)) * (J2+1.d0)* ( & 
  !               sixj(ja,jb,Jtot,jc,jd,J2) * &
  !                ( Vpandya(a,d,i,j,J2,AA,jbas) * Vpandya(i,j,c,b,J2,BB,jbas) - &
  !               Vpandya(a,d,i,j,J2,BB,jbas) * Vpandya(i,j,c,b,J2,AA,jbas)  ) - &
  !               (-1)**((ja+jb-Jtot)/2) * sixj(jb,ja,Jtot,jc,jd,J2) * &
  !               ( Vpandya(b,d,i,j,J2,AA,jbas) * Vpandya(i,j,c,a,J2,BB,jbas) - &
  !               Vpandya(b,d,i,j,J2,BB,jbas) * Vpandya(i,j,c,a,J2,AA,jbas)  ) )
  !       end do
  !    end do
     
  ! end do
  
  scalar_scalar_2body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function EOM_scalar_scalar_2body_comm(AA,BB,a,b,c,d,Jtot,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm,coef9
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

  do i = 1, totorb
     sm = sm + f_elem(a,i,AA,jbas) * pphh_elem( i,b,c,d,Jtot,BB,jbas) &
          + f_elem(b,i,AA,jbas) * pphh_elem( a,i,c,d,Jtot,BB,jbas) &
          - f_elem(i,c,AA,jbas) * pphh_elem( a,b,i,d,Jtot,BB,jbas) &
          - f_elem(i,d,AA,jbas) * pphh_elem( a,b,c,i,Jtot,BB,jbas) 
     
     sm = sm - ph_elem(a,i,BB,jbas) * v_elem( i,b,c,d,Jtot,AA,jbas) &
          - ph_elem(b,i,BB,jbas) * v_elem( a,i,c,d,Jtot,AA,jbas) &
          + ph_elem(i,c,BB,jbas) * v_elem( a,b,i,d,Jtot,AA,jbas) &
          + ph_elem(i,d,BB,jbas) * v_elem( a,b,c,i,Jtot,AA,jbas) 
  
  end do
  

  do i = 1, totorb
     do j = 1, totorb
        
        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,Jtot,AA,jbas)*pphh_elem(i,j,c,d,Jtot,BB,jbas)   &
             - pphh_elem(a,b,i,j,Jtot,BB,jbas)*v_elem(i,j,c,d,Jtot,AA,jbas)) 
     end do
  end do

  do i = 1, totorb
     ji =jbas%jj(i)
     do j = 1,totorb
        jj = jbas%jj(j) 
        
        if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
        do J1 = 0, JTM,2
           do J2 = 0, JTM,2 
              
              sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                   (-1)** ((J1+J2 + ja-jc)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,ja,J2,ji,jb,jc,jd,Jtot) * v_elem(j,a,i,d,J1,AA,jbas) &
                   * pphh_elem(i,b,j,c,J2,BB,jbas) &
                   
                   - (-1)** ((J1+J2 + jb-jc)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,jb,J2,ji,ja,jc,jd,Jtot) * v_elem(j,b,i,d,J1,AA,jbas) &
                   * pphh_elem(i,a,j,c,J2,BB,jbas) *(-1)**((ja+jb-Jtot)/2) &
                   
                   - (-1)** ((J1+J2 + ja-jd)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,ja,J2,ji,jb,jd,jc,Jtot) * v_elem(j,a,i,c,J1,AA,jbas) &
                   * pphh_elem(i,b,j,d,J2,BB,jbas) * (-1)**((jc+jd-Jtot)/2) &
                   
                   + (-1)** ((J1+J2 + jb-jd)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,jb,J2,ji,ja,jd,jc,Jtot) * v_elem(j,b,i,c,J1,AA,jbas) &
                   * pphh_elem(i,a,j,d,J2,BB,jbas) *(-1)**((ja+jb+jc+jd)/2)  )
              
           end do
        end do
     end do
  end do
  
  EOM_scalar_scalar_2body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function scalar_tensor_1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm ,d6ji
  
  rank = BB%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * f_tensor_elem(i,b,BB,jbas) &
          - f_tensor_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  

  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        jj = jbas%jj(j) 
        
        do J1 = 0,JTM,2
           do J2 = 0, JTM,2 
                 
              sm = sm + (jbas%con(i) -jbas%con(j) )* &
                   f_elem(i,j,AA,jbas) * tensor_elem(j,a,i,b,J1,J2,BB,jbas) &
                   * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                   d6ji(J1,J2,rank,jb,ja,ji) 
           
           end do
        end do
     end do
  end do
  
  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        jj = jbas%jj(j) 
        do k = 1, totorb
           jk =jbas%jj(k)
           
           do J1 = 0,JTM,2
              do J2 = 0, JTM,2 
                 
                 sm = sm + 0.5*(jbas%con(k)*jbas%con(j)*(1-jbas%con(i)) + &
                      (1-jbas%con(k))*(1-jbas%con(j))*jbas%con(i) )* &
                      ( v_elem(i,a,j,k,J1,AA,jbas) * tensor_elem(j,k,i,b,J1,J2,BB,jbas) &
                      - tensor_elem(i,a,j,k,J1,J2,BB,jbas) * v_elem(j,k,i,b,J2,AA,jbas) ) &
                      * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                      d6ji(J1,J2,rank,jb,ja,ji) 
           
              end do
           end do
        end do
     end do
  end do
  

  do i =  1, totorb
     ji =jbas%jj(i)
     do j =  1, totorb
        jj = jbas%jj(j) 
        
        do Jtot = 0,JTM,2 

           sm = sm - (jbas%con(i) - jbas%con(j)) * (Jtot+1.d0) * &
                (-1) ** ((ja+jj+Jtot)/2) *d6ji(ji,jj,rank,ja,jb,Jtot) * &
                f_tensor_elem(i,j,BB,jbas) * v_elem(j,a,i,b,Jtot,AA,jbas)
         !  sm = sm + (jbas%con(i)-jbas%con(j)) * (-1) ** ((ja+jb+ji+jj)/2)&
          !      *f_tensor_elem(j,i,BB,jbas)*vpandya(i,j,b,a,rank,AA,jbas) 
           
        end do
     end do
  end do

  scalar_tensor_1body_comm = sm 
  
end function scalar_tensor_1body_comm
!============================================================
!============================================================
real(8) function scalar_tensor_iso1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA
  type(iso_operator) :: BB 
  real(8) :: sm ,d6ji,smx,matml
  
  rank = BB%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
!  print*
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * f_iso_op_elem(i,b,BB,jbas) &
          - f_iso_op_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  
  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        jj = jbas%jj(j) 
        
        do J1 = 0,JTM,2
           do J2 = 0, JTM,2 
                 
              sm = sm + (jbas%con(i) -jbas%con(j) )* &
                   f_elem(i,j,AA,jbas) * iso_op_elem(j,a,i,b,J1,J2,BB,jbas) &
                   * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                   d6ji(J1,J2,rank,jb,ja,ji)
              
           end do
        end do
     end do
  end do

  do J1 = 0,JTM,2
     do J2 = 0, JTM,2 
        
        do i = 1, totorb
           ji = jbas%jj(i)

           matml = 0.d0 
           do j = 1, totorb
              jj = jbas%jj(j) 
              do k = 1, totorb
                 jk =jbas%jj(k)
                 
                 matml = matml + 0.5*(jbas%con(k)*jbas%con(j)*(1-jbas%con(i)) + &
                      (1-jbas%con(k))*(1-jbas%con(j))*jbas%con(i) )* &
                      ( v_elem(i,a,j,k,J1,AA,jbas) * iso_op_elem(j,k,i,b,J1,J2,BB,jbas) &
                      - iso_op_elem(i,a,j,k,J1,J2,BB,jbas) * v_elem(j,k,i,b,J2,AA,jbas) )
 
              end do
           end do

           smx = matml * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                      d6ji(J1,J2,rank,jb,ja,ji) 

           sm = sm + smx

        end do
     end do
  end do
  


  do i =  1, totorb
     ji =jbas%jj(i)
     do j =  1, totorb
        jj = jbas%jj(j) 
        
        do Jtot = 0,JTM,2 

           sm = sm - (jbas%con(i) - jbas%con(j)) * (Jtot+1.d0) * &
                (-1) ** ((ja+jj+Jtot)/2) *d6ji(ji,jj,rank,ja,jb,Jtot) * &
                f_iso_op_elem(i,j,BB,jbas) * v_elem(j,a,i,b,Jtot,AA,jbas)
         !  sm = sm + (jbas%con(i)-jbas%con(j)) * (-1) ** ((ja+jb+ji+jj)/2)&
          !      *f_iso_op_elem(j,i,BB,jbas)*vpandya(i,j,b,a,rank,AA,jbas) 
           
        end do
     end do
  end do

  scalar_tensor_iso1body_comm = sm 
  
end function scalar_tensor_iso1body_comm
!============================================================
!============================================================
real(8) function EOM_scalar_tensor_1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm ,d6ji,sx
  
  rank = BB%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * ph_tensor_elem(i,b,BB,jbas) &
          - ph_tensor_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  

  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        jj = jbas%jj(j) 
        
        do J1 = 0,JTM,2
           do J2 = 0, JTM,2 
                 
              sm = sm + (jbas%con(i) -jbas%con(j) )* &
                   f_elem(i,j,AA,jbas) * pphh_tensor_elem(j,a,i,b,J1,J2,BB,jbas) &
                   * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                   d6ji(J1,J2,rank,jb,ja,ji) 
           
           end do
        end do
     end do
  end do
  
  do i = 1, totorb
     ji = jbas%jj(i)      
     do j = 1, totorb
        jj = jbas%jj(j) 
        do k = 1, totorb
           jk =jbas%jj(k)
           
           do J1 = 0,JTM,2
              do J2 = 0,JTM,2
           
                 sm = sm + 0.5*(jbas%con(k)*jbas%con(j)*(1-jbas%con(i)) + &
                      (1-jbas%con(k))*(1-jbas%con(j))*jbas%con(i) )* &
                      ( v_elem(i,a,j,k,J1,AA,jbas) * pphh_tensor_elem(j,k,i,b,J1,J2,BB,jbas) &
                      - pphh_tensor_elem(i,a,j,k,J1,J2,BB,jbas) * v_elem(j,k,i,b,J2,AA,jbas) ) &
                      * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                      d6ji(J1,J2,rank,jb,ja,ji) 
           
              end do
           end do
        end do
     end do
  end do
  

  do i =  1, totorb
     ji =jbas%jj(i)
     do j =  1, totorb
        jj = jbas%jj(j) 
        
        do Jtot = 0,JTM,2 

           sm = sm - (jbas%con(i) - jbas%con(j)) * (Jtot+1.d0) * &
                (-1) ** ((ja-ji+Jtot)/2) *d6ji(ja,jb,rank,jj,ji,Jtot) * &
                ph_tensor_elem(j,i,BB,jbas) * v_elem(i,a,j,b,Jtot,AA,jbas)
           

           
        end do
     end do
  end do

  
  EOM_scalar_tensor_1body_comm = sm 
  
end function EOM_scalar_tensor_1body_comm
!==================================================================
!==================================================================
real(8) function EOM_scalar_tensor_iso1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA
  type(iso_ladder) :: BB 
  real(8) :: sm ,d6ji,sx
  
  rank = BB%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * f_iso_ladder_elem(i,b,BB,jbas) &
          - f_iso_ladder_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  

  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        jj = jbas%jj(j) 
        
        do J1 = 0,JTM,2
           do J2 = 0, JTM,2 
                 
              sm = sm + (jbas%con(i) -jbas%con(j) )* &
                   f_elem(i,j,AA,jbas) * iso_ladder_elem(j,a,i,b,J1,J2,BB,jbas) &
                   * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                   d6ji(J1,J2,rank,jb,ja,ji) 
           
           end do
        end do
     end do
  end do
  
  do i = 1, totorb
     ji = jbas%jj(i)      
     do j = 1, totorb
        jj = jbas%jj(j) 
        do k = 1, totorb
           jk =jbas%jj(k)
           
           do J1 = 0,JTM,2
              do J2 = 0,JTM,2
           
                 sm = sm + 0.5*(jbas%con(k)*jbas%con(j)*(1-jbas%con(i)) + &
                      (1-jbas%con(k))*(1-jbas%con(j))*jbas%con(i) )* &
                      ( v_elem(i,a,j,k,J1,AA,jbas) * iso_ladder_elem(j,k,i,b,J1,J2,BB,jbas) &
                      - iso_ladder_elem(i,a,j,k,J1,J2,BB,jbas) * v_elem(j,k,i,b,J2,AA,jbas) ) &
                      * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                      d6ji(J1,J2,rank,jb,ja,ji) 
           
              end do
           end do
        end do
     end do
  end do
  

  do i =  1, totorb
     ji =jbas%jj(i)
     do j =  1, totorb
        jj = jbas%jj(j) 
        
        do Jtot = 0,JTM,2 

           sm = sm - (jbas%con(i) - jbas%con(j)) * (Jtot+1.d0) * &
                (-1) ** ((ja-ji+Jtot)/2) *d6ji(ja,jb,rank,jj,ji,Jtot) * &
                f_iso_ladder_elem(j,i,BB,jbas) * v_elem(i,a,j,b,Jtot,AA,jbas)
           

           
        end do
     end do
  end do

  
  EOM_scalar_tensor_iso1body_comm = sm 
  
end function EOM_scalar_tensor_iso1body_comm
!==================================================================
!==================================================================
real(8) function EOM_tensor_prod_1body(AA,BB,p,h,rank_c,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2,ax,ix,p,h,mp,ma,mh,mi,rank_c
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb,jp,jh,jm,holes,parts
  integer :: rank_a,rank_b,mu_a,mu_b,mu_c,bx,jx,mb,mj
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm ,d6ji,sx,sm1,sm2 ,dcgi

  rank_b = BB%rank
  rank_a = AA%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  holes = sum(jbas%con)
  parts = totorb- holes

  jp = jbas%jj(p) 
  jh = jbas%jj(h) 

  do mh = -1*jh,jh,2
     do mp = -1*jp,jp,2
        do mu_c = -1*rank_c,rank_c,2

           do mu_a = -1*rank_a,rank_a,2
              do mu_b = -1*rank_b,rank_b,2
                 sm1 = 0.d0

                 do ax = 1, parts
                    a = jbas%parts(ax)
                    ja = jbas%jj(a) 
                    do ma = -1*ja,ja,2
                       
                       sm1 = sm1 + f_tensor_mscheme(p,mp,a,ma,mu_a,AA,jbas) * &
                            f_tensor_mscheme(a,ma,h,mh,mu_b,BB,jbas)
                    end do
                    
                 end do
                 
                 do ix = 1, holes
                    i = jbas%holes(ix)
                    ji = jbas%jj(i)
                    
                    do mi = -1*ji,ji,2
                       sm1 = sm1 - f_tensor_mscheme(p,mp,i,mi,mu_b,BB,jbas) * &
                            f_tensor_mscheme(i,mi,h,mh,mu_a,AA,jbas)  
                
                    end do
                 end do

                 do ax = 1, parts
                    a = jbas%parts(ax)
                    ja = jbas%jj(a) 
                    do ma = -1*ja,ja,2

                       do ix = 1, holes
                          i = jbas%holes(ix)
                          ji = jbas%jj(i)
                          
                          do mi = -1*ji,ji,2

                             sm1 = sm1 + f_Tensor_mscheme(a,ma,i,mi,mu_b,BB,jbas)*&
                                  tensor_mscheme(p,mp,i,mi,h,mh,a,ma,mu_a,AA,jbas)

                             sm1 = sm1 + f_Tensor_mscheme(i,mi,a,ma,mu_a,AA,jbas)*&
                                  tensor_mscheme(p,mp,a,ma,h,mh,i,mi,mu_b,BB,jbas)
                          end do
                       end do
                    end do
                 end do

                 do ax = 1, parts
                    a = jbas%parts(ax)
                    ja = jbas%jj(a) 
                    do ma = -1*ja,ja,2

                       do bx = 1, parts
                          b = jbas%parts(bx)
                          jb = jbas%jj(b) 
                          do mb = -1*jb,jb,2

                             do ix = 1, holes
                                i = jbas%holes(ix)
                                ji = jbas%jj(i)
                          
                                do mi = -1*ji,ji,2

                                   sm1 = sm1 + tensor_mscheme(i,mi,p,mp,a,ma,b,mb,mu_a,AA,jbas)*&
                                        tensor_mscheme(a,ma,b,mb,i,mi,h,mh,mu_b,BB,jbas) *0.5

                                end do
                             end do
                          end do
                       end do
                    end do
                 end do
                             
                 do ix = 1, holes
                    i = jbas%holes(ix)
                    ji = jbas%jj(i)

                    do mi = -1*ji,ji,2

                       do jx = 1, holes
                          j = jbas%holes(jx)
                          jj = jbas%jj(j)

                          do mj = -1*jj,jj,2

                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a) 
                                do ma = -1*ja,ja,2

                                   sm1 = sm1 - tensor_mscheme(a,ma,p,mp,i,mi,j,mj,mu_b,BB,jbas)*&
                                        tensor_mscheme(i,mi,j,mj,a,ma,h,mh,mu_a,AA,jbas) *0.5

                                end do
                             end do
                          end do
                       end do
                    end do
                 end do
                             
                                   
                 sm = sm + sm1 * dcgi(jh,mh,rank_c,mu_c,jp,mp) &
                            *dcgi(rank_a,mu_a,rank_b,mu_b,rank_c,mu_c)/sqrt(jp+1.d0)  
              end do
           end do
        end do
     end do
  end do
     
     
  
  EOM_tensor_prod_1body = sm 
  
end function EOM_tensor_prod_1body
!==================================================================
!==================================================================
real(8) function EOM_dTz_tensor_prod_1body(AA,BB,p,h,rank_c,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2,ax,ix,p,h,mp,ma,mh,mi,rank_c
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb,jp,jh,jm,holes,parts
  integer :: rank_a,rank_b,mu_a,mu_b,mu_c,bx,jx,mb,mj
  type(spd) :: jbas
  type(sq_op) :: AA
  type(iso_ladder) :: BB 
  real(8) :: sm ,d6ji,sx,sm1,sm2 ,dcgi

  rank_b = BB%rank
  rank_a = AA%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  holes = sum(jbas%con)
  parts = totorb- holes

  jp = jbas%jj(p) 
  jh = jbas%jj(h) 

  do mh = -1*jh,jh,2
     do mp = -1*jp,jp,2
        do mu_c = -1*rank_c,rank_c,2

           do mu_a = -1*rank_a,rank_a,2
              do mu_b = -1*rank_b,rank_b,2
                 sm1 = 0.d0

                 do ax = 1, parts
                    a = jbas%parts(ax)
                    ja = jbas%jj(a) 
                    do ma = -1*ja,ja,2
                       
                       sm1 = sm1 + f_tensor_mscheme(p,mp,a,ma,mu_a,AA,jbas) * &
                            f_iso_ladder_mscheme(a,ma,h,mh,mu_b,BB,jbas)
                    end do
                    
                 end do
                 
                 do ix = 1, holes
                    i = jbas%holes(ix)
                    ji = jbas%jj(i)
                    
                    do mi = -1*ji,ji,2
                       sm1 = sm1 - f_iso_ladder_mscheme(p,mp,i,mi,mu_b,BB,jbas) * &
                            f_tensor_mscheme(i,mi,h,mh,mu_a,AA,jbas)  
                
                    end do
                 end do

                 do ax = 1, parts
                    a = jbas%parts(ax)
                    ja = jbas%jj(a) 
                    do ma = -1*ja,ja,2

                       do ix = 1, holes
                          i = jbas%holes(ix)
                          ji = jbas%jj(i)
                          
                          do mi = -1*ji,ji,2

                             sm1 = sm1 + f_Iso_Ladder_mscheme(a,ma,i,mi,mu_b,BB,jbas)*&
                                  tensor_mscheme(p,mp,i,mi,h,mh,a,ma,mu_a,AA,jbas)

                             sm1 = sm1 + f_Tensor_mscheme(i,mi,a,ma,mu_a,AA,jbas)*&
                                  iso_ladder_mscheme(p,mp,a,ma,h,mh,i,mi,mu_b,BB,jbas)
                          end do
                       end do
                    end do
                 end do

                 do ax = 1, parts
                    a = jbas%parts(ax)
                    ja = jbas%jj(a) 
                    do ma = -1*ja,ja,2

                       do bx = 1, parts
                          b = jbas%parts(bx)
                          jb = jbas%jj(b) 
                          do mb = -1*jb,jb,2

                             do ix = 1, holes
                                i = jbas%holes(ix)
                                ji = jbas%jj(i)
                          
                                do mi = -1*ji,ji,2

                                   sm1 = sm1 + tensor_mscheme(i,mi,p,mp,a,ma,b,mb,mu_a,AA,jbas)*&
                                        iso_ladder_mscheme(a,ma,b,mb,i,mi,h,mh,mu_b,BB,jbas) *0.5

                                end do
                             end do
                          end do
                       end do
                    end do
                 end do
                             
                 do ix = 1, holes
                    i = jbas%holes(ix)
                    ji = jbas%jj(i)

                    do mi = -1*ji,ji,2

                       do jx = 1, holes
                          j = jbas%holes(jx)
                          jj = jbas%jj(j)

                          do mj = -1*jj,jj,2

                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a) 
                                do ma = -1*ja,ja,2

                                   sm1 = sm1 - iso_ladder_mscheme(a,ma,p,mp,i,mi,j,mj,mu_b,BB,jbas)*&
                                        tensor_mscheme(i,mi,j,mj,a,ma,h,mh,mu_a,AA,jbas) *0.5

                                end do
                             end do
                          end do
                       end do
                    end do
                 end do
                             
                                   
                 sm = sm + sm1 * dcgi(jh,mh,rank_c,mu_c,jp,mp) &
                            *dcgi(rank_a,mu_a,rank_b,mu_b,rank_c,mu_c)/sqrt(jp+1.d0)  
              end do
           end do
        end do
     end do
  end do
     
     
  
  EOM_dTz_tensor_prod_1body = sm 
  
end function EOM_dTz_tensor_prod_1body
!==================================================================
!==================================================================
real(8) function EOM_tensor_prod_2body(AA,BB,p1,p2,h1,h2,J1,J2,rank_c,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2,ax,ix,p1,p2,h1,h2,mp1,ma,mh1,mi,rank_c
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb,jp1,jh1,jm,holes,parts,mp2,mh2,M1,M2
  integer :: rank_a,rank_b,mu_a,mu_b,mu_c,bx,jx,mb,mj,jp2,jh2
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm ,d6ji,sx,sm1,sm2,sm3 ,dcgi

  rank_b = BB%rank
  rank_a = AA%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  holes = sum(jbas%con)
  parts = sum(1-jbas%con) 

  jp1 = jbas%jj(p1)
  jp2 = jbas%jj(p2) 
  jh1 = jbas%jj(h1) 
  jh2 = jbas%jj(h2)

  !  !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) REDUCTION(+:sm)
  do M1 = -1*J1,J1,2
     do M2 = -1*J2,J2,2

        do mh2 = -1*jh2,jh2,2
           do mh1 = -1*jh1,jh1,2
              do mp1 = -1*jp1,jp1,2
                 do mp2 = -1*jp2,jp2,2

                    
                    do mu_c = -1*rank_c,rank_c,2                        
                       do mu_a = -1*rank_a,rank_a,2
                          do mu_b = -1*rank_b,rank_b,2


                             sm1 = f_tensor_mscheme(p1,mp1,h1,mh1,mu_a,AA,jbas)*&
                                  f_tensor_mscheme(p2,mp2,h2,mh2,mu_b,BB,jbas)
                             
                             sm1 = sm1 - f_tensor_mscheme(p2,mp2,h1,mh1,mu_a,AA,jbas)*&
                                  f_tensor_mscheme(p1,mp1,h2,mh2,mu_b,BB,jbas)
                             
                             sm1 = sm1 + f_tensor_mscheme(p2,mp2,h2,mh2,mu_a,AA,jbas)*&
                                  f_tensor_mscheme(p1,mp1,h1,mh1,mu_b,BB,jbas)
                             
                             sm1 = sm1 - f_tensor_mscheme(p1,mp1,h2,mh2,mu_a,AA,jbas)*&
                                  f_tensor_mscheme(p2,mp2,h1,mh1,mu_b,BB,jbas)
                             
                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a) 
                                do ma = -1*ja,ja,2

                                   sm1 = sm1 + f_tensor_mscheme(p1,mp1,a,ma,mu_a,AA,jbas) * &
                                        tensor_mscheme(a,ma,p2,mp2,h1,mh1,h2,mh2,mu_b,BB,jbas)

                                   sm1 = sm1 - f_tensor_mscheme(p2,mp2,a,ma,mu_a,AA,jbas) * &
                                        tensor_mscheme(a,ma,p1,mp1,h1,mh1,h2,mh2,mu_b,BB,jbas)

                                end do
                             end do

                             do ix = 1, holes
                                i = jbas%holes(ix)
                                ji = jbas%jj(i) 
                                do mi = -1*ji,ji,2

                                   sm1 = sm1 - f_tensor_mscheme(p1,mp1,i,mi,mu_b,BB,jbas) * &
                                        tensor_mscheme(i,mi,p2,mp2,h1,mh1,h2,mh2,mu_a,AA,jbas)

                                   sm1 = sm1 + f_tensor_mscheme(p2,mp2,i,mi,mu_b,BB,jbas) * &
                                        tensor_mscheme(i,mi,p1,mp1,h1,mh1,h2,mh2,mu_a,AA,jbas)

                                end do
                             end do

                             do ix = 1, holes
                                i = jbas%holes(ix)
                                ji = jbas%jj(i)

                                do mi = -1*ji,ji,2
                                   sm1 = sm1 - tensor_mscheme(p1,mp1,p2,mp2,i,mi,h2,mh2,mu_b,BB,jbas) * &
                                        f_tensor_mscheme(i,mi,h1,mh1,mu_a,AA,jbas)  

                                   sm1 = sm1 + tensor_mscheme(p1,mp1,p2,mp2,i,mi,h1,mh1,mu_b,BB,jbas) * &
                                        f_tensor_mscheme(i,mi,h2,mh2,mu_a,AA,jbas)
                                end do
                             end do

                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a)

                                do ma = -1*ja,ja,2
                                   sm1 = sm1 + tensor_mscheme(p1,mp1,p2,mp2,a,ma,h2,mh2,mu_a,AA,jbas) * &
                                        f_tensor_mscheme(a,ma,h1,mh1,mu_b,BB,jbas)  

                                   sm1 = sm1 - tensor_mscheme(p1,mp1,p2,mp2,a,ma,h1,mh1,mu_a,AA,jbas) * &
                                        f_tensor_mscheme(a,ma,h2,mh2,mu_b,BB,jbas)
                                end do
                             end do




                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a) 
                                do ma = -1*ja,ja,2

                                   do bx = 1, parts
                                      b = jbas%parts(bx)
                                      jb = jbas%jj(b) 
                                      do mb = -1*jb,jb,2


                                         sm1 = sm1 + tensor_mscheme(p1,mp1,p2,mp2,a,ma,b,mb,mu_a,AA,jbas)*&
                                              tensor_mscheme(a,ma,b,mb,h1,mh1,h2,mh2,mu_b,BB,jbas) *0.5

                                      end do
                                   end do
                                end do
                             end do

                             do ix = 1, holes
                                i = jbas%holes(ix)
                                ji = jbas%jj(i)

                                do mi = -1*ji,ji,2

                                   do jx = 1, holes
                                      j = jbas%holes(jx)
                                      jj = jbas%jj(j)

                                      do mj = -1*jj,jj,2


                                         sm1 = sm1 + tensor_mscheme(p1,mp1,p2,mp2,i,mi,j,mj,mu_b,BB,jbas)*&
                                              tensor_mscheme(i,mi,j,mj,h1,mh1,h2,mh2,mu_a,AA,jbas) *0.5

                                      end do
                                   end do
                                end do
                             end do

                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a) 
                                do ma = -1*ja,ja,2

                                   do ix = 1, holes
                                      i = jbas%holes(ix)
                                      ji = jbas%jj(i) 
                                      do mi = -1*ji,ji,2


                                         sm1 = sm1 + tensor_mscheme(p1,mp1,i,mi,h1,mh1,a,ma,mu_a,AA,jbas)*&
                                              tensor_mscheme(p2,mp2,a,ma,h2,mh2,i,mi,mu_b,BB,jbas)

                                         sm1 = sm1 - tensor_mscheme(p2,mp2,i,mi,h1,mh1,a,ma,mu_a,AA,jbas)*&
                                              tensor_mscheme(p1,mp1,a,ma,h2,mh2,i,mi,mu_b,BB,jbas)

                                         sm1 = sm1 + tensor_mscheme(p2,mp2,i,mi,h2,mh2,a,ma,mu_a,AA,jbas)*&
                                              tensor_mscheme(p1,mp1,a,ma,h1,mh1,i,mi,mu_b,BB,jbas)

                                         sm1 = sm1 - tensor_mscheme(p1,mp1,i,mi,h2,mh2,a,ma,mu_a,AA,jbas)*&
                                              tensor_mscheme(p2,mp2,a,ma,h1,mh1,i,mi,mu_b,BB,jbas)


                                      end do
                                   end do
                                end do
                             end do


                             sm = sm + sm1 * dcgi(J2,M2,rank_c,mu_c,J1,M1) * dcgi(jp1,mp1,jp2,mp2,J1,M1)&
                                  *dcgi(jh1,mh1,jh2,mh2,J2,M2) *dcgi(rank_a,mu_a,rank_b,mu_b,rank_c,mu_c)/sqrt(J1+1.d0)  

                          end do
                       end do
                    end do
                 end do

              end do
           end do
        end do
     end do
  end do
!  !$OMP END PARALLEL DO
     
  
  EOM_tensor_prod_2body = sm 
  
end function EOM_tensor_prod_2body
!==================================================================
!==================================================================
real(8) function EOM_dTz_tensor_prod_2body(AA,BB,p1,p2,h1,h2,J1,J2,rank_c,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2,ax,ix,p1,p2,h1,h2,mp1,ma,mh1,mi,rank_c
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb,jp1,jh1,jm,holes,parts,mp2,mh2,M1,M2
  integer :: rank_a,rank_b,mu_a,mu_b,mu_c,bx,jx,mb,mj,jp2,jh2
  type(spd) :: jbas
  type(sq_op) :: AA
  type(iso_ladder) :: BB 
  real(8) :: sm ,d6ji,sx,sm1,sm2 ,dcgi

  rank_b = BB%rank
  rank_a = AA%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  holes = sum(jbas%con)
  parts = sum(1-jbas%con) 

  jp1 = jbas%jj(p1)
  jp2 = jbas%jj(p2) 
  jh1 = jbas%jj(h1) 
  jh2 = jbas%jj(h2)

  !!! ASS 
  !  !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) REDUCTION(+:sm)
  do M1 = -1*J1,J1,2
     do M2 = -1*J2,J2,2

        do mh2 = -1*jh2,jh2,2
           do mh1 = -1*jh1,jh1,2
              do mp1 = -1*jp1,jp1,2
                 do mp2 = -1*jp2,jp2,2


                    do mu_c = -1*rank_c,rank_c,2                        
                       do mu_a = -1*rank_a,rank_a,2
                          do mu_b = -1*rank_b,rank_b,2
                             sm1 = 0.d0

                             sm1 = f_tensor_mscheme(p1,mp1,h1,mh1,mu_a,AA,jbas)*&
                                  f_iso_ladder_mscheme(p2,mp2,h2,mh2,mu_b,BB,jbas)

                             sm1 = sm1 - f_tensor_mscheme(p2,mp2,h1,mh1,mu_a,AA,jbas)*&
                                  f_iso_ladder_mscheme(p1,mp1,h2,mh2,mu_b,BB,jbas)

                             sm1 = sm1 + f_tensor_mscheme(p2,mp2,h2,mh2,mu_a,AA,jbas)*&
                                  f_iso_ladder_mscheme(p1,mp1,h1,mh1,mu_b,BB,jbas)

                             sm1 = sm1 - f_tensor_mscheme(p1,mp1,h2,mh2,mu_a,AA,jbas)*&
                                  f_iso_ladder_mscheme(p2,mp2,h1,mh1,mu_b,BB,jbas)


                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a) 
                                do ma = -1*ja,ja,2

                                   sm1 = sm1 + f_tensor_mscheme(p1,mp1,a,ma,mu_a,AA,jbas) * &
                                        iso_ladder_mscheme(a,ma,p2,mp2,h1,mh1,h2,mh2,mu_b,BB,jbas)

                                   sm1 = sm1 - f_tensor_mscheme(p2,mp2,a,ma,mu_a,AA,jbas) * &
                                        iso_ladder_mscheme(a,ma,p1,mp1,h1,mh1,h2,mh2,mu_b,BB,jbas)

                                end do
                             end do

                             do ix = 1, holes
                                i = jbas%holes(ix)
                                ji = jbas%jj(i) 
                                do mi = -1*ji,ji,2

                                   sm1 = sm1 - f_iso_ladder_mscheme(p1,mp1,i,mi,mu_b,BB,jbas) * &
                                        tensor_mscheme(i,mi,p2,mp2,h1,mh1,h2,mh2,mu_a,AA,jbas)

                                   sm1 = sm1 + f_iso_ladder_mscheme(p2,mp2,i,mi,mu_b,BB,jbas) * &
                                        tensor_mscheme(i,mi,p1,mp1,h1,mh1,h2,mh2,mu_a,AA,jbas)

                                end do
                             end do

                             do ix = 1, holes
                                i = jbas%holes(ix)
                                ji = jbas%jj(i)

                                do mi = -1*ji,ji,2
                                   sm1 = sm1 - iso_ladder_mscheme(p1,mp1,p2,mp2,i,mi,h2,mh2,mu_b,BB,jbas) * &
                                        f_tensor_mscheme(i,mi,h1,mh1,mu_a,AA,jbas)  

                                   sm1 = sm1 + iso_ladder_mscheme(p1,mp1,p2,mp2,i,mi,h1,mh1,mu_b,BB,jbas) * &
                                        f_tensor_mscheme(i,mi,h2,mh2,mu_a,AA,jbas)
                                end do
                             end do

                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a)

                                do ma = -1*ja,ja,2
                                   sm1 = sm1 + tensor_mscheme(p1,mp1,p2,mp2,a,ma,h2,mh2,mu_a,AA,jbas) * &
                                        f_iso_ladder_mscheme(a,ma,h1,mh1,mu_b,BB,jbas)  

                                   sm1 = sm1 - tensor_mscheme(p1,mp1,p2,mp2,a,ma,h1,mh1,mu_a,AA,jbas) * &
                                        f_iso_ladder_mscheme(a,ma,h2,mh2,mu_b,BB,jbas)
                                end do
                             end do




                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a) 
                                do ma = -1*ja,ja,2

                                   do bx = 1, parts
                                      b = jbas%parts(bx)
                                      jb = jbas%jj(b) 
                                      do mb = -1*jb,jb,2


                                         sm1 = sm1 + tensor_mscheme(p1,mp1,p2,mp2,a,ma,b,mb,mu_a,AA,jbas)*&
                                              iso_ladder_mscheme(a,ma,b,mb,h1,mh1,h2,mh2,mu_b,BB,jbas) *0.5

                                      end do
                                   end do
                                end do
                             end do

                             do ix = 1, holes
                                i = jbas%holes(ix)
                                ji = jbas%jj(i)

                                do mi = -1*ji,ji,2

                                   do jx = 1, holes
                                      j = jbas%holes(jx)
                                      jj = jbas%jj(j)

                                      do mj = -1*jj,jj,2


                                         sm1 = sm1 + iso_ladder_mscheme(p1,mp1,p2,mp2,i,mi,j,mj,mu_b,BB,jbas)*&
                                              tensor_mscheme(i,mi,j,mj,h1,mh1,h2,mh2,mu_a,AA,jbas) *0.5

                                      end do
                                   end do
                                end do
                             end do

                             do ax = 1, parts
                                a = jbas%parts(ax)
                                ja = jbas%jj(a) 
                                do ma = -1*ja,ja,2

                                   do ix = 1, holes
                                      i = jbas%holes(ix)
                                      ji = jbas%jj(i) 
                                      do mi = -1*ji,ji,2


                                         sm1 = sm1 + tensor_mscheme(p1,mp1,i,mi,h1,mh1,a,ma,mu_a,AA,jbas)*&
                                              iso_ladder_mscheme(p2,mp2,a,ma,h2,mh2,i,mi,mu_b,BB,jbas)

                                         sm1 = sm1 - tensor_mscheme(p2,mp2,i,mi,h1,mh1,a,ma,mu_a,AA,jbas)*&
                                              iso_ladder_mscheme(p1,mp1,a,ma,h2,mh2,i,mi,mu_b,BB,jbas)

                                         sm1 = sm1 + tensor_mscheme(p2,mp2,i,mi,h2,mh2,a,ma,mu_a,AA,jbas)*&
                                              iso_ladder_mscheme(p1,mp1,a,ma,h1,mh1,i,mi,mu_b,BB,jbas)

                                         sm1 = sm1 - tensor_mscheme(p1,mp1,i,mi,h2,mh2,a,ma,mu_a,AA,jbas)*&
                                              iso_ladder_mscheme(p2,mp2,a,ma,h1,mh1,i,mi,mu_b,BB,jbas)


                                      end do
                                   end do
                                end do
                             end do


                             sm = sm + sm1 * dcgi(J2,M2,rank_c,mu_c,J1,M1) * dcgi(jp1,mp1,jp2,mp2,J1,M1)&
                                  *dcgi(jh1,mh1,jh2,mh2,J2,M2) *dcgi(rank_a,mu_a,rank_b,mu_b,rank_c,mu_c)/sqrt(J1+1.d0)  

                          end do
                       end do
                    end do
                 end do

              end do
           end do
        end do
     end do
  end do
!  !$OMP END PARALLEL DO
     
  
  EOM_dTZ_tensor_prod_2body = sm 
  
end function EOM_dTz_tensor_prod_2body
!==================================================================
!==================================================================
real(8) function scalar_tensor_2body_comm(AA,BB,a,b,c,d,J1,J2,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj,J3,J4,J5,jx
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb,rank
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm,coef9,d6ji,pre,ass,smx,sm1,sm2,sm3,sm4

  rank = BB%rank  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

   do i = 1, totorb
      ji = jbas%jj(i)

     sm = sm + f_elem(a,i,AA,jbas) * tensor_elem( i,b,c,d,J1,J2,BB,jbas) &
          + f_elem(b,i,AA,jbas) * tensor_elem( a,i,c,d,J1,J2,BB,jbas) &
          - f_elem(i,c,AA,jbas) * tensor_elem( a,b,i,d,J1,J2,BB,jbas) &
          - f_elem(i,d,AA,jbas) * tensor_elem( a,b,c,i,J1,J2,BB,jbas) 
     
          
     sm = sm - f_tensor_elem(a,i,BB,jbas) * v_elem( i,b,c,d,J2,AA,jbas) &
          * d6ji(ji,jb,J2,J1,rank,ja) * (-1)**((ja+jb+rank-J2)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
     
          + f_tensor_elem(b,i,BB,jbas) * v_elem( i,a,c,d,J2,AA,jbas) &
          * d6ji(ji,ja,J2,J1,rank,jb) * (-1)**((J1+J2+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &

          - f_tensor_elem(i,c,BB,jbas) * v_elem( a,b,d,i,J1,AA,jbas) &
          *d6ji(ji,jc,rank,J2,J1,jd) * (-1)**((J1+J2+rank)/2) *  &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
         
          
          + f_tensor_elem(i,d,BB,jbas) * v_elem( a,b,c,i,J1,AA,jbas) &
          *d6ji( ji,jd,rank,J2,J1,jc) * (-1)**((jc+jd-J1+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) )
     
   end do
  

  do i = 1, totorb
     do j = 1, totorb
        

        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,J1,AA,jbas)*tensor_elem(i,j,c,d,J1,J2,BB,jbas)   &
             - tensor_elem(a,b,i,j,J1,J2,BB,jbas)*v_elem(i,j,c,d,J2,AA,jbas)) 
     end do
  end do
 

!!$OMP PARALLEL DO PRIVATE( ji,jj,i,j,J3,J4,J5,jx) SHARED(AA,BB) REDUCTION(+:sm)
  ! do i = 1, totorb
  !    ji =jbas%jj(i)
  !    do j = 1,totorb
  !       jj = jbas%jj(j) 
        
  !       if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
  !       do J3 = 0, JTM,2
  !          do J4 = 0, JTM,2 
  !             do J5 = 0,JTM,2
  !                do jx = 1,JTM,2
                    
  !                   sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                    !      (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                    !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)   &
                    ! * coef9(jj,J3,ja,J4,ji,jb,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                    !      d6ji(J1,jx,jd,jc,J2,rank) * v_elem(a,j,d,i,J3,AA,jbas) *&
                    !      tensor_elem(i,b,j,c,J4,J5,BB,jbas)  &
               
                    !      - (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                    !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    ! * coef9(jj,J3,jb,J4,ji,ja,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                    !      d6ji(J1,jx,jd,jc,J2,rank) * v_elem(b,j,d,i,J3,AA,jbas) *&
                    !      tensor_elem(i,a,j,c,J4,J5,BB,jbas) *(-1)**((ja+jb-J1)/2) &
                         
                    !      - (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                    !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    ! * coef9(jj,J3,ja,J4,ji,jb,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                    !      d6ji(J1,jx,jc,jd,J2,rank) * v_elem(a,j,c,i,J3,AA,jbas) *&
                    !      tensor_elem(i,b,j,d,J4,J5,BB,jbas) * (-1)**((jc+jd-J2)/2) &

                    !      + (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                    !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    ! * coef9(jj,J3,jb,J4,ji,ja,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                    !      d6ji(J1,jx,jc,jd,J2,rank) * v_elem(b,j,c,i,J3,AA,jbas) *&
                    !      tensor_elem(i,a,j,d,J4,J5,BB,jbas)  *(-1)**((ja+jb+jc+jd+J1+J2)/2)  )
              
                     
       
!                  end do
!               end do
!            end do
!         end do
!      end do
!   end do
! !!$OMP END PARALLEL DO

  smx = 0.d0 
  do J3 = 0,JTM,2
     do J4 = 0,JTM,2 
        smx = 0.d0 
        sm1=0.d0;sm2=0.d0
        sm3=0.d0;sm4=0.d0
        do i = 1, jbas%total_orbits
           ji = jbas%jj(i)
           do j = 1, jbas%total_orbits
              jj = jbas%jj(j) 
              if (jbas%con(i)-jbas%con(j) == 0) cycle
              
              ! smx = smx- (jbas%con(i)-jbas%con(j))*&
              !      (-1)**((J1+J2+J3+J4)/2) * &
              !      sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
              !      coef9(jb,jd,J3,ja,jc,J4,J1,J2,rank)* &
              !      vcc(b,d,j,i,J3,AA,jbas) * Vgenpandya(i,j,c,a,J3,J4,BB,jbas)
              
              sm1= sm1 +(jbas%con(i)-jbas%con(j))*&
                   (-1)**((ja+jb+J2+J3+J4)/2) * &
                   sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
                   coef9(ja,jd,J3,jb,jc,J4,J1,J2,rank)* &
                   vcc(a,d,j,i,J3,AA,jbas) * Vgenpandya(i,j,c,b,J3,J4,BB,jbas)
              
              ! RAGNAR's expression 'cows'
              ! sm = sm + (jbas%con(i)-jbas%con(j))*&
              !      (-1)**((jb+jd+J2+J4)/2) * &
              !      sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
              !      coef9(ja,jd,J3,jb,jc,J4,J1,J2,rank)* &
              !      Vpandya(a,d,i,j,J3,AA,jbas) * Vgenpandya(i,j,c,b,J3,J4,BB,jbas)

              sm2 =  sm2- (jbas%con(i)-jbas%con(j))*&
                   (-1)**((J1+J2+J3+J4)/2) * &
                   sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
                   coef9(jb,jd,J3,ja,jc,J4,J1,J2,rank)* &
                   vcc(b,d,j,i,J3,AA,jbas) * Vgenpandya(i,j,c,a,J3,J4,BB,jbas)

              sm3 = sm3 + (jbas%con(i)-jbas%con(j))*&
                   (-1)**((jc+jd+J1+J3+J4)/2) * &
                   sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
                   coef9(jb,jc,J3,ja,jd,J4,J1,J2,rank)* &
                   vcc(b,c,j,i,J3,AA,jbas) * Vgenpandya(i,j,d,a,J3,J4,BB,jbas)

              sm4 = sm4 - (jbas%con(i)-jbas%con(j))*&
                   (-1)**((ja+jb+jc+jd+J3+J4)/2) * &
                   sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
                   coef9(ja,jc,J3,jb,jd,J4,J1,J2,rank)* &
                   vcc(a,c,j,i,J3,AA,jbas) * Vgenpandya(i,j,d,b,J3,J4,BB,jbas)
              
             
           end do
        end do
        smx = sm1+sm2+sm3+sm4
        sm = sm + smx

        smx = 0.d0 
     end do
  end do

  scalar_tensor_2body_comm = sm 
  
end function scalar_tensor_2body_comm
!==================================================================
!==================================================================
real(8) function scalar_tensor_iso2body_comm(AA,BB,a,b,c,d,J1,J2,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj,J3,J4,J5,jx
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb,rank
  type(spd) :: jbas
  type(sq_op) :: AA
  type(iso_operator) :: BB 
  real(8) :: sm,coef9,d6ji,pre,ass,smx,sm1,sm2,sm3,sm4
  real(8) :: m1,m2,m3,m4
  
  rank = BB%rank  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

   do i = 1, totorb
      ji = jbas%jj(i)

     sm = sm + f_elem(a,i,AA,jbas) * iso_op_elem( i,b,c,d,J1,J2,BB,jbas) &
          + f_elem(b,i,AA,jbas) * iso_op_elem( a,i,c,d,J1,J2,BB,jbas) &
          - f_elem(i,c,AA,jbas) * iso_op_elem( a,b,i,d,J1,J2,BB,jbas) &
          - f_elem(i,d,AA,jbas) * iso_op_elem( a,b,c,i,J1,J2,BB,jbas) 
     
          
     sm = sm - f_iso_op_elem(a,i,BB,jbas) * v_elem( i,b,c,d,J2,AA,jbas) &
          * d6ji(ji,jb,J2,J1,rank,ja) * (-1)**((ja+jb+rank-J2)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
     
          + f_iso_op_elem(b,i,BB,jbas) * v_elem( i,a,c,d,J2,AA,jbas) &
          * d6ji(ji,ja,J2,J1,rank,jb) * (-1)**((J1+J2+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &

          - f_iso_op_elem(i,c,BB,jbas) * v_elem( a,b,d,i,J1,AA,jbas) &
          *d6ji(ji,jc,rank,J2,J1,jd) * (-1)**((J1+J2+rank)/2) *  &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
         
          
          + f_iso_op_elem(i,d,BB,jbas) * v_elem( a,b,c,i,J1,AA,jbas) &
          *d6ji( ji,jd,rank,J2,J1,jc) * (-1)**((jc+jd-J1+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) )
     
    end do
  

  do i = 1, totorb
     do j = 1, totorb
        
        
        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,J1,AA,jbas)*iso_op_elem(i,j,c,d,J1,J2,BB,jbas)   &
             - iso_op_elem(a,b,i,j,J1,J2,BB,jbas)*v_elem(i,j,c,d,J2,AA,jbas)) 
     end do
  end do
 

!!$OMP PARALLEL DO PRIVATE( ji,jj,i,j,J3,J4,J5,jx) SHARED(AA,BB) REDUCTION(+:sm)
  ! do i = 1, totorb
  !    ji =jbas%jj(i)
  !    do j = 1,totorb
  !       jj = jbas%jj(j) 
        
  !       if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
  !       do J3 = 0, JTM,2
  !          do J4 = 0, JTM,2 
  !             do J5 = 0,JTM,2
  !                do jx = 1,JTM,2
                    
  !                   sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                    !      (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                    !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)   &
                    ! * coef9(jj,J3,ja,J4,ji,jb,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                    !      d6ji(J1,jx,jd,jc,J2,rank) * v_elem(a,j,d,i,J3,AA,jbas) *&
                    !      iso_op_elem(i,b,j,c,J4,J5,BB,jbas)  &
               
                    !      - (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                    !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    ! * coef9(jj,J3,jb,J4,ji,ja,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                    !      d6ji(J1,jx,jd,jc,J2,rank) * v_elem(b,j,d,i,J3,AA,jbas) *&
                    !      iso_op_elem(i,a,j,c,J4,J5,BB,jbas) *(-1)**((ja+jb-J1)/2) &
                         
                    !      - (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                    !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    ! * coef9(jj,J3,ja,J4,ji,jb,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                    !      d6ji(J1,jx,jc,jd,J2,rank) * v_elem(a,j,c,i,J3,AA,jbas) *&
                    !      iso_op_elem(i,b,j,d,J4,J5,BB,jbas) * (-1)**((jc+jd-J2)/2) &

                    !      + (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                    !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    ! * coef9(jj,J3,jb,J4,ji,ja,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                    !      d6ji(J1,jx,jc,jd,J2,rank) * v_elem(b,j,c,i,J3,AA,jbas) *&
                    !      iso_op_elem(i,a,j,d,J4,J5,BB,jbas)  *(-1)**((ja+jb+jc+jd+J1+J2)/2)  )
              
                     
       
!                  end do
!               end do
!            end do
!         end do
!      end do
!   end do
! !!$OMP END PARALLEL DO

  smx = 0.d0 
  do J3 = 0,JTM,2
     do J4 = 0,JTM,2 
        smx = 0.d0 
        sm1=0.d0;sm2=0.d0
        sm3=0.d0;sm4=0.d0
        m1=0.d0;m2=0.d0
        m3=0.d0;m4=0.d0
        
        do i = 1, jbas%total_orbits
           ji = jbas%jj(i)
           do j = 1, jbas%total_orbits
              jj = jbas%jj(j) 
              if (jbas%con(i)-jbas%con(j) == 0) cycle
              
              ! smx = smx- (jbas%con(i)-jbas%con(j))*&
              !      (-1)**((J1+J2+J3+J4)/2) * &
              !      sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
              !      coef9(jb,jd,J3,ja,jc,J4,J1,J2,rank)* &
              !      vcc(b,d,j,i,J3,AA,jbas) * Voppandya(i,j,c,a,J3,J4,BB,jbas)
              
              sm1= sm1 +(jbas%con(i)-jbas%con(j))*&
                   (-1)**((ja+jb+J2+J3+J4)/2) * &
                   sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
                   coef9(ja,jd,J3,jb,jc,J4,J1,J2,rank)* &
                   vcc(a,d,j,i,J3,AA,jbas) * Voppandya(i,j,c,b,J3,J4,BB,jbas)

              m1= m1 -(jbas%con(i)-jbas%con(j))*&                  
                   sqrt((J3+1.d0)*(J4+1.d0))*&
                   vcc(a,d,j,i,J3,AA,jbas) * Voppandya(i,j,c,b,J3,J4,BB,jbas)
              ! if ( (jbas%con(i)-jbas%con(j)) ==1 ) then 
              !    if (abs( vcc(a,d,j,i,J3,AA,jbas) *  Voppandya(i,j,c,b,J3,J4,BB,jbas) ) > 1e-10) then 
              !       print*, i,j, vcc(a,d,j,i,J3,AA,jbas) ,  Voppandya(i,j,c,b,J3,J4,BB,jbas) 
              !    end if
              ! end if              
              ! RAGNAR's expression 'cows'
              ! sm = sm + (jbas%con(i)-jbas%con(j))*&
              !      (-1)**((jb+jd+J2+J4)/2) * &
              !      sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
              !      coef9(ja,jd,J3,jb,jc,J4,J1,J2,rank)* &
              !      Vpandya(a,d,i,j,J3,AA,jbas) * Voppandya(i,j,c,b,J3,J4,BB,jbas)

              sm2 =  sm2- (jbas%con(i)-jbas%con(j))*&
                   (-1)**((J1+J2+J3+J4)/2) * &
                   sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
                   coef9(jb,jd,J3,ja,jc,J4,J1,J2,rank)* &
                   vcc(b,d,j,i,J3,AA,jbas) * Voppandya(i,j,c,a,J3,J4,BB,jbas)

              m2 =  m2- (jbas%con(i)-jbas%con(j))*&
                   sqrt((J3+1.d0)*(J4+1.d0))*&                   
                   vcc(b,d,j,i,J3,AA,jbas) * Voppandya(i,j,c,a,J3,J4,BB,jbas)

              sm3 = sm3 + (jbas%con(i)-jbas%con(j))*&
                   (-1)**((jc+jd+J1+J3+J4)/2) * &
                   sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
                   coef9(jb,jc,J3,ja,jd,J4,J1,J2,rank)* &
                   vcc(b,c,j,i,J3,AA,jbas) * Voppandya(i,j,d,a,J3,J4,BB,jbas)

              m3 = m3 + (jbas%con(i)-jbas%con(j))*&                   
                   sqrt((J3+1.d0)*(J4+1.d0))*&
                   vcc(b,c,j,i,J3,AA,jbas) * Voppandya(i,j,d,a,J3,J4,BB,jbas)
                            
              sm4 = sm4 - (jbas%con(i)-jbas%con(j))*&
                   (-1)**((ja+jb+jc+jd+J3+J4)/2) * &
                   sqrt((J1+1.d0)*(J2+1.d0)*(J3+1.d0)*(J4+1.d0))*&
                   coef9(ja,jc,J3,jb,jd,J4,J1,J2,rank)* &
                   vcc(a,c,j,i,J3,AA,jbas) * Voppandya(i,j,d,b,J3,J4,BB,jbas)

              m4 = m4 - (jbas%con(i)-jbas%con(j))*&                   
                   sqrt((J3+1.d0)*(J4+1.d0))*&
                   vcc(a,c,j,i,J3,AA,jbas) * Voppandya(i,j,d,b,J3,J4,BB,jbas)
                            
             
           end do
        end do
   
 !       if (abs(sm1) > 1e-10 ) print*, '1',sm1,m1, (-1)**((ja-jb+J2+J3+J4)/2) * &
  !                 sqrt((J1+1.d0)*(J2+1.d0))*coef9(ja,jd,J3,jb,jc,J4,J1,J2,rank)
                                      
!        if (abs(sm2) > 1e-10 ) print*, '2',sm2,m2
!        if (abs(sm3) > 1e-10 ) print*, '4',sm3,m3
  !      if (abs(sm4) > 1e-10 ) print*, '3',sm4,m4
        smx = sm1+sm2+sm3+sm4
        sm = sm + smx

        smx = 0.d0 
     end do
  end do

  scalar_tensor_iso2body_comm = sm 
  
end function scalar_tensor_iso2body_comm
!============================================================
!============================================================
real(8) function EOM_scalar_tensor_2body_comm(AA,BB,a,b,c,d,J1,J2,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj,J3,J4,J5,jx
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb,rank
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm,coef9,d6ji

  rank = BB%rank  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

    do i = 1, totorb
       ji = jbas%jj(i)

      sm = sm + f_elem(a,i,AA,jbas) * pphh_tensor_elem( i,b,c,d,J1,J2,BB,jbas) &
           + f_elem(b,i,AA,jbas) * pphh_tensor_elem( a,i,c,d,J1,J2,BB,jbas) &
           - f_elem(i,c,AA,jbas) * pphh_tensor_elem( a,b,i,d,J1,J2,BB,jbas) &
           - f_elem(i,d,AA,jbas) * pphh_tensor_elem( a,b,c,i,J1,J2,BB,jbas) 
     
          
     sm = sm - ph_tensor_elem(a,i,BB,jbas) * v_elem( i,b,c,d,J2,AA,jbas) &
          * d6ji(ji,jb,J2,J1,rank,ja) * (-1)**((ja+jb+rank-J2)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
     
          + ph_tensor_elem(b,i,BB,jbas) * v_elem( i,a,c,d,J2,AA,jbas) &
          * d6ji(ji,ja,J2,J1,rank,jb) * (-1)**((J1+J2+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &

          - ph_tensor_elem(i,c,BB,jbas) * v_elem( a,b,d,i,J1,AA,jbas) &
          *d6ji(ji,jc,rank,J2,J1,jd) * (-1)**((J1+J2+rank)/2) *  &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
         
          
          + ph_tensor_elem(i,d,BB,jbas) * v_elem( a,b,c,i,J1,AA,jbas) &
          *d6ji( ji,jd,rank,J2,J1,jc) * (-1)**((jc+jd-J1+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) )
     
    end do
  

  do i = 1, totorb
     do j = 1, totorb
        

        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,J1,AA,jbas)*pphh_tensor_elem(i,j,c,d,J1,J2,BB,jbas)   &
             - pphh_tensor_elem(a,b,i,j,J1,J2,BB,jbas)*v_elem(i,j,c,d,J2,AA,jbas)) 
     end do
  end do
 
!!$OMP PARALLEL DO PRIVATE( ji,jj,i,j,J3,J4,J5,jx) SHARED(AA,BB) REDUCTION(+:sm)
  do i = 1, totorb
     ji =jbas%jj(i)
     do j = 1,totorb
        jj = jbas%jj(j) 
        
        if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
        do J3 = 0, JTM,2
           do J4 = 0, JTM,2 
              do J5 = 0,JTM,2
                 do jx = 1,JTM,2
                    
                    sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                         (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)   &
                    * coef9(jj,J3,ja,J4,ji,jb,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                         d6ji(J1,jx,jd,jc,J2,rank) * v_elem(a,j,d,i,J3,AA,jbas) *&
                         pphh_tensor_elem(i,b,j,c,J4,J5,BB,jbas) &
                   
                         - (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,jb,J4,ji,ja,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                         d6ji(J1,jx,jd,jc,J2,rank) * v_elem(b,j,d,i,J3,AA,jbas) *&
                         pphh_tensor_elem(i,a,j,c,J4,J5,BB,jbas) *(-1)**((ja+jb-J1)/2) &
                   
                         - (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,ja,J4,ji,jb,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                         d6ji(J1,jx,jc,jd,J2,rank) * v_elem(a,j,c,i,J3,AA,jbas) *&
                         pphh_tensor_elem(i,b,j,d,J4,J5,BB,jbas) * (-1)**((jc+jd-J2)/2) &
                   
                         + (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,jb,J4,ji,ja,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                         d6ji(J1,jx,jc,jd,J2,rank) * v_elem(b,j,c,i,J3,AA,jbas) *&
                         pphh_tensor_elem(i,a,j,d,J4,J5,BB,jbas)  *(-1)**((ja+jb+jc+jd+J1+J2)/2)  &
                         )
                 end do
              end do
           end do
        end do
     end do
  end do
!!$OMP END PARALLEL DO

  EOM_scalar_tensor_2body_comm = sm 
  
end function EOM_scalar_tensor_2body_comm
!============================================================
!============================================================
real(8) function EOM_scalar_tensor_iso2body_comm(AA,BB,a,b,c,d,J1,J2,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj,J3,J4,J5,jx
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb,rank
  type(spd) :: jbas
  type(sq_op) :: AA
  type(iso_ladder) :: BB 
  real(8) :: sm,coef9,d6ji,X1,X2,X3,X4

  rank = BB%rank  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

    do i = 1, totorb
       ji = jbas%jj(i)

      sm = sm + f_elem(a,i,AA,jbas) * iso_ladder_elem( i,b,c,d,J1,J2,BB,jbas) &
           + f_elem(b,i,AA,jbas) * iso_ladder_elem( a,i,c,d,J1,J2,BB,jbas) &
           - f_elem(i,c,AA,jbas) * iso_ladder_elem( a,b,i,d,J1,J2,BB,jbas) &
           - f_elem(i,d,AA,jbas) * iso_ladder_elem( a,b,c,i,J1,J2,BB,jbas) 
     
          
     sm = sm - f_iso_ladder_elem(a,i,BB,jbas) * v_elem( i,b,c,d,J2,AA,jbas) &
          * d6ji(ji,jb,J2,J1,rank,ja) * (-1)**((ja+jb+rank-J2)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
     
          + f_iso_ladder_elem(b,i,BB,jbas) * v_elem( i,a,c,d,J2,AA,jbas) &
          * d6ji(ji,ja,J2,J1,rank,jb) * (-1)**((J1+J2+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &

          - f_iso_ladder_elem(i,c,BB,jbas) * v_elem( a,b,d,i,J1,AA,jbas) &
          *d6ji(ji,jc,rank,J2,J1,jd) * (-1)**((J1+J2+rank)/2) *  &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
         
          
          + f_iso_ladder_elem(i,d,BB,jbas) * v_elem( a,b,c,i,J1,AA,jbas) &
          *d6ji( ji,jd,rank,J2,J1,jc) * (-1)**((jc+jd-J1+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) )
     
    end do
  

  do i = 1, totorb
     do j = 1, totorb
        

        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,J1,AA,jbas)*iso_ladder_elem(i,j,c,d,J1,J2,BB,jbas)   &
             - iso_ladder_elem(a,b,i,j,J1,J2,BB,jbas)*v_elem(i,j,c,d,J2,AA,jbas)) 
     end do
  end do
 
! !!$OMP PARALLEL DO PRIVATE( ji,jj,i,j,J3,J4,J5,jx) SHARED(AA,BB) REDUCTION(+:sm)
!   do i = 1, totorb
!      ji =jbas%jj(i)
!      do j = 1,totorb
!         jj = jbas%jj(j) 
        
!         if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
!         do J3 = 0, JTM,2
!            do J4 = 0, JTM,2 
!               do J5 = 0,JTM,2
!                  do jx = 1,JTM,2
                    
!                     sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
!                     !      (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
!                     !      * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)   &
!                     ! * coef9(jj,J3,ja,J4,ji,jb,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
!                     !      d6ji(J1,jx,jd,jc,J2,rank) * v_elem(a,j,d,i,J3,AA,jbas) *&
!                     !      iso_ladder_elem(i,b,j,c,J4,J5,BB,jbas) &
                         
!                          - (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
!                          * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
!                     * coef9(jj,J3,jb,J4,ji,ja,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
!                          d6ji(J1,jx,jd,jc,J2,rank) * v_elem(b,j,d,i,J3,AA,jbas) *&
!                          iso_ladder_elem(i,a,j,c,J4,J5,BB,jbas) *(-1)**((ja+jb-J1)/2) &
                   
!                          - (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
!                          * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
!                     * coef9(jj,J3,ja,J4,ji,jb,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
!                          d6ji(J1,jx,jc,jd,J2,rank) * v_elem(a,j,c,i,J3,AA,jbas) *&
!                          iso_ladder_elem(i,b,j,d,J4,J5,BB,jbas) * (-1)**((jc+jd-J2)/2) &
                   
!                          + (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
!                          * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
!                     * coef9(jj,J3,jb,J4,ji,ja,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
!                          d6ji(J1,jx,jc,jd,J2,rank) * v_elem(b,j,c,i,J3,AA,jbas) *&
!                          iso_ladder_elem(i,a,j,d,J4,J5,BB,jbas)  *(-1)**((ja+jb+jc+jd+J1+J2)/2)  &
!                          )
!                  end do
!               end do
!            end do
!         end do
!      end do
!   end do
! !!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE( ji,jj,i,j,J3,J4,J5,jx) SHARED(AA,BB) REDUCTION(+:sm)
!  print*
  do J3 = 0, JTM,2
     do J4 = 0, JTM,2 

        X1 = 0.d0
        X2 = 0.d0
        X3 = 0.d0
        X4 = 0.d0
        do i = 1, totorb
           ji =jbas%jj(i)
           do j = 1,totorb
              jj = jbas%jj(j) 

              if ((jbas%con(i)-jbas%con(j)) /= -1) cycle 

              X1 = X1 + (-1)**((J3+J4)/2) * sqrt((J3+1.d0)*(J4+1.d0)) * &
                   VCC(a,d,j,i,J3,AA,jbas) * Visopandya(i,j,c,b,J3,J4,BB,jbas)

              X2 = X2 + (-1)**((J3+J4)/2) * sqrt((J3+1.d0)*(J4+1.d0)) * &
                   VCC(b,d,j,i,J3,AA,jbas) * Visopandya(i,j,c,a,J3,J4,BB,jbas)

              X3 = X3 + (-1)**((J3+J4)/2) * sqrt((J3+1.d0)*(J4+1.d0)) * &
                   VCC(b,c,j,i,J3,AA,jbas) * Visopandya(i,j,d,a,J3,J4,BB,jbas)

              X4 = X4 + (-1)**((J3+J4)/2) * sqrt((J3+1.d0)*(J4+1.d0)) * &
                   VCC(a,c,j,i,J3,AA,jbas) * Visopandya(i,j,d,b,J3,J4,BB,jbas)
                   
           end do
        end do

           sm = sm +  ( &  
             
             - (-1)**((ja+jb+J2)/2) * sqrt((J1+1.d0)*(J2+1.d0)) * &
             coef9(ja,jd,J3,jb,jc,J4,J1,J2,rank) * X1 & 
                   
             +(-1)**((J1+J2)/2) * sqrt((J1+1.d0)*(J2+1.d0)) * &
             coef9(jb,jd,J3,ja,jc,J4,J1,J2,rank) * X2 & 
             
             - (-1)**((jc+jd+J1)/2) * sqrt((J1+1.d0)*(J2+1.d0)) * &
             coef9(jb,jc,J3,ja,jd,J4,J1,J2,rank) * X3 &
             
             + (-1)**((ja+jb+jc+jd)/2) * sqrt((J1+1.d0)*(J2+1.d0)) * &
             coef9(ja,jc,J3,jb,jd,J4,J1,J2,rank)* X4 &
             
             )

     end do
  end do
!!$OMP END PARALLEL DO

  
  EOM_scalar_tensor_iso2body_comm = sm 
  
end function EOM_scalar_tensor_iso2body_comm

end module 

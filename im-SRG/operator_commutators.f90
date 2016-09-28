module operator_commutators
  use cross_coupled
  use isospin_operators
  ! tensor-scalar commutator functions 
  
  ! THE TENSOR MUST BE THE SECOND ARGUMENT
  
contains
!=========================================================
!=========================================================
subroutine operator_commutator_111(L,R,RES,jbas) 
  ! one body part of [L1,R1] 
  implicit none 
  
  type(sq_op) :: L
  type(iso_operator) :: R,RES
  type(spd) :: jbas
  integer :: p,a,b,i,ji,hol,par
  real(8) :: sm 
  real(8),dimension(L%Nsp,L%Nsp) :: Lfock

  hol = L%belowEf
  par = L%Nsp - hol

  do a = 1, L%Nsp 
     do b = 1,L%Nsp
        sm = 0.d0 
        do i = 1, L%Nsp 

           sm = sm + f_elem(a,i,L,jbas)  * f_iso_op_elem(i,b,R,jbas) &
                - f_iso_op_elem(a,i,R,jbas) * f_elem(i,b,L,jbas)  
        end do

        RES%fock(a,b) = sm
     end do
  end do
  
end subroutine
!=========================================================
!=========================================================
subroutine operator_commutator_121(L,R,RES,jbas) 
  ! onebody part of [L1,R2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L
  type(iso_operator) :: R,RES
  integer :: J1,J2,PAR,TZ,ji,ja,jp,jq,a,i,p,q,g
  integer :: ti,ta,tp,tq,li,la,lp,lq,ak,ik,pk,qk,rank
  real(8) :: sm,smx,smy,smx2,smy2,d6ji
  
  rank = R%rank 

  do p = 1,L%Nsp
     
     pk = p
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = 1,L%Nsp 
      
        qk = q 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        ! check if this state is allowed
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
        if (tq .ne. tp-2*RES%dTz) cycle 
        if (.not. (triangle(jq,jp,rank))) cycle
        
        ! internal sum
        sm = 0.d0
        do a = 1,L%Nsp - L%belowEF
           
           ak = jbas%parts(a)
           ja = jbas%jj(ak) 
           ta = jbas%itzp(ak)
           la = jbas%ll(ak)
           
           do i = 1,L%belowEF
              
              ik = jbas%holes(i)
              ji = jbas%jj(ik) 
              ti = jbas%itzp(ik)
              li = jbas%ll(ik)
              
              ! check if this intermediate exists
              if (ji .ne. ja) cycle
              if (li .ne. la) cycle
              if (ti .ne. ta) cycle 
                 
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 
              ! sum over J_total
              do J2 = abs(ji - jq),ji+jq,2
                 do J1 = abs(ja - jp),ja+jp,2 
                    
                    smx = smx + iso_op_elem(ak,pk,ik,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)*&
                         xxxsixj(RES%xindx,J1,J2,rank,jq,jp,ji)             
                 
                 end do
              end do
              
              do J1 = abs(ji - jp),ji+jp,2
                 do J2 = abs(ja - jq),ja+jq,2 
                 
                    smy = smy + iso_op_elem(ik,pk,ak,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)*&
                         xxxsixj(RES%xindx,J1,J2,rank,jq,jp,ji) 
                                     
                 end do
              end do
                            
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) &
              *(-1)**( (rank + jq + ji )/2 ) 
              
           end do 
        end do 
        
        RES%fock(p,q) = RES%fock(p,q) + sm
        
     end do 
  end do 
          
end subroutine             

!====================================================================
!====================================================================
subroutine operator_commutator_211(L,R,RES,jbas) 
  ! onebody part of  - [R1,L2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(iso_operator) :: R,RES
  type(sq_op) :: L
  integer :: J1,J2,PAR,TZ,ji,ja,jp,jq,a,i,p,q,g,JTM
  integer :: ti,ta,tp,tq,li,la,lp,lq,ak,ik,pk,qk,rank
  integer :: rai,rpq,rqp,qx,Ntot
  real(8) :: sm,smx,smy,smx2,smy2,d6ji
  
  rank = R%rank 
  JTM = Jbas%Jtotal_max
  Ntot = R%nsp 
! dfhh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1,R%Nsp
     
     pk = p
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = 1,R%Nsp
      
        qk = q
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        ! check if this state is allowed
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
        if (tq .ne. tp-2*RES%dTz) cycle 
        if (.not. (triangle(jq,jp,rank))) cycle
        
   
        ! internal sum
        sm = 0.d0
        do a = 1,R%Nsp - R%belowEF
           
           ak = jbas%parts(a)
           ja = jbas%jj(ak) 
           ta = jbas%itzp(ak)
           la = jbas%ll(ak)
           
           do i = 1,R%belowEF
              
              ik = jbas%holes(i)
              ji = jbas%jj(ik) 
              ti = jbas%itzp(ik)
              li = jbas%ll(ik)
              
              ! check if this intermediate exists
              if ( mod(li,2) .ne. mod(la+R%dpar/2,2)) cycle
              if (.not. (triangle(ja,ji,rank))) cycle
              Tz = abs(ta -ti)/2 
              if (abs(tp - tq) .ne. Tz*2)  cycle 
              PAR  = mod(la+li,2) 
              if (mod(lp+lq,2) .ne. PAR) cycle                                 
              if (ti  ==  ta-2*RES%dTz) then  

                 sm = sm -  R%fock(ak,ik) * Vpandya(ak,ik,p,q,rank,L,jbas)
              end if

              if (ti  ==  ta+2*RES%dTz) then  

                 sm = sm +  R%fock(ik,ak) * Vpandya(ik,ak,p,q,rank,L,jbas)
              end if
              
           end do 
        end do 
   
        RES%fock(p,q) = RES%fock(p,q) + sm  
        
     end do 
  end do 

end subroutine operator_commutator_211
!==================================================
!==================================================             
subroutine operator_commutator_122(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L
  type(iso_operator) :: R,RES
  integer :: q,IX,JX,nh,np,nb,i,J1,J2
  integer :: a,b,c,d,ja,jb,jc,jd,ji,g_ix,q_sp,i_sp
  integer :: ta,tb,tc,td,ti,la,lb,lc,ld,li,spec
  integer :: jxstart,jxend,ixend,c1,c2,n1,n2
  logical :: square
  real(8) ::  sm
  
  
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1)
     J2 = R%tblck(q)%Jpair(2) 
    
     do g_ix = 1,9 
   
        ! figure out how big the array is
        n1 = size(R%tblck(q)%tgam(g_ix)%X(:,1))
        n2 = size(R%tblck(q)%tgam(g_ix)%X(1,:))

        if ((n1*n2) == 0) cycle 
        
        ! read in information about which 
        ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        square = sqs(g_ix) 
        jxstart = jst(g_ix) 
        
                
     ! main calculation
   
     do IX = 1,n1
        a = R%tblck(q)%tensor_qn(c1,1)%Y(IX,1)
        ja = jbas%jj(a)
        la = jbas%ll(a)
        ta = jbas%itzp(a) 
             
        b = R%tblck(q)%tensor_qn(c1,1)%Y(IX,2)
        jb = jbas%jj(b)
        lb = jbas%ll(b)
        tb = jbas%itzp(b)
 
        do JX = 1,n2
           
           c = R%tblck(q)%tensor_qn(c2,2)%Y(JX,1)
           jc = jbas%jj(c)
           lc = jbas%ll(c)
           tc = jbas%itzp(c)

           d = R%tblck(q)%tensor_qn(c2,2)%Y(JX,2)
           jd = jbas%jj(d)
           ld = jbas%ll(d)
           td = jbas%itzp(d)
                          
           sm = 0.d0 
           
            ! a is replaced
            q_sp = sp_block_index(ja,la,ta,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i)
               
               sm = sm + f_elem(a,i_sp,L,jbas)*iso_op_elem(i_sp,b,c,d,J1,J2,R,jbas)
                  
            end do
           
            
            ! b is replaced
            q_sp = sp_block_index(jb,lb,tb,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               sm = sm + f_elem(b,i_sp,L,jbas)*iso_op_elem(a,i_sp,c,d,J1,J2,R,jbas) 
            end do 
            
            ! c is replaced
            q_sp = sp_block_index(jc,lc,tc,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm - f_elem(i_sp,c,L,jbas)*iso_op_elem(a,b,i_sp,d,J1,J2,R,jbas)
            end do 
            
            ! d is replaced
            q_sp = sp_block_index(jd,ld,td,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm - f_elem(i_sp,d,L,jbas)*iso_op_elem(a,b,c,i_sp,J1,J2,R,jbas)
            end do 
          
            sm = sm / sqrt(1.d0 + kron_del(a,b)) /sqrt(1.d0 + kron_del(c,d)) 

           RES%tblck(q)%tgam(g_ix)%X(IX,JX) = sm 

        end do
     end do
   
     end do 
  end do 

end subroutine operator_commutator_122
!==================================================
!==================================================             
subroutine operator_commutator_212(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L
  type(iso_operator) :: R,RES
  integer :: q,IX,JX,nh,np,nb,i,J1,J2
  integer :: a,b,c,d,ja,jb,jc,jd,ji,g_ix,q_sp,i_sp
  integer :: ta,tb,tc,td,ti,modla,modlb,modlc,modld,modli,spec,rank
  integer :: jxstart,jxend,ixend,c1,c2,n1,n2
  logical :: square
  real(8) ::  sm,sm1,sm2,sm3,sm4,d6ji
  
  rank = R%rank
  
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1)
     J2 = R%tblck(q)%Jpair(2) 
    
     do g_ix = 1,9 
   
        ! figure out how big the array is
        n1 = size(R%tblck(q)%tgam(g_ix)%X(:,1))
        n2 = size(R%tblck(q)%tgam(g_ix)%X(1,:))
        if ((n1*n2) == 0) cycle 
        
        ! read in information about which 
        ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        square = sqs(g_ix) 
        jxstart = jst(g_ix) 
        
                
     ! main calculation
   
     do IX = 1,n1
        a = R%tblck(q)%tensor_qn(c1,1)%Y(IX,1)
        ja = jbas%jj(a)
        modla = mod(jbas%ll(a),2)
        ta = jbas%itzp(a) 
             
        b = R%tblck(q)%tensor_qn(c1,1)%Y(IX,2)
        jb = jbas%jj(b)
        modlb = mod(jbas%ll(b),2)
        tb = jbas%itzp(b)
 
        do JX = 1,n2
           
           c = R%tblck(q)%tensor_qn(c2,2)%Y(JX,1)
           jc = jbas%jj(c)
           modlc = mod(jbas%ll(c),2)
           tc = jbas%itzp(c)

           d = R%tblck(q)%tensor_qn(c2,2)%Y(JX,2)
           jd = jbas%jj(d)
           modld = mod(jbas%ll(d),2)
           td = jbas%itzp(d)
                   
           sm = 0.d0 

           do i = 1,jbas%total_orbits
               
              ji = jbas%jj(i) 
              ti = jbas%itzp(i)
              modli = mod(jbas%ll(i)+R%dpar/2,2) 
              
              sm1 = 0.d0 
              if (ti +2*RES%dTz  == ta) then
                 if (modli == modla ) then 
                    if (triangle(ji,ja,rank)) then  
                       
                       sm1 = sm1 - xxxsixj(RES%xindx,J1,J2,rank,ji,ja,jb)&
                            *f_iso_op_elem(a,i,R,jbas)*v_elem(i,b,c,d,J2,L,jbas)
                     
                    end if
                 end if
              end if
               
              sm1 = sm1*(-1)**((ja + jb-J2)/2) 
               
              
              sm2 = 0.d0 
              if (ti + 2*RES%dTz == tb) then
                 if (modli == modlb ) then 
                    if (triangle(ji,jb,rank)) then  
                       
                       sm2 = sm2 + xxxsixj(RES%xindx,J1,J2,rank,ji,jb,ja)&
                            *f_iso_op_elem(b,i,R,jbas)*v_elem(i,a,c,d,J2,L,jbas)
                       
                    end if
                 end if
              end if
               
              sm2 = sm2*(-1)**((J1+J2)/2) 
              

              sm3 = 0.d0 
              if (ti - 2*RES%dTz == td) then
                 if (modli == modld ) then 
                    if (triangle(ji,jd,rank)) then  
                       
                       sm3 = sm3 + xxxsixj(RES%xindx,J1,J2,rank,jd,ji,jc)&
                            *f_iso_op_elem(i,d,R,jbas)*v_elem(a,b,c,i,J1,L,jbas)
                       
                    end if
                 end if
              end if
               
              sm3 = sm3*(-1)**((jc+jd-J1)/2) 
              

              sm4 = 0.d0 
              if (ti - 2*RES%dTz == tc) then
                 if (modli == modlc ) then 
                    if (triangle(ji,jc,rank)) then  
                       
                       sm4 = sm4 -  xxxsixj(RES%xindx,J1,J2,rank,jc,ji,jd)&
                            *f_iso_op_elem(i,c,R,jbas)*v_elem(a,b,d,i,J1,L,jbas)
                       
                    end if
                 end if
              end if
               
              sm4 = sm4*(-1)**((J2+J1)/2)
              
              sm =  sm + (sm1+sm2+sm3+sm4) 
           end do 
           
           sm = sm * sqrt((J1+1.d0)*(J2+1.d0) / &
              (1.d0 + kron_del(a,b)) /(1.d0 + kron_del(c,d))) * (-1)**(rank/2) 
 
           RES%tblck(q)%tgam(g_ix)%X(IX,JX) = RES%tblck(q)%tgam(g_ix)%X(IX,JX)  +sm 
           
        end do
     end do
   
     end do 
  end do 

end subroutine operator_commutator_212
!===================================================================
!===================================================================
subroutine operator_commutator_222_pp_hh(L,R,RES,jbas) 
  !VERIFIED
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  L
  type(iso_operator) :: R,RES
  integer :: q,qx,q1,q2,J1,J2,Tz,Par,phase,rank,i
  integer :: np1,nb1,nh1,np2,nb2,nh2,pm,jh,ji,ja,jp,ii,jj
  integer :: nhb1,nhb2,ntot1,ntot2,IX,JX,px,hx,p1,h1,kk,ll,a 
  real(8) :: bet_off,al_off,pre1,pre2
  real(8),allocatable,dimension(:,:) :: W,Y
  integer,allocatable,dimension(:,:) :: qn1,qn2
  
  bet_off = 1.d0 
  pm = R%herm*L%herm
  rank = R%rank
!construct temporary matrices
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1) 
     J2 = R%tblck(q)%Jpair(2)
     phase = R%tblck(q)%lam(1)
     par = R%tblck(q)%lam(2) 
     Tz = R%tblck(q)%lam(3)
    
     q1 = block_index(J1,Tz,Par) 
     q2 = block_index(J2,Tz-RES%dTz,mod(Par+R%Dpar/2,2)) 
     qx = iso_ladder_block_index(J2,J1,RANK,Tz-RES%dTz,mod(Par+R%Dpar/2,2))                
     
     nh1 = R%tblck(q)%nhh1
     np1 = R%tblck(q)%npp1
     nb1 = R%tblck(q)%nph1
     nh2 = R%tblck(q)%nhh2
     np2 = R%tblck(q)%npp2
     nb2 = R%tblck(q)%nph2

     nhb1 = nh1 + nb1
     nhb2 = nh2 + nb2
     ntot1 = nhb1 + np1
     ntot2 = nhb2 + np2
     if (ntot1*ntot2==0) cycle
     allocate(W(ntot1,ntot2)) 
     allocate(Y(ntot1,ntot2))
     W = 0.d0
     Y = 0.d0 
     !----------------------------------------------------------------------------
     !         Zpppp 
     !----------------------------------------------------------------------------
     if (np1*np2 .ne. 0)  then 

        al_off = 1.d0 
        call dgemm('N','N', np1,np2,np1,al_off,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%tgam(1)%X,np1,bet_off,W(nhb1+1:ntot1,nhb2+1:ntot2),np1)

        al_off = -1.d0 
        call dgemm('N','N', np1,np2,np2,al_off,R%tblck(q)%tgam(1)%X,np1,&
             L%mat(q2)%gam(1)%X,np2,bet_off,W(nhb1+1:ntot1,nhb2+1:ntot2),np1)


        if(nh1 .ne. 0) then        
           al_off = -1.d0 
           call dgemm('N','N',np1,np2,nh1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(7)%X,nh1,bet_off,Y(nhb1+1:ntot1,nhb2+1:ntot2),np1)
        end if

        if (nh2 .ne. 0 ) then 
           al_off = 1.d0*L%herm        
           call dgemm('N','T',np1,np2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
                L%mat(q2)%gam(3)%X,np2,bet_off,Y(nhb1+1:ntot1,nhb2+1:ntot2),np1)
        end if
     end if

     RES%tblck(q)%tgam(1)%X = RES%tblck(q)%tgam(1)%X  + &
          W(nhb1+1:ntot1,nhb2+1:ntot2) +  Y(nhb1+1:ntot1,nhb2+1:ntot2)   
!----------------------------------------------------------------------------
!         Zhhhh 
!----------------------------------------------------------------------------
     if (nh1*nh2 .ne. 0)  then 
        if(np1 .ne. 0) then        
           al_off = 1.d0 * L%herm  
           call dgemm('T','N', nh1,nh2,np1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(3)%X,np1,bet_off,W(1:nh1,1:nh2),nh1)
        end if

        if (np2 .ne. 0 ) then 
           al_off = -1.d0 
           call dgemm('N','N', nh1,nh2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
                L%mat(q2)%gam(3)%X,np2,bet_off,W(1:nh1,1:nh2),nh1)
        end if

        al_off = -1.d0 
        call dgemm('N','N',nh1,nh2,nh1,al_off,L%mat(q1)%gam(5)%X,nh1,&
             R%tblck(q)%tgam(5)%X,nh1,bet_off,Y(1:nh1,1:nh2),nh1)
                
        al_off = 1.d0        
        call dgemm('N','N',nh1,nh2,nh2,al_off,R%tblck(q)%tgam(5)%X,nh1,&
             L%mat(q2)%gam(5)%X,nh2,bet_off,Y(1:nh1,1:nh2),nh1)
     end if

     RES%tblck(q)%tgam(5)%X = RES%tblck(q)%tgam(5)%X + &
          W(1:nh1,1:nh2) + Y(1:nh1,1:nh2) 
!----------------------------------------------------------------------------
!         Zphph 
!----------------------------------------------------------------------------
     if (nb1*nb2 .ne. 0)  then 

        if (np1 .ne. 0 ) then 
           al_off = 1.d0* L%herm 
           call dgemm('T','N', nb1,nb2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
                R%tblck(q)%tgam(2)%X,np1,bet_off,W(nh1+1:nhb1,nh2+1:nhb2),nb1)
        end if

        if (np2 .ne. 0 ) then 
           al_off = -1.d0 
           call dgemm('N','N', nb1,nb2,np2,al_off,R%tblck(q)%tgam(8)%X,nb1,&
                L%mat(q2)%gam(2)%X,np2,bet_off,W(nh1+1:nhb1,nh2+1:nhb2),nb1)
        end if

        if(nh1 .ne. 0) then        
           al_off = -1.d0 
           call dgemm('N','N',nb1,nb2,nh1,al_off,L%mat(q1)%gam(6)%X,nb1,&
                R%tblck(q)%tgam(9)%X,nh1,bet_off,Y(nh1+1:nhb1,nh2+1:nhb2),nb1)
        end if

        if (nh2 .ne. 0 ) then 
            al_off = 1.d0 * L%herm
            call dgemm('N','T',nb1,nb2,nh2,al_off,R%tblck(q)%tgam(6)%X,nb1,&
                 L%mat(q2)%gam(6)%X,nb2,bet_off,Y(nh1+1:nhb1,nh2+1:nhb2),nb1)
        end if
     end if

     RES%tblck(q)%tgam(4)%X = RES%tblck(q)%tgam(4)%X + &
          W(nh1+1:nhb1,nh2+1:nhb2) + Y(nh1+1:nhb1,nh2+1:nhb2)
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------
     if (np1*nh2 .ne. 0)  then 

        al_off = 1.d0 
        call dgemm('N','N', np1,nh2,np1,al_off,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%tgam(3)%X,np1,bet_off,W(nhb1+1:ntot1,1:nh2),np1)

        if (np2 .ne. 0 ) then         
           al_off = -1.d0 
           call dgemm('N','N', np1,nh2,np2,al_off,R%tblck(q)%tgam(1)%X,np1,&
                L%mat(q2)%gam(3)%X,np2,bet_off,W(nhb1+1:ntot1,1:nh2),np1)
        end if

        if(nh1 .ne. 0) then        
           al_off = -1.d0 
           call dgemm('N','N',np1,nh2,nh1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(5)%X,nh1,bet_off,Y(nhb1+1:ntot1,1:nh2),np1)
        end if

        if (nh2 .ne. 0 ) then 
           al_off = 1.d0        
           call dgemm('N','N',np1,nh2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
                L%mat(q2)%gam(5)%X,nh2,bet_off,Y(nhb1+1:ntot1,1:nh2),np1)
        end if
     end if

     RES%tblck(q)%tgam(3)%X = RES%tblck(q)%tgam(3)%X + &
          W(nhb1+1:ntot1,1:nh2) + Y(nhb1+1:ntot1,1:nh2)
!----------------------------------------------------------------------------
!         Zhhpp 
!----------------------------------------------------------------------------
     if (nh1*np2 .ne. 0)  then 

        if (np1 .ne. 0 ) then 
           al_off = 1.d0 * L%herm  
           call dgemm('T','N', nh1,np2,np1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(1)%X,np1,bet_off,W(1:nh1,nhb2+1:ntot2),nh1)
        end if
        
        al_off = -1.d0 
        call dgemm('N','N', nh1,np2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
             L%mat(q2)%gam(1)%X,np2,bet_off,W(1:nh1,nhb2+1:ntot2),nh1)


        if(nh1 .ne. 0) then        
           al_off = -1.d0 
           call dgemm('N','N',nh1,np2,nh1,al_off,L%mat(q1)%gam(5)%X,nh1,&
                R%tblck(q)%tgam(7)%X,nh1,bet_off,Y(1:nh1,nhb2+1:ntot2),nh1)
        end if

        if (nh2 .ne. 0 ) then 
           al_off = 1.d0  * L%herm       
           call dgemm('N','T',nh1,np2,nh2,al_off,R%tblck(q)%tgam(5)%X,nh1,&
                L%mat(q2)%gam(3)%X,np2,bet_off,Y(1:nh1,nhb2+1:ntot2),nh1)
        end if
     end if

     RES%tblck(q)%tgam(7)%X = RES%tblck(q)%tgam(7)%X  + &
          W(1:nh1,nhb2+1:ntot2) + Y(1:nh1,nhb2+1:ntot2)
!----------------------------------------------------------------------------
!         Zppph 
!----------------------------------------------------------------------------
     if (np1*nb2 .ne. 0)  then 

        al_off = 1.d0
        call dgemm('N','N', np1,nb2,np1,al_off,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%tgam(2)%X,np1,bet_off,W(nhb1+1:ntot1,nh2+1:nhb2),np1)

        if (np2 .ne. 0 ) then 
           al_off = -1.d0 
           call dgemm('N','N', np1,nb2,np2,al_off,R%tblck(q)%tgam(1)%X,np1,&
                L%mat(q2)%gam(2)%X,np2,bet_off,W(nhb1+1:ntot1,nh2+1:nhb2),np1)
        end if

        if(nh1 .ne. 0) then        
           al_off = -1.d0 
           call dgemm('N','N',np1,nb2,nh1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(9)%X,nh1,bet_off,Y(nhb1+1:ntot1,nh2+1:nhb2),np1)
        end if

        if (nh2 .ne. 0 ) then 
           al_off = 1.d0 * L%herm
           call dgemm('N','T',np1,nb2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
                L%mat(q2)%gam(6)%X,nb2,bet_off,Y(nhb1+1:ntot1,nh2+1:nhb2),np1)
        end if

        RES%tblck(q)%tgam(2)%X = RES%tblck(q)%tgam(2)%X + &
             W(nhb1+1:ntot1,nh2+1:nhb2) + Y(nhb1+1:ntot1,nh2+1:nhb2)
     end if

!----------------------------------------------------------------------------
!         Zphpp 
!----------------------------------------------------------------------------
     if (nb1*np2 .ne. 0)  then 

        if (np1 .ne. 0 ) then 
           al_off = 1.d0* L%herm 
           call dgemm('T','N', nb1,np2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
                R%tblck(q)%tgam(1)%X,np1,bet_off,W(nh1+1:nhb1,nhb2+1:ntot2),nb1)
        end if
        
        al_off = -1.d0 
        call dgemm('N','N', nb1,np2,np2,al_off,R%tblck(q)%tgam(8)%X,nb1,&
             L%mat(q2)%gam(1)%X,np2,bet_off,W(nh1+1:nhb1,nhb2+1:ntot2),nb1)


        if(nh1 .ne. 0) then        
           al_off = -1.d0 
           call dgemm('N','N',nb1,np2,nh1,al_off,L%mat(q1)%gam(6)%X,nb1,&
                R%tblck(q)%tgam(7)%X,nh1,bet_off,Y(nh1+1:nhb1,nhb2+1:ntot2),nb1)
        end if

        if (nh2 .ne. 0 ) then 
           al_off = 1.d0 * L%herm
           call dgemm('N','T',nb1,np2,nh2,al_off,R%tblck(q)%tgam(6)%X,nb1,&
                L%mat(q2)%gam(3)%X,np2,bet_off,Y(nh1+1:nhb1,nhb2+1:ntot2),nb1)
        end if
        RES%tblck(q)%tgam(8)%X = RES%tblck(q)%tgam(8)%X + &
             W(nh1+1:nhb1,nhb2+1:ntot2) + Y(nh1+1:nhb1,nhb2+1:ntot2)
     end if
     

!----------------------------------------------------------------------------
!         Zphhh 
!----------------------------------------------------------------------------
     if (nb1*nh2 .ne. 0)  then 

        if (np1 .ne. 0 ) then 
           al_off = 1.d0* L%herm 
           call dgemm('T','N', nb1,nh2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
                R%tblck(q)%tgam(3)%X,np1,bet_off,W(nh1+1:nhb1,1:nh2),nb1)
        end if

        if (np2 .ne. 0 ) then 
           al_off = -1.d0 
           call dgemm('N','N', nb1,nh2,np2,al_off,R%tblck(q)%tgam(8)%X,nb1,&
                L%mat(q2)%gam(3)%X,np2,bet_off,W(nh1+1:nhb1,1:nh2),nb1)
        end if

        if(nh1 .ne. 0) then        
           al_off = -1.d0 
           call dgemm('N','N',nb1,nh2,nh1,al_off,L%mat(q1)%gam(6)%X,nb1,&
                R%tblck(q)%tgam(5)%X,nh1,bet_off,Y(nh1+1:nhb1,1:nh2),nb1)
        end if

        if (nh2 .ne. 0 ) then 
           al_off = 1.d0
           call dgemm('N','N',nb1,nh2,nh2,al_off,R%tblck(q)%tgam(6)%X,nb1,&
                L%mat(q2)%gam(5)%X,nh2,bet_off,Y(nh1+1:nhb1,1:nh2),nb1)
        end if
        RES%tblck(q)%tgam(6)%X = RES%tblck(q)%tgam(6)%X  + &
             W(nh1+1:nhb1,1:nh2) + Y(nh1+1:nhb1,1:nh2)
     end if
          

!----------------------------------------------------------------------------
!         Zhhph 
!----------------------------------------------------------------------------
     if (nh1*nb2 .ne. 0)  then 

        if (np1 .ne. 0 ) then 
           al_off = 1.d0*L%herm
           call dgemm('T','N', nh1,nb2,np1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(2)%X,np1,bet_off,W(1:nh1,nh2+1:nhb2),nh1)
        end if
        
        if (np2 .ne. 0 ) then 
           al_off = -1.d0 
           call dgemm('N','N', nh1,nb2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
                L%mat(q2)%gam(2)%X,np2,bet_off,W(1:nh1,nh2+1:nhb2),nh1)
        end if

        if(nh1 .ne. 0) then        
           al_off = -1.d0 
           call dgemm('N','N',nh1,nb2,nh1,al_off,L%mat(q1)%gam(5)%X,nh1,&
                R%tblck(q)%tgam(9)%X,nh1,bet_off,Y(1:nh1,nh2+1:nhb2),nh1)
        end if

        if (nh2 .ne. 0 ) then 
           al_off = 1.d0 * L%herm
           call dgemm('N','T',nh1,nb2,nh2,al_off,R%tblck(q)%tgam(5)%X,nh1,&
                L%mat(q2)%gam(6)%X,nb2,bet_off,Y(1:nh1,nh2+1:nhb2),nh1)
        end if

        RES%tblck(q)%tgam(9)%X = RES%tblck(q)%tgam(9)%X +&
             W(1:nh1,nh2+1:nhb2) + Y(1:nh1,nh2+1:nhb2)
     end if

     !========!
     !221 part!
     !========!

     W = W * sqrt((J1+1.d0)*(J2+1.d0)) 
     Y = -1*Y * sqrt((J1+1.d0)*(J2+1.d0))
!!!! current state of things:
!!!! I have not declared many of these variables
!!!! There is also an issue with the pp and hh terms, do I have to flip them and add them twice? probably.
          
     allocate(qn1(ntot1,2), qn2(ntot2,2))
     
     ! make one big two particle basis descriptor 

     qn1(1:nh1,:) = L%mat(q1)%qn(3)%Y(:,:)
     qn1(nh1+1:nhb1,:) = L%mat(q1)%qn(2)%Y(:,:)
     qn1(nhb1+1:ntot1,:) = L%mat(q1)%qn(1)%Y(:,:)

     qn2(1:nh2,:) = L%mat(q2)%qn(3)%Y(:,:)
     qn2(nh2+1:nhb2,:) = L%mat(q2)%qn(2)%Y(:,:)
     qn2(nhb2+1:ntot2,:) = L%mat(q2)%qn(1)%Y(:,:)
      


     ! this is the part of the 221 commutator where we sum over a single hole index, and two particles.
     ! the particle indices have already been summed up into W, we loop over all matrix elements of W
     ! and add their contributions to RES where ever possible.
          
     do IX = 1, nhb1
        pre1 = 1.d0
        
        ii = qn1(IX,1) 
        jj = qn1(IX,2)
        
        if (IX .le. nh1) then 

           !!  IX is a hh ket 
           !! this means I have to loop over
           !! all JX configs for each hole state in IX. 
           
           jp = jbas%jj(ii)
           ji = jbas%jj(jj)

           if (ii==jj) pre1 = sqrt(2.d0)
           ! ll is the summing hole
           p1 = ii
           i = jj

           pre1 = pre1 * (-1) ** ((ji-jp-J1)/2)                      
           do JX = 1, nhb2
              
              pre2 = 1.d0

              kk = qn2(JX,1)
              ll = qn2(JX,2)

              if (kk==ll) pre2 = sqrt(2.d0)

              !! figure out which index in the JX config
              !! corresponds the the index being summed (if any) 
              if ( kk == i ) then         
                 jh = jbas%jj(ll)
                 h1 = ll                 
              else if (ll == i )then
                 h1 = kk
                 jh = jbas%jj(kk)
                 pre2 = pre2 * (-1) ** ((ji-jh-J2)/2)
              else
                 cycle
              end if

              ! if one of those indices matches the current "summing hole"
              ! add it to the running sum (RES%fock(px,hx))
              ! despite naming conventions, px and hx are not particle and hole states
              ! this is copied from a different routine where they are. I don't care, get pissed. :) 
              
              px = p1
              hx = h1

              RES%fock(px,hx) = RES%fock(px,hx) + (-1) ** ((jh+ji + rank+J1)/2) &
                   *xxxsixj(RES%xindx,J1,J2,rank,jh,jp,ji)*W(IX,JX)*pre1*pre2

           end do

           
           ! kk is the summing hole
           p1 = jj
           i = ii            

           jp = jbas%jj(jj)
           ji = jbas%jj(ii)

           if (ii==jj) cycle! pre1 = sqrt(2.d0)
           pre1 = 1.d0         

           do JX = 1, nhb2
              
              pre2 = 1.d0

              kk = qn2(JX,1)
              ll = qn2(JX,2)

              if (kk==ll) pre2 = sqrt(2.d0)

              if ( kk == i ) then
                 ! k is a particle
                 jh = jbas%jj(ll)
                 h1 = ll                 
              else if (ll == i )then
                 h1 = kk
                 jh = jbas%jj(kk)
                 pre2 = pre2 * (-1) ** ((ji-jh-J2)/2)
              else                
                 cycle
              end if

              px = p1
              hx = h1
              
              RES%fock(px,hx) = RES%fock(px,hx) + (-1) ** ((jh+ji + rank+J1)/2) &
                   *xxxsixj(RES%xindx,J1,J2,rank,jh,jp,ji)*W(IX,JX)*pre1*pre2

           end do
           
           

        else

           !!  IX is ph config
           !! I only need to loop over JX for the hole state
           
           if ( jbas%con(ii) == 0 ) then
              ! ll is the summing hole
              
              jp = jbas%jj(ii)
              ji = jbas%jj(jj)

              p1 = ii
              i = jj

              pre1 = pre1 * (-1) ** ((ji-jp-J1)/2)

           else    
              ! kk is the summing hole
              p1 = jj
              i = ii              

              jp = jbas%jj(jj)
              ji = jbas%jj(ii)

           end if


           do JX = 1, nhb2
              
              pre2 = 1.d0

              kk = qn2(JX,1)
              ll = qn2(JX,2)

              if (kk==ll) pre2 = sqrt(2.d0)

              if ( kk == i ) then                 
                 jh = jbas%jj(ll)
                 h1 = ll                 
              else if (ll == i )then
                 h1 = kk
                 jh = jbas%jj(kk)
                 pre2 = pre2 * (-1) ** ((ji-jh-J2)/2)
              else
                 cycle
              end if

              px = p1
              hx = h1
              
              RES%fock(px,hx) = RES%fock(px,hx) + (-1) ** ((jh+ji + rank+J1)/2) &
                   *xxxsixj(RES%xindx,J1,J2,rank,jh,jp,ji)*W(IX,JX)*pre1*pre2

           end do

        end if
     end do

         
     do JX = nh2+1, ntot2
        pre2 = 1.d0

        ii = qn2(JX,1) 
        jj = qn2(JX,2)
        
        if ( JX .le. nhb2) then

           !! JX is a ph state. only need one IX loop. GOOD.
           
           if ( jbas%con(ii) == 1 ) then
              ! ii is a hole
              ! jj is the summing particle. 
              jh = jbas%jj(ii)
              ja = jbas%jj(jj)
              h1 = ii
              a = jj
              pre2 = pre2 * (-1) ** ((ja-jh-J2)/2)
           else
              !! ii is the summing particle. 
              h1 = jj
              a = ii
              jh = jbas%jj(jj)
              ja = jbas%jj(ii)
           end if

           do IX = nh1+1, ntot1

              pre1 = 1.d0

              kk = qn1(IX,1) 
              ll = qn1(IX,2)

              if (kk==ll) pre1 = sqrt(2.d0)

              if ( kk == a ) then
                 jp = jbas%jj(ll)
                 p1 = ll                 
              else if (ll == a ) then 
                 p1 = kk
                 jp = jbas%jj(kk)
                 pre1 = pre1 * (-1) ** ((ja-jp-J1)/2)
              else
                 cycle
              end if

              px = p1
              hx = h1

              RES%fock(px,hx) = RES%fock(px,hx) + (-1) ** ((jh+ja + rank+J1)/2) &
                   *xxxsixj(RES%xindx,J1,J2,rank,jh,jp,ja)*Y(IX,JX)*pre1*pre2


           end do

        else

           !! JX is a pp state. Need two IX loops. LAME.
          
           ! ll is the summing particle
           jh = jbas%jj(ii)
           ja = jbas%jj(jj)
           
           if (ii==jj) pre2 = sqrt(2.d0)
           h1 = ii
           a = jj
           pre2 = pre2 * (-1) ** ((ja-jh-J2)/2)
           
           do IX = nh1+1, ntot1

              pre1 = 1.d0

              kk = qn1(IX,1) 
              ll = qn1(IX,2)

              if (kk==ll) pre1 = sqrt(2.d0)

              if ( kk == a ) then
                 jp = jbas%jj(ll)
                 p1 = ll                 
              else if (ll == a ) then 
                 p1 = kk
                 jp = jbas%jj(kk)
                 pre1 = pre1 * (-1) ** ((ja-jp-J1)/2)
              else
                 cycle
              end if

              px = p1
              hx = h1

              RES%fock(px,hx) = RES%fock(px,hx) + (-1) ** ((jh+ja + rank+J1)/2) &
                   *xxxsixj(RES%xindx,J1,J2,rank,jh,jp,ja)*Y(IX,JX)*pre1*pre2

           end do


           !! kk is the summing particle now.
           h1 = jj
           a = ii
           jh = jbas%jj(jj)
           ja = jbas%jj(ii)
           pre2 = 1.d0 
           if (ii==jj) cycle
           do IX = nh1+1, ntot1

              pre1 = 1.d0

              kk = qn1(IX,1) 
              ll = qn1(IX,2)

              if (kk==ll) pre1 = sqrt(2.d0)

              if ( kk == a ) then
                 jp = jbas%jj(ll)
                 p1 = ll                 
              else if (ll == a ) then 
                 p1 = kk
                 jp = jbas%jj(kk)
                 pre1 = pre1 * (-1) ** ((ja-jp-J1)/2)
              else
                 cycle
              end if

              px = p1
              hx = h1

              RES%fock(px,hx) = RES%fock(px,hx) + (-1) ** ((jh+ja + rank+J1)/2) &
                   *xxxsixj(RES%xindx,J1,J2,rank,jh,jp,ja)*Y(IX,JX)*pre1*pre2
           end do

        end if
     end do


     

     
     deallocate(W,Y,qn1,qn2)
  end do

end subroutine operator_commutator_222_pp_hh
!=================================================================
!=================================================================
 subroutine operator_commutator_222_ph(L,R,RES,jbas) 
   ! VERIFIED ph channel 2body operator_commutator. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(iso_operator) :: RES,R
   type(sq_op) :: L
   real(8),allocatable,dimension(:,:) :: PANDYA_A,PANDYA_B,PANDYA_AB
   real(8),allocatable,dimension(:,:) :: PANDYA_revA,PANDYA_revB,PANDYA_revAB
   integer,allocatable,dimension(:,:) :: qn_J3,qn_J4,qn_J3_ph
   integer :: nh,np,nb1,nb2,q,IX,JX,r1,r2,Tz,PAR,JTM,q1,q2,J3,J4,rank,a,b,c,d
   integer :: ta,tb,tc,td,la,lb,lc,ld,n1,n2,c1,c2,jxstart,J4min,J4max,ja,jb,jc,jd
   integer :: J1,J2, Jtot,Ntot,qx,J3min,J3max,ril,rjk,rli,rkj,g_ix,thread,total_threads
   integer :: phase1,phase2,phase3,rik,rki,rjl,rlj,PAR2,J1min,J2min,J1max,J2max
   integer :: phase_34,phase_abcd,phase_ac,phase_bc,n_J3,n_J3_ph, n_J4,PAR_J3,PAR_J4,herm
   integer :: phase_bd,phase_ad,nj_perm,full_int_phase,Tz1_cc,Tz2_cc,omp_get_num_threads
   real(8) :: sm ,pre,pre2,omp_get_wtime ,t1,t2,coef9,factor,sm_ex, nj1,nj2  
   real(8) :: prefac_12,Xelem,Yelem,V,al_off,bet_off
   logical :: square
   
   rank = RES%rank
   ! construct intermediate matrices
   herm = L%herm * R%herm 
   Ntot = RES%Nsp
   JTM = jbas%Jtotal_max*2
   !$OMP PARALLEL
   total_threads = omp_get_num_threads()
   !$OMP END PARALLEL

   ! construct intermediate matrices

   !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(R,L,RES)   
   ! sum over all cross coupled channels for both matrices 

   do thread = 0,total_threads-1

      q = 1 
    
      do J3 = 2*thread, JTM  ,2*total_threads      
         J4min = abs(rank -J3)
         J4max = min(rank + J3,JTM)
         
         do Tz1_cc = -2,2,2 ! this is the cross coupled TZ
            Tz2_cc = abs(Tz1_cc - abs(R%dTz*2))
            
            do PAR_J3 = 0,1
                  
               n_J3 = count_iso_op_configs( J3 ,Tz1_cc, PAR_J3 , jbas, .false. ) !boolean true: ph, false: all configs
               n_J3_ph = count_iso_op_configs( J3 ,Tz1_cc, PAR_J3 , jbas, .true. ) 
               if (n_J3*n_J3_ph == 0) cycle

               allocate(qn_J3(n_J3,2),qn_J3_ph(n_J3_ph,2)) 
               n_J3 = count_iso_op_configs( J3 ,Tz1_cc, PAR_J3 , jbas, .false.,qn_J3) ! this fills the cc- basis descriptor
               n_J3_ph = count_iso_op_configs( J3 ,Tz1_cc, PAR_J3 , jbas, .true.,qn_J3_ph)

               allocate(PANDYA_A(n_J3, n_J3_ph) )
               allocate(PANDYA_revA(n_J3, n_J3_ph) )
               PANDYA_A = 0.d0
               PANDYA_revA = 0.d0

               call fill_rectangle_cc_matrix(J3,PANDYA_A,qn_J3,qn_J3_ph,L,jbas,.true.)  !boolean (true: sum over ph , false: hp) 
               call fill_rectangle_cc_matrix(J3,PANDYA_revA,qn_J3,qn_J3_ph,L,jbas,.false.)
               
               PAR_J4 = mod(PAR_J3 + R%dpar/2,2)

               do J4 = J4min,J4max ,2

                  n_J4 = count_iso_op_configs( J4 ,Tz2_cc, PAR_J4 , jbas, .false. )
                  if (n_J4 == 0) cycle
                  allocate(qn_J4(n_J4,2)) 
                  n_J4 = count_iso_op_configs( J4 ,Tz2_CC, PAR_J4 , jbas, .false.,qn_J4)
                  
                  allocate( PANDYA_B( n_J3_ph, n_J4) )
                  allocate( PANDYA_revB( n_J3_ph, n_J4) )

                  allocate( PANDYA_AB( n_J3, n_J4) )
                  allocate( PANDYA_revAB( n_J3, n_J4) )
                  
                  PANDYA_B = 0.d0
                  PANDYA_revB = 0.d0
                  
                  call fill_generalized_oppandya_matrix(J3,J4,PANDYA_B,qn_J3_ph,qn_J4,R,jbas,.true.)!boolean (true: sum over hp , false: ph)  
                  call fill_generalized_oppandya_matrix(J3,J4,PANDYA_revB,qn_J3_ph,qn_J4,R,jbas,.false.)
                  
                  PANDYA_AB = 0.d0
                  PANDYA_revAB = 0.d0 
                  al_off =  sqrt((J3+1.d0)*(J4+1.d0)) 

                  call dgemm('N','N',n_J3,n_J4,n_J3_ph,al_off,PANDYA_A,n_J3,PANDYA_B,n_J3_ph,bet,PANDYA_AB,n_J3)
                  call dgemm('N','N',n_J3,n_J4,n_J3_ph,al_off,PANDYA_revA,n_J3,PANDYA_revB,n_J3_ph,bet,PANDYA_revAB,n_J3) 

                  do JX = 1,n_J4

                     ! GET KET 
                     c = qn_J4(JX,1)
                     b = qn_J4(JX,2)

                     jc = jbas%jj(c)
                     lc = jbas%ll(c)
                     tc = jbas%itzp(c)

                     jb = jbas%jj(b)
                     lb = jbas%ll(b)
                     tb = jbas%itzp(b)

                     do IX = 1, n_J3 

                        ! GET BRA
                        a = qn_J3(IX,1)
                        d = qn_J3(IX,2)

                        jd = jbas%jj(d)
                        ld = jbas%ll(d)
                        td = jbas%itzp(d)

                        ja = jbas%jj(a)
                        la = jbas%ll(a)
                        ta = jbas%itzp(a)

                        if ( mod(la +lb,2).ne. mod(lc+ld+RES%dpar/2,2) ) cycle
                        Xelem = PANDYA_AB(IX,JX)-PANDYA_revAB(IX,JX)

                        if (abs(Xelem) > 1e-6) then

                        if ( (ta +tb)-2*R%dTz  ==  (tc+td) ) then 

                           
                        ! CALCULATE X CONTRIBUTIONS

                        J1min = abs(ja-jb)
                        J1max = ja+jb
                        
                        J2min = abs(jc-jd) 
                        J2max = jc+jd 
                        
                        ! these are the results of the Matmuls 
                        

                        phase_abcd= (-1)**((ja+jb+jc+jd)/2)
                        
                        

                           do J1 = J1min,J1max,2
                              if ((a==b).and.(mod(J1/2,2)==1)) cycle

                              do J2 = max(J2min,abs(rank-J1)),min(J2max,rank+J1),2 
                                 if ((c==d).and.(mod(J2/2,2)==1)) cycle

                                 nj1 = ninej(RES%xindx,ja,jd,J3,jb,jc,J4,J1,J2,rank)

                                 prefac_12 =  sqrt((J1+1.d0)*(J2+1.d0))
                                       
                                 if (b .ge. a) then 
                                    if (d .ge. c) then 
                                 
                                       ! CALCULATE V^{J1 J2}_{abcd} and V^{J1 J2}_{cdab}
                                       
                                       ! V^{J1 J2}_{abcd} 
                                       V = prefac_12* nj1 * (-1)**((ja+jb+J2+J3+J4)/2) * Xelem
                                       call add_elem_to_iso_op(V,a,b,c,d,J1,J2,RES,jbas) 
                                      
                                    else

                                       V = prefac_12* nj1 * (-1)**((ja-jb+jc+jd+J3+J4)/2) * Xelem 
                                       call add_elem_to_iso_op(V,a,b,d,c,J1,J2,RES,jbas)
                                    end if

                                 else

                                    if (d .ge. c) then
                                       
                                       V = -1*prefac_12* nj1 * (-1)**((J1+J2+J3+J4)/2) * Xelem 
                                       call add_elem_to_iso_op(V,b,a,c,d,J1,J2,RES,jbas)

                                    else
                                       
                                       V = prefac_12* nj1 * (-1)**((jc+jd+J1+J3+J4)/2) * Xelem 
                                       call add_elem_to_iso_op(V,b,a,d,c,J1,J2,RES,jbas)

                                    end if

                                 end if

                              end do
                           end do

                        end if
                                       
                     end if
                     end do
                  end do
                  q = q+ 1 
                  deallocate(PANDYA_B,PANDYA_AB,qn_J4)
                  deallocate(PANDYA_revB,PANDYA_revAB)   
               end do
               deallocate(PANDYA_A,qn_J3,qn_J3_ph)                
               deallocate(PANDYA_revA)                
            end do

         end do

      end do
   end do
   
 end subroutine operator_commutator_222_ph
!=================================================================
!=================================================================
real(8) function operator_commutator_223_single(L,R,ip,iq,ir,is,it,iu,jtot1,jtot2,Jpq,Jst,jbas)
  implicit none 
  
  integer,intent(in) :: ip,iq,ir,is,it,iu,jpq,jst
  integer :: a,b,c,d,Jx,Jy,Jz,J2,J1,J3,phase,rank,ia
  integer :: tp,tq,tr,ts,tt,tu,lp,lq,lr,ls,lt,lu,astart
  integer :: ja,jb,jc,jd,jp,jq,jr,js,jt,ju,jtot1,jtot2
  integer :: j1min,j1max , j2min,j2max,j3min,j3max
  type(sq_op) :: L,R
  type(spd) :: jbas
  real(8) :: sm,sm_sub,multfact,smtot,d6ji,out,otherfact,dsum
  real(8) :: Vs1,Vs2,sj1,sj2,sj3,sj4
  
  smtot = 0.d0 
  rank = R%rank 
  
  jp = jbas%jj(ip)
  jq = jbas%jj(iq)
  jr = jbas%jj(ir)  
  js = jbas%jj(is)
  jt = jbas%jj(it)
  ju = jbas%jj(iu)  

  tp = jbas%itzp(ip)
  tq = jbas%itzp(iq)
  tr = jbas%itzp(ir)  
  ts = jbas%itzp(is)
  tt = jbas%itzp(it)
  tu = jbas%itzp(iu)       

  lp = jbas%ll(ip)
  lq = jbas%ll(iq)
  lr = jbas%ll(ir)  
  ls = jbas%ll(is)
  lt = jbas%ll(it)
  lu = jbas%ll(iu)  

  ! FIRST TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jq+jr+jp-jtot2+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  ! Left-Right 
  sm = 0.d0   
  
  j2min = max(abs(jq - jr),abs(jp-jtot1))
  j2max = min(jq+jr,jp+jtot1)
      
  do a=1,jbas%total_orbits 
     

     if ( (tp+jbas%itzp(a)) .ne. (ts+tt) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 
     
     ja = jbas%jj(a) 
     
     if (.not. triangle(jp,ja,jst) ) cycle
                    
     phase = (-1)**((ja-ju)/2) 
     

     sj1 = v_elem(ip,a,is,it,Jst,L,jbas)*phase

     do J2 = j2min, j2max , 2
        
        sj2 = sj1*sixj(jp,jq,Jpq,jr,jtot1,J2)*sqrt(J2+1.d0) 
        
        j3min = min(abs(ja - ju),abs(jp-jtot2),abs(J2-rank))
        j3max = max(ja+ju,jp+jtot2,J2+rank) 
     
        do J3 = j3min,j3max,2

           sm = sm - (-1)**(J3/2) * sqrt(J3 + 1.d0) &            
            * sj2 * sixj(jp,ja,Jst,ju,jtot2,J3) *  &
             xxxsixj(R%xindx,J2,J3,rank,jtot2,jtot1,jp) * &
             tensor_elem(iq,ir,a,iu,J2,J3,R,jbas)
           
        end do
     end do
  end do
 
  smtot = smtot + sm*multfact
  
  !Right-left  

  multfact = (-1)**((jq+jr+jp-jtot2+rank)/2) *sqrt((Jpq+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) ) 
  sm = 0.d0 
  do a=1,jbas%total_orbits 
     
     
     if ( (tu+jbas%itzp(a)) .ne. (tq+tr) ) cycle
     if ( mod(lu+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 
     
     ja = jbas%jj(a) 
          
     j1min = min(abs(jp - ja),abs(Jst-rank))
     j1max = max(jp + ja,Jst+rank) 
     
     j3min = max( abs(jq - jr) , abs(ja - ju) , abs(jp-jtot1)) 
     j3max = min( jq+jr , ja+ju , jp+jtot1) 
     
     phase = (-1)**((ja - jp)/2)
             
     do J1 = j1min, j1max , 2      

        Vs1 = tensor_elem(ip,a,is,it,J1,Jst,R,jbas)*(-1)**(J1/2) 
        sj1 = sqrt(J1+1.d0)*xxxsixj(R%xindx,J1,Jst,rank,jtot2,jtot1,ju) 

        do J3 = j3min,j3max,2

           sm = sm +  phase *sj1*(J3 + 1.d0) &
                * sixj(jp,jq,Jpq,jr,jtot1,J3) * sixj(jp,ja,J1,ju,jtot1,J3) &
                * Vs1 * v_elem(iq,ir,a,iu,J3,L,jbas)

        end do
     end do
  end do 
  
  
  smtot = smtot + sm*multfact 
  
  !SECOND TERM
  !Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jp+jq+jr+jt+ju+jtot2+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 
  ! added a minus sign
  do a=1,jbas%total_orbits
     
     
     if ( (tp+jbas%itzp(a)) .ne. (tt+tu) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 
     
     
     ja = jbas%jj(a)

     j2min = max(abs(jp - ja) , abs(jt - ju) ,abs(js-jtot2)) 
     j2max = min(jp+ja,jt+ju,js+jtot2) 
     
     j1min = max(abs(ja - js),abs(jp-jtot2)) 
     j1max = min(ja+js,jp+jtot2) 
     
     phase = (-1) ** ((ja - js)/2) 
     
     do J1 = j1min,j1max,2
        
        sj1 = sqrt(J1+1.d0)*phase  
       
        do J2 = j2min,j2max,2
           
           sj2 = sj1 * (J2+1.d0) * (-1)**(J2/2) * sixj(js,jt,Jst,ju,jtot2,J2) * &
                sixj(jp,ja,J2,js,jtot2,J1) * v_elem(ip,a,it,iu,J2,L,jbas) 
     
           j3min = min(abs(jq - jr),abs(J1-rank)) 
           j3max = max(jq+jr,J1+rank) 
     
           do J3 = j3min,j3max,2
           
              sm = sm - sqrt(J3+1.d0) * sj2 * sixj(jp,jq,Jpq,jr,jtot1,J3) &
                   * xxxsixj(R%xindx,J1,J3,rank,jtot1,jtot2,jp) * (-1)**(J1/2) * &
                     tensor_elem(iq,ir,a,is,J3,J1,R,jbas)
           end do
        end do         
        
     end do
 
  end do


  do a=1,jbas%total_orbits
     
     
     if ( (ts+jbas%itzp(a)) .ne. (tq+tr) ) cycle
     if ( mod(ls+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 
     
     ja = jbas%jj(a)
     
     j1min = min(abs(jp - ja),abs(js-jtot1))
     j1max = max(jp+ja ,js+jtot1) 
     
     j3min = max( abs(jq - jr) , abs(ja - js) ,abs(jp-jtot1)) 
     j3max = min( jq+jr,ja+js,jp+jtot1)
     
     phase = (-1) ** ((ja - jp)/2) 
     
     do J1 = j1min,j1max,2
        
        sj1 = sqrt(J1+1.d0) * phase *(-1)**(J1/2)
        
        j2min = min(abs(jt - ju),abs(J1-rank)) 
        j2max = max(jt+ju ,J1+rank) 
        
        do J2 = j2min,j2max,2 
           
           sj2 = sj1 * sqrt(J2+1.d0) * sixj(js,jt,Jst,ju,jtot2,J2) * (-1)**(J2/2) *&
                xxxsixj(R%xindx,J2,J1,rank,jtot1,jtot2,js) * tensor_elem(ip,a,it,iu,J1,J2,R,jbas) 
           
           do J3 = j3min,j3max,2
              
              sm = sm + sj2 * (J3+1.d0) * sixj(jp,jq,Jpq,jr,jtot1,J3) &
                   * sixj(jp,ja,J1,js,jtot1,J3) * v_elem(iq,ir,a,is,J3,L,jbas)
           end do 
        end do 
     end do
 
  end do
     
  smtot = smtot + sm*multfact

  ! THIRD TERM    
  !Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jq+jr+jtot2+js+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
          
     if ( (tp+jbas%itzp(a)) .ne. (ts+tu) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(ls+lu,2) ) cycle 
     
     ja = jbas%jj(a)     
     
     j2min = max( abs(jp - ja) , abs(js - ju),abs(jt-jtot2) ) 
     j2max = min( jp+ja , js+ju,jt+jtot2) 
     
     j3min = max(abs(jq - jr) ,abs(jp-jtot1))
     j3max = min(jq+jr,jp+jtot1) 
   
     
     phase = (-1) ** ((ja + jp)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        
        sj1 = phase*sixj(jp,jq,Jpq,jr,jtot1,J3)*sqrt(J3+1.d0)  
        
        j1min = max(abs(ja - jt),abs(rank-J3))
        j1max = min(ja+jt,rank+J3) 
             
        do J1 = j1min,j1max,2 
           
           sj2 =  sj1*(-1)**(J1/2)*tensor_elem(iq,ir,a,it,J3,J1,R,jbas)*&
                xxxsixj(R%xindx,J1,J3,rank,jtot1,jtot2,jp)*sqrt(J1+1.d0)

           do J2 = j2min,j2max,2
              sm = sm - (J2+1.d0) *sj2* sixj(js,jt,Jst,jtot2,ju,J2) &
                   * sixj(jp,ja,J2,jt,jtot2,J1) * v_elem(ip,a,iu,is,J2,L,jbas)
           end do
        end do
     end do
 
  end do

  ! right-left

  do a=1,jbas%total_orbits
     
     
     if ( (tt+jbas%itzp(a)) .ne. (tq+tr) ) cycle
     if ( mod(lt+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 

     ja = jbas%jj(a)     
     
     j2min =  max(abs(js - ju),abs(jt-jtot2))
     j2max =  min(js+ju,jt+jtot2)
     
     j3min = max( abs(jq - jr) , abs(ja - jt) ,abs(jp-jtot1)) 
     j3max = min( jq+jr , ja+jt,jp+jtot1)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J2 = j2min,j2max,2
        
        sj1 = phase* sixj(js,jt,Jst,jtot2,ju,J2)*sqrt(J2+1.d0) 
    
        j1min = max(abs(jp - ja),abs(rank-J2))
        j1max = min(jp+ja ,rank+J2)
    
        do J1 = j1min,j1max,2 
           sj2 = sj1*(-1)**(J1/2)*xxxsixj(R%xindx,J1,J2,rank,jtot2,jtot1,jt) *sqrt(J1+1.d0) *&
                tensor_elem(ip,a,iu,is,J1,J2,R,jbas)      

           do J3 = j3min,j3max,2
              sm = sm + (J3+1.d0) *sj2* sixj(jp,jq,Jpq,jr,jtot1,J3) &
                   * sixj(jp,ja,J1,jt,jtot1,J3) * v_elem(iq,ir,a,it,J3,L,jbas)
           end do
        
        end do
        
     end do
 
  end do
     
  smtot = smtot + sm*multfact



!  FOURTH TERM    
!  Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jp+jtot2+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     
     if ( (tq+jbas%itzp(a)) .ne. (ts+tt) ) cycle
     if ( mod(lq+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 

     ja = jbas%jj(a)     
     
     
     j1min = max(abs(jp - jr),abs(jtot1-jq))
     j1max = min(jp+jr ,jtot1+jq) 
        
     phase = (-1) ** ((ja - ju)/2) ! changed to ja+js rather than ja-js 
     
     sj1 = v_elem(iq,a,is,it,Jst,L,jbas)*phase
     
     do J1 = j1min,j1max,2
        
        sj2 = sj1*sixj(jp,jq,Jpq,jtot1,jr,J1)*sqrt(J1+1.d0)*(-1)**(J1/2)
        
        j2min = max(abs(ja - ju),abs(rank-J1),abs(jq-jtot2))
        j2max = min(ja+ju,rank+J1,jq+jtot2)
        
        do J2 = j2min,j2max,2 
           
           sm = sm -  sj2*sqrt(J2+1.d0)*(-1)**(J2/2)*tensor_elem(ir,ip,a,iu,J1,J2,R,jbas)*&
                xxxsixj(R%xindx,J1,J2,rank,jtot2,jtot1,jq)*sixj(ja,jq,Jst,jtot2,ju,J2)

        end do
     end do

  end do

  smtot = smtot + sm*multfact

  ! right-left
      
  sm = 0.d0 
  
  multfact = (-1)**((jp-jtot2+Jpq+rank)/2) *sqrt((Jpq+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     

     if ( (tu+jbas%itzp(a)) .ne. (tr+tp) ) cycle
     if ( mod(lu+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

     ja = jbas%jj(a)     
          
     j1min = max(abs(jq - ja),abs(rank-Jst))
     j1max = min(jq+ja,rank+Jst) 
   
     j2min = max(abs(ja - ju),abs(jr-jp),abs(jq-jtot1)) 
     j2max = min(ja+ju,jr+jp,jq+jtot1) 
     
     phase = (-1) ** ((ja + jq)/2) ! changed to ja+js rather than ja-js 
          
     do J1 = j1min,j1max,2
        
        sj1 = phase*sqrt(J1+1.d0)* &
             tensor_elem(iq,a,is,it,J1,Jst,R,jbas)*(-1)**(J1/2)
        
        do J2 = j2min,j2max,2 
           
           sm = sm +  sj1*(J2+1.d0)*(-1)**(J2/2)*v_elem(ir,ip,a,iu,J2,L,jbas) *&
                xxxsixj(R%xindx,J1,Jst,rank,jtot2,jtot1,ju)*sixj(ja,jq,J1,jtot1,ju,J2) &
                *sixj(jp,jq,Jpq,jtot1,jr,J2)
        end do
     end do

  end do

  smtot = smtot + sm*multfact
 
  ! FIFTH TERM    
  !Left-Right
  sm = 0.d0 

  multfact = (-1)**((ju+jt+jtot2+js+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     
     
     if ( (tq+jbas%itzp(a)) .ne. (tt+tu) ) cycle
     if ( mod(lq+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 

     ja = jbas%jj(a)     
     
     j1min = max( abs(jq - ja) , abs(jt - ju) ) 
     j1max = min( jq+ja , jt+ju) 
     
     j3min = max(abs(ja - js),abs(jq-jtot2))
     j3max = min(ja+js,jq+jtot2) 
     
     phase = (-1) ** ((ja + jp)/2) ! changed to ja+js rather than ja-js 
     
     do J1 = j1min,j1max,2
        
        sj1 = phase*sixj(js,jt,Jst,ju,jtot2,J1)*(J1+1.d0)*(-1)**(J1/2)&
             *v_elem(iq,a,it,iu,J1,L,jbas) 
        
        do J3 = j3min,j3max,2 
           
           sj2 =  sj1*(-1)**(J3/2)*&
                sixj(jq,ja,J1,js,jtot2,J3)*sqrt(J3+1.d0)

           j2min = max(abs(jp - jr),abs(rank-J3),abs(jq-jtot1))
           j2max = min(jp+jr,rank+J3,jq+jtot1) 
   
           do J2 = j2min,j2max,2
              sm = sm - sqrt(J2+1.d0)*(-1)**(J2/2)*sj2* sixj(jp,jq,Jpq,jtot1,jr,J2) &
                   * xxxsixj(R%xindx,J3,J2,rank,jtot1,jtot2,jq) * tensor_elem(ir,ip,a,is,J2,J3,R,jbas)
           end do
        end do
     end do
 
  end do
  smtot = smtot + sm*multfact

  sm = 0.d0 
  ! right-left
  multfact = (-1)**((jp+jq+jtot2+ju+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     if ( (ts+jbas%itzp(a)) .ne. (tr+tp) ) cycle
     if ( mod(ls+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

     ja = jbas%jj(a)     
    
     
     j2min =  max(abs(jt - ju),abs(js-jtot2))
     j2max =  max(jt+ju,js+jtot2)
     
     j3min = max( abs(jp - jr) , abs(ja - js) ) 
     j3max = min( jp+jr , ja+js)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        
        sj1 = phase*(-1)**(J3/2)* sixj(jp,jq,Jpq,jtot1,jr,J3)*(J3+1.d0)*&
             v_elem(ir,ip,a,is,J3,L,jbas) 
        
        do J2 = j2min,j2max,2 
           sj2 = sj1*(-1)**(J2/2)*sixj(js,jt,Jst,ju,jtot2,J2)*sqrt(J2+1.d0)

           j1min = max(abs(jq - ja),abs(rank-J2),abs(js-jtot1)) 
           j1max = min(jq+ja,rank+J2,js+jtot1)
                    
           do J1 = j1min,j1max,2
              sm = sm + sqrt(J1+1.d0) *(-1)**(J1/2) *sj2* sixj(jq,ja,J1,js,jtot1,J3) &
                   * xxxsixj(R%xindx,J2,J1,rank,jtot1,jtot2,js) * tensor_elem(iq,a,it,iu,J1,J2,R,jbas)
           end do
        
        end do
        
     end do
 
  end do
     
  smtot = smtot + sm*multfact 

  ! SIXTH TERM    
  !Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jtot2+js+Jpq+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
        
     
          
     if ( (tq+jbas%itzp(a)) .ne. (tu+ts) ) cycle
     if ( mod(lq+jbas%ll(a),2) .ne. mod(lu+ls,2) ) cycle 

     ja = jbas%jj(a)     
     
     j1min = max( abs(jq - ja) , abs(js - ju) ,abs(jt-jtot2) ) 
     j1max = min( jq+ja , js+ju , jt+jtot2 ) 
     
     j2min = max(abs(jp - jr),abs(jq-jtot1))
     j2max = min(jp+jr ,jq+jtot1) 
       
     phase = (-1) ** ((ja - jp)/2) ! changed to ja+js rather than ja-js 
     
     do J1 = j1min,j1max,2
        
        sj1 = phase*sixj(js,jt,Jst,jtot2,ju,J1)*(J1+1.d0) &
             *v_elem(iq,a,iu,is,J1,L,jbas) 
        
        do J2 = j2min,j2max,2 
           
           sj2 =  sj1*(-1)**(J2/2)*&
                sixj(jp,jq,Jpq,jtot1,jr,J2)*sqrt(J2+1.d0)
           
           j3min = max(abs(ja - jt),abs(rank-J2)) 
           j3max = min(ja+jt,rank+J2)
     
           do J3 = j3min,j3max,2
              sm = sm - sqrt(J3+1.d0)*(-1)**(J3/2)*sj2* sixj(jq,ja,J1,jt,jtot2,J3) &
                   * xxxsixj(R%xindx,J3,J2,rank,jtot1,jtot2,jq) * tensor_elem(ir,ip,a,it,J2,J3,R,jbas)
           end do
        end do
     end do
 
  end do
  
  smtot = smtot + sm*multfact

  sm = 0.d0 
  ! right-left
  multfact = (-1)**((jp+jq+js+jtot2+Jpq+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     if ( (tt+jbas%itzp(a)) .ne. (tr+tp) ) cycle
     if ( mod(lt+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

     ja = jbas%jj(a)     
    
     j2min = max(abs(js - ju),abs(jt-jtot2)) 
     j2max = min(js+ju ,jt+jtot2) 
     
     j3min = max( abs(jp - jr) , abs(ja - jt) ) 
     j3max = min( jp+jr , ja+jt)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        
        sj1 = phase*(-1)**(J3/2) * sixj(jp,jq,Jpq,jtot1,jr,J3)*(J3+1.d0)*&
             v_elem(ir,ip,a,it,J3,L,jbas) 
        
        do J2 = j2min,j2max,2 

           sj2 = sj1*sixj(js,jt,Jst,jtot2,ju,J2)*sqrt(J2+1.d0)
                    
           j1min = max(abs(jq - ja),abs(rank-J2),abs(jt-jtot1))
           j1max = min(jq+ja,rank+J2,jt+jtot1)
           
           do J1 = j1min,j1max,2
              sm = sm + sqrt(J1+1.d0) *(-1)**(J1/2) *sj2* sixj(jq,ja,J1,jt,jtot1,J3) &
                   * xxxsixj(R%xindx,J2,J1,rank,jtot1,jtot2,jt) * tensor_elem(iq,a,iu,is,J1,J2,R,jbas)
           end do
        
        end do
        
     end do
 
  end do
  smtot = smtot + sm*multfact


  ! SEVENTH TERM
  
  ! Left-right
  sm = 0.d0
  multfact = (-1)**((Jpq+rank+jr-jtot2)/2) *sqrt((Jst+1.d0)*(jtot1+1.d0)*&
       (jtot2+1.d0)) 

  do a=1,jbas%total_orbits

     
     
     if ( (tr+jbas%itzp(a)) .ne. (ts+tt) ) cycle
     if ( mod(lr+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 
     
     ja = jbas%jj(a)

     if (.not. triangle(jr,ja,jst) ) cycle
     
     j3min = max(abs(ja-ju),abs(Jpq-rank),abs(jr-jtot2))
     j3max = min(ja+ju,Jpq+rank,jr+jtot2)  
     
     sj1 = v_elem(ir,a,is,it,Jst,L,jbas) * (-1)**((ja+ju)/2) 

     do J3=j3min,j3max,2 
        sm = sm - sj1* sixj(jr,ja,Jst,ju,jtot2,J3) * (-1)**(J3/2) * sqrt(J3+1.d0) &
             * xxxsixj(R%xindx,J3,Jpq,rank,jtot1,jtot2,jr) * tensor_elem(ip,iq,a,iu,Jpq,J3,R,jbas)
     end do 

  end do 
  
  smtot = smtot + sm*multfact

  ! right-left
  sm = 0.d0
  multfact = (-1)**((Jpq+rank)/2) *sqrt((Jpq+1.d0)*(jtot1+1.d0)*&
       (jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     if ( (tu+jbas%itzp(a)) .ne. (tp+tq) ) cycle
     if ( mod(lu+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 

     ja = jbas%jj(a)

     if (.not. triangle(ju,ja,Jpq) ) cycle
     
     j3min = max(abs(ja-jr),abs(rank-Jst),abs(ju-jtot1))
     j3max = min(ja+jr,rank+Jst,ju+jtot1)  
     
     sj1 = v_elem(ip,iq,a,iu,Jpq,L,jbas) * (-1)**((ja+jtot2)/2) 

     do J3=j3min,j3max,2 

        sm = sm + sj1* sixj(jr,ja,J3,ju,jtot1,Jpq) * sqrt(J3+1.d0) *(-1)**(J3/2) &
             * xxxsixj(R%xindx,J3,Jst,rank,jtot2,jtot1,ju) * tensor_elem(ir,a,is,it,J3,Jst,R,jbas)

     end do 

  end do 
  
  smtot = smtot + sm*multfact


  ! EIGHTH TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jt+jr+js+jtot2+Jpq+rank)/2) *sqrt((Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  ! Left-Right 
  sm = 0.d0   
  
  do a=1,jbas%total_orbits 

     
          
     if ( (tr+jbas%itzp(a)) .ne. (tt+tu) ) cycle
     if ( mod(lr+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 

     ja = jbas%jj(a) 
     
     j2min = max(abs(ja - js),abs(rank-Jpq),abs(jtot2-jr))
     j2max = min(ja+js ,rank+Jpq,jtot2+jr)
          
     j1min = max(abs(ja - jr),abs(jt-ju)) 
     j1max = min(ja+jr,jt+ju)  
          
     phase = (-1)**((ja+ju)/2) 
     
     do J1 = j1min, j1max , 2
        sj1 = phase*(-1)**(J1/2)*(J1+1.d0) * sixj(js,jt,Jst,ju,jtot2,J1)&
             * v_elem(ir,a,it,iu,J1,L,jbas) 
     
        do J2 = j2min,j2max,2
        
           sm = sm - (-1)**(J2/2) * sqrt(J2 + 1.d0) &            
            * sj1 * sixj(ja,js,J2,jtot2,jr,J1) *  &
             xxxsixj(R%xindx,Jpq,J2,rank,jtot2,jtot1,jr) * &
             tensor_elem(ip,iq,a,is,Jpq,J2,R,jbas)
        end do
     end do
  end do
 
  smtot = smtot + sm*multfact
  
  ! !Right-left  
  sm = 0.d0 
  multfact = (-1)**((jt+jtot2+Jpq+rank)/2) *sqrt((Jst+1.d0)*(Jpq+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  
  do a=1,jbas%total_orbits 
     
     if ( (ts+jbas%itzp(a)) .ne. (tp+tq) ) cycle
     if ( mod(ls+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 

     ja = jbas%jj(a) 
          
     if (.not. triangle(ja,js,Jpq)) cycle 
     
     j1min = max(abs(jr - ja),abs(js-jtot1))
     j1max = min(jr + ja  ,js+jtot1) 
     
     
     phase = (-1)**((ja - ju)/2)
           
     sj1 = phase * v_elem(ip,iq,a,is,Jpq,L,jbas)
     do J1 = j1min, j1max , 2      

        sj2 = sqrt(J1+1.d0)*(-1)**(J1/2)*sixj(jr,ja,J1,js,jtot1,Jpq) * sj1

        j2min = max(abs(jt - ju),abs(J1-rank),abs(js-jtot2))
        j2max = min(jt+ju,J1+rank,js+jtot2)
     
        do J2 = j2min,j2max,2

           sm = sm + sj2* sqrt(J2 + 1.d0) * (-1)**(J2/2) &
                * sixj(js,jt,Jst,ju,jtot2,J2) * xxxsixj(R%xindx,J1,J2,rank,jtot2,jtot1,js) &
                * tensor_elem(ir,a,it,iu,J1,J2,R,jbas)

        end do
     end do
  end do 
  
  
  smtot = smtot + sm*multfact
  
! NINTH TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jr+jtot2+Jpq+Jst+rank)/2) *sqrt((Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  ! Left-Right 
  sm = 0.d0   
  
  do a=1,jbas%total_orbits 
     
     
     if ( (tr+jbas%itzp(a)) .ne. (tu+ts) ) cycle
     if ( mod(lr+jbas%ll(a),2) .ne. mod(lu+ls,2) ) cycle 

     ja = jbas%jj(a) 
     
     j2min = max(abs(ja - jt),abs(rank-Jpq))
     j2max = min(ja+jt,rank+Jpq) 
          
     j1min = max(abs(ja - jr),abs(js-ju),abs(jt-jtot2)) 
     j1max = min(ja+jr,js+ju,jt+jtot2)  
          
     phase = (-1)**((ja-js)/2) 
     
     do J1 = j1min, j1max , 2
        sj1 = phase*(J1+1.d0) * sixj(js,jt,Jst,jtot2,ju,J1)&
             * v_elem(ir,a,iu,is,J1,L,jbas) 
        
        do J2 = j2min,j2max,2
        
           sm = sm - (-1)**(J2/2) * sqrt(J2 + 1.d0) &            
            * sj1 * sixj(ja,jr,J1,jtot2,jt,J2) *  &
             xxxsixj(R%xindx,Jpq,J2,rank,jtot2,jtot1,jr) * &
             tensor_elem(ip,iq,a,it,Jpq,J2,R,jbas)
        end do
        
     end do
  end do
 
  smtot = smtot + sm*multfact
  
  ! !Right-left  
  sm = 0.d0 
  multfact = (-1)**((jt-jtot2+Jpq+Jst+rank)/2) *sqrt((Jst+1.d0)*(Jpq+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  
  do a=1,jbas%total_orbits 
     
     if ( (tt+jbas%itzp(a)) .ne. (tp+tq) ) cycle
     if ( mod(lt+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 

     ja = jbas%jj(a) 
          
     if (.not. triangle(ja,jt,Jpq)) cycle 
     
     j1min = max(abs(jr - ja),abs(jt-jtot1))
     j1max = max(jr + ja,jt+jtot1) 
     
     
     phase = (-1)**((ja + js)/2)
           
     sj1 = phase * v_elem(ip,iq,a,it,Jpq,L,jbas)
     do J1 = j1min, j1max , 2      

        sj2 = sqrt(J1+1.d0)*(-1)**(J1/2)*sixj(jr,ja,J1,jt,jtot1,Jpq) * sj1

        j2min = max(abs(js - ju),abs(rank-J1),abs(jt-jtot2))
        j2max = min(js+ju,rank+J1,jt+jtot2)
     
        do J2 = j2min,j2max,2

           sm = sm + sj2* sqrt(J2 + 1.d0) &
                * sixj(js,jt,Jst,jtot2,ju,J2) * xxxsixj(R%xindx,J1,J2,rank,jtot2,jtot1,jt) &
                * tensor_elem(ir,a,iu,is,J1,J2,R,jbas)

        end do
     end do
  end do 
  
  
  smtot = smtot + sm*multfact
    

     
  operator_commutator_223_single = smtot

end function operator_commutator_223_single
!=====================================================================================
!=====================================================================================
end module 
  
  
  

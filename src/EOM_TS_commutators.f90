 module EOM_TS_commutators
  use cross_coupled
  ! tensor-scalar commutator functions 
  
  ! THE TENSOR MUST BE THE SECOND ARGUMENT
  
contains
!=========================================================
!=========================================================
subroutine EOM_TS_commutator_111(L,R,RES,jbas) 
  ! one body part of [L1,R1] 
  implicit none 
  
  type(sq_op) :: L,R,RES
  type(spd) :: jbas
  integer :: p,a,i,ji,hol,par
  real(8) :: sm 
  real(8),dimension(L%belowEF,L%belowEF) :: th1,th2
  real(8),dimension(L%Nsp-L%belowEF,L%belowEF) :: tb1,tb2
  real(8),dimension(L%Nsp-L%belowEF,L%Nsp-L%belowEF) :: tp1,tp2 
 
  hol = L%belowEf
  par = L%Nsp - hol

!dfph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  call dgemm('N','N',par,hol,par,al,L%fpp,par,R%fph,par,bet,tb2,par) 
  call dgemm('N','N',par,hol,hol,al,R%fph,par,L%fhh,hol,bet,tb1,par) 
  
  RES%fph = tb2 - tb1 
    
end subroutine
!=========================================================
!=========================================================
subroutine EOM_TS_commutator_121(L,R,RES,jbas) 
  ! onebody part of [L1,R2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: J1,J2,PAR,TZ,ji,ja,jp,jq,a,i,p,q,g
  integer :: ti,ta,tp,tq,li,la,lp,lq,ak,ik,pk,qk,rank
  real(8) :: sm,smx,smy,smx2,smy2,d6ji
  
  rank = R%rank 
 !dfph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1,L%nsp-L%belowEF
     
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
    
     do q = 1,L%belowEF
        
        qk = jbas%holes(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
       
        ! check if this state is allowed
        
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
        if (tq .ne. tp) cycle 
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
              ! sum over J_total
              do J2 = abs(ji - jq),ji+jq,2
                 do J1 = abs(ja - jp),ja+jp,2 

                    smx = smx + tensor_elem(ak,pk,ik,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)  &
                          * xxxsixj(R%xindx,J1,J2,rank,jq,jp,ji)
                    
                 end do
              end do
             
              sm = sm + L%fph(a,i) * L%herm*smx &
              *(-1)**( (rank + jq + ji )/2 )
       
           end do 
        end do 
        
        RES%fph(p,q) = RES%fph(p,q) + sm  
        
     end do 
  end do 
 
end subroutine             

!====================================================================
!====================================================================
subroutine  EOM_TS_commutator_211(LCC,R,RES,jbas) 
  ! onebody part of  - [R1,L2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: R,RES
  type(ex_cc_mat) :: LCC 
  integer :: J1,J2,PAR,TZ,ji,ja,jp,jq,a,i,p,q,g,JTM
  integer :: ti,ta,tp,tq,li,la,lp,lq,ak,ik,pk,qk,rank
  integer :: rai,rpq,rqp,qx,Ntot
  real(8) :: sm,smx,smy,smx2,smy2,d6ji,LXXX
  
  rank = R%rank 
  JTM = Jbas%Jtotal_max
  Ntot = R%nsp 

! dfph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1,R%nsp-R%belowEF
     
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = 1,R%belowEF
      
        qk = jbas%holes(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        ! check if this state is allowed
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
        if (tq .ne. tp) cycle 
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
              if (ti .ne. ta) cycle 
              if (.not. (triangle(ja,ji,rank))) cycle
               
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 

              Tz = abs(ta -ti)/2 
              if (abs(tp - tq) .ne. Tz*2)  cycle 
              PAR = mod(la+li,2) 
              if (mod(lp+lq,2) .ne. PAR) cycle                                 

              qx = rank/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
                
              rai = fetch_rval(ak,ik,R%belowEF,qx,LCC)
              rpq = fetch_rval(pk,qk,R%belowEF,qx,LCC)

              sm = sm - (-1)**(( jp + jq + rank)/2) * R%fph(a,i) &               
              * LCC%CCX(qx)%X(rpq,rai) / sqrt(rank + 1.d0 ) 

              ! the last (rank + 1) is divided out because
              ! the CC matrix elements are scaled by that, 
              ! which is wrong here. 
         
           end do 
        end do 
        
        RES%fph(p,q) = RES%fph(p,q) + sm  
        
     end do 
  end do 

end subroutine
!==================================================
!==================================================             
subroutine EOM_TS_commutator_122(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: q,IX,JX,nh,np,nb,i,J1,J2
  integer :: a,b,c,d,ja,jb,jc,jd,ji,g_ix,q_sp,i_sp
  integer :: ta,tb,tc,td,ti,la,lb,lc,ld,li,spec
  integer :: jxstart,jxend,ixend,c1,c2,n1,n2
  logical :: square
  real(8) ::  sm
  
  
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1)
     J2 = R%tblck(q)%Jpair(2) 
    
     do g_ix = 3,7,4 
   
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
               if (jbas%con(i_sp) .ne. jbas%con(a) ) cycle
               sm = sm + f_elem(a,i_sp,L,jbas)*tensor_elem(i_sp,b,c,d,J1,J2,R,jbas)

            end do 
              
            ! b is replaced
            q_sp = sp_block_index(jb,lb,tb,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               if (jbas%con(i_sp) .ne. jbas%con(b) ) cycle               
               sm = sm + f_elem(b,i_sp,L,jbas)*tensor_elem(a,i_sp,c,d,J1,J2,R,jbas) 
            end do 
            
            ! c is replaced
            q_sp = sp_block_index(jc,lc,tc,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               if (jbas%con(i_sp) .ne. jbas%con(c) ) cycle 
               sm = sm - f_elem(i_sp,c,L,jbas)*tensor_elem(a,b,i_sp,d,J1,J2,R,jbas)
            end do 
            
            ! d is replaced
            q_sp = sp_block_index(jd,ld,td,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               if (jbas%con(i_sp) .ne. jbas%con(d) ) cycle 
               sm = sm - f_elem(i_sp,d,L,jbas)*tensor_elem(a,b,c,i_sp,J1,J2,R,jbas)
            end do 
          
            sm = sm / sqrt(1.d0 + kron_del(a,b)) /sqrt(1.d0 + kron_del(c,d)) 
            
            if (g_ix ==7) sm = -1*sm*L%herm 
            RES%tblck(q)%tgam(g_ix)%X(IX,JX) = sm 


        end do
     end do
   
     end do 
  end do 

end subroutine           
!==================================================
!==================================================             
subroutine EOM_TS_commutator_212(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: q,IX,JX,nh,np,nb,i,J1,J2
  integer :: a,b,c,d,ja,jb,jc,jd,ji,g_ix,q_sp,i_sp
  integer :: ta,tb,tc,td,ti,modla,modlb,modlc,modld,modli,spec,rank
  integer :: jxstart,jxend,ixend,c1,c2,n1,n2,corrector
  logical :: square
  real(8) ::  sm,sm1,sm2,sm3,sm4,d6ji
  
  rank = R%rank
  
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1)
     J2 = R%tblck(q)%Jpair(2) 
    
     do g_ix = 3,7,4 
   
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
              
              sm1=0.d0
              sm2=0.d0
              if (jbas%con(i) .ne. jbas%con(a) ) then             
                 if (ti == ta) then
                    if (modli == modla ) then 
                       if (triangle(ji,ja,rank)) then  
                          
                          sm1 = sm1 - xxxsixj(R%xindx,J1,J2,rank,ji,ja,jb)&
                               *f_tensor_elem(a,i,R,jbas)*v_elem(i,b,c,d,J2,L,jbas)
                       end if
                    end if
                 end if
               
                 sm1 = sm1*(-1)**((ja + jb-J2)/2) 
               
                 if (ti == tb) then
                    if (modli == modlb ) then 
                       if (triangle(ji,jb,rank)) then  
                          
                          sm2 = sm2 + xxxsixj(R%xindx,J1,J2,rank,ji,jb,ja)&
                               *f_tensor_elem(b,i,R,jbas)*v_elem(i,a,c,d,J2,L,jbas)
                       
                       end if
                    end if
                 end if
               
                 sm2 = sm2*(-1)**((J1+J2)/2) 
              end if

              sm3=0.d0
              sm4=0.d0
    
              if (jbas%con(i) .ne. jbas%con(c) ) then 
                 if (ti == td) then
                    if (modli == modld ) then 
                       if (triangle(ji,jd,rank)) then  
                       
                          sm3 = sm3 + xxxsixj(R%xindx,J1,J2,rank,jd,ji,jc)&
                               *f_tensor_elem(i,d,R,jbas)*v_elem(a,b,c,i,J1,L,jbas)
                       
                       end if
                    end if
                 end if
               
                 sm3 = sm3*(-1)**((jc+jd-J1)/2) 
                 
                 if (ti == tc) then
                    if (modli == modlc ) then 
                       if (triangle(ji,jc,rank)) then  
                          
                          sm4 = sm4 -  xxxsixj(R%xindx,J1,J2,rank,jc,ji,jd)&
                               *f_tensor_elem(i,c,R,jbas)*v_elem(a,b,d,i,J1,L,jbas)
                       
                       end if
                    end if
                 end if
                 
                 sm4 = sm4*(-1)**((J2+J1)/2)
              end if 
              
              sm =  sm + (sm1+sm2+sm3+sm4) 
           end do 
           
           sm = sm * sqrt((J1+1.d0)*(J2+1.d0) / &
              (1.d0 + kron_del(a,b)) /(1.d0 + kron_del(c,d))) * (-1)**(rank/2) 

           if (g_ix==7) sm = sm *(-1)**((J1+J2+2)/2)*L%herm
           
           RES%tblck(q)%tgam(g_ix)%X(IX,JX) = RES%tblck(q)%tgam(g_ix)%X(IX,JX)  + sm 
           
        end do
     end do
   
     end do 
  end do 

end subroutine           
!===================================================================
!===================================================================
subroutine EOM_TS_commutator_221(w1,w2,pm,RES,jbas) 
  ! verified
  ! THIS NEEDS TO BE RUN AFTER 222_pp_hh 
  ! 222_pp_hh sets up the intermediary matrices (w1,w2) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES,w1,w2
  integer :: i,j,q,Abody,Ntot,nh,np,nb,a,c
  integer :: ik,jk,ck,ji,jj,ti,tj,li,lj,jc,J1,J2
  integer,intent(in) :: pm
  real(8) :: sm,sm1,sm2,d6ji
  
  Abody = w1%belowEF
  Ntot = w1%Nsp
 
  ! fph
  do ik = 1 , Ntot-Abody
     i = jbas%parts(ik) 
     ji = jbas%jj(i) 
     li = jbas%ll(i) 
     ti = jbas%itzp(i) 
    
     do jk = 1 , Abody

        j = jbas%holes(jk)       
        jj = jbas%jj(j) 
        if (.not. (triangle(jj,ji,w1%rank))) cycle
        lj = jbas%ll(j) 
        if (mod(li,2) .ne. mod(lj+w1%dpar/2,2))  cycle
        tj = jbas%itzp(j)
        if (tj .ne. ti) cycle 
              
        sm = 0.d0 
      
        do ck = 1, Abody
           c = jbas%holes(ck) 
           jc = jbas%jj(c)
           ! w1 matrix results from multiplying the pp channel
           sm1 = 0.d0 
           do J1 = abs(jc - ji),jc+ji,2
            
              ! NOTE: 
              ! THESE SUMS HAVE TO BE BROKEN UP SO the J on the left side is 
              ! smaller. I don't have the other matrix multiplication.
              do J2 = abs(jc - jj),min(J1-2,jc+jj),2
                 
                ! use w1, because it sums over the pp indices
                sm1 = sm1 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(RES%xindx,J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w1,jbas)*(-1)**(J1/2)
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
             
                ! use w1, because it sums over the pp indices
                sm1 = sm1 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(RES%xindx,J1,J2,w1%rank,jj,ji,jc) &
                     *tensor_elem(c,i,c,j,J1,J2,w1,jbas)*(-1)**(J1/2)
                
             end do

          end do
          sm = sm + sm1*(-1)**((jc+jj)/2)
        end do 
        

        do ck = 1, Ntot - Abody
           c = jbas%parts(ck) 
           jc = jbas%jj(c)
           sm2 = 0.d0
           do J1 = abs(jc - ji),jc+ji,2
             do J2 = abs(jc - jj),min(J1-2,jc+jj),2

                ! use w1, because it sums over the pp indices
                sm2 = sm2 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(RES%xindx,J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w2,jbas)*(-1)**(J1/2)
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
           
                ! use w1, because it sums over the pp indices
                sm2 = sm2 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(RES%xindx,J1,J2,w1%rank,jj,ji,jc) &
                     *tensor_elem(c,i,c,j,J1,J2,w2,jbas)*(-1)**(J1/2)
                
             end do
            
           end do 
           sm = sm + sm2 *  (-1)**((jc+jj)/2)
        end do 
     
        RES%fph(ik,jk) = RES%fph(ik,jk) + sm * (-1)**(w1%rank/2)  

     end do 
  end do       

end subroutine
!===================================================================
!===================================================================
subroutine EOM_TS_commutator_222_pp_hh(L,R,RES,w1,w2,jbas) 
  !VERIFIED
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  L,R,RES,w1,w2
  integer :: q,q1,q2,J1,J2,Tz,Par,phase,rank
  integer :: np1,nb1,nh1,np2,nb2,nh2,pm,i
  real(8) :: bet_off,al_off
  
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
     q2 = block_index(J2,Tz,mod(Par+R%Dpar/2,2)) 
     
     nh1 = R%tblck(q)%nhh1
     np1 = R%tblck(q)%npp1
     nb1 = R%tblck(q)%nph1
     nh2 = R%tblck(q)%nhh2
     np2 = R%tblck(q)%npp2
     nb2 = R%tblck(q)%nph2
       
     do i = 1, 9 
        if (allocated(w1%tblck(q)%tgam(i)%X)) then 
           w1%tblck(q)%tgam(i)%X=0.d0
           w2%tblck(q)%tgam(i)%X=0.d0
        end if 
     end do
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------
     if (np1*nh2 .ne. 0)  then 
        
        !w1pphh = Apppp.Bpphh 
 
        call dgemm('N','N',np1,nh2,np1,al,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%tgam(3)%X,np1,bet,w1%tblck(q)%tgam(3)%X,np1)
        
        !w2pphh = -Bpphh.Ahhhh 
        al_off = -1 
        call dgemm('N','N',np1,nh2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
             L%mat(q2)%gam(5)%X,nh2,bet,w2%tblck(q)%tgam(3)%X,np1)
             
     end if 
       
     RES%tblck(q)%tgam(3)%X = RES%tblck(q)%tgam(3)%X + &
          w1%tblck(q)%tgam(3)%X - w2%tblck(q)%tgam(3)%X

!----------------------------------------------------------------------------
!         Zhhpp 
!----------------------------------------------------------------------------
     if (np2*nh1 .ne. 0)  then 
     
        al_off = -1*L%herm
        !w1hhpp = -Bhhpp.Apppp 
        call dgemm('N','N',nh1,np2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
             L%mat(q2)%gam(1)%X,np2,bet,w1%tblck(q)%tgam(7)%X,nh1)
        
        
        al_off = L%herm 
        !w1pphh = Ahhhh.Bhhpp 
        call dgemm('N','N',nh1,np2,nh1,al_off,L%mat(q1)%gam(5)%X,nh1,&
             R%tblck(q)%tgam(7)%X,nh1,bet,w2%tblck(q)%tgam(7)%X,nh1)

     end if 
       
     RES%tblck(q)%tgam(7)%X = RES%tblck(q)%tgam(7)%X - &
          w1%tblck(q)%tgam(7)%X + w2%tblck(q)%tgam(7)%X

!----------------------------------------------------------------------------
!         Zppph 
!----------------------------------------------------------------------------
     if (np1*nb2 .ne. 0)  then 
     
     
        if (nh2 .ne. 0) then 
           !w2ppph = -Bpphh.Ahhph
           al_off = -1*L%herm
           call dgemm('N','T',np1,nb2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
                L%mat(q2)%gam(6)%X,nb2,bet,w2%tblck(q)%tgam(2)%X,np1)
        end if 
        
     end if 
       
!----------------------------------------------------------------------------
!         Zphpp 
!----------------------------------------------------------------------------
     if (nb1*np2 .ne. 0)  then 
        
        if (nh1 .ne. 0) then 
           al_off = -1*L%herm
           !w2phpp = Aphhh.Bhhpp
           call dgemm('N','N',nb1,np2,nh1,al_off,L%mat(q1)%gam(6)%X,nb1,&
                R%tblck(q)%tgam(7)%X,nh1,bet,w2%tblck(q)%tgam(8)%X,nb1)
        end if 
        
     end if 
       

!----------------------------------------------------------------------------
!         Zphhh 
!----------------------------------------------------------------------------
     if (nb1*nh2 .ne. 0)  then 
     
        if (np1 .ne. 0) then 
           !w1phhh = Aphpp.Bpphh
           al_off = L%herm
           call dgemm('T','N',nb1,nh2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
                R%tblck(q)%tgam(3)%X,np1,bet,w1%tblck(q)%tgam(6)%X,nb1)
        end if 
 
     end if 

!----------------------------------------------------------------------------
!         Zhhph 
!----------------------------------------------------------------------------
     if (nh1*nb2 .ne. 0)  then 
     
        if (np2 .ne. 0)  then 
           al_off = L%herm
           !w1hhph = -Bhhpp.Appph 
           call dgemm('N','N',nh1,nb2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
                L%mat(q2)%gam(2)%X,np2,bet,w1%tblck(q)%tgam(9)%X,nh1)
        end if
     end if 
  end do
  
end subroutine 
!=================================================================
!=================================================================
 subroutine EOM_TS_commutator_222_ph(LCC,RCC,RES,WCC,jbas) 
   ! VERIFIED ph channel 2body EOM_TS_commutator. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(sq_op) :: RES
   type(ex_pandya_mat) :: RCC,WCC
   type(ex_cc_mat) :: LCC
   integer :: nh,np,nb1,nb2,q,IX,JX,i,j,k,l,r1,r2,Tz,PAR,JTM,q1,q2,J3,J4,rank
   integer :: ji,jj,jk,jl,ti,tj,tk,tl,li,lj,lk,ll,n1,n2,c1,c2,jxstart,J4min,J4max
   integer :: J1,J2, Jtot,Ntot,qx,J3min,J3max,ril,rjk,rli,rkj,g_ix,thread,total_threads
   integer :: phase1,phase2,phase3,rik,rki,rjl,rlj,PAR2,f1,f2,Atot
   real(8) :: sm ,pre,pre2,omp_get_wtime ,t1,t2,coef9,factor,sm_ex
   logical :: square
   
  rank = RES%rank
  Ntot = RES%Nsp
  JTM = jbas%Jtotal_max
  Atot = RES%belowEF
  total_threads = size(RES%direct_omp) - 1
   ! construct intermediate matrices

  do q = 1,RCC%nblocks
     if (RCC%jval2(q) > jbas%jtotal_max*2) cycle

     nb2 = size(RCC%CCX(q)%X(1,:))
     nb1 = size(RCC%CCR(q)%X(:,1))
     r1 = size(RCC%CCX(q)%X(:,1))
     r2 = size(RCC%CCR(q)%X(1,:))      

     if (r1 * r2 == 0) cycle
      
     PAR = mod(q-1,2)
     Tz = mod((q-1)/2,2) 
         
     if (nb1 .ne. 0 ) then 
        J1 = RCC%Jval(q) 
        factor = 1.d0/sqrt(J1+1.d0)
        q1 = J1/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
        
        call dgemm('N','N',r1,r2,nb1,factor,LCC%CCX(q1)%X,r1,&
             RCC%CCR(q)%X,nb1,bet,WCC%CCX(q)%X,r1) 
        ! J1<= J2 
     end if
         
     if (nb2 .ne. 0 ) then 
        
        J2 = RCC%Jval2(q) 
        PAR2 = mod(PAR+RCC%dpar/2,2) 
        q2 = J2/2+1 + Tz*(JTM+1) + 2*PAR2*(JTM+1)
        factor = 1.d0/sqrt(J2+1.d0)
     
        call dgemm('N','T',r1,r2,nb2,factor,RCC%CCX(q)%X,r1,&
             LCC%CCX(q2)%X,r2,bet,WCC%CCR(q)%X,r1) 
      !  J1 >= J2 
     
     end if
     
  end do
  
!!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE), SHARED(RES,WCC)  
 ! do thread = 1, total_threads
  !   do q = 1+RES%direct_omp(thread),RES%direct_omp(thread+1) 

   do q = 1,RES%nblocks
      J1 = RES%tblck(q)%jpair(1)
      J2 = RES%tblck(q)%jpair(2)
               
       g_ix = 3
       
       ! figure out how big the array is
       n1 = size(RES%tblck(q)%tgam(g_ix)%X(:,1))
       n2 = size(RES%tblck(q)%tgam(g_ix)%X(1,:))
       if ((n1*n2) .ne. 0) then  
        
       ! read in information about which 
       ! array we are using from public arrays
       c1 = sea1(g_ix) 
       c2 = sea2(g_ix) 
       square = sqs(g_ix) 
       jxstart = jst(g_ix) 
        
       do  IX =  1, n1 
          pre = 1.d0 

          i = RES%tblck(q)%tensor_qn(c1,1)%Y(IX,1)
          j = RES%tblck(q)%tensor_qn(c1,1)%Y(IX,2)
          
          if (i == j )  pre  = .70710678118d0
          ji = jbas%jj(i) 
          jj = jbas%jj(j) 
          li = jbas%ll(i) 
          lj = jbas%ll(j)
          ti = jbas%itzp(i) 
          tj = jbas%itzp(j)
          
          do JX =1,n2
             pre2 = 1.d0 
             k = RES%tblck(q)%tensor_qn(c2,2)%Y(JX,1)
             l = RES%tblck(q)%tensor_qn(c2,2)%Y(JX,2)
             
             if (k == l )  pre2 = .70710678118d0
             jk = jbas%jj(k) 
             jl = jbas%jj(l) 
             lk = jbas%ll(k) 
             ll = jbas%ll(l)
             tk = jbas%itzp(k) 
             tl = jbas%itzp(l)
             
             phase1 = (-1) ** (( ji + jj + jk + jl )/2) 
             
             sm = 0.d0 
             sm_ex = 0.d0 
            
             J3min = abs(ji - jl) 
             J3max = ji + jl
            
             J4min = abs(jj - jk)
             J4max = jj + jk 
            
            
             Tz = abs(ti -tl)/2                         
             PAR = mod(li+ll,2) 
             
             if (mod(lk+lj+RCC%dpar/2,2) == PAR) then 
                if (abs(tk - tj) == Tz*2)  then 
             
                   do J3 = J3min,J3max,2
                      
                      q1 = block_index(J3,Tz,PAR)
                      
                      if (jbas%con(i)-jbas%con(l) > 0) then 
                         rli = fetch_rval(l,i,Atot,q1,RCC)
                      else
                         rli = fetch_rval(i,l,Atot,q1,RCC)
                      end if 
                      
                      do J4 = max( J3 , J4min ) , J4max,2 
                         
                         if (.not. (triangle(J3,J4,rank))) cycle
                  
                         PAR2 = mod(PAR + RCC%dpar/2,2) 
                         q2 = block_index(J4,Tz,PAR2)
               
                         if (jbas%con(j)-jbas%con(k) > 0) then 
                            rkj = fetch_rval(k,j,Atot,q2,RCC)
                         else
                            rkj = fetch_rval(j,k,Atot,q2,RCC)
                         end if

              
                         qx = CCtensor_block_index(J3,J4,rank,Tz,PAR)
                         sm = sm + sqrt((J3+1.d0)*(J4+1.d0))* &
                              ninej(RES%xindx,ji,jl,J3,jj,jk,J4,J1,J2,rank)  * ( &
                              (-1)**((J3+J4+2)/2) * WCC%CCX(qx)%X(rli,rkj) &
                              + (-1)**((rank+2)/2) * WCC%CCR(qx)%X(rli,rkj) &
                              )

                      end do 

                      do J4 = J4min , min(J4max,J3-2),2 
                         if (.not. (triangle(J3,J4,rank))) cycle


                         PAR2 = mod(PAR + RCC%dpar/2,2)                  
                         q2 = block_index(J4,Tz,PAR2)
                         
                         if (jbas%con(j)-jbas%con(k) > 0) then 
                            rkj = fetch_rval(k,j,Atot,q2,RCC)
                         else
                            rkj = fetch_rval(j,k,Atot,q2,RCC)
                         end if

                         qx = CCtensor_block_index(J4,J3,rank,Tz,PAR2)

                         sm = sm + sqrt((J3+1.d0)*(J4+1.d0))* &
                              ninej(RES%xindx,ji,jl,J3,jj,jk,J4,J1,J2,rank)  * (  &
                              (-1)**((J3+J4+2)/2) * WCC%CCR(qx)%X(rkj,rli) &
                              + (-1)**((rank+2)/2) * WCC%CCX(qx)%X(rkj,rli) &
                              )

                      end do


                   end do 

                end if 
             end if 

             ! exchange of 1 set of indeces
             J3min = abs(jj - jl) 
             J3max = jj + jl

             J4min = abs(ji - jk)
             J4max = ji + jk 

             Tz = abs(tl -tj)/2 
             PAR = mod(ll+lj,2) 

             if (mod(li+lk+RCC%dpar/2,2) == PAR) then 
                if (abs(ti - tk) == Tz*2)  then 

                   do J3 = J3min,J3max,2
                      q1 = block_index(J3,Tz,PAR)
                      
                      if (jbas%con(j)-jbas%con(l) > 0) then 
                         rlj = fetch_rval(l,j,Atot,q1,RCC)
                      else
                         rlj = fetch_rval(j,l,Atot,q1,RCC)
                      end if
                      
                      do J4 = max( J3 , J4min ) , J4max,2 

                         if (.not. (triangle(J3,J4,rank))) cycle


                         PAR2 = mod(PAR + RCC%dpar/2,2)
                         q2 = block_index(J4,Tz,PAR2)
                         
                         if (jbas%con(i)-jbas%con(k) > 0) then 
                            rki = fetch_rval(k,i,Atot,q2,RCC)
                         else
                            rki = fetch_rval(i,k,Atot,q2,RCC)
                         end if

                         qx = CCtensor_block_index(J3,J4,rank,Tz,PAR)

                         sm_ex = sm_ex - sqrt((J3+1.d0)*(J4+1.d0))* &
                              ninej(RES%xindx,jj,jl,J3,ji,jk,J4,J1,J2,rank) * ( &
                              (-1)**((J3+J4+2)/2) * WCC%CCX(qx)%X(rlj,rki) &
                              + (-1)**((rank+2)/2) * WCC%CCR(qx)%X(rlj,rki) &       
                              )

                      end do

                      do J4 = J4min , min(J4max,J3-2),2 
                         if (.not. (triangle(J3,J4,rank))) cycle

                         PAR2 = mod(PAR + RCC%dpar/2,2)
                         q2 = block_index(J4,Tz,PAR2)

                         
                         if (jbas%con(i)-jbas%con(k) > 0) then 
                            rki = fetch_rval(k,i,Atot,q2,RCC)
                         else
                            rki = fetch_rval(i,k,Atot,q2,RCC)
                         end if

                         qx = CCtensor_block_index(J4,J3,rank,Tz,PAR2)

                         sm_ex = sm_ex - sqrt((J3+1.d0)*(J4+1.d0))* &
                              ninej(RES%xindx,jj,jl,J3,ji,jk,J4,J1,J2,rank) * ( &
                              (-1)**((J3+J4+2)/2) * WCC%CCR(qx)%X(rki,rlj) &
                              + (-1)**((rank+2)/2) * WCC%CCX(qx)%X(rki,rlj) &                  
                              ) 

                      end do


                   end do
               end if 
            end if 
           
          
          RES%tblck(q)%tgam(g_ix)%X(IX,JX) = RES%tblck(q)%tgam(g_ix)%X(IX,JX) + ( sm * &
               (-1) ** ((ji+jj+J2)/2) + sm_ex * (-1)**((J1+J2)/2) )&
               *   pre * pre2 *   sqrt((J1+1.d0)*(J2+1.d0))
   
         end do 
      end do

      end if
      g_ix = 7
        
      ! figure out how big the array is
      n1 = size(RES%tblck(q)%tgam(g_ix)%X(:,1))
      n2 = size(RES%tblck(q)%tgam(g_ix)%X(1,:))
      if ((n1*n2) == 0) cycle 
        
      ! read in information about which 
      ! array we are using from public arrays
      c1 = sea1(g_ix) 
      c2 = sea2(g_ix) 
      square = sqs(g_ix) 
      jxstart = jst(g_ix) 
       
      do  IX =  1, n1 
         pre = 1.d0 
         
         k = RES%tblck(q)%tensor_qn(c1,1)%Y(IX,1)
         l = RES%tblck(q)%tensor_qn(c1,1)%Y(IX,2)
          
          if (k == l )  pre  = .70710678118d0
          jk = jbas%jj(k) 
          jl = jbas%jj(l) 
          lk = jbas%ll(k) 
          ll = jbas%ll(l)
          tk = jbas%itzp(k) 
          tl = jbas%itzp(l)
          
          do JX =1,n2
             pre2 = 1.d0 
             i = RES%tblck(q)%tensor_qn(c2,2)%Y(JX,1)
             j = RES%tblck(q)%tensor_qn(c2,2)%Y(JX,2)
             
             if (i == j )  pre2 = .70710678118d0
             ji = jbas%jj(i) 
             jj = jbas%jj(j) 
             li = jbas%ll(i) 
             lj = jbas%ll(j)
             ti = jbas%itzp(i) 
             tj = jbas%itzp(j)
             
             phase1 = (-1) ** (( ji + jj + jk + jl )/2) 

             sm = 0.d0 
             sm_ex = 0.d0 
            
             J3min = abs(ji - jl) 
             J3max = ji + jl
            
             J4min = abs(jj - jk)
             J4max = jj + jk 
            
            
             Tz = abs(ti -tl)/2                         
             PAR = mod(li+ll,2) 

             
             if (mod(lk+lj+RCC%dpar/2,2) == PAR) then 
                if (abs(tk - tj) == Tz*2)  then 

                   do J3 = J3min,J3max,2
                      
                      q1 = block_index(J3,Tz,PAR)
                      
                      if (jbas%con(i)-jbas%con(l) > 0) then 
                         rli = fetch_rval(l,i,Atot,q1,RCC)
                      else
                         rli = fetch_rval(i,l,Atot,q1,RCC)
                      end if 
                      
                      do J4 = max( J3 , J4min ) , J4max,2 
                         
                         if (.not. (triangle(J3,J4,rank))) cycle
                  
                         PAR2 = mod(PAR + RCC%dpar/2,2) 
                         q2 = block_index(J4,Tz,PAR2)

                         if (jbas%con(j)-jbas%con(k) > 0) then 
                            rkj = fetch_rval(k,j,Atot,q2,RCC)
                         else
                            rkj = fetch_rval(j,k,Atot,q2,RCC)
                         end if
              
                         qx = CCtensor_block_index(J3,J4,rank,Tz,PAR)
                         sm = sm + sqrt((J3+1.d0)*(J4+1.d0))* &
                              ninej(RES%xindx,ji,jl,J3,jj,jk,J4,J2,J1,rank)  * ( &
                              (-1)**((J3+J4+2)/2) * WCC%CCX(qx)%X(rli,rkj) &
                              + (-1)**((rank+2)/2) * WCC%CCR(qx)%X(rli,rkj) &
                              )

                      end do 

                      do J4 = J4min , min(J4max,J3-2),2 
                         if (.not. (triangle(J3,J4,rank))) cycle


                         PAR2 = mod(PAR + RCC%dpar/2,2)                  
                         q2 = block_index(J4,Tz,PAR2)

                         if (jbas%con(j)-jbas%con(k) > 0) then 
                            rkj = fetch_rval(k,j,Atot,q2,RCC)
                         else
                            rkj = fetch_rval(j,k,Atot,q2,RCC)
                         end if

                         qx = CCtensor_block_index(J4,J3,rank,Tz,PAR2)

                         sm = sm + sqrt((J3+1.d0)*(J4+1.d0))* &
                              ninej(RES%xindx,ji,jl,J3,jj,jk,J4,J2,J1,rank)  * (  &
                              (-1)**((J3+J4+2)/2) * WCC%CCR(qx)%X(rkj,rli) &
                              + (-1)**((rank+2)/2) * WCC%CCX(qx)%X(rkj,rli) &
                              )
                       
                      end do


                   end do 

                end if 
             end if 

             ! exchange of 1 set of indeces
             J3min = abs(jj - jl) 
             J3max = jj + jl

             J4min = abs(ji - jk)
             J4max = ji + jk 

             Tz = abs(tl -tj)/2 
             PAR = mod(ll+lj,2) 

             if (mod(li+lk+RCC%dpar/2,2) == PAR) then 
                if (abs(ti - tk) == Tz*2)  then 

                   do J3 = J3min,J3max,2
                      q1 = block_index(J3,Tz,PAR)
                      
                      if (jbas%con(j)-jbas%con(l) > 0) then 
                         rlj = fetch_rval(l,j,Atot,q1,RCC)
                      else
                         rlj = fetch_rval(j,l,Atot,q1,RCC)
                      end if
                      
                      do J4 = max( J3 , J4min ) , J4max,2 

                         if (.not. (triangle(J3,J4,rank))) cycle


                         PAR2 = mod(PAR + RCC%dpar/2,2)
                         q2 = block_index(J4,Tz,PAR2)

                         if (jbas%con(i)-jbas%con(k) > 0) then 
                            rki = fetch_rval(k,i,Atot,q2,RCC)
                         else
                            rki = fetch_rval(i,k,Atot,q2,RCC)
                         end if

                         qx = CCtensor_block_index(J3,J4,rank,Tz,PAR)

                         sm_ex = sm_ex - sqrt((J3+1.d0)*(J4+1.d0))* &
                              ninej(RES%xindx,jj,jl,J3,ji,jk,J4,J2,J1,rank) * ( &
                              (-1)**((J3+J4+2)/2) * WCC%CCX(qx)%X(rlj,rki) &
                              + (-1)**((rank+2)/2) * WCC%CCR(qx)%X(rlj,rki) &       
                              )
                        
                      end do

                      do J4 = J4min , min(J4max,J3-2),2 
                         if (.not. (triangle(J3,J4,rank))) cycle

                         PAR2 = mod(PAR + RCC%dpar/2,2)
                         q2 = block_index(J4,Tz,PAR2)


                         if (jbas%con(i)-jbas%con(k) > 0) then 
                            rki = fetch_rval(k,i,Atot,q2,RCC)
                         else
                            rki = fetch_rval(i,k,Atot,q2,RCC)
                         end if                         
                         qx = CCtensor_block_index(J4,J3,rank,Tz,PAR2)

                         sm_ex = sm_ex - sqrt((J3+1.d0)*(J4+1.d0))* &
                              ninej(RES%xindx,jj,jl,J3,ji,jk,J4,J2,J1,rank) * ( &
                              (-1)**((J3+J4+2)/2) * WCC%CCR(qx)%X(rki,rlj) &
                              + (-1)**((rank+2)/2) * WCC%CCX(qx)%X(rki,rlj) &                  
                              ) 
                          
                      end do


                   end do
               end if 
            end if 
          
          RES%tblck(q)%tgam(g_ix)%X(IX,JX) = RES%tblck(q)%tgam(g_ix)%X(IX,JX) + ( sm * &
               (-1) ** ((ji+jj+J1)/2) + sm_ex * (-1)**((J1+J2)/2) )&
               *   pre * pre2 *   sqrt((J1+1.d0)*(J2+1.d0))
   
         end do 
      end do


!      end do 
 !  end do

   end do 
!!$OMP END PARALLEL DO 
   
end subroutine 
!=====================================================
!=====================================================      
real(8) function EOM_TS_commutator_223_single(L,R,ip,iq,ir,is,it,iu,jtot1,jtot2,Jpq,Jst,jbas)
  implicit none 
  
  integer,intent(in) :: ip,iq,ir,is,it,iu,jpq,jst
  integer :: a,b,c,d,Jx,Jy,Jz,J2,J1,J3,phase,rank,ia
  integer :: tp,tq,tr,ts,tt,tu,lp,lq,lr,ls,lt,lu,astart
  integer :: ja,jb,jc,jd,jp,jq,jr,js,jt,ju,jtot1,jtot2
  integer :: j1min,j1max , j2min,j2max,j3min,j3max
  type(sq_op) :: L,R
  type(spd) :: jbas
  real(8) :: sm,sm_sub,multfact,smtot,d6ji,out,otherfact,dsum
  real(8) :: Vs1,Vs2,sj1,sj2,sj3,sj4,stemp
  
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
      
  do ia=1,L%belowEF 
     a = jbas%holes(ia)

     if ( (tp+jbas%itzp(a)) .ne. (ts+tt) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 
     
     ja = jbas%jj(a) 
     
     if (.not. triangle(jp,ja,jst) ) cycle
                    
     phase = (-1)**((ja-ju)/2) 
     

     sj1 = v_elem(ip,a,is,it,Jst,L,jbas)*phase
     
     if (abs(sj1) < 1e-8) cycle
     
     do J2 = j2min, j2max , 2
        if ((iq == ir) .and. (mod(J2/2,2)==1)) cycle
        sj2 = sj1*sixj(jp,jq,Jpq,jr,jtot1,J2)*sqrt(J2+1.d0) 
        if (abs(sj2) < 1e-8) cycle        
        
        j3min = max(abs(ja - ju),abs(jp-jtot2),abs(J2-rank))
        j3max = min(ja+ju,jp+jtot2,J2+rank) 
     
        do J3 = j3min,j3max,2
           if ((a == iu) .and. (mod(J3/2,2)==1)) cycle
            
           sm = sm - (-1)**(J3/2) * sqrt(J3 + 1.d0) &            
            * sj2 * sixj(jp,ja,Jst,ju,jtot2,J3) *  &
             xxxsixj(R%xindx,J2,J3,rank,jtot2,jtot1,jp) * &
             tensor_elem(iq,ir,a,iu,J2,J3,R,jbas)
              
        end do
     end do
  end do
 
  smtot = smtot + sm*multfact
  
  !Right-left  
  if (.not.((is == it) .and. (mod(Jst/2,2)==1))) then 
     multfact = (-1)**((jq+jr+jp-jtot2+rank)/2) *sqrt((Jpq+1.d0) &
          *(jtot1+1.d0)*(jtot2+1.d0) ) 
     sm = 0.d0 
     do ia=1,L%Nsp-L%belowEF 
        a = jbas%parts(ia)

        if ( (tu+jbas%itzp(a)) .ne. (tq+tr) ) cycle
        if ( mod(lu+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 

        ja = jbas%jj(a) 

        j1min = max(abs(jp - ja),abs(Jst-rank))
        j1max = min(jp + ja,Jst+rank) 

        j3min = max( abs(jq - jr) , abs(ja - ju) , abs(jp-jtot1)) 
        j3max = min( jq+jr , ja+ju , jp+jtot1) 

        phase = (-1)**((ja - jp)/2)

        do J1 = j1min, j1max , 2      
           if ((ip == a) .and. (mod(J1/2,2)==1)) cycle     

           Vs1 = tensor_elem(ip,a,is,it,J1,Jst,R,jbas)*(-1)**(J1/2) 
           sj1 = sqrt(J1+1.d0)*xxxsixj(R%xindx,J1,Jst,rank,jtot2,jtot1,ju) 
           if (abs(sj1) < 1e-8) cycle
           do J3 = j3min,j3max,2
              if ((iq == ir) .and. (mod(J3/2,2)==1)) cycle
              if ((a == iu) .and. (mod(J3/2,2)==1)) cycle
              sm = sm +  phase *sj1*(J3 + 1.d0) &
                   * sixj(jp,jq,Jpq,jr,jtot1,J3) * sixj(jp,ja,J1,ju,jtot1,J3) &
                   * Vs1 * v_elem(iq,ir,a,iu,J3,L,jbas)

           end do
        end do
     end do


     smtot = smtot + sm*multfact 
  end if
  !SECOND TERM
  !Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jp+jq+jr+jt+ju+jtot2+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 
  ! added a minus sign
  do ia=1,L%belowEF
     a = jbas%holes(ia)
     
     if ( (tp+jbas%itzp(a)) .ne. (tt+tu) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 
     
     
     ja = jbas%jj(a)

     j2min = max(abs(jp - ja) , abs(jt - ju) ,abs(js-jtot2)) 
     j2max = min(jp+ja,jt+ju,js+jtot2) 
     
     j1min = max(abs(ja - js),abs(jp-jtot2)) 
     j1max = min(ja+js,jp+jtot2) 
     
     phase = (-1) ** ((ja - js)/2) 
     
     do J1 = j1min,j1max,2
        if ((a == is) .and. (mod(J1/2,2)==1)) cycle
        sj1 = sqrt(J1+1.d0)*phase  
        
        do J2 = j2min,j2max,2
           if ((ip == a) .and. (mod(J2/2,2)==1)) cycle
           if ((it == iu) .and. (mod(J2/2,2)==1)) cycle

           sj2 = sj1 * (J2+1.d0) * (-1)**(J2/2) * sixj(js,jt,Jst,ju,jtot2,J2) * &
                sixj(jp,ja,J2,js,jtot2,J1) * v_elem(ip,a,it,iu,J2,L,jbas) 
           if (abs(sj2) < 1e-8) cycle
           j3min = max(abs(jq - jr),abs(J1-rank)) 
           j3max = min(jq+jr,J1+rank) 
     
           do J3 = j3min,j3max,2
              if ((iq == ir) .and. (mod(J3/2,2)==1)) cycle
              sm = sm - sqrt(J3+1.d0) * sj2 * sixj(jp,jq,Jpq,jr,jtot1,J3) &
                   * xxxsixj(R%xindx,J1,J3,rank,jtot1,jtot2,jp) * (-1)**(J1/2) * &
                     tensor_elem(iq,ir,a,is,J3,J1,R,jbas)
           end do
        end do         
        
     end do
 
  end do


  do ia=1,L%Nsp-L%belowEF
     a = jbas%parts(ia)
     
     if ( (ts+jbas%itzp(a)) .ne. (tq+tr) ) cycle
     if ( mod(ls+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 
     
     ja = jbas%jj(a)
     
     j1min = max(abs(jp - ja),abs(js-jtot1))
     j1max = min(jp+ja ,js+jtot1) 
     
     j3min = max( abs(jq - jr) , abs(ja - js) ,abs(jp-jtot1)) 
     j3max = min( jq+jr,ja+js,jp+jtot1)
     
     phase = (-1) ** ((ja - jp)/2) 
     
     do J1 = j1min,j1max,2
        if ((ip == a ) .and. (mod(J1/2,2)==1)) cycle
        sj1 = sqrt(J1+1.d0) * phase *(-1)**(J1/2)
        
        j2min = max(abs(jt - ju),abs(J1-rank)) 
        j2max = min(jt+ju ,J1+rank) 
        
        do J2 = j2min,j2max,2 
           if ((it == iu) .and. (mod(J2/2,2)==1)) cycle
           sj2 = sj1 * sqrt(J2+1.d0) * sixj(js,jt,Jst,ju,jtot2,J2) * (-1)**(J2/2) *&
                xxxsixj(R%xindx,J2,J1,rank,jtot1,jtot2,js) * tensor_elem(ip,a,it,iu,J1,J2,R,jbas) 
           if (abs(sj2) < 1e-8) cycle
           do J3 = j3min,j3max,2
              if ((iq == ir) .and. (mod(J3/2,2)==1)) cycle
              if ((a == is) .and. (mod(J3/2,2)==1)) cycle
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

  do ia=1,L%belowEF
     a = jbas%holes(ia)
          
     if ( (tp+jbas%itzp(a)) .ne. (ts+tu) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(ls+lu,2) ) cycle 
     
     ja = jbas%jj(a)     
     
     j2min = max( abs(jp - ja) , abs(js - ju),abs(jt-jtot2) ) 
     j2max = min( jp+ja , js+ju,jt+jtot2) 
     
     j3min = max(abs(jq - jr) ,abs(jp-jtot1))
     j3max = min(jq+jr,jp+jtot1) 
   
     
     phase = (-1) ** ((ja + jp)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        if ((iq == ir) .and. (mod(J3/2,2)==1)) cycle
        sj1 = phase*sixj(jp,jq,Jpq,jr,jtot1,J3)*sqrt(J3+1.d0)  

        j1min = max(abs(ja - jt),abs(rank-J3))
        j1max = min(ja+jt,rank+J3) 
             
        do J1 = j1min,j1max,2 
           if ((a == it) .and. (mod(J1/2,2)==1)) cycle
           sj2 =  sj1*(-1)**(J1/2)*tensor_elem(iq,ir,a,it,J3,J1,R,jbas)*&
                xxxsixj(R%xindx,J1,J3,rank,jtot1,jtot2,jp)*sqrt(J1+1.d0)
           if (abs(sj2) < 1e-8) cycle
           do J2 = j2min,j2max,2
              if ((ip == a) .and. (mod(J2/2,2)==1)) cycle
              if ((iu == is) .and. (mod(J2/2,2)==1)) cycle
              sm = sm - (J2+1.d0) *sj2* sixj(js,jt,Jst,jtot2,ju,J2) &
                   * sixj(jp,ja,J2,jt,jtot2,J1) * v_elem(ip,a,iu,is,J2,L,jbas)
           end do
        end do
     end do
 
  end do

  ! right-left

  do ia=1,L%Nsp-L%belowEF
     a = jbas%parts(ia)
     
     if ( (tt+jbas%itzp(a)) .ne. (tq+tr) ) cycle
     if ( mod(lt+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 

     ja = jbas%jj(a)     
     
     j2min =  max(abs(js - ju),abs(jt-jtot2))
     j2max =  min(js+ju,jt+jtot2)
     
     j3min = max( abs(jq - jr) , abs(ja - jt) ,abs(jp-jtot1)) 
     j3max = min( jq+jr , ja+jt,jp+jtot1)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J2 = j2min,j2max,2
        if ((iu == is) .and. (mod(J2/2,2)==1)) cycle        
        sj1 = phase* sixj(js,jt,Jst,jtot2,ju,J2)*sqrt(J2+1.d0) 
        
        j1min = max(abs(jp - ja),abs(rank-J2))
        j1max = min(jp+ja ,rank+J2)
    
        do J1 = j1min,j1max,2
           if ((ip == a) .and. (mod(J1/2,2)==1)) cycle

           sj2 = sj1*(-1)**(J1/2)*xxxsixj(R%xindx,J1,J2,rank,jtot2,jtot1,jt) *sqrt(J1+1.d0) *&
                tensor_elem(ip,a,iu,is,J1,J2,R,jbas)      
           if (abs(sj2) < 1e-8) cycle
           do J3 = j3min,j3max,2
              if ((iq == ir) .and. (mod(J3/2,2)==1)) cycle
              if ((a == it) .and. (mod(J3/2,2)==1)) cycle
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

  if (.not.((is==it).and.(mod(Jst/2,2)==1))) then 
     do ia=1,L%belowEF
        a = jbas%holes(ia)
        if ((a == iq) .and. (mod(Jst/2,2)==1)) cycle
        if ( (tq+jbas%itzp(a)) .ne. (ts+tt) ) cycle
        if ( mod(lq+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 

        ja = jbas%jj(a)     


        j1min = max(abs(jp - jr),abs(jtot1-jq))
        j1max = min(jp+jr ,jtot1+jq) 

        phase = (-1) ** ((ja - ju)/2) ! changed to ja+js rather than ja-js 

        sj1 = v_elem(iq,a,is,it,Jst,L,jbas)*phase
        if (abs(sj1) < 1e-8) cycle
        do J1 = j1min,j1max,2
           if ((ip == ir) .and. (mod(J1/2,2)==1)) cycle          
           sj2 = sj1*sixj(jp,jq,Jpq,jtot1,jr,J1)*sqrt(J1+1.d0)*(-1)**(J1/2)
           if (abs(sj2) < 1e-8) cycle
           j2min = max(abs(ja - ju),abs(rank-J1),abs(jq-jtot2))
           j2max = min(ja+ju,rank+J1,jq+jtot2)

           do J2 = j2min,j2max,2 
              if ((a == iu) .and. (mod(J2/2,2)==1)) cycle

              sm = sm -  sj2*sqrt(J2+1.d0)*(-1)**(J2/2)*tensor_elem(ir,ip,a,iu,J1,J2,R,jbas)*&
                   xxxsixj(R%xindx,J1,J2,rank,jtot2,jtot1,jq)*sixj(ja,jq,Jst,jtot2,ju,J2)

           end do
        end do

     end do

     smtot = smtot + sm*multfact
  end if
  ! right-left
      
  sm = 0.d0 
  
  multfact = (-1)**((jp-jtot2+Jpq+rank)/2) *sqrt((Jpq+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 
  if (.not.((is==it).and.(mod(Jst/2,2)==1))) then 
     do ia=1,L%Nsp-L%belowEF
        a = jbas%parts(ia)

        if ( (tu+jbas%itzp(a)) .ne. (tr+tp) ) cycle
        if ( mod(lu+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

        ja = jbas%jj(a)     

        j1min = max(abs(jq - ja),abs(rank-Jst))
        j1max = min(jq+ja,rank+Jst) 

        j2min = max(abs(ja - ju),abs(jr-jp),abs(jq-jtot1)) 
        j2max = min(ja+ju,jr+jp,jq+jtot1) 

        phase = (-1) ** ((ja + jq)/2) ! changed to ja+js rather than ja-js 

        do J1 = j1min,j1max,2
           if ((a == iq) .and. (mod(J1/2,2)==1)) cycle
           sj1 = phase*sqrt(J1+1.d0)* &
                tensor_elem(iq,a,is,it,J1,Jst,R,jbas)*(-1)**(J1/2)
           if (abs(sj1) < 1e-8) cycle
           do J2 = j2min,j2max,2 
              if ((a == iu) .and. (mod(J2/2,2)==1)) cycle
              if ((ir == ip) .and. (mod(J2/2,2)==1)) cycle 
              sm = sm +  sj1*(J2+1.d0)*(-1)**(J2/2)*v_elem(ir,ip,a,iu,J2,L,jbas) *&
                   xxxsixj(R%xindx,J1,Jst,rank,jtot2,jtot1,ju)*sixj(ja,jq,J1,jtot1,ju,J2) &
                   *sixj(jp,jq,Jpq,jtot1,jr,J2)
           end do
        end do

     end do

     smtot = smtot + sm*multfact
  end if
  ! FIFTH TERM    
  !Left-Right
  sm = 0.d0 

  multfact = (-1)**((ju+jt+jtot2+js+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do ia=1,L%belowEF
     
     a = jbas%holes(ia)
     
     if ( (tq+jbas%itzp(a)) .ne. (tt+tu) ) cycle
     if ( mod(lq+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 

     ja = jbas%jj(a)     
     
     j1min = max( abs(jq - ja) , abs(jt - ju) ) 
     j1max = min( jq+ja , jt+ju) 
     
     j3min = max(abs(ja - js),abs(jq-jtot2))
     j3max = min(ja+js,jq+jtot2) 
     
     phase = (-1) ** ((ja + jp)/2) ! changed to ja+js rather than ja-js 
     
     do J1 = j1min,j1max,2
        if ((iq == a) .and. (mod(J1/2,2)==1)) cycle 
        if ((it == iu) .and. (mod(J1/2,2)==1)) cycle 
        sj1 = phase*sixj(js,jt,Jst,ju,jtot2,J1)*(J1+1.d0)*(-1)**(J1/2)&
             *v_elem(iq,a,it,iu,J1,L,jbas) 
        if (abs(sj1) < 1e-8) cycle
        
        do J3 = j3min,j3max,2 
           if ((a == is) .and. (mod(J3/2,2)==1)) cycle 
           sj2 =  sj1*(-1)**(J3/2)*&
                sixj(jq,ja,J1,js,jtot2,J3)*sqrt(J3+1.d0)
           
           j2min = max(abs(jp - jr),abs(rank-J3),abs(jq-jtot1))
           j2max = min(jp+jr,rank+J3,jq+jtot1) 
   
           do J2 = j2min,j2max,2
              if ((ir == ip) .and. (mod(J2/2,2)==1)) cycle 
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

  do ia=1,L%Nsp-L%belowEF
     a = jbas%parts(ia)
     if ( (ts+jbas%itzp(a)) .ne. (tr+tp) ) cycle
     if ( mod(ls+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

     ja = jbas%jj(a)     
    
     
     j2min =  max(abs(jt - ju),abs(js-jtot2))
     j2max =  min(jt+ju,js+jtot2)
     
     j3min = max( abs(jp - jr) , abs(ja - js) ) 
     j3max = min( jp+jr , ja+js)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        if ((ir == ip) .and. (mod(J3/2,2)==1)) cycle
        if ((a == is) .and. (mod(J3/2,2)==1)) cycle 
        sj1 = phase*(-1)**(J3/2)* sixj(jp,jq,Jpq,jtot1,jr,J3)*(J3+1.d0)*&
             v_elem(ir,ip,a,is,J3,L,jbas) 
        if (abs(sj1) < 1e-8) cycle
        do J2 = j2min,j2max,2
           if ((it == iu) .and. (mod(J2/2,2)==1)) cycle 
           sj2 = sj1*(-1)**(J2/2)*sixj(js,jt,Jst,ju,jtot2,J2)*sqrt(J2+1.d0)

           j1min = max(abs(jq - ja),abs(rank-J2),abs(js-jtot1)) 
           j1max = min(jq+ja,rank+J2,js+jtot1)
                    
           do J1 = j1min,j1max,2
              if ((a == iq) .and. (mod(J1/2,2)==1)) cycle 
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

  do ia=1,L%belowEF
        
     a = jbas%holes(ia)
          
     if ( (tq+jbas%itzp(a)) .ne. (tu+ts) ) cycle
     if ( mod(lq+jbas%ll(a),2) .ne. mod(lu+ls,2) ) cycle 

     ja = jbas%jj(a)     
     
     j1min = max( abs(jq - ja) , abs(js - ju) ,abs(jt-jtot2) ) 
     j1max = min( jq+ja , js+ju , jt+jtot2 ) 
     
     j2min = max(abs(jp - jr),abs(jq-jtot1))
     j2max = min(jp+jr ,jq+jtot1) 
       
     phase = (-1) ** ((ja - jp)/2) ! changed to ja+js rather than ja-js 
     
     do J1 = j1min,j1max,2
        if ((iq == a) .and. (mod(J1/2,2)==1)) cycle
        if ((iu == is) .and. (mod(J1/2,2)==1)) cycle 
        sj1 = phase*sixj(js,jt,Jst,jtot2,ju,J1)*(J1+1.d0) &
             *v_elem(iq,a,iu,is,J1,L,jbas) 
        if (abs(sj1) < 1e-8) cycle
        do J2 = j2min,j2max,2 
           if ((ir == ip) .and. (mod(J2/2,2)==1)) cycle 
           sj2 =  sj1*(-1)**(J2/2)*&
                sixj(jp,jq,Jpq,jtot1,jr,J2)*sqrt(J2+1.d0)
           
           j3min = max(abs(ja - jt),abs(rank-J2)) 
           j3max = min(ja+jt,rank+J2)
     
           do J3 = j3min,j3max,2
              if ((it == a) .and. (mod(J3/2,2)==1)) cycle 
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

  do ia=1,L%Nsp-L%belowEF
     a = jbas%parts(ia)
     if ( (tt+jbas%itzp(a)) .ne. (tr+tp) ) cycle
     if ( mod(lt+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

     ja = jbas%jj(a)     
    
     j2min = max(abs(js - ju),abs(jt-jtot2)) 
     j2max = min(js+ju ,jt+jtot2) 
     
     j3min = max( abs(jp - jr) , abs(ja - jt) ) 
     j3max = min( jp+jr , ja+jt)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        if ((ir == ip) .and. (mod(J3/2,2)==1)) cycle
        if ((it == a) .and. (mod(J3/2,2)==1)) cycle 
        sj1 = phase*(-1)**(J3/2) * sixj(jp,jq,Jpq,jtot1,jr,J3)*(J3+1.d0)*&
             v_elem(ir,ip,a,it,J3,L,jbas) 
        if (abs(sj1) < 1e-8) cycle
        do J2 = j2min,j2max,2 
           if ((iu == is) .and. (mod(J2/2,2)==1)) cycle 
           sj2 = sj1*sixj(js,jt,Jst,jtot2,ju,J2)*sqrt(J2+1.d0)
                    
           j1min = max(abs(jq - ja),abs(rank-J2),abs(jt-jtot1))
           j1max = min(jq+ja,rank+J2,jt+jtot1)
           
           do J1 = j1min,j1max,2
              if ((iq == a) .and. (mod(J1/2,2)==1)) cycle 
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
  if (.not.((is==it).and.(mod(Jst/2,2)==1))) then
     if (.not.((ip==iq).and.(mod(Jpq/2,2)==1))) then 
        multfact = (-1)**((Jpq+rank+jr-jtot2)/2) *sqrt((Jst+1.d0)*(jtot1+1.d0)*&
             (jtot2+1.d0)) 

        do ia=1,L%belowEF
           if ((ir == a) .and. (mod(Jst/2,2)==1)) cycle 
           a = jbas%holes(ia)

           if ( (tr+jbas%itzp(a)) .ne. (ts+tt) ) cycle
           if ( mod(lr+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 

           ja = jbas%jj(a)

           if (.not. triangle(jr,ja,jst) ) cycle

           j3min = max(abs(ja-ju),abs(Jpq-rank),abs(jr-jtot2))
           j3max = min(ja+ju,Jpq+rank,jr+jtot2)  

           sj1 = v_elem(ir,a,is,it,Jst,L,jbas) * (-1)**((ja+ju)/2) 
           if (abs(sj1) < 1e-8) cycle
           do J3=j3min,j3max,2 
              if ((iu == a) .and. (mod(J3/2,2)==1)) cycle 
              sm = sm - sj1* sixj(jr,ja,Jst,ju,jtot2,J3) * (-1)**(J3/2) * sqrt(J3+1.d0) &
                   * xxxsixj(R%xindx,J3,Jpq,rank,jtot1,jtot2,jr) * tensor_elem(ip,iq,a,iu,Jpq,J3,R,jbas)
           end do

        end do

        smtot = smtot + sm*multfact
     end if
  end if
  ! right-left
  sm = 0.d0
  if (.not.((is==it).and.(mod(Jst/2,2)==1))) then
     if (.not.((ip==iq).and.(mod(Jpq/2,2)==1))) then 
        
        multfact = (-1)**((Jpq+rank)/2) *sqrt((Jpq+1.d0)*(jtot1+1.d0)*&
             (jtot2+1.d0)) 

        do ia=1,L%Nsp-L%belowEF
           a = jbas%parts(ia)

           if ( (tu+jbas%itzp(a)) .ne. (tp+tq) ) cycle
           if ( mod(lu+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 
           if ((iu == a) .and. (mod(Jpq/2,2)==1)) cycle
           
           ja = jbas%jj(a)

           if (.not. triangle(ju,ja,Jpq) ) cycle

           j3min = max(abs(ja-jr),abs(rank-Jst),abs(ju-jtot1))
           j3max = min(ja+jr,rank+Jst,ju+jtot1)  

           sj1 = v_elem(ip,iq,a,iu,Jpq,L,jbas) * (-1)**((ja+jtot2)/2) 
           if (abs(sj1) < 1e-8) cycle
           do J3=j3min,j3max,2 
              if ((ir == a) .and. (mod(J3/2,2)==1)) cycle 
              sm = sm + sj1* sixj(jr,ja,J3,ju,jtot1,Jpq) * sqrt(J3+1.d0) *(-1)**(J3/2) &
                   * xxxsixj(R%xindx,J3,Jst,rank,jtot2,jtot1,ju) * tensor_elem(ir,a,is,it,J3,Jst,R,jbas)

           end do

        end do

        smtot = smtot + sm*multfact
     end if
  end if

  ! EIGHTH TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jt+jr+js+jtot2+Jpq+rank)/2) *sqrt((Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  ! Left-Right 
  sm = 0.d0   
  if (.not.((ip==iq).and.(mod(Jpq/2,2)==1))) then 
     do ia=1,L%belowEF 

        a = jbas%holes(ia)

        if ( (tr+jbas%itzp(a)) .ne. (tt+tu) ) cycle
        if ( mod(lr+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 

        ja = jbas%jj(a) 

        j2min = max(abs(ja - js),abs(rank-Jpq),abs(jtot2-jr))
        j2max = min(ja+js ,rank+Jpq,jtot2+jr)

        j1min = max(abs(ja - jr),abs(jt-ju)) 
        j1max = min(ja+jr,jt+ju)  

        phase = (-1)**((ja+ju)/2) 

        do J1 = j1min, j1max , 2
           if ((ir == a) .and. (mod(J1/2,2)==1)) cycle 
           if ((it == iu) .and. (mod(J1/2,2)==1)) cycle 
           sj1 = phase*(-1)**(J1/2)*(J1+1.d0) * sixj(js,jt,Jst,ju,jtot2,J1)&
                * v_elem(ir,a,it,iu,J1,L,jbas) 
           if (abs(sj1) < 1e-8) cycle
           do J2 = j2min,j2max,2
              if ((is == a) .and. (mod(J2/2,2)==1)) cycle 
              sm = sm - (-1)**(J2/2) * sqrt(J2 + 1.d0) &            
                   * sj1 * sixj(ja,js,J2,jtot2,jr,J1) *  &
                   xxxsixj(R%xindx,Jpq,J2,rank,jtot2,jtot1,jr) * &
                   tensor_elem(ip,iq,a,is,Jpq,J2,R,jbas)
           end do
        end do
     end do

     smtot = smtot + sm*multfact
  end if
  ! !Right-left  
  sm = 0.d0
  if (.not.((ip==iq).and.(mod(Jpq/2,2)==1))) then 
     multfact = (-1)**((jt+jtot2+Jpq+rank)/2) *sqrt((Jst+1.d0)*(Jpq+1.d0) &
          *(jtot1+1.d0)*(jtot2+1.d0) )  

     do ia=1,L%Nsp-L%belowEF 
        a = jbas%parts(ia)
        if ( (ts+jbas%itzp(a)) .ne. (tp+tq) ) cycle
        if ( mod(ls+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 
        if ((a==is).and.(mod(Jpq/2,2)==1)) cycle 
        ja = jbas%jj(a) 

        if (.not. triangle(ja,js,Jpq)) cycle 

        j1min = max(abs(jr - ja),abs(js-jtot1))
        j1max = min(jr + ja  ,js+jtot1) 


        phase = (-1)**((ja - ju)/2)

        sj1 = phase * v_elem(ip,iq,a,is,Jpq,L,jbas)
        if (abs(sj1) < 1e-8) cycle
        do J1 = j1min, j1max , 2      
           if ((ir==a).and.(mod(J1/2,2)==1)) cycle
           sj2 = sqrt(J1+1.d0)*(-1)**(J1/2)*sixj(jr,ja,J1,js,jtot1,Jpq) * sj1

           j2min = max(abs(jt - ju),abs(J1-rank),abs(js-jtot2))
           j2max = min(jt+ju,J1+rank,js+jtot2)

           do J2 = j2min,j2max,2
              if ((it==iu).and.(mod(J2/2,2)==1)) cycle
              sm = sm + sj2* sqrt(J2 + 1.d0) * (-1)**(J2/2) &
                   * sixj(js,jt,Jst,ju,jtot2,J2) * xxxsixj(R%xindx,J1,J2,rank,jtot2,jtot1,js) &
                   * tensor_elem(ir,a,it,iu,J1,J2,R,jbas)

           end do
        end do
     end do


     smtot = smtot + sm*multfact
  end if
! NINTH TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jr+jtot2+Jpq+Jst+rank)/2) *sqrt((Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  ! Left-Right
  
  sm = 0.d0   
  if (.not.((ip==iq).and.(mod(Jpq/2,2)==1))) then 
     do ia=1,L%belowEF 
        a = jbas%holes(ia)

        if ( (tr+jbas%itzp(a)) .ne. (tu+ts) ) cycle
        if ( mod(lr+jbas%ll(a),2) .ne. mod(lu+ls,2) ) cycle 

        ja = jbas%jj(a) 

        j2min = max(abs(ja - jt),abs(rank-Jpq))
        j2max = min(ja+jt,rank+Jpq) 

        j1min = max(abs(ja - jr),abs(js-ju),abs(jt-jtot2)) 
        j1max = min(ja+jr,js+ju,jt+jtot2)  

        phase = (-1)**((ja-js)/2) 

        do J1 = j1min, j1max , 2
           if ((ir==a).and.(mod(J1/2,2)==1)) cycle
           if ((iu==is).and.(mod(J1/2,2)==1)) cycle
           sj1 = phase*(J1+1.d0) * sixj(js,jt,Jst,jtot2,ju,J1)&
                * v_elem(ir,a,iu,is,J1,L,jbas) 
           if (abs(sj1) < 1e-8) cycle
           do J2 = j2min,j2max,2
              if ((it==a).and.(mod(J2/2,2)==1)) cycle
              sm = sm - (-1)**(J2/2) * sqrt(J2 + 1.d0) &            
                   * sj1 * sixj(ja,jr,J1,jtot2,jt,J2) *  &
                   xxxsixj(R%xindx,Jpq,J2,rank,jtot2,jtot1,jr) * &
                   tensor_elem(ip,iq,a,it,Jpq,J2,R,jbas)
           end do

        end do
     end do

     smtot = smtot + sm*multfact
  end if

  if (.not.((ip==iq).and.(mod(Jpq/2,2)==1))) then   
  ! !Right-left  
     sm = 0.d0 
     multfact = (-1)**((jt-jtot2+Jpq+Jst+rank)/2) *sqrt((Jst+1.d0)*(Jpq+1.d0) &
          *(jtot1+1.d0)*(jtot2+1.d0) )  

     do ia=1,L%Nsp-L%belowEF 
        a = jbas%parts(ia)
        if ( (tt+jbas%itzp(a)) .ne. (tp+tq) ) cycle
        if ( mod(lt+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 
        if ((it==a).and.(mod(Jpq/2,2)==1)) cycle
        ja = jbas%jj(a) 
        
        if (.not. triangle(ja,jt,Jpq)) cycle 

        j1min = max(abs(jr - ja),abs(jt-jtot1))
        j1max = min(jr + ja,jt+jtot1) 


        phase = (-1)**((ja + js)/2)

        sj1 = phase * v_elem(ip,iq,a,it,Jpq,L,jbas)
        if (abs(sj1) < 1e-8) cycle
        do J1 = j1min, j1max , 2      
           if ((ir==a).and.(mod(J1/2,2)==1)) cycle
           sj2 = sqrt(J1+1.d0)*(-1)**(J1/2)*sixj(jr,ja,J1,jt,jtot1,Jpq) * sj1

           j2min = max(abs(js - ju),abs(rank-J1),abs(jt-jtot2))
           j2max = min(js+ju,rank+J1,jt+jtot2)

           do J2 = j2min,j2max,2
              if ((iu==is).and.(mod(J2/2,2)==1)) cycle
              sm = sm + sj2* sqrt(J2 + 1.d0) &
                   * sixj(js,jt,Jst,jtot2,ju,J2) * xxxsixj(R%xindx,J1,J2,rank,jtot2,jtot1,jt) &
                   * tensor_elem(ir,a,iu,is,J1,J2,R,jbas)

           end do
        end do
     end do


     smtot = smtot + sm*multfact
  end if

     
  EOM_TS_commutator_223_single = smtot

end function EOM_TS_commutator_223_single

end module 
  
  
  

module EOM_dTZ_commutators
  use isospin_operators 
  use cross_coupled
  ! tensor-scalar commutator functions 

  ! THE TENSOR MUST BE THE SECOND ARGUMENT
  
contains
!=========================================================
!=========================================================
subroutine EOM_dTZ_commutator_111(L,R,RES,jbas) 
  ! one body part of [L1,R1] 
  implicit none 
  
  type(sq_op) :: L
  type(iso_ladder) :: R,RES
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
subroutine EOM_dTZ_commutator_121(L,R,RES,jbas) 
  ! onebody part of [L1,R2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L
  type(iso_ladder) :: R,RES
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
        if (tq .ne. tp - RES%dTz*2) cycle 
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

                    smx = smx + iso_ladder_elem(ak,pk,ik,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)  &
                          * xxxsixj(RES%xindx,J1,J2,rank,jq,jp,ji)
                    
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
subroutine  EOM_dTZ_commutator_211(L,R,RES,jbas) 
  ! onebody part of  - [R1,L2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(iso_ladder) :: R,RES
  type(sq_op) :: L
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
        if (tq .ne. tp-RES%dTz*2) cycle 
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
              if (ti .ne. ta-R%dTz*2) cycle 
              if (.not. (triangle(ja,ji,rank))) cycle
               
              Tz = abs(ta -ti)/2 
              if (abs(tp - tq) .ne. Tz*2)  cycle 
              PAR = mod(la+li,2) 
              if (mod(lp+lq,2) .ne. PAR) cycle                                 

              sm = sm + R%fph(a,i) * Vpandya(pk,qk,ak,ik,rank,L,jbas)
         
           end do 
        end do 
        
        RES%fph(p,q) = RES%fph(p,q) + sm  
        
     end do 
  end do 

end subroutine
!==================================================
!==================================================             
subroutine EOM_dTZ_commutator_122(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L
  type(iso_ladder) :: R,RES
  integer :: q,IX,JX,nh,np,nb,i,J1,J2
  integer :: a,b,c,d,ja,jb,jc,jd,ji,g_ix,q_sp,i_sp
  integer :: ta,tb,tc,td,ti,la,lb,lc,ld,li,spec
  integer :: jxstart,jxend,ixend,c1,c2,n1,n2
  logical :: square
  real(8) ::  sm
  
  
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1)
     J2 = R%tblck(q)%Jpair(2) 
       
     ! figure out how big the array is
     n1 = R%tblck(q)%npp1
     n2 = R%tblck(q)%nhh2
     if ((n1*n2) == 0) cycle 
     
     ! main calculation
     do IX = 1,n1
        a = R%tblck(q)%qn(1)%Y(IX,1)
        ja = jbas%jj(a)
        la = jbas%ll(a)
        ta = jbas%itzp(a) 
             
        b = R%tblck(q)%qn(1)%Y(IX,2)
        jb = jbas%jj(b)
        lb = jbas%ll(b)
        tb = jbas%itzp(b)
 
        do JX = 1,n2
           
           c = R%tblck(q)%qn(2)%Y(JX,1)
           jc = jbas%jj(c)
           lc = jbas%ll(c)
           tc = jbas%itzp(c)

           d = R%tblck(q)%qn(2)%Y(JX,2)
           jd = jbas%jj(d)
           ld = jbas%ll(d)
           td = jbas%itzp(d)
                   
           sm = 0.d0 

            ! a is replaced
            q_sp = sp_block_index(ja,la,ta,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               if (jbas%con(i_sp) .ne. jbas%con(a) ) cycle
               sm = sm + f_elem(a,i_sp,L,jbas)*iso_ladder_elem(i_sp,b,c,d,J1,J2,R,jbas)

            end do 
              
            ! b is replaced
            q_sp = sp_block_index(jb,lb,tb,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               if (jbas%con(i_sp) .ne. jbas%con(b) ) cycle               
               sm = sm + f_elem(b,i_sp,L,jbas)*iso_ladder_elem(a,i_sp,c,d,J1,J2,R,jbas) 
            end do 
            
            ! c is replaced
            q_sp = sp_block_index(jc,lc,tc,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               if (jbas%con(i_sp) .ne. jbas%con(c) ) cycle 
               sm = sm - f_elem(i_sp,c,L,jbas)*iso_ladder_elem(a,b,i_sp,d,J1,J2,R,jbas)
            end do 
            
            ! d is replaced
            q_sp = sp_block_index(jd,ld,td,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               if (jbas%con(i_sp) .ne. jbas%con(d) ) cycle 
               sm = sm - f_elem(i_sp,d,L,jbas)*iso_ladder_elem(a,b,c,i_sp,J1,J2,R,jbas)
            end do 
          
            sm = sm / sqrt(1.d0 + kron_del(a,b)) /sqrt(1.d0 + kron_del(c,d)) 
            
            RES%tblck(q)%Xpphh(IX,JX) = sm 


        end do
     end do
  end do 

end subroutine EOM_dTZ_commutator_122
!==================================================
!==================================================             
subroutine EOM_dTZ_commutator_212(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L
  type(iso_ladder) :: R,RES
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

     n1 = R%tblck(q)%npp1
     n2 = R%tblck(q)%nhh2
     
     ! figure out how big the array is
     if ((n1*n2) == 0) cycle 
     
     ! read in information about which 
     ! array we are using from public arrays
     
     ! main calculation
   
     do IX = 1,n1
        a = R%tblck(q)%qn(1)%Y(IX,1)
        ja = jbas%jj(a)
        modla = mod(jbas%ll(a),2)
        ta = jbas%itzp(a) 
             
        b = R%tblck(q)%qn(1)%Y(IX,2)
        jb = jbas%jj(b)
        modlb = mod(jbas%ll(b),2)
        tb = jbas%itzp(b)
 
        do JX = 1,n2
           
           c = R%tblck(q)%qn(2)%Y(JX,1)
           jc = jbas%jj(c)
           modlc = mod(jbas%ll(c),2)
           tc = jbas%itzp(c)

           d = R%tblck(q)%qn(2)%Y(JX,2)
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
                 if (ti == ta-R%dTz*2) then
                    if (modli == modla ) then 
                       if (triangle(ji,ja,rank)) then  
                          
                          sm1 = sm1 - xxxsixj(RES%xindx,J1,J2,rank,ji,ja,jb)&
                               *f_iso_ladder_elem(a,i,R,jbas)*v_elem(i,b,c,d,J2,L,jbas)
                       end if
                    end if
                 end if
               
                 sm1 = sm1*(-1)**((ja + jb-J2)/2) 
               
                 if (ti == tb-R%dTz*2) then
                    if (modli == modlb ) then 
                       if (triangle(ji,jb,rank)) then  
                          
                          sm2 = sm2 + xxxsixj(RES%xindx,J1,J2,rank,ji,jb,ja)&
                               *f_iso_ladder_elem(b,i,R,jbas)*v_elem(i,a,c,d,J2,L,jbas)
                       
                       end if
                    end if
                 end if
               
                 sm2 = sm2*(-1)**((J1+J2)/2) 
              end if

              sm3=0.d0
              sm4=0.d0
    
              if (jbas%con(i) .ne. jbas%con(c) ) then 
                 if (ti-R%dTz*2 == td) then
                    if (modli == modld ) then 
                       if (triangle(ji,jd,rank)) then  
                       
                          sm3 = sm3 + xxxsixj(RES%xindx,J1,J2,rank,jd,ji,jc)&
                               *f_iso_ladder_elem(i,d,R,jbas)*v_elem(a,b,c,i,J1,L,jbas)
                       
                       end if
                    end if
                 end if
               
                 sm3 = sm3*(-1)**((jc+jd-J1)/2) 
                 
                 if (ti-R%dTz*2 == tc) then
                    if (modli == modlc ) then 
                       if (triangle(ji,jc,rank)) then  
                          
                          sm4 = sm4 -  xxxsixj(RES%xindx,J1,J2,rank,jc,ji,jd)&
                               *f_iso_ladder_elem(i,c,R,jbas)*v_elem(a,b,d,i,J1,L,jbas)
                       
                       end if
                    end if
                 end if
                 
                 sm4 = sm4*(-1)**((J2+J1)/2)
              end if 
              
              sm =  sm + (sm1+sm2+sm3+sm4) 
           end do 
           
           sm = sm * sqrt((J1+1.d0)*(J2+1.d0) / &
              (1.d0 + kron_del(a,b)) /(1.d0 + kron_del(c,d))) * (-1)**(rank/2) 
           
           RES%tblck(q)%Xpphh(IX,JX) = RES%tblck(q)%Xpphh(IX,JX)  + sm 
           
        end do
     end do
   
  end do 

end subroutine EOM_dTZ_commutator_212
!===================================================================
!===================================================================
subroutine EOM_dTZ_commutator_222_pp_hh(L,R,RES,jbas) 
  !VERIFIED
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  L
  type(iso_ladder) :: R,RES
  integer :: q,q1,q2,J1,J2,Tz1,Tz2,Par,phase,rank,a,ja ,IX,JX,hx,px
  integer :: np1,nb1,nh1,np2,nb2,nh2,pm,i,p1,h1,kk,ll,jp,jh,ji
  real(8) :: bet_off,al_off,pre1,pre2,d6ji
  real(8),allocatable,dimension(:,:) :: WINT 
  
  pm = R%herm*L%herm
  rank = R%rank
  bet_off = 1.d0 
!construct temporary matrices
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1) 
     J2 = R%tblck(q)%Jpair(2)
     phase = R%tblck(q)%lam(1)
     par = R%tblck(q)%lam(2) 
     Tz1 = R%tblck(q)%lam(3)
     Tz2 = R%tblck(q)%lam(4)

     if (abs(Tz2) > 1 ) cycle
     if (abs(Tz1) > 1 ) cycle
     if (J2 > jbas%jtotal_max*2 ) cycle
     
     q1 = block_index(J1,Tz1,Par) 
     q2 = block_index(J2,Tz2,mod(Par+R%Dpar/2,2)) 


     
     np1 = R%tblck(q)%npp1
     nh2 = R%tblck(q)%nhh2

     nb1 = L%mat(q1)%nph
     nb2 = L%mat(q2)%nph 
       
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------
     if (np1*nh2 .ne. 0)  then 
        
        !Cpphh = Apppp.Bpphh 
    
        call dgemm('N','N',np1,nh2,np1,al,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%Xpphh,np1,bet_off,RES%tblck(q)%Xpphh,np1)
        
        !Cpphh = Bpphh.Ahhhh 

        call dgemm('N','N',np1,nh2,nh2,al,R%tblck(q)%Xpphh,np1,&
             L%mat(q2)%gam(5)%X,nh2,bet_off,RES%tblck(q)%Xpphh,np1)
             
     end if

     if (np1*nh2*nb1 .ne. 0)  then 

        allocate(WINT(nb1,nh2))
        WINT= 0.d0 
        
        al_off = L%herm * sqrt((J1+1.d0)*(J2+1.d0))
        
        call dgemm('T','N',nb1,nh2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
             R%tblck(q)%Xpphh,np1,bet,WINT,nb1)

        do IX = 1, nb1
           pre1 = 1.d0

           kk = L%mat(q1)%qn(2)%Y(IX,1)
           ll = L%mat(q1)%qn(2)%Y(IX,2)            

           if ( jbas%con(kk) == 0 ) then
              ! k is a particle

              jp = jbas%jj(kk)
              ji = jbas%jj(ll)

              p1 = kk
              i = ll
              
              pre1 = pre1 * (-1) ** ((ji-jp-J1)/2)

           else

              p1 = ll
              i = kk              

              jp = jbas%jj(ll)
              ji = jbas%jj(kk)

           end if

           
           do JX = 1, nh2

              pre2 = 1.d0

              kk = L%mat(q2)%qn(3)%Y(JX,1)
              ll = L%mat(q2)%qn(3)%Y(JX,2)            

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

              px = p1 - sum(jbas%con(1:p1-1))
              hx = h1 - sum(1-jbas%con(1:h1-1))

              RES%fph(px,hx) = RES%fph(px,hx) + (-1) ** ((jh+ji + rank+J1)/2) &
                   *xxxsixj(RES%xindx,J1,J2,rank,jh,jp,ji)*WINT(IX,JX)*pre1*pre2

           end do
        end do
        deallocate(WINT)
     end if

     if (np1*nh2*nb2 .ne. 0)  then 

        allocate(WINT(np1,nb2))
        WINT = 0.d0 
        al_off = -1*L%herm * sqrt((J1+1.d0)*(J2+1.d0))
        
        call dgemm('N','T',np1,nb2,nh2,al_off,R%tblck(q)%Xpphh,np1,&
             L%mat(q2)%gam(6)%X,nb2,bet,WINT,np1)

        do JX = 1, nb2
           pre2 = 1.d0

           kk = L%mat(q2)%qn(2)%Y(JX,1)
           ll = L%mat(q2)%qn(2)%Y(JX,2)            

           if ( jbas%con(kk) == 1 ) then
              ! k is a hole
              jh = jbas%jj(kk)
              ja = jbas%jj(ll)
              h1 = kk
              a = ll
              pre2 = pre2 * (-1) ** ((ja-jh-J2)/2)
           else
              h1 = ll
              a = kk
              jh = jbas%jj(ll)
              ja = jbas%jj(kk)
           end if
           
           do IX = 1, np1

              pre1 = 1.d0

              kk = L%mat(q1)%qn(1)%Y(IX,1)
              ll = L%mat(q1)%qn(1)%Y(IX,2)            

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

              px = p1 - sum(jbas%con(1:p1-1))
              hx = h1 - sum(1-jbas%con(1:h1-1))

              RES%fph(px,hx) = RES%fph(px,hx) + (-1) ** ((jh+ja + rank+J1)/2) &
                   *xxxsixj(RES%xindx,J1,J2,rank,jh,jp,ja)*WINT(IX,JX)*pre1*pre2

           end do
        end do
        deallocate(WINT)
     end if

  end do
     
end subroutine EOM_dTZ_commutator_222_pp_hh
!=================================================================
!=================================================================
 subroutine EOM_dTz_commutator_222_ph(L,R,RES,jbas) 
   ! VERIFIED ph channel 2body EOM_dTz_commutator. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(iso_ladder) :: R,RES
   real(8),allocatable,dimension(:,:) :: PANDYA_A,PANDYA_B,PANDYA_AB
   type(sq_op) :: L
   integer :: nh,np,nb1,nb2,q,IX,JX,i,r1,r2,Tz1_cc,Tz2_CC,PAR_J3,JTM,q1,q2,J3,J4,a,p1,p2,h1,h2
   integer :: ji,jh1,jh2,ti,th1,th2,lh1,lh2,li,n1,n2,c1,c2,jxstart,J4min,J4max,ja,jp1,jp2
   integer :: J1,J2, Jtot,Ntot,qx,J3min,J3max,ril,rjk,rli,rkj,g_ix,thread,total_threads
   integer :: phase1,phase2,phase3,rik,rki,rjl,rlj,PAR2,J1min,J2min,J1max,J2max,n_J3,n_J4,n_J5
   integer :: phase_34,phase_abcd,phase_ac,phase_bc,j5,PAR_J4,PAR_J5,rank,omp_get_num_threads
   integer :: phase_bd,phase_ad,nj_perm,full_int_phase,J5min,J5max,tp1,tp2,lp1,lp2
   integer,allocatable,dimension(:,:) :: qn_J3,qn_J4,qn_J5
   real(8) :: sm ,pre,pre2,omp_get_wtime ,t1,t2,coef9,factor,sm_ex,nj1,nj2  
   real(8) :: prefac_12,prefac_1,nj,Xelem,Yelem,V,al_off,d6ji
   logical :: square
   
   Ntot = RES%Nsp
   JTM = jbas%Jtotal_max*2
   !$OMP PARALLEL
   total_threads = omp_get_num_threads()
   !$OMP END PARALLEL
   rank = R%rank

   ! construct intermediate matrices

   !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(R,L,RES)   
   ! sum over all cross coupled channels for both matrices 
   do thread = 0,total_threads-1
      
   
   
   do J3 = 2*thread, JTM  ,2*total_threads      
      J4min = abs(rank -J3)
      J4max = min(rank + J3,JTM)

      do Tz1_cc = 0,2,2 ! this is the cross coupled TZ
         Tz2_cc = abs(Tz1_cc - abs(R%dTz*2))
         
         do PAR_J3 = 0,1

            n_J3 = count_dTz_configs( J3 ,Tz1_cc, PAR_J3 , jbas, .true. ) !boolean true: ph, false: hp configs
            if (n_J3 == 0) cycle
            allocate(qn_J3(n_J3,2)) 
            n_J3 = count_dTz_configs( J3 ,Tz1_cc, PAR_J3 , jbas, .true.,qn_J3) ! this fills the cc- basis descriptor

            allocate(PANDYA_A(n_J3, n_J3) )
            PANDYA_A = 0.d0 
            call fill_cc_matrix(J3,PANDYA_A,qn_J3,L,jbas)
            
            PAR_J4 = mod(PAR_J3 + R%dpar/2,2)
            
            do J4 = J4min,J4max ,2

               n_J4 = count_dTz_configs( J4 ,Tz2_cc, PAR_J4 , jbas, .true. )
               if (n_J4 == 0) cycle
               allocate(qn_J4(n_J4,2)) 
               n_J4 = count_dTz_configs( J4 ,Tz2_CC, PAR_J4 , jbas, .true.,qn_J4)

               allocate( PANDYA_B( n_J3, n_J4) )
               allocate( PANDYA_AB( n_J3, n_J4) )
               PANDYA_B = 0.d0
               
               call fill_generalized_isopandya_matrix(J3,J4,PANDYA_B,qn_J3,qn_J4,R,jbas)
              
               PANDYA_AB = 0.d0 
               
               al_off =  (-1)** ((J4+J3)/2) *  sqrt((J3+1.d0)*(J4+1.d0)) 

               call dgemm('N','N',n_J3,n_J4,n_J3,al_off,PANDYA_A,n_J3,PANDYA_B,n_J3,bet,PANDYA_AB,n_J3) 

               do JX = 1,n_J4

                  ! GET KET 
                  p2 = qn_J4(JX,1)
                  h1 = qn_J4(JX,2)

                  jp2 = jbas%jj(p2)
                  lp2 = jbas%ll(p2)
                  tp2 = jbas%itzp(p2)
                  jh1 = jbas%jj(h1)
                  lh1 = jbas%ll(h1)
                  th1 = jbas%itzp(h1)

                  do IX = 1, n_J3 
                     
                     ! GET BRA
                     p1 = qn_J3(IX,1)
                     h2 = qn_J3(IX,2)

                     jp1 = jbas%jj(p1)
                     jh2 = jbas%jj(h2)
                     lp1 = jbas%ll(p1)
                     tp1 = jbas%itzp(p1)
                     lh2 = jbas%ll(h2)
                     th2 = jbas%itzp(h2)

                     if ( (tp1 +tp2)-2*R%dTz .ne. (th1+th2) ) cycle
                     if ( mod(lp1 +lp2,2) .ne. mod(lh1+lh2+RES%dpar/2,2) ) cycle
                        
                     ! CALCULATE X CONTRIBUTIONS

                     J1min = abs(jp1-jp2)
                     J1max = jp2+jp1
                     
                     ! these are the results of the Matmuls 

                     Xelem = PANDYA_AB(IX,JX)
                     if (abs(Xelem) <1e-6) cycle

                     do J1 = J1min,J1max,2

                        if ((p1==p2).and.(mod(J1/2,2)==1)) cycle                           
                        J2min = max(abs(jh1-jh2),abs(rank-J1))
                        J2max = min(jh1+jh2 ,rank+J1) 
                        prefac_1 = sqrt(J1+1.d0)
                        
                        do J2 = J2min,J2max,2

                           if ((h1==h2).and.(mod(J2/2,2)==1)) cycle                           
                           nj = ninej(RES%xindx,jp1,jh2,J3,jp2,jh1,J4,J1,J2,rank)
                           prefac_12 = prefac_1 *sqrt(J2+1.d0)
                           
                           
                           if ( h2 .ge. h1 ) then
                              if ( p2 .ge. p1 ) then

                                 V = prefac_12*nj*(-1)**((jp1-jp2+J2)/2)*Xelem

                                 call add_elem_to_ladder(V,p1,p2,h1,h2,J1,J2,RES,jbas)                                                                             
                              else

                                 V = prefac_12*nj*(-1)**((J1+J2)/2)*Xelem

                                 call add_elem_to_ladder(V,p2,p1,h1,h2,J1,J2,RES,jbas)                                                                                                                
                              end if
                           else
                              if (p1 > p2) then                                    
                                 V = prefac_12*nj*(-1)**((jh1-jh2+J1)/2)*Xelem                                    

                                 call add_elem_to_ladder(V,p2,p1,h2,h1,J1,J2,RES,jbas)                                         
                              else
                                 
                                 V = prefac_12*nj*(-1)**((jp1+jp2+jh1+Jh2)/2)*Xelem

                                 call add_elem_to_ladder(V,p1,p2,h2,h1,J1,J2,RES,jbas)                                         

                              end if
                           end if                           

                        end do
                     end do
                  end do
               end do


               deallocate(PANDYA_B,PANDYA_AB,qn_J4) 
                                 
            end do
            deallocate(qn_J3,PANDYA_A)
         end do
      end do
   end do
   end do
   !$OMP END PARALLEL DO
 end subroutine EOM_dTz_commutator_222_ph


 end module 
  
  
  

module interaction_IO
  use basic_IMSRG
  implicit none
  
  contains

subroutine get_me2j_spfile()
  ! This constructs the sps file for
  ! me2j's format, after it's been converted to a pn-basis 
  implicit none 
  
  integer :: e,eMax,l,jj,tz,q,n,lmax,nmax
  character(2) :: eMaxchr
  character(13) :: fmt
  
  print*, 'eMax' 
  read*, emax
  print*, 'lMax' 
  read*, lmax
  print*, 'nmax' 
  read*, nmax
  open(unit=24,file='hk12_lmax10.sps')
  
  q = 1
  
  fmt = '(5(I5),e17.7)'
  do  e = 0,eMax
     
     do l = mod(e,2),min(e,lmax),2
        
        n = (e-l)/2
        
        if (n > nmax) cycle 
        
        do jj = abs(2*l-1),2*l+1,2
           
           do tz = -1,1,2
              write(*,fmt) q, n , l , jj , tz, float(e) 
              write(24,fmt) q, n , l , jj , tz, float(e) 
              q = q+1
           end do 
        end do 
     end do 
  end do 
  
  close(24) 

  STOP
end subroutine
!=====================================================================================
!=====================================================================================
subroutine read_me2j_interaction(H,jbas,jbx,htype)
   ! reads three-body force from darmstadt .me2j.gz
  use gzipmod
  implicit none 
  
  integer :: nlj1,nlj2,nnlj1,nnlj2,j,T,Mt,nljMax,endpoint,j_min,j_max,htype,lMax,eMaxfile
  integer :: l1,l2,ll1,ll2,j1,j2,jj1,jj2,Ntot,i,q,bospairs,qx,ta,tb,tc,td,n1,n2,nn1,nn2
  integer :: eMax,iMax,jmax,jmin,JT,a,b,c,d,C1,C2,i1,i2,pre,COM,x,PAR,endsz,ax,bx,cx,dx,Ntot_file
  integer,allocatable,dimension(:) :: indx 
  real(8),allocatable,dimension(:) :: ME,MEpp,MErr,me_fromfile,ppff,rrff
  real(8) :: V,g1,g2,g3,hw,pre2
  type(spd) :: jbas,jbx
  type(sq_op) :: H
  logical :: pp_calc,rr_calc
  character(1) :: rem
  character(2) :: eMaxchr
  type(c_ptr) :: buf,buf2,buf3
  integer(c_int) :: hndle,hndle2,hndle3,sz,sz2,sz3,rx
  character(kind=C_CHAR,len=200) :: buffer,buffer2,buffer3
  
  hw = H%hospace
  COM = 0
  if (htype == 1) COM = 1

  lMax = H%lmax
  eMax = H%eMax
  Ntot_file = jbx%total_orbits/2
  Ntot = jbas%total_orbits/2
  
  i = 0
  q = 0

  ! move in increments of two, because I use a pn basis,
  !  heiko's is isospin coupled (half the states) 

  eMax = maxval(jbas%e)   
  
  ! counting the states, and labeling them
  do nlj1=1, 2*Ntot_file,2 
     l1= jbx%ll(nlj1)
     j1= jbx%jj(nlj1)
     
     do nlj2 = 1, nlj1,2
        l2= jbx%ll(nlj2)
        j2= jbx%jj(nlj2)
        
        do nnlj1 = 1 ,nlj1 , 2
           ll1= jbx%ll(nnlj1)
           jj1= jbx%jj(nnlj1)
           
           endpoint = nnlj1 
           if ( nlj1==nnlj1 ) endpoint = nlj2
         
           do nnlj2 = 1, endpoint , 2 
              ll2= jbx%ll(nnlj2)
              jj2= jbx%jj(nnlj2)
              
              if (mod(l1+l2,2) .ne. mod(ll1+ll2,2)) cycle
              jmin = max( abs(j1-j2) , abs(jj1-jj2) ) 
              jmax = min( j1+j2  , jj1+jj2) 
              
              if (jmin > jmax) cycle 
      
              do JT = jmin,jmax,2
                 i = i + 4
              end do 
                    
           end do
        end do 
     end do 
  end do 
  iMax = i 
  
  allocate(me(iMax)) 
  allocate(me_fromfile(10)) 
  

  ! using zlib c library, which is bound with fortran in file "gzipmod.f90" 
  
  ! I don't know why you have to tack on those //achars(0) but it seems nessecary 
  hndle=gzOpen(trim(TBME_DIR)//trim(adjustl(intfile))//achar(0),"r"//achar(0)) 
  
  sz=200
  
  buf=gzGets(hndle,buffer,sz) 
  endpoint = 10 
  write(rem,'(I1)') endpoint-1
  endsz = 130 
  
  do i = 1,iMax,10
  
     if (i+10 > iMax) then 
        deallocate(me_fromfile)
        allocate(me_fromfile( iMax - i + 1) ) 
        endpoint = iMax-i + 1
        endsz = 13+(endpoint-1)*13 
        write(rem,'(I1)') endpoint-1
     end if
  
     buf = gzGets(hndle,buffer(1:sz),sz)

     
  
     read(buffer(1:endsz),'(f12.7,'//rem//'(f13.7))') me_fromfile 
        
     do j = 1,endpoint 
        ME(i+j-1) = me_fromfile(j)   
     end do

  end do

  rx = gzClose(hndle) 
  
  deallocate(me_fromfile)
  allocate(me_fromfile(4))
  
  ! redo this loop to put everything in pn basis
  
  i=0
  do nlj1=1, 2*Ntot_file,2 
     l1= jbx%ll(nlj1)
     j1= jbx%jj(nlj1)
     n1= jbx%nn(nlj1) 
     
     do nlj2 = 1, nlj1,2
        l2= jbx%ll(nlj2)
        j2= jbx%jj(nlj2)
        n2= jbx%nn(nlj2) 
        
        do nnlj1 = 1 ,nlj1 , 2
           ll1= jbx%ll(nnlj1)
           jj1= jbx%jj(nnlj1)
           nn1= jbx%nn(nnlj1) 
   
           endpoint = nnlj1 
           if ( nlj1==nnlj1 ) endpoint = nlj2
         
           do nnlj2 = 1, endpoint , 2 
              ll2= jbx%ll(nnlj2)
              jj2= jbx%jj(nnlj2)
              nn2= jbx%nn(nnlj2) 
              
              if (mod(l1+l2,2) .ne. mod(ll1+ll2,2)) cycle
              jmin = max( abs(j1-j2) , abs(jj1-jj2) ) 
              jmax = min( j1+j2  , jj1+jj2) 
              PAR = mod(l1+l2,2) 
           
              if (jmin > jmax) cycle 
            
              do JT = jmin,jmax,2
              
                 me_fromfile=ME(i+1:i+4)
                 i = i + 4 ! four different TMt qnums
                 
                 if ((l1 > lmax).or.(l2 > lmax).or.&
                      (ll1 > lmax).or.(ll2 > lmax))  cycle 
                 
                 if ((2*n1+l1 > eMax ).or.(2*n2+l2 > eMax )) cycle 
                 if ((2*nn1+ll1 > eMax ).or.(2*nn2+ll2 > eMax )) cycle 
                 

                 
                 !sum over all isospin configs
                 do ax = nlj1,nlj1+1
                    a = jbx%con(ax) !converting from full sp basis to restricted one

                    do bx = nlj2,nlj2+1
                       b = jbx%con(bx)

                       do cx = nnlj1,nnlj1+1
                          c = jbx%con(cx)

                          do dx = nnlj2,nnlj2+1  
                             d = jbx%con(dx) 

                             ! conversion factor to mT scheme 
                             pre2 = 1.d0 
                             if ( a == b ) pre2 = pre2*sqrt(0.5d0) 
                             if ( c == d ) pre2 = pre2*sqrt(0.5d0) 

                             ! heikos labeling is backwards
                             ta = -jbas%itzp(a)
                             tb = -jbas%itzp(b)
                             tc = -jbas%itzp(c)
                             td = -jbas%itzp(d)

                             T = ta+tb
                             if (tc+td .ne. T) cycle
                             T = -T/2
                             q = block_index(JT,T,Par)

                             ! convert to pn matrix element       
                             V =  0.125d0*(ta-tb)*(tc-td)*me_fromfile(1)+&   ! 00 clebsch
                                  kron_del(ta+tb,-2)*kron_del(tc+td,-2)*me_fromfile(2)+& ! 1-1 
                                  kron_del(ta+tb,2)*kron_del(tc+td,2)*me_fromfile(4)+& !11 
                                  0.125d0*abs((ta-tb)*(tc-td))*me_fromfile(3) !10 

                             ! pipj ,rirj
                             g2 = r1_r2( a, b, c, d, JT ,jbas ) 
                             g3 = p1_p2( a, b, c, d, JT ,jbas )
                             g1 = (g3 + H%com_hw**2 /hw**2 *g2)*H%lawson_beta 

                             ! center of mass subtraction
                             V = V*pre2 +(g1-g3*COM)*hw/(H%Aneut+H%Aprot)

                             C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
                             C2 = jbas%con(c)+jbas%con(d) + 1

                             qx = C1*C2
                             qx = qx + adjust_index(qx)   !Vpppp nature  

                             ! get the indeces in the correct order
                             pre = 1
                             if ( a > b )  then 
                                x = bosonic_tp_index(b,a,Ntot*2) 
                                j_min = jbas%xmap(x)%Z(1)  
                                i1 = jbas%xmap(x)%Z( (JT-j_min)/2 + 2) 
                                pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -JT)/2 ) 
                             else
                                ! if (a == b) pre = pre * sqrt( 2.d0 )

                                x = bosonic_tp_index(a,b,Ntot*2) 
                                j_min = jbas%xmap(x)%Z(1)  
                                i1 = jbas%xmap(x)%Z( (JT-j_min)/2 + 2) 

                             end if

                             if (c > d)  then     

                                x = bosonic_tp_index(d,c,Ntot*2) 
                                j_min = jbas%xmap(x)%Z(1)  
                                i2 = jbas%xmap(x)%Z( (JT-j_min)/2 + 2) 

                                pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -JT)/2 ) 
                             else 
                                ! if (c == d) pre = pre * sqrt( 2.d0 )

                                x = bosonic_tp_index(c,d,Ntot*2) 
                                j_min = jbas%xmap(x)%Z(1)  
                                i2 = jbas%xmap(x)%Z( (JT-j_min)/2 + 2) 

                             end if

                             ! kets/bras are pre-scaled by sqrt(2) if they 
                             ! have two particles in the same sp-shell

                             ! get the units right. I hope 

                             if ((qx == 1) .or. (qx == 5) .or. (qx == 4)) then 
                                H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
                                H%mat(q)%gam(qx)%X(i1,i2)  = V *pre

                             else if (C1>C2) then
                                H%mat(q)%gam(qx)%X(i2,i1)  = V *pre

                             else
                                H%mat(q)%gam(qx)%X(i1,i2) = V * pre

                             end if
                             ! I shouldn't have to worry about hermiticity here, input is assumed to be hermitian
                          end do
                       end do
                    end do
                 end do
                 
              end do ! end sum over j 
              
           
           end do  !end sums over Tcoupled lables
        end do
     end do
  end do
  close(7)
end subroutine read_me2j_interaction
!=====================================================================================
!=====================================================================================
subroutine read_me3j(store_3b,jbas,jbx,eMax,lmax)
  ! reads three-body force from darmstadt .me3j.gz or .me3j.bin
  use three_body_routines
  implicit none 
  
  type(spd) :: jbas,jbx
  type(three_body_force) :: store_3b
  integer :: nlj1,nlj2,nlj3,nnlj1,nnlj2,nnlj3,aux,aux1,aux2,aux3,aux4
  integer :: nnlj2_end,nnlj3_end,twoTMin,twoTmax,twoJCMin,twoJCMax
  integer :: twoJCMindown,twoJCMaxup,twoJCMindownket,twoJCMindownbra
  integer :: twoJCMaxupket,twoJCMaxupbra,la,lb,lc,ld,le,lf
  integer :: nsp,ja,jb,jc,jd,je,jf,iblock,Jab,JJab,Tab,TTab,ttot,jtot
  integer :: ea,eb,ec,ed,ef,ee,e1max,E3max,JabMax,JabMin,JJabMax,JJabMin
  integer :: Tij_indx,Tlm_indx,blocksize,endpoint,endsz,endsz2
  integer :: lmax3,jtot_max,jtot_max_1,jtot_max_2,jtot_min,jtot_min_1
  integer :: jtot_min_2,i,II,JJ,Jab_max,Jab_min,jc_max,jc_min
  integer :: Jde_min,Jde_max,x1,x2,q,NN,nsp_iso,tc_min,tc_max,j
  integer :: spot_in_gz_file,iMax,PAR,Tab_indx,TTab_indx,r,w,eMax,lmax
  integer :: E3Max_file,lmax_file,eMax_file,a,b,c,d,e,f
  real(8) :: szofblock,V,elem
  character(1)::rem
  real(8),allocatable,dimension(:) :: xblock,me_fromfile
  real,allocatable,dimension(:) :: ME
  type(c_ptr) :: buf,buf2,buf3
  integer(c_int) :: hndle,hndle2,hndle3,sz,sz2,sz3,rx
  character(kind=C_CHAR,len=200) :: buffer,buffer2,buffer3
  logical :: autozero ,thing
  character(255) :: header 
 
  iMax = store_3b%num_elems 
  allocate(ME(iMax)) 
  nsp = jbx%total_orbits
  nsp_iso = nsp/2
  !open file 

  E3max = store_3b%e3max

  E3Max_file = jbas%E3Max_3file 
  lMax_file = jbas%lMax_3file 
  eMax_file = jbas%eMax_3file 

  if( threebody_file(len(trim(threebody_file))-6:len(trim(threebody_file))) == 'me3j.gz') then 
     
     allocate(me_fromfile(10)) 

     hndle=gzOpen(trim(TBME_DIR)//trim(adjustl(threebody_file))//achar(0),"r"//achar(0))   
  
     sz=200
  
     buf=gzGets(hndle,buffer,sz) !first line 
     endpoint = 10 
     write(rem,'(I1)') endpoint-1
     endsz = 130 
     
     do i = 1,iMax,10
        
        if (i+10 > iMax) then 
           deallocate(me_fromfile)
           allocate(me_fromfile( iMax - i + 1) ) 
           endpoint = iMax-i + 1
           endsz = 13+(endpoint-1)*13 
           write(rem,'(I1)') endpoint-1
        end if
  
        buf = gzGets(hndle,buffer(1:sz),sz)

     
  
        read(buffer(1:endsz),'(f12.7,'//rem//'(f13.7))') me_fromfile 
        
        do j = 1,endpoint 
           ME(i+j-1) = me_fromfile(j)   
        end do

     end do

     rx = gzClose(hndle) 
  
     deallocate(me_fromfile)
  
  else
     
     open(unit=77,file=trim(TBME_DIR)//trim(adjustl(threebody_file)), form = 'unformatted', &
          access = 'stream') 
     
     read(77) header,ME
     close(77) 
  
  end if 

  i = 1
  r = 0 
  w = 0
  spot_in_gz_file = 0 
  do nlj1 = 1, nsp_iso 
     la = jbx%ll(2*nlj1)
     ja = jbx%jj(2*nlj1)
     ea = 2*jbx%nn(2*nlj1)+la 
     if (ea > E3Max_file) exit
     
     do nlj2 = 1,nlj1
        lb = jbx%ll(2*nlj2)
        jb = jbx%jj(2*nlj2)
        eb = 2*jbx%nn(2*nlj2)+lb 
        if (ea + eb > E3Max_file) exit
        
        do nlj3 = 1,nlj2
           lc = jbx%ll(2*nlj3)
           jc = jbx%jj(2*nlj3)
           ec = 2*jbx%nn(2*nlj3)+lc 
           if (ea + eb + ec > E3Max_file) exit
           
     do nnlj1 = 1, nlj1 
        ld = jbx%ll(2*nnlj1)
        jd = jbx%jj(2*nnlj1)
        ed = 2*jbx%nn(2*nnlj1)+ld 
        if (ed > E3Max_file) exit
        
        if ( nlj1 == nnlj1 ) then
           nnlj2_end = nlj2
        else
           nnlj2_end = nnlj1
        end if
                  
        do nnlj2 = 1, nnlj2_end            
           le = jbx%ll(2*nnlj2)
           je = jbx%jj(2*nnlj2)
           ee = 2*jbx%nn(2*nnlj2)+le 
           if (ed + ee > E3Max_file) exit
        
           if ( (nlj1 == nnlj1).and.(nlj2==nnlj2)) then 
              nnlj3_end = nlj3
           else
              nnlj3_end = nnlj2
           end if 
           
           do nnlj3 = 1,nnlj3_end
              lf = jbx%ll(2*nnlj3)
              jf = jbx%jj(2*nnlj3)
              ef = 2*jbx%nn(2*nnlj3)+lf 
              if (ed + ee + ef > E3Max_file) exit

              
              ! check parity 
              PAR = mod(la+lb+lc,2) 
              if ( mod(ld+le+lf,2) .ne. PAR ) cycle
              
              JabMax =  ja+jb 
              JabMin = abs(ja-jb) 
              
              JJabMax = jd+je
              JJabMin = abs(jd-je) 
              
              !determine roughly the two*J bounds
              
              if (abs(ja-jb) > jc) then 
                 twoJCMindownbra = abs(ja-jb)-jc
              else if ( jc < ja+jb) then 
                 twoJCMindownbra = 1
              else
                 twoJCMindownbra = jc - ja -jb
              end if 
              
              if (abs(jd-je) > jf) then 
                 twoJCMindownket = abs(jd-je)-jf
              else if ( jf < jd+je) then 
                 twoJCMindownket = 1
              else
                 twoJCMindownket = jf - jd -je
              end if 
              
              twoJCMaxupbra = ja+jb+jc
              twoJCMaxupket = jd+je+jf
              
              twoJCMindown = max(twoJCMindownket,twoJCMindownbra) 
              twoJCMaxup = min(twoJCMaxupket,twoJCMaxupbra)
              
              if (twoJCMindown > twoJCMaxup) cycle
              
              do Jab = JabMin,JabMax,2
                 do JJab=JJabMin,JJabMax,2
                    
                    twoJCMin = max(abs(Jab-jc),abs(JJab-jf))
                    twoJCMax = min(Jab+jc,JJab+jf) 
                    
                    do jtot = twoJCMin,twoJCMax,2
                       
                       do Tab = 0,1
                          do TTab = 0,1
                             
                             twoTMin = max(abs(2*Tab-1),abs(2*TTab-1))
                             twoTMax = min(2*Tab+1,2*TTab+1) 
                             
                             do ttot = twoTMin,twoTmax,2
                                
                                q = block_index_3b(jtot,ttot,PAR)                                 
                                elem = ME(i)
                               
                                i = i + 1 
                                 
                                if ((ea > eMax).or.(eb > eMax).or.(ec > eMax).or.(ed > eMax)&
                                     .or.(ee > eMax).or.(ef > eMax)) cycle 
                                
                                if ((la > lMax).or.(lb > lMax).or.(lc > lMax).or.(ld > lMax)&
                                     .or.(le > lMax).or.(lf > lMax)) cycle 
                                      
                                if ( ea+eb+ec > E3Max ) cycle 
                                if ( ed+ee+ef > E3Max ) cycle
                                                                

                                a = jbx%con(2*nlj1) 
                                b = jbx%con(2*nlj2)
                                c = jbx%con(2*nlj3)
                                d = jbx%con(2*nnlj1)
                                e = jbx%con(2*nnlj2)
                                f = jbx%con(2*nnlj3)

                                call SetME(Jab,JJab,jtot,2*Tab,2*TTab,ttot,a,b,c,d,e,f,elem,store_3b,jbas)
                              
                             end do !ttot
                          end do !TTab
                       end do ! Tab
                    end do !jtot
                 end do !JJab
              end do !Jab
           end do !nnlj3
        end do !nnlj2
     end do !nnlj1 
        end do !nlj3
     end do !nlj2
  end do !nlj1
               
end subroutine read_me3j


end module




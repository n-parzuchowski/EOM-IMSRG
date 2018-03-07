module interaction_IO
  use basic_IMSRG
  implicit none
  
  contains
!==================================================================  
!==================================================================
subroutine read_interaction(H,jbas,htype) 
  ! read interaction from ASCII file produced by VRenormalize
  implicit none
  
  type(sq_op) :: H
  type(spd) :: jbas
  integer :: ist,J,Tz,Par,a,b,c,d,q,qx,N,j_min,x
  real(8) :: V,Vcm,g1,g2,g3,pre,hw,V1
  integer :: C1,C2,int1,int2,i1,i2,htype,COM
  
  hw = H%hospace
  
  open(unit=39,file = trim(TBME_DIR)//trim(adjustl(intfile))) 
  
  read(39,*);read(39,*);read(39,*);read(39,*)
  read(39,*);read(39,*);read(39,*);read(39,*) !skip all of the garbage 
  
  COM = 0
  if (htype == 1) COM = 1 ! remove center of mass hamiltonian? 
  
  N = jbas%total_orbits
  
  do 
     read(39,*,iostat=ist) Tz,Par,J,a,b,c,d,V!,g1,g2,g3
     !read(39,*) Tz,Par,J,a,b,c,d,V,g1,g2,g3

     ! g1 is COM expectation value, NOT CALCULATED WITH SCOTT'S CODE
     ! g2 is r1*r2 ME 
     ! g3 is p1*p2 ME 
     
     if (ist > 0) STOP 'interaction file error' 
     if (ist < 0) exit

     g2 = r1_r2( a, b, c, d, J ,jbas ) ! morten and Koshiroh are inconsistent with their definitions
     g3 = p1_p2( a, b, c, d, J ,jbas ) ! morten and Koshiroh are inconsistent with their definitions
     ! so I am doing it explicitly. 
     
     g1 = (g3 + H%com_hw**2 /hw**2 *g2)*H%lawson_beta ! lawson term    

     V = V + (g1 - g3*COM) *hw/(H%Aneut + H%Aprot) ! center of mass correction
 
     
     q = block_index(J,Tz,Par) 
     
     C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
     C2 = jbas%con(c)+jbas%con(d) + 1
    
     qx = C1*C2
     qx = qx + adjust_index(qx)   !Vpppp nature  

     ! get the indeces in the correct order
     pre = 1

     if ( a > b )  then 
        
        x = bosonic_tp_index(b,a,N) 
        j_min = jbas%xmap(x)%Z(1)  
        i1 = jbas%xmap(x)%Z( (J-j_min)/2 + 2) 
        pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J)/2 ) 
     else
       ! if (a == b) pre = pre * sqrt( 2.d0 )
       
        x = bosonic_tp_index(a,b,N) 
        j_min = jbas%xmap(x)%Z(1)  
        i1 = jbas%xmap(x)%Z( (J-j_min)/2 + 2) 
     end if
  
     if (c > d)  then     
        
        x = bosonic_tp_index(d,c,N) 
        j_min = jbas%xmap(x)%Z(1)  
        i2 = jbas%xmap(x)%Z( (J-j_min)/2 + 2) 
        
        pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J)/2 ) 
     else 
       ! if (c == d) pre = pre * sqrt( 2.d0 )
      
        x = bosonic_tp_index(c,d,N) 
        j_min = jbas%xmap(x)%Z(1)  
        i2 = jbas%xmap(x)%Z( (J-j_min)/2 + 2) 
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
        H%mat(q)%gam(qx)%X(i1,i2) =  V *pre 
     end if 
     ! I shouldn't have to worry about hermiticity here, input is assumed to be hermitian
     
  end do
  close(39)
      
end subroutine read_interaction
!==================================================================  
!==================================================================
subroutine read_gz(H,jbas,htype) 
  ! read interaction from gz file produced by VRenormalize 
  implicit none
  
  type(sq_op) :: H
  type(spd) :: jbas
  integer :: ist,J,Tz,Par,a,b,c,d,q,qx,N,j_min,x
  real(8) :: V,Vcm,g1,g2,g3,pre,hw
  integer :: C1,C2,int1,int2,i1,i2,htype,COM
  integer :: ntot,npp,npn,nnn,count,i
  type(c_ptr) :: filehandle,rx
  character(200) :: string_in, fixed
  
  hw = H%hospace 
  filehandle = gzOpen(trim(TBME_DIR)//trim(adjustl(intfile))//achar(0),"r"//achar(0))  
  
  COM = 0
  if (htype == 1) COM = 1 ! remove center of mass hamiltonian? 
  
  N = jbas%total_orbits
  
  ! the first line is the number of matrix elements... I hope. 
  string_in = read_morten_gz(filehandle) 
  
  do i = 3, 20
     if (string_in(i:i+3) == 'XXXX') exit
  end do 
  fixed = string_in(1:i-2) 
  read(fixed,'(I49)') ntot  
  write(*,'(A11,I20)') '# of TBME: ',ntot

  do count = 1, ntot

     string_in = read_morten_gz(filehandle) 
     string_in = adjustl(string_in)

     do i = 20, 80
        if (string_in(i:i+3) == 'XXXX') exit
     end do
     if (string_in(1:1) == '-') then 
        fixed = '  '//string_in(1:i-2) 
     else 
        fixed = '   '//string_in(1:i-2) 
     end if 

     read(fixed(1:4),'(I4)') Tz
     read(fixed(5:8),'(I4)') Par
     read(fixed(9:12),'(I4)') J
     read(fixed(13:16),'(I4)') a
     read(fixed(17:20),'(I4)') b
     read(fixed(21:24),'(I4)') c
     read(fixed(25:28),'(I4)') d     
     read(fixed(29:49),'(e20.6)') V

     g2 = r1_r2( a, b, c, d, J ,jbas ) ! morten and Koshiroh are inconsistent with their definitions
     g3 = p1_p2( a, b, c, d, J ,jbas ) ! morten and Koshiroh are inconsistent with their definitions
     
     ! g1 is COM expectation value, NOT CALCULATED WITH SCOTT'S CODE
     ! g2 is r1*r2 ME 
     ! g3 is p1*p2 ME 
     g1 = (g3 + H%com_hw**2 /hw**2 *g2)*H%lawson_beta ! lawson term    

     V = V + (g1 - g3*COM) *hw/(H%Aneut + H%Aprot) ! center of mass correction
          
     q = block_index(J,Tz,Par) 
     
     C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
     C2 = jbas%con(c)+jbas%con(d) + 1
    
     qx = C1*C2
     qx = qx + adjust_index(qx)   !Vpppp nature  

     ! get the indeces in the correct order
     pre = 1
     if ( a > b )  then 
        
        x = bosonic_tp_index(b,a,N) 
        j_min = jbas%xmap(x)%Z(1)  
        i1 = jbas%xmap(x)%Z( (J-j_min)/2 + 2) 
        pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J)/2 ) 
     else
       ! if (a == b) pre = pre * sqrt( 2.d0 )
       
        x = bosonic_tp_index(a,b,N) 
        j_min = jbas%xmap(x)%Z(1)  
        i1 = jbas%xmap(x)%Z( (J-j_min)/2 + 2) 
     end if
  
     if (c > d)  then     
        
        x = bosonic_tp_index(d,c,N) 
        j_min = jbas%xmap(x)%Z(1)  
        i2 = jbas%xmap(x)%Z( (J-j_min)/2 + 2) 
        
        pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J)/2 ) 
     else 
       ! if (c == d) pre = pre * sqrt( 2.d0 )
      
        x = bosonic_tp_index(c,d,N) 
        j_min = jbas%xmap(x)%Z(1)  
        i2 = jbas%xmap(x)%Z( (J-j_min)/2 + 2) 
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
  rx = gzclose(filehandle)
      
end subroutine read_gz
!=============================================================================
!=============================================================================
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
!  read(eMaxchr,'(I2)') eMax
!  eMaxchr = adjustl(eMaxchr) 
  !open(unit=24,file=trim(SP_DIRECTORY_LIST(1))//'hk'//trim(eMaxchr)//'_lmax10.sps')
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
!  print*, trim(SP_DIRECTORY_LIST(1))//'hk'//trim(eMaxchr)//'Lmax10.sps'
  STOP
end subroutine
!=====================================================================================
!=====================================================================================
subroutine read_me2j_interaction(H,jbas,jbx,htype) 
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
  type(c_ptr) :: buf,buf2,buf3,hndle,hndle2,hndle3,rx
  integer(c_int) :: sz,sz2,sz3
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
  ! allocate array to store positions of matrix elements

  ! move in increments of two, because I use a pn basis,
  ! heiko's is isospin coupled (half the states) 

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


     ! pipj 
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
end do !end sums over isospin  
end do ! end sum over j 


end do  !end sums over Tcoupled lables
end do
end do
end do
close(7)

end subroutine read_me2j_interaction

subroutine read_me2b_interaction(H,jbas,htype,rr,pp,Lawson) 
  use gzipmod
  implicit none 
  
  integer :: nlj1,nlj2,nnlj1,nnlj2,j,T,Mt,nljMax,endpoint,j_min,j_max,htype,Lmax
  integer :: l1,l2,ll1,ll2,j1,j2,jj1,jj2,Ntot,i,q,bospairs,qx,ta,tb,tc,td,bMax
  integer :: eMax,iMax,jmax,jmin,JT,a,b,c,d,C1,C2,i1,i2,COM,x,PAR,endsz,aMax
  integer :: t1,t2,lj1,lj2,n1,n2,Pi,Tz,AA,BB,qq,iq,jq,a_hh,a_ph,a_pp
  integer,allocatable,dimension(:) :: indx , nMax_lj
  real(8),allocatable,dimension(:) :: ME,MEpp,MErr,me_fromfile,ppff,rrff
  real(8) :: V,g1,g2,g3,hw,pre2,pre,sm
  type(spd) :: jbas 
  type(sq_op) :: H,stors
  type(sq_op),optional :: pp,rr
  logical :: pp_calc,rr_calc,file_there,noteffedup
  character(1),optional :: Lawson
  character(2) :: eMaxchr
  character(200) :: me1bfile
  integer :: lj,twol,twoj,ljMax,idx,idxx
  integer,allocatable,dimension(:,:) :: SPBljs 
  type(c_ptr) :: buf,buf2,buf3,hndle,hndle2,hndle3,rx
  integer(c_int) :: sz,sz2,sz3
  character(kind=C_CHAR,len=200) :: buffer,buffer2,buffer3

  
  hw = jbas%hw
  rr_calc = .false.
  pp_calc = .false. 
  Ntot = jbas%total_orbits
  Lmax = maxval(jbas%ll) 
  eMax = maxval(jbas%e)
! populate lj array
  lj = 0
  do twol = 0, 2 * Lmax , 2
     do  twoj = abs(twol - 1) , twol+1 , 2
        lj=lj+1
     end do 
  end do 
  ljMax = lj 
  allocate(SPBljs(lj,2)) 
  allocate(nMax_lj(lj))
  
  lj = 0
  do twol = 0, 2 * Lmax , 2
     do  twoj = abs(twol - 1) , twol+1 , 2
        lj=lj+1
        SPBljs(lj,1) = twol
        sPBljs(lj,2) = twoj
        nMax_lj(lj) = (eMax - twol/2)/2
     end do
  end do
  
  allocate(stors%mat(H%nblocks))


 
  ! using zlib c library, which is bound with fortran in file "gzipmod.f90" 
  
  ! I don't know why you have to tack on those //achars(0) but it seems nessecary 
  if ( present(Lawson) ) then 
     hndle=gzOpen(trim(TBME_DIR)//"O16_Hcm_eMax10_hwHO020.ham0.me2b.gz"//achar(0),"r"//achar(0)) 
  else 
     hndle=gzOpen(trim(TBME_DIR)//trim(adjustl(intfile))//achar(0),"r"//achar(0)) 
  end if 
  
! here is where we start dealing with the two body piece
  
  sz=200
  
  buf=gzGets(hndle,buffer,sz) 
  buf=gzGets(hndle,buffer,sz) 
  buf=gzGets(hndle,buffer,sz) 
  
  read(buffer(6:9),'(I4)',iostat=qq) bMax 
  if (qq .ne. 0 ) then 
     read(buffer(6:8),'(I3)',iostat=qq) bMax 
     if (qq .ne. 0 ) then 
        read(buffer(6:7),'(I2)',iostat=qq) bMax
     end if 
  end if

 q = 0
! heiko's code calls protons 1 and neutrons 0

do Tz = 1 , -1, -1  
  do Pi = 0,1
     do JT = 0, 2*jbas%Jtotal_max,2 
        if ((Lmax == eMax) .and. (JT == 2*jbas%Jtotal_max)&
             .and. (Abs(Tz)==1)) cycle
        if ((JT == 2*jbas%Jtotal_max) .and. (Pi==1)) cycle
        q = q+1
     
        stors%mat(q)%lam(1) = JT
        stors%mat(q)%lam(2) = Pi
        stors%mat(q)%lam(3) = Tz
        
        select case ( Tz)
           case ( -1 ) 
              t1 = -1
              t2 = -1
           case ( 0 ) 
              t1 = 1 
              t2 = -1
           case ( 1 ) 
              t1 = 1
              t2 = 1 
        end select
                 
        a = 0
        a_hh = 0
        a_ph = 0 
        a_pp = 0
        do lj1 = 1, ljMax
           do lj2 = 1, ljMax

              j1 = SPBljs(lj1,2) 
              j2 = SPBljs(lj2,2)
              l1 = SPBljs(lj1,1)/2
              l2 = SPBljs(lj2,1)/2
           
              if ( ( JT < abs(j1-j2) ) .or. (JT > j1 + j2) ) cycle
              if ( mod(l1 + l2 ,2 ) .ne.Pi ) cycle 

              
              do n1 = 0,nMax_lj(lj1)
                 idx = (lj1-1) * (nMax_lj(1) +1 ) +n1 
                 do n2 = 0,nMax_lj(lj2) 
                    idxx = (lj2-1) * (nMax_lj(1) +1 ) +n2                 
                 
                    if ( (Tz .ne. 0) .and. (idx > idxx) ) cycle
                    if ( (mod(JT/2,2) == 1) .and. (lj1==lj2) .and. &
                         (n1==n2) .and. (Tz .ne. 0) ) cycle
                  
                    ! now search for sp labels
                    do i = 1, jbas%total_orbits 
                       if ( jbas%jj(i) .ne. j1 ) cycle
                       if ( jbas%nn(i) .ne. n1 ) cycle
                       if ( jbas%ll(i) .ne. l1 ) cycle
                       if ( jbas%itzp(i) .ne. t1 ) cycle                     
                       exit
                    end do
                  
                    do j = 1, jbas%total_orbits 
                       if ( jbas%jj(j) .ne. j2 ) cycle
                       if ( jbas%nn(j) .ne. n2 ) cycle
                       if ( jbas%ll(j) .ne. l2 ) cycle
                       if ( jbas%itzp(j) .ne. t2 ) cycle                     
                       exit
                    end do
                  
                    a = a + 1
                    select case(jbas%con(i) + jbas%con(j))
                       case(0)
                          a_pp = a_pp + 1
                       case(1)
                          a_ph = a_ph + 1
                       case(2)
                          a_hh = a_hh + 1
                    end select

                 end do
              end do
           end do
        end do
        
    
        stors%mat(q)%npp = a_pp 
        stors%mat(q)%nph = a_ph
        stors%mat(q)%nhh = a_hh
        stors%mat(q)%ntot = a 

    
     end do
  end do
end do

sz = 20
 q = 0
! heiko's code calls protons 1 and neutrons 0

do Tz = 1 , -1, -1  
  do Pi = 0,1
     do JT = 0, 2*jbas%Jtotal_max,2 
        if ((Lmax == eMax) .and. (JT == 2*jbas%Jtotal_max)&
             .and. (Abs(Tz)==1)) cycle
        if ((JT == 2*jbas%Jtotal_max) .and. (Pi==1)) cycle
        q = q+1
     
        buf=gzGets(hndle,buffer,sz) 
  
        read(buffer(10:16),'(I6)') aMax 
        
 
        allocate(stors%mat(q)%qn(1)%Y( aMax+1, 2) ) 

        select case ( Tz)
           case ( -1 ) 
              t1 = -1
              t2 = -1
           case ( 0 ) 
              t1 = 1 
              t2 = -1
           case ( 1 ) 
              t1 = 1
              t2 = 1 
        end select
                 
        a = 0
        a_hh = 0
        a_ph = 0 
        a_pp = 0
        do lj1 = 1, ljMax
           do lj2 = 1, ljMax

              j1 = SPBljs(lj1,2) 
              j2 = SPBljs(lj2,2)
              l1 = SPBljs(lj1,1)/2
              l2 = SPBljs(lj2,1)/2
           
              if ( ( JT < abs(j1-j2) ) .or. (JT > j1 + j2) ) cycle
              if ( mod(l1 + l2 ,2 ) .ne.Pi ) cycle 

              
              do n1 = 0,nMax_lj(lj1)
                 idx = (lj1-1) * (nMax_lj(1) +1 ) +n1 
                 do n2 = 0,nMax_lj(lj2) 
                    idxx = (lj2-1) * (nMax_lj(1) +1 ) +n2                 
                 
                    if ( (Tz .ne. 0) .and. (idx > idxx) ) cycle
                    if ( (mod(JT/2,2) == 1) .and. (lj1==lj2) .and. &
                         (n1==n2) .and. (Tz .ne. 0) ) cycle
                  
                    ! now search for sp labels
                    do i = 1, jbas%total_orbits 
                       if ( jbas%jj(i) .ne. j1 ) cycle
                       if ( jbas%nn(i) .ne. n1 ) cycle
                       if ( jbas%ll(i) .ne. l1 ) cycle
                       if ( jbas%itzp(i) .ne. t1 ) cycle                     
                       exit
                    end do
                  
                    do j = 1, jbas%total_orbits 
                       if ( jbas%jj(j) .ne. j2 ) cycle
                       if ( jbas%nn(j) .ne. n2 ) cycle
                       if ( jbas%ll(j) .ne. l2 ) cycle
                       if ( jbas%itzp(j) .ne. t2 ) cycle                     
                       exit
                    end do
                  
                    
                    select case(jbas%con(i) + jbas%con(j))
                       case(0)
                          a_pp = a_pp + 1   
                          a = stors%mat(q)%nhh + &
                               stors%mat(q)%nph + a_pp
                       case(1)
                          a_ph = a_ph + 1
                          a = stors%mat(q)%nhh + a_ph
                       case(2)
                          a_hh = a_hh + 1
                          a = a_hh 
                    end select
                    
                    stors%mat(q)%qn(1)%Y(a,1) = i
                    stors%mat(q)%qn(1)%Y(a,2) = j
                    
                      
                 end do
              end do
           end do
        end do
            
     end do
  end do
end do


! okay for now on there is a space, then a line that specifies the block
! and aMax for the block. We already know that stuff so we will just ignore
! it and read in the matrix elements
sz = 200
qq = 0
noteffedup = .true. 
do Tz = 1 , -1, -1  
  do Pi = 0,1
     do JT = 0, 2*jbas%Jtotal_max,2 

        if ((Lmax == eMax) .and. (JT == 2*jbas%Jtotal_max)&
             .and. (Abs(Tz)==1)) cycle
        if ((JT == 2*jbas%Jtotal_max) .and. (Pi==1)) cycle
        qq = qq+1 ! heikos block index

        q = block_index(JT,Tz,Pi) ! my block index
        
        ! space then label
        buf=gzGets(hndle,buffer,sz) 
        
        if (noteffedup) then  
           buf=gzGets(hndle,buffer,sz)
        end if
        noteffedup=.true. 
        ! ignore
        sz = 200
        

      do 
         
      buf=gzGets(hndle,buffer,sz)
      
      read(buffer(1:1),'(I1)',iostat=iq) AA
      if (iq .ne. 0 ) then 
         ! because it's damned near impossible to 
         ! figure out where the block ends
         noteffedup=.false. 
         exit
      end if
         ! figure out where the spaces are that separate things 
      i = 1
      ! first space
       do 
         if ( buffer(i:i) == ' ' ) then
            read(buffer(1:i-1),'(I5)')  AA 
            i = i + 2
            j = i 
            exit
         end if
         i = i + 1
      end do
      ! second space
      do 
         if ( buffer(i:i) == ' ' ) then 
            read(buffer(j:i-1),'(I5)') BB 
            i = i + 2
            exit
         end if
         i = i + 1
      end do
      

      AA = AA + 1
      BB = BB + 1

      
      if ( buffer(i:i) == '-' ) then 
         ! negative number
         read(buffer(i:i+10), '( f11.8 )' )  V 
      else 
         ! positive
         read(buffer(i:i+9), '( f10.8 )' )  V 
      end if 
      

      ! oTay should have the matrix element now. 
     
      ! indeces     
      a = stors%mat(qq)%qn(1)%Y(AA,1)
      b = stors%mat(qq)%qn(1)%Y(AA,2)      
      c = stors%mat(qq)%qn(1)%Y(BB,1)
      d = stors%mat(qq)%qn(1)%Y(BB,2)      
      
      ! i think the scaling and COM subtraction have already been done

!!!!=========================================================================
      ! start the classical method of sorting these into my arrays now
!!!!=========================================================================     
      C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
      C2 = jbas%con(c)+jbas%con(d) + 1
    
      qx = C1*C2
      qx = qx + adjust_index(qx)   !Vpppp nature  
      
      ! get the indeces in the correct order
      pre = 1.d0
      if ( a > b )  then 

         x = bosonic_tp_index(b,a,Ntot) 
         j_min = jbas%xmap(x)%Z(1)  
         i1 = jbas%xmap(x)%Z( (JT-j_min)/2 + 2) 
         pre = pre * (-1.)**( 1 + (jbas%jj(a) + jbas%jj(b) -JT)/2 ) 

      else
         if (a == b) pre = pre / sqrt( 2.d0 )

         x = bosonic_tp_index(a,b,Ntot) 

         j_min = jbas%xmap(x)%Z(1)  
         i1 = jbas%xmap(x)%Z( (JT-j_min)/2 + 2) 
      end if

      if (c > d)  then     

         x = bosonic_tp_index(d,c,Ntot)

         j_min = jbas%xmap(x)%Z(1)  

         i2 = jbas%xmap(x)%Z( (JT-j_min)/2 + 2) 

         pre = pre * (-1.)**( 1 + (jbas%jj(c) + jbas%jj(d) -JT)/2 ) 

      else 
         if (c == d) pre = pre / sqrt( 2.d0 )

         x = bosonic_tp_index(c,d,Ntot) 
         j_min = jbas%xmap(x)%Z(1)  
         i2 = jbas%xmap(x)%Z( (JT-j_min)/2 + 2) 

      end if

      ! kets/bras are pre-scaled by sqrt(2) if they 
      ! have two particles in the same sp-shell

      ! get the units right. I hope 

      if ((qx == 1) .or. (qx == 5) .or. (qx == 4)) then 
         H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
         H%mat(q)%gam(qx)%X(i1,i2)  = V *pre

         if (rr_calc) then 
            STOP 'Darn!!, this is not implemented yet' 
            rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot)
            rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot)
         end if

         if (pp_calc) then 
            STOP 'Darn!!, this is not implemented yet' 
            pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot)
            pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot)
         end if


      else if (C1>C2) then
         H%mat(q)%gam(qx)%X(i2,i1)  = V *pre

         if (rr_calc) then 
            STOP 'Darn!!, this is not implemented yet' 
            rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot) 
         end if

         if (pp_calc) then 
            STOP 'Darn!!, this is not implemented yet' 
            pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot) 
         end if

      else
         H%mat(q)%gam(qx)%X(i1,i2) = V * pre

         if (rr_calc) then 
            STOP 'Darn!!, this is not implemented yet' 
            rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot) 
         end if

         if (pp_calc) then 
            STOP 'Darn!!, this is not implemented yet' 
            pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot) 
         end if

      end if
      ! I shouldn't have to worry about hermiticity here, input is assumed to be hermitian

      if (AA == stors%mat(qq)%ntot) then 
         if (BB == stors%mat(qq)%ntot) then 
            exit
         end if
      end if

   end do 
end do   ! end loops over conserved quantities
end do 
end do 



! i guess we are done with the two body piece

14 me1bfile = intfile(1:len(trim(intfile))-5)//'1b.gz'
  
  if ( present(Lawson) ) then 
     hndle=gzOpen(trim(TBME_DIR)//"O16_Hcm_eMax10_hwHO020.ham0.me1b.gz"//achar(0),"r"//achar(0)) 
  else 
     hndle=gzOpen(trim(TBME_DIR)//trim(adjustl(me1bfile))//achar(0),"r"//achar(0)) 
  end if 

  sz=200

  ! read verion line, and then some integer
  buf=gzGets(hndle,buffer,sz) 
  ! the integer probably has to do with the file size
  buf=gzGets(hndle,buffer,sz) 
 
  
  read(buffer(1:4),'(I4)',iostat=qq) aMax 
  if (qq .ne. 0 ) then 
     read(buffer(1:3),'(I3)',iostat=qq) aMax 
     if (qq .ne. 0 ) then 
        read(buffer(1:2),'(I2)',iostat=qq) aMax
     end if 
  end if
  sz = 20
  ! I assume this is the zero body piece right here
  buf=gzGets(hndle,buffer,sz) 
 
  ! lets get it right since it can have up to 4 digits before the decimal
  if (buffer(4:4) == '.') then 
     read(buffer(1:10),'(f10.6)') H%E0 
  else if (buffer(5:5) == '.') then 
     read(buffer(1:11),'(f11.6)') H%E0 
  else if (buffer(6:6) == '.') then 
     read(buffer(1:12),'(f12.6)') H%E0
  else 
     print*, 'what the darn is going on with the me1b file? ' 
  end if

  ! now lets read the 1 body piece


  sz=200

  do a= 1, aMax


     buf=gzGets(hndle,buffer,sz) 

     read(buffer(2:2),'(I1)') t1
     read(buffer(4:5),'(I2)') lj
     read(buffer(8:9),'(I2)') n1
     read(buffer(11:12),'(I2)') n2

     t1 = (-2*t1+ 1)  

     lj = lj + 1
     l1 = SPBljs(lj,1)/2
     j1 = SPBljs(lj,2)


     read(buffer(15:28),'(f14.10)') V


     ! V is now the one body matrix element
     ! now search for sp labels
     do i = 1, jbas%total_orbits 
        if ( jbas%jj(i) .ne. j1 ) cycle
        if ( jbas%nn(i) .ne. n1 ) cycle
        if ( jbas%ll(i) .ne. l1 ) cycle
        if ( jbas%itzp(i) .ne. t1 ) cycle                     
        exit
     end do

     do j = 1, jbas%total_orbits 
        if ( jbas%jj(j) .ne. j1 ) cycle
        if ( jbas%nn(j) .ne. n2 ) cycle
        if ( jbas%ll(j) .ne. l1 ) cycle
        if ( jbas%itzp(j) .ne. t1 ) cycle                     
        exit
     end do

     ! okay now I have my indeces 


     if (jbas%con(i) + jbas%con(j) == 2 ) then 

        !fhh 

        H%fhh( hb4(i)+1 , hb4(j)+1 ) = V 

     else if (jbas%con(i) + jbas%con(j) == 0 ) then 

        !fpp 

        H%fpp( pb4(i)+1 , pb4(j)+1 ) = V 

     else if ((jbas%con(i)==0) .and. (jbas%con(j) == 1) ) then 
        !fph
        H%fph( pb4(i)+1 , hb4(j)+1 ) = V
     end if


  end do

  rx = gzClose(hndle)

 end subroutine read_me2b_interaction
!==========================================================
subroutine read_me3j(store_3b,jbas,jbx,eMax,lmax) 
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
  type(c_ptr) :: buf,buf2,buf3,hndle,hndle2,hndle3,rx
  integer(c_int) :: sz,sz2,sz3
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

     hndle=gzOpen(trim(THREE_DIR)//trim(adjustl(threebody_file))//achar(0),"r"//achar(0))   
  
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
     
     open(unit=77,file=trim(THREE_DIR)//trim(adjustl(threebody_file)), form = 'unformatted', &
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
 
subroutine write_onebody_human(Op,jbas,lab)    
  ! write Op to readable file  
  implicit none 

  type(sq_op) :: Op
  type(spd) :: jbas
  integer :: n1,l1,jj1,tt1,  n2,l2,jj2,tt2, AA,BB,ii,jj,JT
  real(8) :: xx,delE,delF ,yy
  character(4) :: lab
  
  
  open(unit=78,file=trim(OUTPUT_DIR)//trim(prefix)//"_unnorm_"//&
       trim(adjustl(lab))//"_1b.dat")
  open(unit=79,file=trim(OUTPUT_DIR)//trim(prefix)//"_"//&
       trim(adjustl(lab))//"_1b.dat")
  
  write(78,"(A)")  "####  ZERO-BODY OFFSET  ####"
  write(79,"(A)")  "####  ZERO-BODY VALUE  ####"

  delE=Op%E0
  do ii = 1,jbas%total_orbits
     if (jbas%con(ii) == 0 ) cycle
     delE = delE - f_elem(ii,ii,Op,jbas)* (jbas%jj(ii)+1.d0)
  end do
  
  do ii = 1,jbas%total_orbits
     if (jbas%con(ii) == 0 ) cycle
     do jj = 1,jbas%total_orbits
        if (jbas%con(jj) == 0 ) cycle
        
        do JT = abs(jbas%jj(ii)-jbas%jj(jj)) ,jbas%jj(ii)+jbas%jj(jj) ,2 
           delE = delE + 0.5*v_elem(ii,jj,ii,jj,JT,Op,jbas)*(JT+1.d0)
        end do
  
     end do
  end do

  write(78,"(f20.10)") delE
  write(79,"(f20.10)") Op%E0

  write(78,"(A)")  "####  ONE-BODY MATRIX ELEMENTS  ####"
  write(78,"(A)")  "####  T(2,1) == T(1,2) "
  write(78,"(A)")  "# n1  l1  2*j1  2*tz1     n2  l2  2*j2  2*tz2      T(1,2)"           

  write(79,"(A)")  "####  ONE-BODY MATRIX ELEMENTS  ####"
  write(79,"(A)")  "####  f(2,1) == f(1,2) "
  write(79,"(A)")  "# n1  l1  2*j1  2*tz1     n2  l2  2*j2  2*tz2      f(1,2)"           

  
  do AA = 1,jbas%total_orbits

     n1 = jbas%nn(AA)
     l1 = jbas%ll(AA)
     jj1 = jbas%jj(AA)
     tt1 = jbas%itzp(AA) 

     do BB = AA,jbas%total_orbits

        n2 = jbas%nn(BB)
        l2 = jbas%ll(BB)
        jj2 = jbas%jj(BB)
        tt2 = jbas%itzp(BB) 

        if (jj1 .ne. jj2) cycle         
        if (tt1 .ne. tt2) cycle 

        delF  = 0.d0         
        do ii = 1,jbas%total_orbits
           if (jbas%con(ii) == 0 ) cycle
           do JT = abs(jbas%jj(ii) - jj2),jbas%jj(ii)+jj2,2  
              delF = delF + v_elem(AA,ii,BB,ii,JT,Op,jbas)*(JT+1.d0)/(jj2+1.d0)
           end do
        end do

        yy = f_elem(AA,BB,Op,jbas)
        xx = yy - delF

        if (abs(yy) > 1e-10) then
           write(79,"(4(I5),I8,3(I5),f20.7)") n1,l1,jj1,tt1,n2,l2,jj2,tt2, yy
        end if

        if (abs(xx) > 1e-10) then
           write(78,"(4(I5),I8,3(I5),f20.7)") n1,l1,jj1,tt1,n2,l2,jj2,tt2, xx
        end if

     end do
  end do

  close(78)
  close(79)
end subroutine write_onebody_human 


subroutine write_onebody_tensor_human(Op,jbas,lab)    
  ! write Op to readable file  
  implicit none 

  type(sq_op) :: Op
  type(spd) :: jbas
  integer :: n1,l1,jj1,tt1,  n2,l2,jj2,tt2, AA,BB,  ii,jj,Jtot1,Jtot2
  real(8) :: xx,delE ,d6ji,delF,yy
  character(4) :: lab
  
  
  open(unit=78,file=trim(OUTPUT_DIR)//trim(prefix)//"_unnorm_"//&
       trim(adjustl(lab))//"_1b.dat")
  open(unit=79,file=trim(OUTPUT_DIR)//trim(prefix)//"_"//&
       trim(adjustl(lab))//"_1b.dat")

  
  write(78,"(A)")  "####  ZERO-BODY OFFSET  ####"
  write(79,"(A)")  "####  ZERO-BODY VALUE  ####"
  write(79,"(f20.10)") Op%E0
  write(78,"(f20.10)") Op%E0

  write(78,"(A)")  "####  ONE-BODY REDUCED MATRIX ELEMENTS  ####"
  write(78,"(A)")  "####  T(2,1) == T(1,2)  * (-1) ** (j1-j2) "
  write(78,"(A)")  "# n1  l1  2*j1  2*tz1     n2  l2  2*j2  2*tz2      T(1,2)"           

  write(79,"(A)")  "####  ONE-BODY REDUCED MATRIX ELEMENTS  ####"
  write(79,"(A)")  "####  f(2,1) == f(1,2)  * (-1) ** (j1-j2) "
  write(79,"(A)")  "# n1  l1  2*j1  2*tz1     n2  l2  2*j2  2*tz2      f(1,2)"           

  
  do AA = 1,jbas%total_orbits

     n1 = jbas%nn(AA)
     l1 = jbas%ll(AA)
     jj1 = jbas%jj(AA)
     tt1 = jbas%itzp(AA) 

     do BB = AA,jbas%total_orbits

        n2 = jbas%nn(BB)
        l2 = jbas%ll(BB)
        jj2 = jbas%jj(BB)
        tt2 = jbas%itzp(BB) 

        if (tt1 .ne. tt2) cycle 

        delF  = 0.d0         
        do ii = 1,jbas%total_orbits
           if (jbas%con(ii) == 0 ) cycle
           do Jtot1 = abs(jbas%jj(ii) - jj1),jbas%jj(ii)+jj1,2  
              do Jtot2 = abs(jbas%jj(ii) - jj2),jbas%jj(ii)+jj2,2
                 delF = delF + tensor_elem(AA,ii,BB,ii,Jtot1,Jtot2,Op,jbas)&
                      *sqrt((Jtot1+1.d0)*(Jtot2+1.d0)) *(-1)**((jbas%jj(ii)+jj1+Jtot2+Op%Rank)/2)&
                      *d6ji(Jtot1,Jtot2,Op%rank,jj2,jj1,jbas%jj(ii)) 
              end do
           end do
        end do

        yy = f_tensor_elem(AA,BB,Op,jbas)
        xx = yy - delF

        if (abs(yy) > 1e-10) then
           write(79,"(4(I5),I8,3(I5),f20.7)") n1,l1,jj1,tt1,n2,l2,jj2,tt2, yy
        end if
        
        if (abs(xx) > 1e-10) then
           write(78,"(4(I5),I8,3(I5),f20.7)") n1,l1,jj1,tt1,n2,l2,jj2,tt2, xx
        end if

     end do
  end do

  close(78)
  close(79)
end subroutine write_onebody_tensor_human 



subroutine write_twobody_human(Op,jbas,lab)    
  ! write Op to readable file  
  implicit none 

  type(sq_op) :: Op
  type(spd) :: jbas
  integer :: n1,l1,jj1,tt1,  n2,l2,jj2,tt2, AA,BB,CC,DD,JTOT
  integer :: n3,l3,jj3,tt3,  n4,l4,jj4,tt4,TTot, PAR,jmin,jmax
  real(8) :: xx
  character(4) :: lab
  
  
  open(unit=78,file=trim(OUTPUT_DIR)//trim(prefix)//"_unnorm_"//&
       trim(adjustl(lab))//"_2b.dat")

  write(78,"(A)")  "####  SINGLE-PARTICLE LABELS  ####"
  write(78,"(A)")  "### index    n    l    2*j    2*tz      occ. "

  do AA=1,jbas%total_orbits
     write(78,"(6(I5))") AA,  jbas%nn(AA) , jbas%ll(AA),  jbas%jj(AA) , jbas%itzp(AA)  ,jbas%con(AA)
  end do
     
  write(78,"(A)")  "####  TWO-BODY UNNORMALIZED MATRIX ELEMENTS  ####"
  write(78,"(A)")  "####  GAMMA((ab)J (cd)J) ==  GAMMA((cd)J (ab)J) "
  write(78,"(A)")  "# a   b   c   d    2*J    GAMMA((ab)J (cd)J)"           
  
  do AA = 1,jbas%total_orbits

     n1 = jbas%nn(AA)
     l1 = jbas%ll(AA)
     jj1 = jbas%jj(AA)
     tt1 = jbas%itzp(AA) 

     do BB = AA,jbas%total_orbits

        n2 = jbas%nn(BB)
        l2 = jbas%ll(BB)
        jj2 = jbas%jj(BB)
        tt2 = jbas%itzp(BB)

        PAR = mod(l1+l2,2)
        Ttot = tt1+tt2 
        
        do CC = 1,jbas%total_orbits
             
           n3 = jbas%nn(CC)
           l3 = jbas%ll(CC)
           jj3 = jbas%jj(CC)
           tt3 = jbas%itzp(CC) 
           
           do DD = CC,jbas%total_orbits
              
              n4 = jbas%nn(DD)
              l4 = jbas%ll(DD)
              jj4 = jbas%jj(DD)
              tt4 = jbas%itzp(DD) 

              if ((tt3+tt4) .ne. Ttot) cycle
              if (mod(l3+l4,2) .ne. PAR) cycle

              jmin = max(abs(jj1-jj2),abs(jj3-jj4))
              jmax = min(jj1+jj2,jj3+jj4) 


              do JTot= jmin,jmax,2 
                 xx = v_elem(AA,BB,CC,DD,Jtot,Op,jbas)   
             
                 if (abs(xx) > 1e-10) then
                    write(78,"(5(I5),f20.7)") aa,bb,cc,dd,Jtot,xx
                 end if
              end do
              
           end do
        end do
     end do
  end do

  close(78)

end subroutine write_twobody_human



subroutine write_twobody_tensor_human(Op,jbas,lab)    
  ! write Op to readable file  
  implicit none 

  type(sq_op) :: Op
  type(spd) :: jbas
  integer :: n1,l1,jj1,tt1,  n2,l2,jj2,tt2, AA,BB,CC,DD,JTOT1,JTOT2
  integer :: n3,l3,jj3,tt3,  n4,l4,jj4,tt4,TTot, PAR,jmin,jmax
  real(8) :: xx
  character(4) :: lab
  
  
  open(unit=78,file=trim(OUTPUT_DIR)//trim(prefix)//"_unnorm_"//&
       trim(adjustl(lab))//"_2b.dat")

  write(78,"(A)")  "####  SINGLE-PARTICLE LABELS  ####"
  write(78,"(A)")  "### index    n    l    2*j    2*tz      occ. "

  do AA=1,jbas%total_orbits
     write(78,"(6(I5))") AA,  jbas%nn(AA) , jbas%ll(AA),  jbas%jj(AA) , jbas%itzp(AA)  ,jbas%con(AA)
  end do
     
  write(78,"(A)")  "####  TWO-BODY UNNORMALIZED REDUCED MATRIX ELEMENTS  ####"
  write(78,"(A)")  "####  GAMMA((ab)J1 (cd)J2) ==  GAMMA((cd)J2 (ab)J1)*(-1)**(J1-J2) "
  write(78,"(A)")  "# a   b   c   d    2*J1  2*J2    GAMMA((ab)J1 (cd)J2)"           
  
  do AA = 1,jbas%total_orbits

     n1 = jbas%nn(AA)
     l1 = jbas%ll(AA)
     jj1 = jbas%jj(AA)
     tt1 = jbas%itzp(AA) 

     do BB = AA,jbas%total_orbits

        n2 = jbas%nn(BB)
        l2 = jbas%ll(BB)
        jj2 = jbas%jj(BB)
        tt2 = jbas%itzp(BB)

        PAR = mod(l1+l2,2)
        Ttot = tt1+tt2 
        
        do CC = 1,jbas%total_orbits
             
           n3 = jbas%nn(CC)
           l3 = jbas%ll(CC)
           jj3 = jbas%jj(CC)
           tt3 = jbas%itzp(CC) 
           
           do DD = CC,jbas%total_orbits
              
              n4 = jbas%nn(DD)
              l4 = jbas%ll(DD)
              jj4 = jbas%jj(DD)
              tt4 = jbas%itzp(DD) 

              if ((tt3+tt4) .ne. Ttot) cycle
              if (mod(l3+l4,2) .ne. PAR) cycle

              jmin = max(abs(jj1-jj2),abs(jj3-jj4))
              jmax = min(jj1+jj2,jj3+jj4) 


              do JTot1= abs(jj1-jj2),jj1+jj2,2
                 do JTot2= abs(jj3-jj4),jj3+jj4,2
                    xx = tensor_elem(AA,BB,CC,DD,Jtot1,Jtot2,Op,jbas)   
             
                    if (abs(xx) > 1e-10) then
                       write(78,"(6(I5),f20.7)") aa,bb,cc,dd,Jtot1,Jtot2,xx
                    end if
                 end do
              end do
           end do
        end do
     end do
  end do

  close(78)

end subroutine write_twobody_tensor_human


end module interaction_IO




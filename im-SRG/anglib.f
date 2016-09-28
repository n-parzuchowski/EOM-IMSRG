c     horrifying but fast and functional ang. mom. routines
c     abandon hope ye who enter here

      function  dcg(aj1,am1,aj2,am2,aj3,am3c)
c *****
c ***    revised  by  t. otsuka  in  order  to  fit
c ***      facom    computer
c *****
      implicit  real*8  (a-h,o-z)
      parameter  (nm=52)
      logical  ibc
      common /bcodcm/q(nm,nm), p(nm,nm)
cc      data  im / -2147483647 /
      data  ibc / .false. /
cc      ie(k) = iand(im,k)
c   ** q(i,m) (i.ge.m) binomial coef. c(i-1,m-1)
c   ** q(m,i) (i.ge.m) q(i,m)**(-1/2)
c      p(i,j) = q(j,i)
c** **
c]    write(6,660)  aj1,aj2,aj3,am1,am2,am3c
c]660 format(/t6,'input j & m :',6f8.3/)
      j1= aj1 + aj1 + 0.1d0
      j2 = aj2 + aj2 + 0.1d0
      j3 = aj3 + aj3 + 0.1d0
      m1 = am1 + am1 + dsign(0.1d0, am1)
      m2 = am2 + am2 + dsign(0.1d0, am2)
      m3 = - ( am3c + am3c + dsign(0.1d0, am3c) )
c]    write(6,640) j1,j2,j3,m1,m2,m3
c]640 format(/t6,'transfmed  integers :',6i5/)
10    if( ibc )  goto 12
 
      call bcodin
 
      ibc = .true.
12    jt=j1+j2+j3
      l1 = jt-2*j1
      l2 = jt-2*j2
      l3 = jt-2*j3
      j1m = j1-m1
      j2m = j2 + m2
      j3m = j3 + m3
      if(m1+m2+m3.ne.0)  goto  990
c** ** triangular and j-m condition check
      if( (l1.lt.0).or.(l2.lt.0).or.(l3.lt.0).or.
     &    (j1m.lt.0).or.(j2m.lt.0).or.(j3m.lt.0).or.
     &    (j1+j1.lt.0).or.(j2-m2.lt.0).or.(j3-m3.lt.0) ) goto 990
cc      if((ie(l1) + ie(l2) + ie(l3) + ie(j1m) + ie(j2m) + ie(j1+m1)
cc     (  + ie(j2-m2) + ie(j3-m3) + ie(j3m)) .lt. 0 )  goto  990
      jt = jt / 2
      j1m = j1m / 2
      j2m = j2m / 2
      j3m = j3m / 2
      l1 = l1 / 2
      l2 = l2 / 2
      l3 = l3 / 2
      sum = 0
      kf = max0(0,j1m-l2,j2m-l1)+1
      kl = min0(l3,j1m,j2m)+1
      if(kl.gt.nm)go to 900
      do 1 k = kf,kl
c       sum = -sum+q(l3+1,k)*q(l2+1,j1m-k+2)*q(l1+1,j2m-k+2)
        sum = -sum+p(k,l3+1)*p(j1m-k+2,l2+1)*p(j2m-k+2,l1+1)
1     continue
      dcg = sum*q(j1m+1,j1+1)*q(j2m+1,j2+1)*q(j3m+1,j3+1)
     &      *q(l3+1,jt+1)*q(2,jt+2)/(q(l3+1,j2+1)*q(l3+1,j1+1))
      jm=kl+l2-j3m+1
      if(jm-2*(jm/2).ne.0)  dcg  =  - dcg
      dcg = dcg/q(2,j3+2)
      jm=(j1-j2-m3)/2
      if(jm-2*(jm/2).ne.0)  dcg  =  - dcg
c]    write(6,650)  dcg
c]650 format(t6,'cal  cg  in  dcg :',f9.5/)
      return
  900 mm3=-m3
      write(6,602) j1,m1,j2,m2,j3,mm3
  990 dcg = 0
      return
  602 format('dcg/dcgi  out of range  ;(',
     &       3(i4,'/2,'),i4,'/2:',i4,'/2',i4,'/2);')
c --------------------------------------------------------------------
      entry dcgi1(j1i,m1i,j2i,m2i,j3i,m3ci)
c** ** clebsch-gordan coefficient program  (integer  input)
      j1 = j1i
      j2 = j2i
      j3 = j3i
      m1 = m1i
      m2 = m2i
      m3=-m3ci
      go to 10
      end
 
 
 
      subroutine bcodin
      implicit real*8(a-h,o-z)
      parameter  (nm=52, nm2=nm*nm)
      common /bcodcm/q(nm2), p(nm,nm)
      dimension      qx(nm,nm)
      equivalence   (q, qx)
 
c   ** main routine de hitsuyo na common   coded by m.oka
c   **  common /bcodcm/q(nm,nm)
c   ** q(m,i) (i.ge.m) binomial coef. c(i-1,m-1)
c   ** q(i,m) (i.ge.m) q(m,i)**(-1/2)
c   ** p(i,m) = q(m,i)
      data  nmax / 0 /
c %%%    original  form  in  the  kohchan  library .
c %%%      fitted  to  fortran  77  ( june  1982 )
c %%% na(i,j)=i+(j-1)*nm
      na(i,j) = j + (i-1)*nm
c** **
      if(nm.lt.3) go to 24
      if(nmax.eq.nm) return
      if(nmax.ne.0) go to 23
      nmax=nm
      q(1) = 1
      q(2) = 1
      q(nm+1) = 1
      q(nm+2) = 1
      do 10 i=3,nm
      iq = na(i,1)
      iq1 = na(i,i)
      iq2 = na(1,i)
      q( iq )=1
      q(iq1)=1
      q(iq2)=1
      do 10 j=2,i-1
      iq = na(j,i)
      iq1 = na(j-1,i-1)
      iq2 = na(j,i-1)
      q(iq) = q( iq1 ) + q( iq2 )
   10 continue
      do 20 j=3,nm
      do 20 i=1,j-1
      iq = na(j,i)
      iq1 = na(i,j)
      q( iq ) = 1 / dsqrt( q( iq1 ) )
   20 continue
 
      do  240  ka = 1, nm
        do  242  kb = 1, nm
          p( ka, kb ) = qx( kb, ka )
242     continue
240   continue
 
c     write(6,600)  nm
  600 format(/t12,'***  bcodin  called  (nm=',i2,')  ***'/)
      return
   23 write(6,101)  nmax,nm
      go to 25
   24 write(6,100)  nm
   25 return
  100 format(t15,'nmax.lt.3 in (bcodin)  nmax=',i4)
  101 format(t15,'/bcodcm/ length mismatched  old nmax=',i4,
     &       ', new nmax =',i4)
      end
 
 
 
      function dcgi00()
      implicit real*8 (a-h,o-z)
      real*8 qfactl(0:100),w0,w1
      save
 
      qfactl(0)=0.0d0
      do 10 i=1,100
        qfactl(i)=qfactl(i-1)+dlog(1.0d0*i)
10    continue
      dcgi00=0.0d0
      return
 
      entry dcgi(j1,m1,j2,m2,j,m)
c     this ludicrous thing seems to give the CG coefs, 
c     but you need to supply twice the values, and only use half in j1,j2
c     which you then mult by 2.
**** condition check ****
      if ( (j1.lt.0).or.(iabs(m1).gt.j1).or.
     &     (j2.lt.0).or.(iabs(m2).gt.j2).or.
     &     (j.lt.0).or.(iabs(m).gt.j) ) then
        dcgi=0.0d0
        return
      end if
 
      if ( (j.lt.iabs(j1-j2)).or.(j.gt.j1+j2) ) then
        dcgi=0.0d0
        return
      end if
 
      if (m1+m2.ne.m) then
        dcgi=0.0d0
        return
      end if
 
**** calc ****
      w0=0.5d0*(
     &   dlog(1.0d0*(j+1))
     &  +qfactl((j1+j2-j)/2)+qfactl((j1-j2+j)/2)
     &  +qfactl((-j1+j2+j)/2)
     &  -qfactl((j1+j2+j)/2+1)
     &  +qfactl((j1+m1)/2)+qfactl((j1-m1)/2)
     &  +qfactl((j2+m2)/2)+qfactl((j2-m2)/2)
     &  +qfactl((j+m)/2)+qfactl((j-m)/2) )
      w1=0.0d0
      k1=(j1+j2-j)/2
      k2=(j1-m1)/2
      k3=(j2+m2)/2
      k4=(-j+j2-m1)/2
      k5=(-j+j1+m2)/2
      imax=min0(k1,k2,k3)
      imin=max0(0,k4,k5)
      do 100 i=imin,imax
        sign=(-1.0d0)**i
        w1=w1+sign*dexp(w0
     &        -qfactl(i)-qfactl(k1-i)-qfactl(k2-i)-qfactl(k3-i)
     &        -qfactl(i-k4)-qfactl(i-k5))
100   continue
 
      dcgi=w1
      return
      end



      function arand(ir)
****
****   generate uniform (0,1) pesudo-random numbers
****
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(ml=1664501,ia=1229,ic=351750)
      parameter(anorm=1.0d0/ml)

      ir=mod(ia*ir+ic,ml)
      arand=ir*anorm

      return
      end



c
c        for check
c
c      implicit real*8 (a-h,o-z)
c
c      s=0.0d0
c      do 100 j1=0,6,2
c      do 100 j2=0,6,2
c      do 100 j3=0,6,2
c      do 100 l1=4,4
c      do 100 l2=0,6,2
c      do 100 l3=0,6,2
c      s1=d6ji(j1,j2,j3,l1,l2,l3)
c      s2=an6j(j1/2,j2/2,j3/2,l1/2,l2/2,l3/2)
c      s=s+dabs(s1-s2)
c100   continue
c      write(*,*) 's=',s
c      end



      subroutine dfact0()
      implicit real*8 (a-h,o-z)

      common /factor/ dfactl(0:1000)
      save

c      write(*,*) 'dfact0 is called'

c      if (n1.gt.1000) then
c        write(*,*) '[ dfact0 ] too large n1'
c        stop
c      end if
      dfactl(0)=0.0d0
c      do 100 i=1,n1
      do 100 i=1,1000
        dfactl(i)=dfactl(i-1)+dlog(dfloat(i))
100   continue

      return
      end


      function d6ji(j1,j2,j3,l1,l2,l3)
****
****  six j calculation 
****            ( ref. edmonds p.90 )
****
****       ( 2*j1  2*j2  2*j3 )
****       ( 2*l1  2*l2  2*l3 )
****
****   when you use, call dfact0 before!
****
      implicit real*8 (a-h,o-z)
c      data iflag / 0 /
c      save iflag

c      if (iflag.eq.0) then
c        call dfact0(1000)
c        iflag=1
c      end if

**** condition check ****
      if ( (j1.lt.0).or.(j2.lt.0).or.(j3.lt.0).or.
     &     (l1.lt.0).or.(l2.lt.0).or.(l3.lt.0) ) then
        d6ji=0.0d0
        return
c       write(*,*) '[ d6ji ] error. negative j!!'
c       stop
      end if

      if ( (j1.lt.iabs(j2-j3)).or.(j1.gt.j2+j3).or.
     &     (j3.lt.iabs(l1-l2)).or.(j3.gt.l1+l2).or.
     &     (j2.lt.iabs(l1-l3)).or.(j2.gt.l1+l3).or.
     &     (j1.lt.iabs(l2-l3)).or.(j1.gt.l2+l3) ) then
        d6ji=0.0d0
        return
c       write(*,*) '[ d6ji ] error. triangle condition'
c       stop
      end if

      if ( (iand(j1+j2+j3,1).eq.1).or.
     &     (iand(j3+l1+l2,1).eq.1).or.
     &     (iand(j2+l1+l3,1).eq.1).or.
     &     (iand(j1+l2+l3,1).eq.1) ) then
        d6ji=0.0d0
        return
      end if

**** calc ****
      df=d(j1,j2,j3)+d(j1,l2,l3)+d(l1,j2,l3)+d(l1,l2,j3)
      d6ji=wwwwww(df,j1,j2,j3,l1,l2,l3)

      return
      end


      function d(j,k,l)
**** triangle function ****
      implicit real*8 (a-h,o-z)
      common /factor/ dfactl(0:1000)
      save

      d=dfactl((j+k-l)/2)+dfactl((j-k+l)/2)
     &  +dfactl((-j+k+l)/2)-dfactl((j+k+l)/2+1)
      d=d*0.5d0
      return
      end

      function wwwwww(df,j1,j2,j3,l1,l2,l3)
**** summation function ****
      implicit real*8 (a-h,o-z)
      common /factor/ dfactl(0:1000)
      save

      m1=(j1+j2+j3)/2
      m2=(j1+l2+l3)/2
      m3=(l1+j2+l3)/2
      m4=(l1+l2+j3)/2
      m5=(j1+j2+l1+l2)/2
      m6=(j2+j3+l2+l3)/2
      m7=(j3+j1+l3+l1)/2
      imin=max0(m1,m2,m3,m4)
      imax=min0(m5,m6,m7)
      w0=0.0d0
      do 100 i=imin,imax
        sign=1.0d0
        if (iand(i,1).eq.1) sign=-sign
        w0=w0+sign
     &    *dexp(df+dfactl(i+1)
     &    -dfactl(i-m1)-dfactl(i-m2)-dfactl(i-m3)-dfactl(i-m4)
     &    -dfactl(m5-i)-dfactl(m6-i)-dfactl(m7-i))
100   continue
      wwwwww=w0
      return
      end




      function an6j(j1,j2,j3,l1,l2,l3)
****
****       six j calculation for special case of l1 = 2
****          all j and l are integers
****          half integer is not supported
****
      implicit real*8 (a-h,o-z)
**** triangle condition ****
      if ( ( iabs(j1-j2).gt.j3 ).or.( j1+j2.lt.j3 ).or.
     &     ( iabs(l1-l2).gt.j3 ).or.( l1+l2.lt.j3 ).or.
     &     ( iabs(j1-l2).gt.l3 ).or.( j1+l2.lt.l3 ).or.
     &     ( iabs(l1-j2).gt.l3 ).or.( l1+j2.lt.l3 ) ) then
c        write(*,*) '[ an6j ] error. triangle condition'
c        stop
      an6j=0.0d0
      return
      end if

**** argument arrangement ****
      if ( l1.eq.2 ) then
        ia=j1
        ia2=j1*2
        if ( iabs(j2-l3).le.iabs(j3-l2) ) then
          if ( l2.le.j3 ) then
            ib=j2
            ic=j3
            ib2=j2*2
            ibd2=l3*2
            ic2=j3*2
            icd2=l2*2
          else
            ib=l3
            ic=l2
            ib2=l3*2
            ibd2=j2*2
            ic2=l2*2
            icd2=j3*2
          end if
        else
          if ( l3.le.j2 ) then
            ib=j3
            ic=j2
            ib2=j3*2
            ibd2=l2*2
            ic2=j2*2
            icd2=l3*2
         else
            ib=l2
            ic=l3
            ib2=l2*2
            ibd2=j3*2
            ic2=l3*2
            icd2=j2*2
          end if
        end if
      else
        write(*,*) '[ an6j ] argument error.'
        stop
      end if

      is=ia+ib+ic
      if ( iand(is,1).eq.1 ) then
        sign=-1.0d0
      else
        sign=1.0d0
      end if

      if ( (icd2.eq.ic2-4).and.(ibd2.eq.ib2-4) ) then
c        write(*,*) 'type 1'
        w=dsqrt(
     &    dfloat( (is-2)*(is-1) )*dfloat( (is)*(is+1) )*
     &    dfloat( (is-ia2-3)*(is-ia2-2) )*
     &    dfloat( (is-ia2-1)*(is-ia2) ) /
     &    dfloat( (ib2-3)*(ib2-2) )/dfloat( (ib2-1)*(ib2)*(ib2+1) ) /
     &    dfloat( (ic2-3)*(ic2-2) )/dfloat( (ic2-1)*(ic2)*(ic2+1) ) )
      else if ( (icd2.eq.ic2-4).and.(ibd2.eq.ib2-2) ) then
c        write(*,*) 'type 2'
        w=dsqrt(
     &    dfloat( 4*(is-1)*(is) )*dfloat( (is+1)*(is-ia2-2) )*
     &    dfloat( (is-ia2-1)*(is-ia2) )*
     &    dfloat( (is-ib2)*(is-ic2+1) ) /
     &    dfloat( (ib2-2)*(ib2-1)*(ib2) )/dfloat( (ib2+1)*(ib2+2) ) /
     &    dfloat( (ic2-3)*(ic2-2)*(ic2-1) )/dfloat( (ic2)*(ic2+1) ) )
      else if ( (icd2.eq.ic2-4).and.(ibd2.eq.ib2) ) then
c        write(*,*) 'type 3'
        w=dsqrt(
     &    dfloat( 6*(is)*(is+1) )*dfloat( (is-ia2-1)*(is-ia2) )*
     &    dfloat( (is-ib2)*(is-ib2-1) )*
     &    dfloat( (is-ic2+1)*(is-ic2+2) ) /
     &    dfloat( (ib2-1)*(ib2)*(ib2+1) )/dfloat( (ib2+2)*(ib2+3) ) /
     &    dfloat( (ic2-3)*(ic2-2)*(ic2-1) )/dfloat( (ic2)*(ic2+1) ) )
      else if ( (icd2.eq.ic2-4).and.(ibd2.eq.ib2+2) ) then
c        write(*,*) 'type 4'
        w=dsqrt(
     &    dfloat( 4*(is+1)*(is-ia2) )*dfloat( (is-ib2-2)*(is-ib2-1) )*
     &    dfloat( (is-ib2)*(is-ic2+1) )*
     &    dfloat( (is-ic2+2)*(is-ic2+3) ) /
     &    dfloat( (ib2)*(ib2+1)*(ib2+2) )/dfloat( (ib2+3)*(ib2+4) ) /
     &    dfloat( (ic2-3)*(ic2-2)*(ic2-1) )/dfloat( (ic2)*(ic2+1) ) )
      else if ( (icd2.eq.ic2-4).and.(ibd2.eq.ib2+4) ) then
c        write(*,*) 'type 5'
        w=dsqrt(
     &    dfloat( (is-ib2-3)*(is-ib2-2) )*dfloat( (is-ib2-1)*(is-ib2) )*
     &    dfloat( (is-ic2+1)*(is-ic2+2) )*
     &    dfloat( (is-ic2+3)*(is-ic2+4) ) /
     &    dfloat( (ib2+1)*(ib2+2)*(ib2+3) )/dfloat( (ib2+4)*(ib2+5) ) /
     &    dfloat( (ic2-3)*(ic2-2)*(ic2-1) )/dfloat( (ic2)*(ic2+1) ) )
      else if ( (icd2.eq.ic2-2).and.(ibd2.eq.ib2-2) ) then
c        write(*,*) 'type 6'
        w=dfloat( 4*( (ia+ib)*(ia-ib+1)-(ic-1)*(ic-ib+1) ) )*
     &    dsqrt(
     &    dfloat( is*(is+1) )*dfloat( (is-ia2-1)*(is-ia2) ) /
     &    dfloat( (ib2-2)*(ib2-1) )/dfloat( (ib2)*(ib2+1)*(ib2+2) ) /
     &    dfloat( (ic2-2)*(ic2-1) )/dfloat( (ic2)*(ic2+1)*(ic2+2) ) )
      else if ( (icd2.eq.ic2-2).and.(ibd2.eq.ib2) ) then
c        write(*,*) 'type 7'
        w=dfloat( 2*( (ia+ib+1)*(ia-ib)-ic*ic+1 ) )*
     &    dsqrt(
     &    dfloat( 6*(is+1)*(is-ia2) )*dfloat( (is-ib2)*(is-ic2+1) ) /
     &    dfloat( (ib2-1)*(ib2)*(ib2+1) )/dfloat( (ib2+2)*(ib2+3) ) /
     &    dfloat( (ic2-2)*(ic2-1)*(ic2) )/dfloat( (ic2+1)*(ic2+2) ) )
      else if ( (icd2.eq.ic2-2).and.(ibd2.eq.ib2+2) ) then
c        write(*,*) 'type 8'
        w=dfloat( 4*( (ia+ib+2)*(ia-ib-1)-(ic-1)*(ib+ic+2) ) )*
     &    dsqrt(
     &    dfloat( (is-ib2-1)*(is-ib2) )*
     &    dfloat( (is-ic2+1)*(is-ic2+2) ) /
     &    dfloat( (ib2)*(ib2+1)*(ib2+2) )/dfloat( (ib2+3)*(ib2+4) ) /
     &    dfloat( (ic2-2)*(ic2-1)*(ic2) )/dfloat( (ic2+1)*(ic2+2) ) )
      else if ( (icd2.eq.ic2).and.(ibd2.eq.ib2) ) then
c        write(*,*) 'type 9'
        w=dfloat( 2*( 3*( ib*(ib+1)+ic*(ic+1)-ia*(ia+1) )*
     &                  ( ib*(ib+1)+ic*(ic+1)-ia*(ia+1)-1 )
     &                  -4*ib*(ib+1)*ic*(ic+1) ) ) /
     &    dsqrt(
     &    dfloat( (ib2-1)*(ib2)*(ib2+1) )*dfloat( (ib2+2)*(ib2+3) )*
     &    dfloat( (ic2-1)*(ic2)*(ic2+1) )*dfloat( (ic2+2)*(ic2+3) ) )
      else
c        write(*,*) '[ an6j ] error. argument'
        stop
      end if

      an6j=sign*w
      return
      end




      function wcoeff(ia,ib,ic,id,ie,if)
****
****    racah oceff. W(ia,ib,ic,id ; ie,if)
****
      implicit real*8 (a-h,o-z)

      wcoeff=d6ji(ia,ib,ie,id,ic,if)
      if (iand((ia+ib+ic+id)/2,1).eq.1) wcoeff=-wcoeff

      return
      end




      function coef9(j1,j2,j3,j4,j5,j6,j7,j8,j9)
****
****  to calc 9j           1       2       3
****                       4       5       6
****                       7       8       9
****
      implicit real*8 (a-h,o-z)

      coef9=0.0d0
      kmin=max0(iabs(j1-j9),iabs(j4-j8),iabs(j2-j6))
      kmax=min0(j1+j9,j4+j8,j2+j6)
      if (kmin.gt.kmax) return
      do 100 k=kmin,kmax,2
      coef9=coef9+(k+1)
     &  *wcoeff(j1,j9,j4,j8,k,j7)
     &  *wcoeff(j2,j6,j8,j4,k,j5)
     &  *wcoeff(j1,j9,j2,j6,k,j3)
100   continue

      return
      end


      function dfun(nj2t,nm2t,nk2t,beta)
c      function dfun(nj2,nm2,nk2,beta)
c      function dfun(a,b,c,cost)
c
c     d_function d^a_b,c(cos(t))
c      dfun=<a,b| exp(-i*t*Jy) |a,c>
c      a,b,c : real*8
c     okabe san no hon p.353
c
      implicit real*8 (a-h,o-z)
      parameter(lmax=100,ld=lmax*2+4)
      dimension cm(-1:1),ck(-1:1),wk(-ld:ld,0:1),df(-2:2,-2:2)
      data eps /1.d-14/
      dfun=0.d0
      cost = dcos(beta)
c      if (a.lt.0.d0) return
      if (abs(cost).gt.1.d0+eps) return
      nj2 = nj2t
      nm2 = nm2t
      nk2 = nk2t
c      nj2=nint(a+a)
c      nm2=nint(b+b)
c      nk2=nint(c+c)
      if (abs(nm2).gt.nj2) return
      if (abs(nk2).gt.nj2) return
      if (mod(nj2+nm2,2).eq.1) return
      if (mod(nj2+nk2,2).eq.1) return
      if (nj2.eq.0) then
         dfun=1.d0
         return
      endif
      ch=sqrt((1.d0+cost)/2.d0)
      sh=sqrt((1.d0-cost)/2.d0)
      df(-1,-1)= ch
      df(-1, 1)= sh
      df( 1,-1)=-sh
      df( 1, 1)= ch
      if (nj2.eq.1) then
         dfun=df(nm2,nk2)
         return
      endif
      z=cost
      if (nm2.ge.0) then
        if (nk2.ge.0) then
           f=1.d0
        else
           nk2=-nk2
           z=-z
           f=(-1)**((nj2+nm2)/2)
        endif
      else
        if (nk2.gt.0) then
           nm2=-nm2
           z=-z
           f=(-1)**((nj2-nk2)/2)
        else
           nm2=-nm2
           nk2=-nk2
           f=(-1)**((nm2-nk2)/2)
        endif
      endif
      if (nk2.ge.nm2 ) then
         m2=nm2
         k2=nk2
         d=f
      else
         m2=nk2
         k2=nm2
         d=f*(-1)**((nm2-nk2)/2)
      endif
      ch=sqrt((1.d0+z)/2.d0)
      sh=sqrt((1.d0-z)/2.d0)
      ch2=ch*ch
      sh2=sh*sh
      cs2=sqrt(2.d0)*ch*sh
      df(0,-2)=-cs2
      df(0, 0)=-sh2+ch2
      df(0, 2)= cs2
      df(2,-2)= sh2
      df(2, 0)=-cs2
      df(2, 2)= ch2
      ji=mod(nj2,2)+2
      do 4 k=-3,3
 4    wk(k,0)=0.d0
      if (ji.eq.2) then
         wk( 0,0)=1.d0
      else
         wk(-1,0)=-sh
         wk( 1,0)= ch
      endif
      kc=0
      kd=1
      mm=ji-2
      do 1 j2=ji,nj2,2
         ml=mm
         mm=min(j2,m2)
         m1=(mm-ml)/2
         xj=j2/2.d0-1.d0
         xm=mm/2.d0
         s2=sqrt((2*xj+1)*(2*xj+2))
         cm(0)=sqrt((xj-xm+1)*(xj+xm+1)*2)/s2
         cm(1)=sqrt((xj+xm)*(xj+xm+1))/s2
         ki=max(-j2,k2-nj2+j2)
         kf=min(j2,k2+nj2-j2)
         do 2 kk=ki,kf,2
            xk=kk/2.d0
            ck(-1)=sqrt((xj-xk)*(xj-xk+1))/s2
            ck( 0)=sqrt((xj-xk+1)*(xj+xk+1)*2)/s2
            ck( 1)=sqrt((xj+xk)*(xj+xk+1))/s2
            s=0.d0
            do 3 k1=-1,1
 3          s=s+ck(k1)*df(m1*2,k1*2)*wk(kk-k1*2,kc)
            wk(kk,kd)=s/cm(m1)
 2       continue
         wk(ki-2,kd)=0.d0
         wk(ki-4,kd)=0.d0
         wk(kf+2,kd)=0.d0
         wk(kf+4,kd)=0.d0
         ke=kd
         kd=kc
         kc=ke
 1    continue
      dfun=wk(k2,kc)*d
      return
      end



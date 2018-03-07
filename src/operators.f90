module operators
  use cross_coupled
  use isospin_operators
  implicit none 
  
  interface initialize_transition_operator
     module procedure initialize_EM_operator,initialize_beta_operator 
  end interface initialize_transition_operator

contains
!==================================================================  
!==================================================================
subroutine initialize_scalar_operator(name,Op,jbas)
  implicit none

  type(sq_op) :: Op
  type(spd) :: jbas
  character(10) :: name

  ! you can add new scalar operators here. 
  if (trim(name)=="rho21") then     
     call  initialize_rho21_zerorange(Op,jbas)
  else
     write(*,*) "Operator named "//trim(name)//" not found. FAIL."
     STOP
  end if
  
end subroutine initialize_scalar_operator
!==================================================================  
!==================================================================
subroutine write_scalar_operator(name,Op,jbas)
  implicit none

  type(sq_op) :: Op
  type(spd) :: jbas
  character(10) :: name

  ! you can add new scalar operators here. 
  if (trim(name)=="rho21") then
     print*, 'SRC DENSITY:', Op%E0
     
     open(unit=81,file=trim(OUTPUT_DIR)//trim(prefix)//"_SRC_density.dat")
     write(81,*) nint(Op%hospace), Op%eMax, Op%E0
     close(81)
     open(unit=81,file=trim(OUTPUT_DIR)//trim(prefix)//"_energy.dat")
     write(81,*) nint(Op%hospace), Op%eMax, Op%E0
     close(81)
  else
     write(*,*) "Operator named "//trim(name)//" not found. FAIL."
     STOP "How did this not get caught in initialize_scalar_operator?"
  end if
  
end subroutine write_scalar_operator
!==================================================================  
!==================================================================
subroutine calculate_pipj(pp,jbas) 
  implicit none
  
  type(sq_op) :: pp
  type(spd) :: jbas
  integer :: ist,J,Tz,Par,a,b,c,d,q,qx,N,j_min,x,II,JJ,mass
  real(8) :: V,Vcm,g1,g2,g3,pre,hw,V1
  integer :: C1,C2,int1,int2,i1,i2,htype,COM,n1,n2,g_ix  
  
  N = jbas%total_orbits
  
  mass = 0.d0 
  do a = 1, N
     mass = mass + jbas%con(a) *(jbas%jj(a)+1) 
  end do

  call calculate_h0_harm_osc(jbas,pp,4)

  do q = 1, pp%nblocks
     
     J = pp%mat(q)%lam(1) 
  
     do g_ix = 1,6 
   
        ! figure out how big the array is
        n1 = size(pp%mat(q)%gam(g_ix)%X(:,1))
        n2 = size(pp%mat(q)%gam(g_ix)%X(1,:))
        if ((n1*n2) == 0) cycle 
        
        ! read in information about which 
        ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        
        do II = 1, n1
           a = pp%mat(q)%qn(c1)%Y(II,1)
           b = pp%mat(q)%qn(c1)%Y(II,2)
           do JJ = 1, n2 
              c = pp%mat(q)%qn(c2)%Y(JJ,1)
              d = pp%mat(q)%qn(c2)%Y(JJ,2)        

              pp%mat(q)%gam(g_ix)%X(II,JJ) = p1_p2( a, b, c, d, J ,jbas )* &
                   pp%hospace/mass

           end do 
        end do 

     end do 
  end do   

      
end subroutine calculate_pipj
!==================================================================  
!==================================================================
subroutine calculate_rirj(rr,jbas) 
  implicit none
  
  type(sq_op) :: rr
  type(spd) :: jbas
  integer :: ist,J,Tz,Par,a,b,c,d,q,qx,N,j_min,x,II,JJ,mass
  real(8) :: V,Vcm,g1,g2,g3,pre,hw,V1
  integer :: C1,C2,int1,int2,i1,i2,htype,COM,n1,n2,g_ix
  logical :: rr_calc,pp_calc
  
  
  N = jbas%total_orbits
  
  mass = 0.d0 
  do a = 1, N
     mass = mass + jbas%con(a) *(jbas%jj(a)+1) 
  end do
  ! check if we are concerned with other operators

  call calculate_h0_harm_osc(jbas,rr,5)
     
  do q = 1, rr%nblocks
     
     J = rr%mat(q)%lam(1) 
  
     do g_ix = 1,6 
   
        ! figure out how big the array is
        n1 = size(rr%mat(q)%gam(g_ix)%X(:,1))
        n2 = size(rr%mat(q)%gam(g_ix)%X(1,:))
        if ((n1*n2) == 0) cycle 
        
        ! read in information about which 
        ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        
        do II = 1, n1
           a = rr%mat(q)%qn(c1)%Y(II,1)
           b = rr%mat(q)%qn(c1)%Y(II,2)
           do JJ = 1, n2 
              c = rr%mat(q)%qn(c2)%Y(JJ,1)
              d = rr%mat(q)%qn(c2)%Y(JJ,2)        

              rr%mat(q)%gam(g_ix)%X(II,JJ) = r1_r2( a, b, c, d, J ,jbas )* &
                   rr%hospace/mass

           end do 
        end do 

     end do 
  end do   

end subroutine calculate_rirj
!=====================================================================================
!=====================================================================================
subroutine initialize_EM_operator(trs_type,rank,Op,zr,jbas,tcalc)
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: Op,zr
  character(1) :: trs_type,ranklab
  integer :: rank 
  logical :: tcalc
 
  Select case (trs_type) 
     case ('E') ! ELECTRIC TRANSITION
         Op%rank = 2*rank
         Op%herm = 1
         Op%dpar = abs((-1)**rank-1) 
         write(ranklab,'(I1)') rank 
         Op%trans_label = trs_type//ranklab         
         if (rank==0) then 
            print*, 'not implemented'
         else            
            if (allocated(phase_pp)) then
               deallocate(phase_hh,phase_pp)
            end if
            call allocate_tensor(jbas,Op,zr)
            Op%hospace=zr%hospace
            call calculate_EX(Op,jbas)
         end if 
     case ('M') ! MAGNETIC TRANSITION
        Op%rank = 2*rank
        Op%herm = 1
        Op%dpar = abs((-1)**(rank-1)-1) 
        write(ranklab,'(I1)') rank 
        Op%trans_label = trs_type//ranklab
        if (rank==0) then 
           print*, 'not implemented'
        else
            if (allocated(phase_pp)) then
               deallocate(phase_hh,phase_pp)
            end if
           call allocate_tensor(jbas,Op,zr)
           Op%hospace=zr%hospace
           call calculate_MX(Op,jbas)
        end if
     case ('F') ! FERMI TRANSITION
        stop 'not implemented'
     case ('G') ! GAMOW-TELLER TRANSITION
        stop 'Beta-decay operators need to use type(iso_operator) structure'        
     case default
        tcalc = .false. 
     end select
end subroutine initialize_EM_operator
!=====================================================================================
!=====================================================================================
subroutine initialize_beta_operator(trs_type,dTz,Op,zr,jbas,tcalc)
  implicit none
  
  type(spd) :: jbas
  type(iso_operator) :: Op
  type(sq_op) :: zr
  character(1) :: trs_type,ranklab
  character(1):: dTz 
  logical :: tcalc
  
  Select case (trs_type) 
     case ('E') ! ELECTRIC TRANSITION
        STOP 'Electromagnetic operators need to use type(sq_op) structure'        
     case ('M') ! MAGNETIC TRANSITION
        STOP 'Electromagnetic operators need to use type(sq_op) structure'
     case ('F') ! FERMI TRANSITION
        stop 'not implemented'
     case ('G') ! GAMOW-TELLER TRANSITION
        OP%rank = 2
        Op%herm = 1
        Op%dpar = 0        
        Op%trans_label = 'GT'
        if (dTZ == '+') then
           Op%dTz = 1
        else if (dTZ == '-') then
           Op%dTz = -1
        else
           stop 'not implemented'
        end if
        
        call allocate_isospin_operator(jbas,Op,zr)
        Op%hospace=zr%hospace
        call calculate_GT(Op,jbas)               
     case default
        tcalc = .false. 
  end select
end subroutine initialize_beta_operator
!===================================================================
!===================================================================
subroutine initialize_TDA(TDA,jbas,Jtarget,PARtarget,cut) 
  ! figure out how big the TDA matrix has to be
  ! allocate it
  ! uses full_sp_block_mat, just because it's got the right shape
  implicit none 
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: TDA
  integer :: JT,PI,Jmax,q,i,a,r,nh,np,ix,ax
  integer :: Jtarget, PARtarget,cut
 
  Jmax = jbas%Jtotal_max
  nh = sum(jbas%con) 
  np = jbas%total_orbits - nh 
  TDA%blocks = 1   
  
  if (allocated(TDA%blkM)) then 
     deallocate(TDA%blkM,TDA%map) 
  end if 
  
  allocate(TDA%blkM(1))
  allocate(TDA%map(1)) 
  
  q = 1
  
  !do PI = 0,1
   !  do JT = 0,2*Jmax,2
    
  PI = PARtarget 
  JT = Jtarget 
      
        r = 0
        ! count the states
        do ix = 1,nh
           do ax = 1,np 
              i =jbas%holes(ix) 
              a =jbas%parts(ax)
              if (a > cut) cycle 
              if (i < 10) cycle
              if (jbas%itzp(i) .ne. jbas%itzp(a) ) cycle ! cannot change from p to n
            
              if (.not. (triangle(jbas%jj(i),jbas%jj(a),JT))) cycle 
              if (mod(jbas%ll(i)+jbas%ll(a),2) .ne. PI )  cycle ! technically l_a - l_i  
              
              r = r + 1 
           end do 
        end do 
              
        ! allocate
        TDA%map(q) = r
        allocate(TDA%blkM(q)%matrix(r,r)) 
        allocate(TDA%blkM(q)%extra(10*r))
        allocate(TDA%blkM(q)%eigval(r)) 
        allocate(TDA%blkM(q)%labels(r,2)) 
        allocate(TDA%blkM(q)%states(r))
        TDA%blkM(q)%states = 0.d0
        
        ! fill the states
        r = 0
        do ix = 1,nh
           do ax = 1,np 
              i =jbas%holes(ix) 
              a =jbas%parts(ax)              
              if (a > cut) cycle 
              if (i < 10) cycle
              if (jbas%itzp(i) .ne. jbas%itzp(a) ) cycle ! cannot change from p to n
            
              if (.not. (triangle(jbas%jj(i),jbas%jj(a),JT))) cycle 
              if (mod(jbas%ll(i)+jbas%ll(a),2) .ne. PI )  cycle ! technically l_a - l_i 
              r = r + 1 
              TDA%blkM(q)%labels(r,1) = i 
              TDA%blkM(q)%labels(r,2) = a 
           end do 
        end do 
        TDA%blkM(q)%lmda(1) = JT
        TDA%blkM(q)%lmda(2) = PI 
        TDA%blkM(q)%lmda(3) = 0 
 !       q = q + 1
 !    end do 
 ! end do 

       print*, TDA%blkM(1)%labels(:,1)
       print*, TDA%blkM(1)%labels(:,2)
end subroutine 
!==========================================
!==========================================
subroutine calc_TDA(TDA,HS,HSCC,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: TDA
  type(sq_op) :: HS 
  type(cc_mat) :: HSCC
  integer :: q,JT,r1,r2,a,b,i,j,x,g,NBindx,Rindx,q1,Tz,PAR,JTM
  integer :: ji,jj,ja,jb,phase
  
  JTM = jbas%Jtotal_max
  do q = 1, TDA%blocks
     JT = TDA%blkM(q)%lmda(1) 
     Tz = 0
     PAR = TDA%blkM(q)%lmda(2)
  
     q1 = JT/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
     
     do r1 = 1, TDA%map(q) 
        i = TDA%blkM(q)%labels(r1,1)
        a = TDA%blkM(q)%labels(r1,2)
        
        ji = jbas%jj(i)
        ja = jbas%jj(a) 
        ! get CC index for this pair
        x = CCindex(a,i,HS%Nsp)
        g = 1
        do while (HSCC%qmap(x)%Z(g) .ne. q1 )
           g = g + 1
        end do

        NBindx = HSCC%nbmap(x)%Z(g) 
           
        do r2 = r1, TDA%map(q)
           j = TDA%blkM(q)%labels(r2,1)
           b = TDA%blkM(q)%labels(r2,2)
           
           jj = jbas%jj(j)
           jb = jbas%jj(b)
           phase = (-1)**((jj+jb)/2)
           ! get CCindex for this pair
           x = CCindex(b,j,HS%Nsp) 
           g = 1
           do while (HSCC%qmap(x)%Z(g) .ne. q1 )
              g = g + 1
           end do
              
           Rindx = HSCC%rmap(x)%Z(g)  
           
           ! one body piece
           TDA%blkM(q)%matrix(r1,r2) = f_elem(a,b,HS,jbas) * &
                kron_del(i,j) - f_elem(j,i,HS,jbas) * kron_del(a,b)
           
           ! two body piece
!           TDA%blkM(q)%matrix(r1,r2) = TDA%blkM(q)%matrix(r1,r2) - &
 !               Vcc(i,a,b,j,JT,HS,jbas)*(-1)**((ja+ji)/2) * &
  !              (-1) ** ((JT)/2) /sqrt(JT+1.d0) 
           
           TDA%blkM(q)%matrix(r1,r2) = TDA%blkM(q)%matrix(r1,r2) - &
                HSCC%CCX(q1)%X(Rindx,NBindx) * HSCC%herm * & !need scaling.
                (-1) ** (JT/2) /sqrt(JT+1.d0) * phase  
         
           ! I have to account for the fact that the CCME are scaled by 
           ! Sqrt[2J+1] * (-1)^(j_i + j_a) 
        
           ! hermiticity 
           TDA%blkM(q)%matrix(r2,r1) = TDA%blkM(q)%matrix(r1,r2) 
        end do 
     end do 
  end do 
        
end subroutine 
!===================================================
!===================================================
subroutine TDA_expectation_value(TDA_HAM,TDA_OP) 
  ! takes expectation values of TDA_OP for the eigenvectors defined in TDA_HAM
  implicit none 
  
  type(full_sp_block_mat) :: TDA_HAM, TDA_OP,AUX  
  integer :: q,dm,i 
  
  call duplicate_sp_mat(TDA_HAM,AUX) 
 
  do q = 1, TDA_OP%blocks
        
     dm = TDA_OP%map(q)
     if (dm == 0)  cycle
     
     ! transform operator into basis defined by TDA vectors
     call dgemm('N','N',dm,dm,dm,al,TDA_OP%blkM(q)%matrix&
          ,dm,TDA_HAM%blkM(q)%matrix,dm,bet,AUX%blkM(q)%matrix,dm) 
     call dgemm('T','N',dm,dm,dm,al,TDA_HAM%blkM(q)%matrix&
          ,dm,AUX%blkM(q)%matrix,dm,bet,TDA_OP%blkM(q)%matrix,dm)
     
     ! diagonal elements are the expectation values
     do i = 1, dm 
        TDA_OP%blkM(q)%eigval(i) = TDA_OP%blkM(q)%matrix(i,i) 
     end do    
     
  end do

end subroutine 
!=========================================================
subroutine initialize_rms_radius(rms,rr,jbas)
  implicit none 
  
  type(sq_op) :: rr, rms 
  type(spd) :: jbas
  real(8) :: mass_factor 
  integer :: q,i
  
  mass_factor = 1.d0-1.d0/dfloat(rr%Aprot + rr%Aneut) 

  jbas%hw = 1
  call calculate_h0_harm_osc(jbas,rms,5)
  jbas%hw = rr%hospace
  
  ! multiply by scale factors to make it into r^2 instead of u_ho 
  rms%fhh = rms%fhh * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  rms%fpp = rms%fpp * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  rms%fph = rms%fph * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  
  
  do q = 1, rms%nblocks
     do i = 1,6
        
        rms%mat(q)%gam(i)%X = -2.d0 * hbarc2_over_mc2 &
             / rr%hospace**2 / dfloat(rr%Aneut + rr%Aprot) * &
             rr%mat(q)%gam(i)%X 
        
     end do 
  end do 


end subroutine 
!=========================================================
subroutine initialize_CM_radius(rms,rr,jbas)
  implicit none 
  
  type(sq_op) :: rr, rms 
  type(spd) :: jbas
  real(8) :: mass_factor 
  integer :: q,i
  
  mass_factor = 1.d0/dfloat(rr%Aprot + rr%Aneut) 

  jbas%hw=1
  call calculate_h0_harm_osc(jbas,rms,5)
  jbas%hw=rr%hospace
  
  ! multiply by scale factors to make it into r^2 instead of u_ho 
  rms%fhh = rms%fhh * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  rms%fpp = rms%fpp * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  rms%fph = rms%fph * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  
  
  do q = 1, rms%nblocks
     do i = 1,6
        
        rms%mat(q)%gam(i)%X = 2.d0 * hbarc2_over_mc2 &
             / rr%hospace**2 * mass_factor * &
             rr%mat(q)%gam(i)%X 
        
     end do 
  end do 


end subroutine
!====================================================================
!=========================================================
subroutine initialize_CM_radius_onebody(rms,rr,jbas)
  implicit none 
  
  type(sq_op) :: rr, rms 
  type(spd) :: jbas
  real(8) :: mass_factor, elem
  integer :: q,a,b,ak,bk,na,nb,la,lb
  
  mass_factor = 1.d0/dfloat(rr%Aprot + rr%Aneut) 
  
  
  ! multiply by scale factors to make it into r^2 instead of u_ho 
  
  do a = 1, rr%belowEF
     
     ak = jbas%holes(a)
     na = jbas%nn(ak)
     la = jbas%ll(ak) 
     
     do b = a, rr%belowEF
     
        bk = jbas%holes(b)
        nb = jbas%nn(bk)
        lb = jbas%ll(bk)

  
        if (la == lb - 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 0.5d0)
           
           else if (na == nb + 1)  then 
              
              elem = -1.d0*sqrt(nb + 1.d0)

           else 
              elem = 0.d0 
              
           end if 
        
        else if (la == lb + 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 1.5d0)
           
           else if (na == nb - 1) then 
              
              elem = -1.d0*sqrt(float(nb)) 
           
           else 
              
              elem = 0.d0 
           end if 
        else 
           elem = 0.d0 
        end if 
  
        rms%fhh(a,b) = sqrt(hbarc2_over_mc2 / rms%hospace) * elem * mass_factor
        rms%fhh(b,a) = rms%fhh(a,b) 
     end do  
  end do 

!fph  
  do a = 1, rr%Nsp - rr%belowEF
     
     ak = jbas%parts(a)
     na = jbas%nn(ak)
     la = jbas%ll(ak) 
     
     do b = 1, rr%belowEF
     
        bk = jbas%holes(b)
        nb = jbas%nn(bk)
        lb = jbas%ll(bk)

  
        if (la == lb - 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 0.5d0)
           
           else if (na == nb + 1)  then 
              
              elem = -1.d0*sqrt(nb + 1.d0)

           else 
              elem = 0.d0 
              
           end if 
        
        else if (la == lb + 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 1.5d0)
           
           else if (na == nb - 1) then 
              
              elem = -1.d0*sqrt(float(nb)) 
           
           else 
              
              elem = 0.d0
           end if 
        else 
           elem = 0.d0 
        end if 
  
        rms%fph(a,b) = sqrt(hbarc2_over_mc2 / rms%hospace) * elem * mass_factor
       
     end do  
  end do 
 
!fpp
   do a = 1,rr%Nsp-rr%belowEF
     
     ak = jbas%parts(a)
     na = jbas%nn(ak)
     la = jbas%ll(ak) 
     
     do b = a, rr%Nsp - rr%belowEF
     
        bk = jbas%parts(b)
        nb = jbas%nn(bk)
        lb = jbas%ll(bk)

  
        if (la == lb - 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 0.5d0)
           
           else if (na == nb + 1)  then 
              
              elem = -1.d0*sqrt(nb + 1.d0)

           else 
              elem = 0.d0 
              
           end if 
        
        else if (la == lb + 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 1.5d0)
           
           else if (na == nb - 1) then 
              
              elem = -1.d0*sqrt(float(nb)) 
           
           else 
              
              elem = 0.d0 
           end if
        else 
           elem = 0.d0 
        end if 
  
        rms%fpp(a,b) = sqrt(hbarc2_over_mc2 / rms%hospace) * elem * mass_factor
        rms%fpp(b,a) = rms%fpp(a,b) 
     end do  
  end do 
 
  

end subroutine initialize_CM_radius_onebody
!==================================================================== 
subroutine build_Hcm(pp,rr,Hcm,jbas)
  implicit none 

  type(spd) :: jbas
  type(sq_op) :: pp,rr,Htemp,Hcm
  real(8) :: hw_tilde,hw
  integer :: i

  hw = pp%hospace
  hw_tilde = pp%com_hw 
  if (allocated(phase_pp)) then
     deallocate(phase_hh,phase_pp)
  end if
  call allocate_tensor(jbas,Hcm,pp)
  Hcm%hospace=pp%hospace
            
  call duplicate_sq_op(pp,Htemp)
  call add_sq_op(pp,1.d0,rr,hw_tilde**2/hw**2,Htemp)  
  Hcm%E0 = Htemp%E0 - 1.5d0 * hw_tilde

  call copy_rank0_to_tensor_format(Htemp,Hcm,jbas) 
end subroutine build_Hcm
!==================================================================== 
subroutine calculate_CM_energy(pp,rr) 
  implicit none 
  
  type(sq_op) :: pp,rr,Hcm
  real(8) :: hw,wTs(2),Ecm(3) 
  integer :: i

  hw = pp%hospace
  call duplicate_sq_op(pp,Hcm)
  call add_sq_op(pp,1.d0,rr,1.d0,Hcm)
  Ecm(1) = Hcm%E0 - 1.5d0*hw ! store Ecm for this Hcm frequency 
    
  ! new frequencies
  wTs = optimum_omega_for_CM_hamiltonian(hw,Ecm(1)) 
  
  do i = 1, 2
     call add_sq_op(pp,1.d0,rr,wTs(i)**2/hw**2,Hcm)
     Ecm(i+1) = Hcm%E0 - 1.5d0 * wTs(i) 
  end do   

  !write results to file
  open(unit=42,file=trim(OUTPUT_DIR)//&
         trim(adjustl(prefix))//'_Ecm.dat')
  write(42,'(6(e14.6))') hw, wTs, Ecm 
  close(42)
  
end subroutine
!=====================================================
!==================================================================== 
subroutine calculate_CM_energy_TDA(TDA,rr,pp,ppTDA,rrTDA,hw) 
  implicit none 
  
  type(sq_op) :: rr,pp,Hcm 
  type(full_sp_block_mat) :: TDA,HcmTDA,ppTDA,rrTDA
  real(8) :: hw,wTs(2)
  real(8),allocatable,dimension(:) :: energies,omegas 
  integer :: i,q
  character(3) :: args
 
  
  call duplicate_sp_mat(TDA,HcmTDA) 
  call duplicate_sq_op(rr,Hcm)
  allocate(energies(3*ppTDA%map(1))) 
  allocate(omegas(2*TDA%map(1)))
  
  call add_sp_mat(ppTDA,1.d0,rrTDA,1.d0,HcmTDA)
  call add_sq_op(pp,1.d0,rr,1.d0,Hcm) 
  call TDA_expectation_value(TDA,HcmTDA) 

  energies(1:TDA%map(1)) = HcmTDA%blkM(1)%eigval + Hcm%E0 - 1.5d0*hw
  
  do q = 1, TDA%map(1) 
     ! new frequencies
     wTs = optimum_omega_for_CM_hamiltonian(hw,energies(q)) 
     omegas(q) = wTs(1) 
     omegas(q+TDA%map(1)) = wTs(2)
     
     do i = 1, 2
        call add_sq_op(pp,1.d0,rr,wTs(i)**2/hw**2,Hcm) 
        call add_sp_mat(ppTDA,1.d0,rrTDA,wTs(i)**2/hw**2,HcmTDA)
        call TDA_expectation_value(TDA,HcmTDA) 
        energies(i*TDA%map(1)+q) = HcmTDA%blkM(1)%eigval(q) + Hcm%E0 - 1.5d0*wTs(i) 
     end do
  end do 
  
  ! for formatting
   i = 1+5*TDA%map(1) 
   write(args,'(I3)') i 
   args = adjustl(args) 

  !write results to file
  open(unit=42,file=trim(OUTPUT_DIR)//&
         trim(adjustl(prefix))//'_Ecm_excited.dat')
  write(42,'('//trim(args)//'(e14.6))') hw, omegas, energies 
  close(42)
  
end subroutine calculate_CM_energy_TDA

subroutine calculate_EX(op,jbas) 
  ! calculates electromagnetic transition operator 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: op   
  integer :: a,b,ax,bx,rank,na,nb
  integer :: ta,tb,tc,td,la,lb,lc,ld,ja,jb,jc,jd
  real(8) :: pol,charge(2),dcgi,dcgi00,x,dcg,hw,holength

  pol = 0.0d0  ! between 0 and 1 
  ! if (op%rank == 2 )then  !!! electricdipole transition
  !    pol = -0.5 ! N=Z 
  ! end if
  
  hw = op%hospace
  holength = sqrt(hbarc2_over_mc2/hw) 
  x = dcgi00()
  ! right now the units of this operator are: e fm^(X)
  ! so the charge is not multiplied by the fundamental charge. 
  charge(1) = (1+pol) 
  charge(2) = pol 
  
  !!! Brown sd-Shell values (2008)
 ! charge(1) = 1.36
 ! charge(2) = 0.45  

  !!! Richard Trippel values
  ! charge(1) = 1  
  ! charge(2)=1

  ! access the correct charge with: 
  ! charge( (jbas%tz(i) + 1)/2 + 1 )  
  
  rank = op%rank
  
  do ax = 1, op%belowEF 
     do bx = 1, op%belowEF 
 
        a = jbas%holes(ax)
        b = jbas%holes(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
       
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 
       
        if ( mod(la+lb+rank/2,2) == 1) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
               
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 

        op%fhh(ax,bx) = charge((ta + 1)/2+1)/sqrt(4.d0*Pi_const)* &
             (-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2)*holength**(rank/2)  
     end do 
  end do 
        
  do ax = 1, op%nsp-op%belowEF 
     do bx = 1, op%nsp-op%belowEF 
  
        a = jbas%parts(ax)
        b = jbas%parts(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
        
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 
        
        if ( mod(la+lb+rank/2,2) == 1) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
        
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 
        
        op%fpp(ax,bx) = charge((ta + 1)/2+1)/sqrt(4.d0*Pi_const)* &
             (-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2)*holength**(rank/2)  
     end do 
  end do     

  do ax = 1, op%nsp-op%belowEF 
     do bx = 1, op%belowEF 
  
        a = jbas%parts(ax)
        b = jbas%holes(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
        
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 
        
        if ( mod(la+lb+rank/2,2) == 1) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
        
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 
        
        op%fph(ax,bx) = charge((ta + 1)/2+1)/sqrt(4.d0*Pi_const)* &
             (-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2)*holength**(rank/2)  
     end do
  end do     
        
end subroutine calculate_EX
!==================================================================================
!==================================================================================
subroutine calculate_MX(op,jbas) 
  ! calculates electromagnetic transition operator 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: op   
  integer :: a,b,ax,bx,rank,na,nb,kappa
  integer :: ta,tb,tc,td,la,lb,lc,ld,ja,jb,jc,jd
  real(8) :: pol,charge(2),dcgi,dcgi00,x,dcg,hw,holength
  real(8) :: mu_N_overC,gl(2),gs(2)
  

  pol = 0.0d0   ! between -1 and 0 
  hw = op%hospace
  holength = sqrt(hbarc2_over_mc2/hw) 
  mu_N_overc = 1.d0!sqrt(hbarc2_over_mc2/4.d0) 
  x = dcgi00()
  ! right now the units of this operator are: e fm^(X)
  ! so the charge is not multiplied by the fundamental charge. 
  charge(1) = (1+pol) 
  charge(2) = pol 
  gl(1) = 1
  gl(2) = 0
  gs(1) = 5.586
  gs(2) = -3.826
  ! access the correct charge with: 
  ! charge( (jbas%tz(i) + 1)/2 + 1 )  
  
  rank = op%rank
  
  do ax = 1, op%belowEF 
     do bx = 1, op%belowEF 
 
        a = jbas%holes(ax)
        b = jbas%holes(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
       
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 

        if ( mod(la+lb+rank/2,2) == 0) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
               
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 
        
        kappa = (-1) ** ((ja+1)/2+la)*(ja+1) &
             +(-1) ** ((jb+1)/2+lb)*(jb+1) 
        
        op%fhh(ax,bx) = mu_N_overC*charge(1) &
             /sqrt(4.d0*Pi_const)*(-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2-1)*holength**(rank/2-1)*&
             (rank-kappa)/2*(gl((ta+1)/2+1)*(1+kappa/(rank+2.d0))-gs((ta+1)/2+1)/2.d0) 
     end do 
  end do 
        
  do ax = 1, op%nsp-op%belowEF 
     do bx = 1, op%nsp-op%belowEF 
  
        a = jbas%parts(ax)
        b = jbas%parts(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
        
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 

        if ( mod(la+lb+rank/2,2) == 0) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
        
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 

        kappa = (-1) ** ((ja+1)/2+la)*(ja+1) &
             +(-1) ** ((jb+1)/2+lb)*(jb+1) 
        
        op%fpp(ax,bx) = mu_N_overC*charge(1) &
             /sqrt(4.d0*Pi_const)*(-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2-1)*holength**(rank/2-1)*&
             (rank-kappa)/2*(gl((ta+1)/2+1)*(1+kappa/(rank+2.d0))-gs((ta+1)/2+1)/2.d0) 
       
     end do 
  end do     

  do ax = 1, op%nsp-op%belowEF 
     do bx = 1, op%belowEF 
  
        a = jbas%parts(ax)
        b = jbas%holes(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
        
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 
!        if (la .ne. lb) cycle
        if ( mod(la+lb+rank/2,2) == 0) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
        
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 

        kappa = (-1) ** ((ja+1)/2+la)*(ja+1) &
             +(-1) ** ((jb+1)/2+lb)*(jb+1) 

        op%fph(ax,bx) = mu_N_overC*charge(1) &
             /sqrt(4.d0*Pi_const)*(-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2-1)*holength**(rank/2-1)*&
             (rank-kappa)/2*(gl((ta+1)/2+1)*(1+kappa/(rank+2.d0))-gs((ta+1)/2+1)/2.d0) 
        
     end do
  end do     
        
end subroutine calculate_MX
!==================================================================================
!==================================================================================
subroutine calculate_GT(op,jbas) 
  ! calculates Gamow-Teller transition operator 
  implicit none 
  
  type(spd) :: jbas
  type(iso_operator) :: op   
  integer :: a,b,ax,bx,rank,na,nb,kappa,N
  integer :: ta,tb,tc,td,la,lb,lc,ld,ja,jb,jc,jd
  real(8) :: pol,charge(2),dcgi,dcgi00,x,dcg,hw,holength
  real(8) :: mu_N_overC,gl(2),gs(2),d6ji
  

  pol = 0.0d0   ! between -1 and 0 
  hw = op%hospace
  holength = sqrt(hbarc2_over_mc2/hw) 
  mu_N_overc = 1.d0!sqrt(hbarc2_over_mc2/4.d0) 
  x = dcgi00()
  N = jbas%total_orbits
  ! access the correct charge with: 
  ! charge( (jbas%tz(i) + 1)/2 + 1 )  
  
  rank = op%rank
  
  do a = 1, N  
     do b = 1,N
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
       
        if (ta - 2*OP%dTz .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 

        if (la .ne. lb) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 

        if (na .ne. nb) cycle
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
               
        if (.not.(triangle(ja,jb,rank))) cycle

        op%fock(a,b) = (-1) ** ((ja+3)/2 + la ) * sqrt(2.d0) * d6ji(1,1,2,jb,ja,la*2)*&
             sqrt((ja+1.d0)*(jb+1.d0))

        op%fock(b,a) = op%fock(a,b) ** ((ja-jb)/2)
     end do
  end do 
             
        
end subroutine calculate_GT
!====================================================================
!====================================================================
real(8) function iso_transition_to_ground_ME( Trans_op , Qdag,jbas ) 
  implicit none
  
  type(spd) :: jbas
  type(iso_operator) :: Trans_op
  type(iso_ladder) :: Qdag 
  integer :: ja,jb,ji,jj,rank,Nsp,Abody
  integer :: a,b,i,j,ax,ix,jx,bx,J1,J2
  real(8) :: sm , phase
  
  sm = 0.d0 
  
  Nsp = jbas%total_orbits 
  Abody = sum(jbas%con) 
   
  rank = Trans_op%rank 
  
  if (rank .ne. Qdag%rank)  then 
     iso_transition_to_ground_ME = 0.d0 
     return
  end if 

  
  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do ix = 1,Abody 
        i = jbas%holes(ix)
        ji = jbas%jj(i) 
        
        
        phase = (-1) ** ((ja-ji)/2)

        sm = sm + f_iso_op_elem(i,a,Trans_op,jbas)*&
             f_iso_ladder_elem(a,i,Qdag,jbas)*phase 
     end do
  end do 

  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do bx = 1,Nsp-Abody
        b = jbas%parts(bx)
        jb = jbas%jj(b) 
        
        do ix = 1,Abody 
           i = jbas%holes(ix)
           ji = jbas%jj(i) 
           
           do jx = 1,Abody 
              j = jbas%holes(jx)
              jj = jbas%jj(j) 
  
              do J1 = abs(ji-jj),ji+jj,2
                 do J2 = abs(ja-jb),ja+jb,2
                    
                    phase = (-1) **((J1+J2)/2) 
                    sm = sm + phase*0.25d0*&
                         iso_op_elem(i,j,a,b,J1,J2,Trans_op,jbas)*&
                         iso_ladder_elem(a,b,i,j,J2,J1,Qdag,jbas) 
                    
                 end do 
              end do 
           end do
        end do
     end do
  end do

  iso_transition_to_ground_ME = sm * (-1.d0)**(rank/2)
  !  BY SUHONEN'S DEFINITION, I SHOULD BE DEVIDING BY Sqrt(2J+1) 
  !  BUT THE LANCZOS ALGORITHM DID THAT FOR US ALREADY. 
  
end function iso_transition_to_ground_ME
!==============================================================
!==============================================================
real(8) function transition_ME_Tz_EM_Tz( Xout,Trans_op ,Xin,jbas ) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: Trans_op
  type(iso_ladder) :: Xout,Xin,product
  integer :: ja,jb,ji,jj,rank_out,rank_in,rank_op,Nsp,Abody,M
  integer :: a,b,i,j,ax,ix,jx,bx,J1,J2,dpar_in,dpar_out,dpar_op
  real(8) :: sm , phase,dcgi,mult,aaa,bbb

  sm = 0.d0
  Nsp = jbas%total_orbits 
  Abody = sum(jbas%con) 
   
  rank_op = Trans_op%rank 
  rank_in = Xin%rank 
  rank_out = Xout%rank

  if (.not. triangle(rank_in,rank_op,rank_out))  then 
     transition_ME_Tz_EM_Tz = 0.d0 
     return
  end if 

  if (Xin%dTz .ne. Xout%dTz) then
     transition_ME_Tz_EM_Tz = 0.d0
     return
  end if
  
  dpar_in = Xin%dpar/2
  dpar_out = Xout%dpar/2
  dpar_op = Trans_op%dpar/2
  
  if (mod(dpar_in+dpar_op+dpar_out,2).ne.0) then
     transition_ME_Tz_EM_Tz = 0.d0
  end if

  call duplicate_isospin_ladder(Xout,product)
  call tensor_dTz_tensor_product(Trans_op,Xin,product,jbas)
  
  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do ix = 1,Abody 
        i = jbas%holes(ix)
        ji = jbas%jj(i) 
        
        sm = sm + f_iso_ladder_elem(a,i,Xout,jbas)*&
             f_iso_ladder_elem(a,i,product,jbas)
     end do 
  end do 

  
  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do bx = 1,Nsp-Abody
        b = jbas%parts(bx)
        jb = jbas%jj(b) 
        
        do ix = 1,Abody 
           i = jbas%holes(ix)
           ji = jbas%jj(i) 
           
           do jx = 1,Abody 
              j = jbas%holes(jx)
              jj = jbas%jj(j) 
  
              do J1 = abs(ji-jj),ji+jj,2
                 do J2 = abs(ja-jb),ja+jb,2

                    if (.not. triangle(J1,J2,rank_out)) cycle

                    sm = sm + 0.25d0*&
                        iso_ladder_elem(a,b,i,j,J2,J1,product,jbas)*&
                        iso_ladder_elem(a,b,i,j,J2,J1,Xout,jbas) 
                    
                 end do
              end do 
           end do
        end do
     end do
  end do
  
  transition_ME_Tz_EM_Tz = sm * sqrt(rank_in+1.d0) &
       * ( -1 ) ** ((rank_in+rank_out+rank_op)/2)  
  !  BY SUHONEN'S DEFINITION, I SHOULD BE DEVIDING BY Sqrt(2J+1) 
  !  BUT THE LANCZOS ALGORITHM DID THAT FOR US ALREADY. 
 
end function transition_ME_Tz_EM_Tz
!====================================================================
!====================================================================
real(8) function transition_ME( Xout,Trans_op ,Xin,jbas ) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: Trans_op,Xout,Xin,product
  integer :: ja,jb,ji,jj,rank_out,rank_in,rank_op,Nsp,Abody,M
  integer :: a,b,i,j,ax,ix,jx,bx,J1,J2,dpar_in,dpar_out,dpar_op
  real(8) :: sm , phase,dcgi,mult
  
  sm = 0.d0 
  
  Nsp = jbas%total_orbits 
  Abody = sum(jbas%con) 
   
  rank_op = Trans_op%rank 
  rank_in = Xin%rank 
  rank_out = Xout%rank
  
  if (.not. triangle(rank_in,rank_op,rank_out))  then 
     transition_ME = 0.d0 
     return
  end if 

  dpar_in = Xin%dpar/2
  dpar_out = Xout%dpar/2
  dpar_op = Trans_op%dpar/2
  
  if (mod(dpar_in+dpar_op+dpar_out,2).ne.0) then
     transition_ME = 0.d0
  end if

  call duplicate_sq_op(Xout,product)
  call tensor_product(Trans_op,Xin,product,jbas)
  
  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do ix = 1,Abody 
        i = jbas%holes(ix)
        ji = jbas%jj(i) 
        
        sm = sm + f_tensor_elem(a,i,Xout,jbas)*&
             f_tensor_elem(a,i,product,jbas)
     end do 
  end do 

  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do bx = 1,Nsp-Abody
        b = jbas%parts(bx)
        jb = jbas%jj(b) 
        
        do ix = 1,Abody 
           i = jbas%holes(ix)
           ji = jbas%jj(i) 
           
           do jx = 1,Abody 
              j = jbas%holes(jx)
              jj = jbas%jj(j) 
  
              do J1 = abs(ji-jj),ji+jj,2
                 do J2 = abs(ja-jb),ja+jb,2
                    
                    sm = sm + 0.25d0*&
                         tensor_elem(a,b,i,j,J2,J1,product,jbas)*&
                         tensor_elem(a,b,i,j,J2,J1,Xout,jbas) 
                    
                 end do 
              end do 
           end do
        end do
     end do
  end do

  transition_ME = sm * sqrt(rank_in+1.d0)  * ( -1 ) ** ((rank_in+rank_out+rank_op)/2)  
  !  BY SUHONEN'S DEFINITION, I SHOULD BE DEVIDING BY Sqrt(2J+1) 
  !  BUT THE LANCZOS ALGORITHM DID THAT FOR US ALREADY. 
  
end function transition_ME
!====================================================================
!====================================================================
real(8) function transition_to_ground_ME( Trans_op , Qdag,jbas ) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: Trans_op,Qdag 
  integer :: ja,jb,ji,jj,rank,Nsp,Abody
  integer :: a,b,i,j,ax,ix,jx,bx,J1,J2
  real(8) :: sm , phase, sm1,sm2
  
  sm = 0.d0 
  
  Nsp = jbas%total_orbits 
  Abody = sum(jbas%con) 
   
  rank = Trans_op%rank 
  
  if (rank .ne. Qdag%rank)  then 
     transition_to_ground_ME = 0.d0 
     return
  end if 

  
  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do ix = 1,Abody 
        i = jbas%holes(ix)
        ji = jbas%jj(i) 
        
        phase = (-1) ** ((ja-ji)/2)
        sm = sm + f_tensor_elem(i,a,Trans_op,jbas)*&
             f_tensor_elem(a,i,Qdag,jbas)*phase
     end do
  end do 
  sm1 = sm
  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do bx = 1,Nsp-Abody
        b = jbas%parts(bx)
        jb = jbas%jj(b) 
        
        do ix = 1,Abody 
           i = jbas%holes(ix)
           ji = jbas%jj(i) 
           
           do jx = 1,Abody 
              j = jbas%holes(jx)
              jj = jbas%jj(j) 
  
              do J1 = abs(ji-jj),ji+jj,2
                 do J2 = abs(ja-jb),ja+jb,2
                    
                    phase = (-1) **((J1+J2)/2) 
                    sm = sm + phase*0.25d0*&
                         tensor_elem(i,j,a,b,J1,J2,Trans_op,jbas)*&
                         tensor_elem(a,b,i,j,J2,J1,Qdag,jbas) 
                    
                 end do 
              end do 
           end do
        end do
     end do
  end do
  sm2 = sm - sm1
  
  transition_to_ground_ME = sm * (-1.d0)**(rank/2)
  !  BY SUHONEN'S DEFINITION, I SHOULD BE DEVIDING BY Sqrt(2J+1) 
  !  BUT THE LANCZOS ALGORITHM DID THAT FOR US ALREADY. 
  
end function transition_to_ground_ME
!==========================================================================
!==========================================================================
subroutine construct_number_operator(op,H,jbas) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: op,H
  integer :: i,ix,ji
  
  call duplicate_sq_op(H,op) 
  
  do ix = 1, H%belowEF
     i= jbas%holes(ix) 
     ji = jbas%jj(i)
     op%fhh(ix,ix) = (ji+1.d0)
  end do 
  
  do ix = 1,jbas%total_orbits-H%belowEF
     i= jbas%parts(ix) 
     ji = jbas%jj(i)
     op%fpp(ix,ix) = (ji+1.d0)
  end do 

  op%E0 = sum(op%fhh)
end subroutine

real(8) function RabLAM(na,la,nb,lb,LAM) 
  ! radial integrals in suhonen chapter 6. 
  implicit none 
  
  integer,intent(in) :: na,la,nb,lb,LAM
  integer :: tau_A , tau_b,sig,sig_min,sig_max
  real(8) :: coef,sm,xxx
  
  if (mod(la+lb+LAM,2).ne.0) then 
     RabLAM = 0.d0 
     return 
  end if

  tau_a = max((lb - la + LAM)/2,0)
  tau_b = max((la - lb + LAM)/2,0)
  sig_min = max(0,na-tau_a,nb-tau_b)
  sig_max = min(na,nb) 

  coef = sqrt(factorial(na)/half_int_gamma(na+la+1))
  coef = coef*sqrt(factorial(nb)/half_int_gamma(nb+lb+1))
  coef = coef * ( -1) ** (na+nb)*factorial(tau_a)*factorial(tau_b) 

  sm = 0.d0
  do sig = sig_min,sig_max      
     xxx = (half_int_gamma((la+lb+lam)/2+sig+1)/factorial(sig))
     xxx = xxx / factorial(na-sig)
     xxx = xxx / factorial(nb-sig)
     xxx = xxx / factorial(sig+tau_a-na)
     xxx = xxx / factorial(sig+tau_b-nb)
     sm = sm + xxx
  end do

  RabLAM = coef*sm
end function

! These functions are often an overflow threat, but 
! here I don't think I need to worry, because emax= 2n+l doesn't get much
! bigger than 15, and these functions definitely work well into the 20s. 

real(8) function factorial(N) 
  ! N!
  implicit none
  integer,intent(in) :: N
  integer :: i
  real(8) :: sm 
  
  sm = 1
  
  do i = 1,N
     sm = sm * i 
  end do 
  
  factorial = sm
end function factorial

real(8) function double_factorial(N) 
  ! N!! 
  implicit none
  integer,intent(in) :: N
  integer :: i
  real(8) :: sm 
  
  sm = 1
  
  do i = N,1,-2
     sm = sm * i 
  end do 
  
  double_factorial = sm
end function double_factorial

real(8) function fac_over_doublefac(N,M)
  !  N!/M!!
  implicit none 
  integer,intent(in) :: N,M
  integer :: i,j,k 
  real(8) :: sm

  sm = 1.d0

  i = N
  j = M 
  do while ( (i > 1) .and. (j > 1))  
     sm = sm * (dfloat(i)/dfloat(j))
     i = i -1
     j = j -2
  end do

  if ( i.le.1) then
     do k = j,1,-2
        sm = sm/k
     end do
  else if (j.le.1) then
     do k = i,1,-1
        sm = sm*k
     end do
  end if

  fac_over_doublefac= sm
  
end function fac_over_doublefac
     
     
real(8) function half_int_gamma(N)
  !GAMMA(N+1/2) 
  implicit none
  
  integer,intent(in) :: N
  integer :: i
  real(8) :: sm

  sm = sqrt(acos(-1.d0))
  
  do i = 1,N 
     sm = (i-0.5)*sm
  end do 

  half_int_gamma = sm
end function   
!=====================================================================================
!=====================================================================================
subroutine tensor_product(AA,BB,CC,jbas) 
  !TENSOR PRODUCT A . B = C 
  ! mind you that B and C must be of pphh form 
  use tensor_products
  implicit none 
  
  type(spd) :: jbas 
  type(sq_op) :: AA,BB,CC 
  integer :: a,b,c,d,i,j,k,l,p1,p2,h1,h2 
  integer :: p1x,p2x,h1x,h2x,holes,parts
  integer :: jp1,jp2,jh1,jh2,ja,jb,ji,jj
  integer :: J1,J2,J1x,J2x,rank_a,rank_b,rank_c,g
  integer :: par_a,par_b,par_c,J1max,J1min,J5min,J5max
  integer :: ax,bx,ix,jx,lh1,lh2,lp1,lp2,c1,c2
  integer :: J2min,J2max,J3min,J3max,J3,q,I_BIG,J_BIG
  integer :: tp1,tp2,th1,th2,n1,n2,J4,J5,J4min,J4max
  real(8) :: sm1,sm2,sm3,sm4,sm,d6ji,coef9,pre,presum
  
  rank_a = AA%rank
  rank_b = BB%rank
  rank_c = CC%rank 

  par_a = AA%dpar/2
  par_b = BB%dpar/2
  par_c = CC%dpar/2 
  
  if (.not. triangle(rank_a,rank_b,rank_c)) then
     STOP 'non-triangular tensor product' 
  end if

  if ( mod(par_a + par_b,2) .ne. (mod(par_c,2)) ) then
     STOP 'parity violating tensor product'
  end if

  holes = sum(jbas%con)
  parts = sum(1-jbas%con) 
  
  do p1x = 1, parts
     p1 = jbas%parts(p1x)

     jp1 = jbas%jj(p1) 
     lp1 = jbas%ll(p1)
     
     do h1x = 1,holes
        h1 = jbas%holes(h1x)
        jh1 = jbas%jj(h1) 
        lh1 = jbas%ll(h1) 

        if (jbas%itzp(h1) .ne. jbas%itzp(p1)) cycle 
        if (.not. triangle(jp1,jh1,rank_c) ) cycle
        if ( mod(lh1+lp1+par_c,2) == 1 ) cycle
        sm = 0.d0 
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !   1 + 1 -> 1
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        sm1 = 0.d0 

        do ax = 1, parts
           a = jbas%parts(ax)
           ja = jbas%jj(a)
           if ( mod(jbas%ll(a)+jbas%ll(p1)+AA%dpar/2,2) .ne. 0) cycle
           if ( jbas%itzp(a).ne.jbas%itzp(p1)) cycle           

           sm1 = sm1 + f_tensor_elem(p1,a,AA,jbas) * f_tensor_elem(a,h1,BB,jbas) &
                * (-1) ** ((rank_c + jp1 + jh1)/2) * d6ji(rank_a,rank_b,rank_c,jh1,jp1,ja)

        end do

        sm2 = 0.d0 

        do ix = 1, holes
           i = jbas%holes(ix)
           ji = jbas%jj(i)
           if ( mod(jbas%ll(i)+jbas%ll(p1)+BB%dpar/2,2) .ne. 0) cycle
           if ( jbas%itzp(i).ne.jbas%itzp(p1)) cycle           
           sm2 = sm2 - f_tensor_elem(p1,i,BB,jbas) * f_tensor_elem(i,h1,AA,jbas) &
                * (-1) ** ((rank_a + rank_b + jp1 + jh1)/2) * d6ji(rank_a,rank_b,rank_c,jp1,jh1,ji)

        end do

        sm = sm + sqrt(rank_c+1.d0) * (sm1+sm2) 
        ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! !   1 + 2 -> 1
        ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        sm1 = 0.d0
        sm2 = 0.d0 
        do ix = 1,holes
           i = jbas%holes(ix)
           do ax = 1,parts
              a = jbas%parts(ax)
             
              sm1 = sm1 +  f_tensor_elem(i,a,AA,jbas) * Vgenpandya(p1,h1,i,a,rank_c,rank_a,BB,jbas)
              sm2 = sm2 +  Vgenpandya(p1,h1,a,i,rank_c,rank_b,AA,jbas) * f_tensor_elem(a,i,BB,jbas)  
           end do
        end do

        sm = sm + 1/sqrt(rank_a+1.d0) * sm1 + (-1)**((rank_a+rank_b+rank_c)/2)/sqrt(rank_b+1.d0)*sm2

        ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! !   2 + 2 -> 1
        ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        do ix = 1, holes
           i = jbas%holes(ix)
           ji = jbas%jj(i)

           J1min = abs(jp1-ji)
           J1max = jp1+ji

           do ax = 1, parts
              a = jbas%parts(ax)
              ja = jbas%jj(a)

              do bx = ax,parts
                 b = jbas%parts(bx)
                 jb = jbas%jj(b)

                 pre = 1.d0
                 if (a==b) then
                    pre = 0.5d0
                 end if 

                 do J1 = J1min,J1max,2
                    J3min = max(abs(jh1-ji),abs(rank_c-J1)) 
                    J3max = min(jh1+ji,rank_c+J1)
                    do J3=J3min,J3max,2
                       J2min = max(abs(ja-jb),abs(J1-rank_a),abs(J3-rank_b))
                       J2max = max(ja+jb,J1+rank_a,J3+rank_b)
                       do J2 = J2min,J2max,2                 
                          
                          sm = sm - d6ji(J3,J1,rank_c,jp1,jh1,ji) &
                               * d6ji(rank_a,rank_b,rank_c,J3,J1,J2) &
                               * tensor_elem(i,p1,a,b,J1,J2,AA,jbas) &
                               * tensor_elem(a,b,h1,i,J2,J3,BB,jbas) &
                               * sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0))*pre
                          
                       end do
                    end do
                 end do
              end do
           end do
        end do

        sm2 = 0.d0 
        do ax = 1, parts
           a = jbas%parts(ax)
           ja = jbas%jj(a)

           J1min = abs(jp1-ja)
           J1max = jp1+ja

           do ix = 1, holes
              i = jbas%holes(ix)
              ji = jbas%jj(i)

              do jx = ix,holes
                 j = jbas%holes(jx)
                 jj = jbas%jj(j)

                 pre =1.d0

                 if (i ==j ) pre = 0.5d0 

                 do J1 = J1min,J1max,2
                    J3min = max(abs(jh1-ja),abs(J1-rank_c))
                    J3max = min(jh1+ja,J1+rank_c)
           
                    do J3 = J3min,J3max,2
                       J2min = max(abs(ji-jj),abs(J1-rank_b),abs(J3-rank_a))
                       J2max = max(ji+jj,J1+rank_b,J3+rank_a)                
                       do J2=J2min,J2max,2

                          sm2 = sm2 + d6ji(J3,J1,rank_c,jp1,jh1,ja) &
                               * d6ji(rank_b,rank_a,rank_c,J3,J1,J2) &
                               * tensor_elem(a,p1,i,j,J1,J2,BB,jbas) &
                               * tensor_elem(i,j,h1,a,J2,J3,AA,jbas) &
                               * sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0))*pre
                       end do
                    end do
                 end do
              end do

           end do
        end do
                 

        sm = sm + sm2*(-1)**((rank_a+rank_b+rank_c)/2) 

        CC%fph(p1x,h1x) = sm
        
     end do
  end do

  do q = 1, CC%nblocks
     J1x = CC%tblck(q)%Jpair(1)
     J2x = CC%tblck(q)%Jpair(2)
     do g=3,7,4

        n1 = size(CC%tblck(q)%tgam(g)%X(:,1))
        n2 = size(CC%tblck(q)%tgam(g)%X(1,:))
        ! main calculation

        do I_BIG = 1,n1

           if (g == 3)  then 
              p1 = CC%tblck(q)%tensor_qn(1,1)%Y(I_BIG,1)
              jp1 = jbas%jj(p1)           
              lp1 = jbas%ll(p1)
              tp1 = jbas%itzp(p1)

              p2 = CC%tblck(q)%tensor_qn(1,1)%Y(I_BIG,2)
              jp2 = jbas%jj(p2)
              lp2 = jbas%ll(p2)
              tp2 = jbas%itzp(p2)
              J1 = J1x
              J2 = J2x             
           else
              h1 = CC%tblck(q)%tensor_qn(3,1)%Y(I_BIG,1)
              jh1 = jbas%jj(h1)
              lh1 = jbas%ll(h1)
              th1 = jbas%itzp(h1)

              h2 = CC%tblck(q)%tensor_qn(3,1)%Y(I_BIG,2)
              jh2 = jbas%jj(h2)
              lh2 = jbas%ll(h2)
              th2 = jbas%itzp(h2)
              J2 = J1x
              J1 = J2x 
           end if

           do J_BIG = 1,n2

              if (g==3) then 
                 h1 = CC%tblck(q)%tensor_qn(3,2)%Y(J_BIG,1)
                 jh1 = jbas%jj(h1)                 
                 lh1 = jbas%ll(h1)
                 th1 = jbas%itzp(h1)

                 h2 = CC%tblck(q)%tensor_qn(3,2)%Y(J_BIG,2)
                 jh2 = jbas%jj(h2)
                 lh2 = jbas%ll(h2)
                 th2 = jbas%itzp(h2)

              else
                 p1 = CC%tblck(q)%tensor_qn(1,2)%Y(J_BIG,1)
                 jp1 = jbas%jj(p1)                 
                 lp1 = jbas%ll(p1)
                 tp1 = jbas%itzp(p1)

                 p2 = CC%tblck(q)%tensor_qn(1,2)%Y(J_BIG,2)
                 jp2 = jbas%jj(p2)
                 lp2 = jbas%ll(p2)
                 tp2 = jbas%itzp(p2)
              end if

              pre = 1.d0 
              if (p1==p2) pre = sqrt(0.5d0)
              if (h1==h2) pre = pre*sqrt(0.5d0)

!!! DISCONNECTED DIAGRAMS ==============================================================
              sm = sqrt((J1+1.d0)*(J2+1.d0)*(rank_c+1.d0)) * ( &
                   f_tensor_elem(p1,h1,AA,jbas) * f_tensor_elem(p2,h2,BB,jbas) * &
                   coef9(jp1,jh1,rank_a,jp2,jh2,rank_b,J1,J2,rank_c)&
                   
                   - f_tensor_elem(p2,h1,AA,jbas) * f_tensor_elem(p1,h2,BB,jbas) * &
                   coef9(jp2,jh1,rank_a,jp1,jh2,rank_b,J1,J2,rank_c) * (-1)**((jp1+jp2+J1)/2)&

                   + f_tensor_elem(p2,h2,AA,jbas) * f_tensor_elem(p1,h1,BB,jbas) * &
                   coef9(jp2,jh2,rank_a,jp1,jh1,rank_b,J1,J2,rank_c) * (-1)**((jp1+jp2+jh1+jh2+J1+J2)/2)&

                   - f_tensor_elem(p1,h2,AA,jbas) * f_tensor_elem(p2,h1,BB,jbas) * &
                   coef9(jp1,jh2,rank_a,jp2,jh1,rank_b,J1,J2,rank_c) * (-1)**((jh1+jh2+J2)/2) ) 
!!! ==================================================================================
              do ax = 1, parts
                 a = jbas%parts(ax)
                 ja = jbas%jj(a)

                 sm1 = 0.d0
                 if (jbas%itzp(a)== tp1 ) then
                    if (mod(jbas%ll(a)+par_a+lp1,2)==0) then 
                       if (triangle(ja,jp1,rank_a) )then

                          J3min = max(abs(ja - jp2),abs(rank_a-J1),abs(rank_b-J2))
                          J3max =  min(ja + jp2,rank_a+J1,rank_b+J2)

                          do J3 = J3min,J3max,2 

                             sm1 = sm1 + (-1) ** (J3/2) * d6ji(J3,J1,rank_a,jp1,ja,jp2) &
                                  * d6ji(rank_a,rank_b,rank_c,J2,J1,J3) * f_tensor_elem(p1,a,AA,jbas) &
                                  * tensor_elem(a,p2,h1,h2,J3,J2,BB,jbas)* sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 
                          end do

                          sm = sm + (-1) ** ((jp1+jp2+J1+J2+rank_a+rank_c)/2) * sm1 

                       end if
                    end if
                 end if

                 sm1 = 0.d0
                 if (jbas%itzp(a)== tp2 ) then
                    if (mod(jbas%ll(a)+par_a+lp2,2)==0) then 
                       if (triangle(ja,jp2,rank_a) )then

                          J3min = max(abs(ja - jp1),abs(rank_a-J1),abs(rank_b-J2))
                          J3max = min(ja + jp1,rank_a+J1,rank_b+J2)

                          do J3 = J3min,J3max,2 

                             sm1 = sm1 - (-1) ** (J3/2) * d6ji(J3,J1,rank_a,jp2,ja,jp1) &
                                  * d6ji(rank_a,rank_b,rank_c,J2,J1,J3) * f_tensor_elem(p2,a,AA,jbas) &
                                  * tensor_elem(a,p1,h1,h2,J3,J2,BB,jbas)*sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0))  

                          end do
                          sm = sm + (-1) ** ((J2+rank_a+rank_c)/2) * sm1 
                       end if
                    end if
                 end if

                 sm1 = 0.d0 
                 if (jbas%itzp(a)== th1 ) then
                    if (mod(jbas%ll(a)+par_b+lh1,2)==0) then 
                       if (triangle(ja,jh1,rank_b) )then

                          J3min = max(abs(ja - jh2),abs(rank_a-J1),abs(rank_b-J2))
                          J3max = min(ja + jh2,rank_a+J1,rank_b+J2)

                          do J3 = J3min,J3max,2
                             sm1 = sm1 - (-1)**(J3/2) * d6ji(J3,J2,rank_b,jh1,ja,jh2)&
                                  * d6ji(rank_a,rank_b,rank_c,J2,J1,J3) * f_tensor_elem(a,h1,BB,jbas)&
                                  * tensor_elem(p1,p2,h2,a,J1,J3,AA,jbas)*sqrt((J2+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 
                          end do
                          sm = sm + sm1 * (-1)**((J1+rank_b+rank_c)/2) 
                       end if
                    end if
                 end if

                 sm1 = 0.d0 
                 if (jbas%itzp(a)== th2 ) then
                    if (mod(jbas%ll(a)+par_b+lh2,2)==0) then 
                       if (triangle(ja,jh2,rank_b) )then

                          J3min = max(abs(ja - jh1),abs(rank_a-J1),abs(rank_b-J2))
                          J3max = min(ja + jh1,rank_a+J1,rank_b+J2)

                          do J3 = J3min,J3max,2
                             sm1 = sm1 + (-1)**(J3/2) * d6ji(J3,J2,rank_b,jh2,ja,jh1)&
                                  * d6ji(rank_a,rank_b,rank_c,J2,J1,J3) * f_tensor_elem(a,h2,BB,jbas)&
                                  * tensor_elem(p1,p2,h1,a,J1,J3,AA,jbas)*sqrt((J2+1.d0)*(J3+1.d0)*(rank_c+1.d0))  
                          end do
                          sm = sm + sm1 * (-1)**((jh1+jh2+J2+J1+rank_b+rank_c)/2) 
                       end if
                    end if
                 end if
              end do


              do ix = 1, holes
                 i = jbas%holes(ix)
                 ji = jbas%jj(i)

                 sm1 = 0.d0
                 if (jbas%itzp(i)== th1 ) then
                    if (mod(jbas%ll(i)+par_a+lh1,2)==0) then 
                       if (triangle(ji,jh1,rank_a) )then

                          J3min = max(abs(ji - jh2),abs(rank_a-J2),abs(rank_b-J1))
                          J3max = min(ji + jh2,rank_a+J2,rank_b+J1)

                          do J3 = J3min,J3max,2 

                             sm1 = sm1 + (-1) ** (J3/2) * d6ji(J3,J2,rank_a,jh1,ji,jh2) &
                                  * d6ji(rank_a,rank_b,rank_c,J1,J2,J3) * f_tensor_elem(i,h1,AA,jbas) &
                                  * tensor_elem(p1,p2,h2,i,J1,J3,BB,jbas)*sqrt((J2+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 

                          end do
                          sm = sm + (-1) ** ((J1+rank_b)/2) * sm1 
                       end if
                    end if
                 end if

                 sm1 = 0.d0
                 if (jbas%itzp(i)== th2 ) then
                    if (mod(jbas%ll(i)+par_a+lh2,2)==0) then 
                       if (triangle(ji,jh2,rank_a) )then

                          J3min = max(abs(ji - jh1),abs(rank_a-J2),abs(rank_b-J1))
                          J3max = min(ji + jh1,rank_a+J2,rank_b+J1)

                          do J3 = J3min,J3max,2 

                             sm1 = sm1 - (-1) ** (J3/2) * d6ji(J3,J2,rank_a,jh2,ji,jh1) &
                                  * d6ji(rank_a,rank_b,rank_c,J1,J2,J3) * f_tensor_elem(i,h2,AA,jbas) &
                                  * tensor_elem(p1,p2,h1,i,J1,J3,BB,jbas)*sqrt((J2+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 

                          end do
                          sm = sm + (-1) ** ((jh1+jh2+J1+J2+rank_b)/2) * sm1 
                       end if
                    end if
                 end if

                 sm1 = 0.d0 
                 if (jbas%itzp(i)== tp1 ) then
                    if (mod(jbas%ll(i)+par_b+lp1,2)==0) then 
                       if (triangle(ji,jp1,rank_b) )then

                          J3min = max(abs(ji - jp2),abs(rank_a-J2),abs(rank_b-J1))
                          J3max = min(ji + jp2,rank_a+J2,rank_b+J1)

                          do J3 = J3min,J3max,2
                             sm1 = sm1 - (-1)**(J3/2) * d6ji(J3,J1,rank_b,jp1,ji,jp2)&
                                  * d6ji(rank_a,rank_b,rank_c,J1,J2,J3) * f_tensor_elem(p1,i,BB,jbas)&
                                  * tensor_elem(i,p2,h1,h2,J3,J2,AA,jbas)*sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 
                          end do
                          sm = sm + sm1 * (-1)**((jp1+jp2+J2+J1+rank_a)/2) 
                       end if
                    end if
                 end if

                 sm1 = 0.d0 
                 if (jbas%itzp(i)== tp2 ) then
                    if (mod(jbas%ll(i)+par_b+lp2,2)==0) then 
                       if (triangle(ji,jp2,rank_b) )then

                          J3min = max(abs(ji - jp1),abs(rank_a-J2),abs(rank_b-J1))
                          J3max = min(ji + jp1,rank_a+J2,rank_b+J1)

                          do J3 = J3min,J3max,2
                             sm1 = sm1 + (-1)**(J3/2) * d6ji(J3,J1,rank_b,jp2,ji,jp1)&
                                  * d6ji(rank_a,rank_b,rank_c,J1,J2,J3) * f_tensor_elem(p2,i,BB,jbas)&
                                  * tensor_elem(i,p1,h1,h2,J3,J2,AA,jbas)*sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0))  
                          end do
                          sm = sm + sm1 * (-1)**((J2+rank_a)/2) 
                       end if
                    end if
                 end if

              end do
              CC%tblck(q)%tgam(g)%X(I_BIG,J_BIG) = sm*pre
           end do
        end do
     end do
  end do

   call tensor_product_222_pp_hh(AA,BB,CC,jbas)
   call tensor_product_222_ph(AA,BB,CC,jbas) 
  ! do nothing
end subroutine tensor_product
!=====================================================================================
!=====================================================================================
subroutine tensor_dTz_tensor_product(AA,BB,CC,jbas) 
  !TENSOR PRODUCT A . B = C 
  ! mind you that B and C must be of pphh form 
  use tensor_products
  implicit none 
  
  type(spd) :: jbas 
  type(sq_op) :: AA
  type(iso_ladder) :: BB,CC 
  integer :: a,b,c,d,i,j,k,l,p1,p2,h1,h2 
  integer :: p1x,p2x,h1x,h2x,holes,parts
  integer :: jp1,jp2,jh1,jh2,ja,jb,ji,jj
  integer :: J1,J2,J1x,J2x,rank_a,rank_b,rank_c,g
  integer :: par_a,par_b,par_c,J1max,J1min,J5min,J5max
  integer :: ax,bx,ix,jx,lh1,lh2,lp1,lp2,c1,c2
  integer :: J2min,J2max,J3min,J3max,J3,q,I_BIG,J_BIG
  integer :: tp1,tp2,th1,th2,n1,n2,J4,J5,J4min,J4max
  real(8) :: sm1,sm2,sm3,sm4,sm,d6ji,coef9,pre,presum
  
  rank_a = AA%rank
  rank_b = BB%rank
  rank_c = CC%rank 

  par_a = AA%dpar/2
  par_b = BB%dpar/2
  par_c = CC%dpar/2 
  
  if (.not. triangle(rank_a,rank_b,rank_c)) then
     STOP 'non-triangular tensor product' 
  end if

  if ( mod(par_a + par_b,2) .ne. (mod(par_c,2)) ) then
     STOP 'parity violating tensor product'
  end if

  holes = sum(jbas%con)
  parts = sum(1-jbas%con) 
  
  do p1x = 1, parts
     p1 = jbas%parts(p1x)

     jp1 = jbas%jj(p1) 
     lp1 = jbas%ll(p1)
     
     do h1x = 1,holes
        h1 = jbas%holes(h1x)
        jh1 = jbas%jj(h1) 
        lh1 = jbas%ll(h1) 

        if (jbas%itzp(h1)  .ne. jbas%itzp(p1) - 2*CC%dTz ) cycle 
        if (.not. triangle(jp1,jh1,rank_c) ) cycle
        if ( mod(lh1+lp1+par_c,2) == 1 ) cycle
        sm = 0.d0 
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !   1 + 1 -> 1
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        sm1 = 0.d0 

        do ax = 1, parts
           a = jbas%parts(ax)
           ja = jbas%jj(a)
           if ( mod(jbas%ll(a)+jbas%ll(p1)+AA%dpar/2,2) .ne. 0) cycle
           if ( jbas%itzp(a).ne.jbas%itzp(p1)) cycle           

           sm1 = sm1 + f_tensor_elem(p1,a,AA,jbas) * f_iso_ladder_elem(a,h1,BB,jbas) &
                * (-1) ** ((rank_c + jp1 + jh1)/2) * d6ji(rank_a,rank_b,rank_c,jh1,jp1,ja)

        end do

        sm2 = 0.d0 

        do ix = 1, holes
           i = jbas%holes(ix)
           ji = jbas%jj(i)
           if ( mod(jbas%ll(i)+jbas%ll(p1)+BB%dpar/2,2) .ne. 0) cycle
           if ( jbas%itzp(i).ne.jbas%itzp(p1)-BB%dTz*2) cycle           
           sm2 = sm2 - f_iso_ladder_elem(p1,i,BB,jbas) * f_tensor_elem(i,h1,AA,jbas) &
                * (-1) ** ((rank_a + rank_b + jp1 + jh1)/2) * d6ji(rank_a,rank_b,rank_c,jp1,jh1,ji)

        end do

        sm = sm + sqrt(rank_c+1.d0) * (sm1+sm2) 
        ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! !   1 + 2 -> 1
        ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        sm1 = 0.d0
        sm2 = 0.d0 
        do ix = 1,holes
           i = jbas%holes(ix)
           do ax = 1,parts
              a = jbas%parts(ax)
             
              sm1 = sm1 +  f_tensor_elem(i,a,AA,jbas) * Visopandya(p1,h1,i,a,rank_c,rank_a,BB,jbas)
              sm2 = sm2 +  Vgenpandya(p1,h1,a,i,rank_c,rank_b,AA,jbas) * f_iso_ladder_elem(a,i,BB,jbas)  
           end do
        end do

        sm = sm + 1/sqrt(rank_a+1.d0) * sm1 + (-1)**((rank_a+rank_b+rank_c)/2)/sqrt(rank_b+1.d0)*sm2

        ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! !   2 + 2 -> 1
        ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        do ix = 1, holes
           i = jbas%holes(ix)
           ji = jbas%jj(i)

           J1min = abs(jp1-ji)
           J1max = jp1+ji

           do ax = 1, parts
              a = jbas%parts(ax)
              ja = jbas%jj(a)

              do bx = ax,parts
                 b = jbas%parts(bx)
                 jb = jbas%jj(b)

                 pre = 1.d0
                 if (a==b) then
                    pre = 0.5d0
                 end if 

                 do J1 = J1min,J1max,2
                    J3min = max(abs(jh1-ji),abs(rank_c-J1)) 
                    J3max = min(jh1+ji,rank_c+J1)
                    do J3=J3min,J3max,2
                       J2min = max(abs(ja-jb),abs(J1-rank_a),abs(J3-rank_b))
                       J2max = max(ja+jb,J1+rank_a,J3+rank_b)
                       do J2 = J2min,J2max,2                 
                          
                          sm = sm - d6ji(J3,J1,rank_c,jp1,jh1,ji) &
                               * d6ji(rank_a,rank_b,rank_c,J3,J1,J2) &
                               * tensor_elem(i,p1,a,b,J1,J2,AA,jbas) &
                               * iso_ladder_elem(a,b,h1,i,J2,J3,BB,jbas) &
                               * sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0))*pre
                          
                       end do
                    end do
                 end do
              end do
           end do
        end do

        sm2 = 0.d0 
        do ax = 1, parts
           a = jbas%parts(ax)
           ja = jbas%jj(a)

           J1min = abs(jp1-ja)
           J1max = jp1+ja

           do ix = 1, holes
              i = jbas%holes(ix)
              ji = jbas%jj(i)

              do jx = ix,holes
                 j = jbas%holes(jx)
                 jj = jbas%jj(j)

                 pre =1.d0

                 if (i ==j ) pre = 0.5d0 

                 do J1 = J1min,J1max,2
                    J3min = max(abs(jh1-ja),abs(J1-rank_c))
                    J3max = min(jh1+ja,J1+rank_c)
           
                    do J3 = J3min,J3max,2
                       J2min = max(abs(ji-jj),abs(J1-rank_b),abs(J3-rank_a))
                       J2max = max(ji+jj,J1+rank_b,J3+rank_a)                
                       do J2=J2min,J2max,2

                          sm2 = sm2 + d6ji(J3,J1,rank_c,jp1,jh1,ja) &
                               * d6ji(rank_b,rank_a,rank_c,J3,J1,J2) &
                               * iso_ladder_elem(a,p1,i,j,J1,J2,BB,jbas) &
                               * tensor_elem(i,j,h1,a,J2,J3,AA,jbas) &
                               * sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0))*pre
                       end do
                    end do
                 end do
              end do

           end do
        end do
                 

       sm = sm + sm2*(-1)**((rank_a+rank_b+rank_c)/2) 

        CC%fph(p1x,h1x) = sm
        
     end do
  end do

  do q = 1, CC%nblocks
     J1x = CC%tblck(q)%Jpair(1)
     J2x = CC%tblck(q)%Jpair(2)
     
     n1 = CC%tblck(q)%npp1
     n2 = CC%tblck(q)%nhh2
     ! main calculation
     
     do I_BIG = 1,n1
       
        p1 = CC%tblck(q)%qn(1)%Y(I_BIG,1)
        jp1 = jbas%jj(p1)           
        lp1 = jbas%ll(p1)
        tp1 = jbas%itzp(p1)
        
        p2 = CC%tblck(q)%qn(1)%Y(I_BIG,2)
        jp2 = jbas%jj(p2)
        lp2 = jbas%ll(p2)
        tp2 = jbas%itzp(p2)
        J1 = J1x
        J2 = J2x             
   
        do J_BIG = 1,n2


           h1 = CC%tblck(q)%qn(2)%Y(J_BIG,1)
           jh1 = jbas%jj(h1)                 
           lh1 = jbas%ll(h1)
           th1 = jbas%itzp(h1)

           h2 = CC%tblck(q)%qn(2)%Y(J_BIG,2)
           jh2 = jbas%jj(h2)
           lh2 = jbas%ll(h2)
           th2 = jbas%itzp(h2)

           pre = 1.d0 
           if (p1==p2) pre = sqrt(0.5d0)
           if (h1==h2) pre = pre*sqrt(0.5d0)
           sm = 0.d0 
!!! DISCONNECTED DIAGRAMS ==============================================================
              sm = sqrt((J1+1.d0)*(J2+1.d0)*(rank_c+1.d0)) * ( &
                   f_tensor_elem(p1,h1,AA,jbas) * f_iso_ladder_elem(p2,h2,BB,jbas) * &
                   coef9(jp1,jh1,rank_a,jp2,jh2,rank_b,J1,J2,rank_c)&
                   
                   - f_tensor_elem(p2,h1,AA,jbas) * f_iso_ladder_elem(p1,h2,BB,jbas) * &
                   coef9(jp2,jh1,rank_a,jp1,jh2,rank_b,J1,J2,rank_c) * (-1)**((jp1+jp2+J1)/2)&

                   + f_tensor_elem(p2,h2,AA,jbas) * f_iso_ladder_elem(p1,h1,BB,jbas) * &
                   coef9(jp2,jh2,rank_a,jp1,jh1,rank_b,J1,J2,rank_c) * (-1)**((jp1+jp2+jh1+jh2+J1+J2)/2)&

                   - f_tensor_elem(p1,h2,AA,jbas) * f_iso_ladder_elem(p2,h1,BB,jbas) * &
                   coef9(jp1,jh2,rank_a,jp2,jh1,rank_b,J1,J2,rank_c) * (-1)**((jh1+jh2+J2)/2) )
              
!!! ==================================================================================
           
           do ax = 1, parts
              a = jbas%parts(ax)
              ja = jbas%jj(a)

              sm1 = 0.d0
              if (jbas%itzp(a)== tp1 ) then
                 if (mod(jbas%ll(a)+par_a+lp1,2)==0) then 
                    if (triangle(ja,jp1,rank_a) )then
                       
                       J3min = max(abs(ja - jp2),abs(rank_a-J1),abs(rank_b-J2))
                       J3max =  min(ja + jp2,rank_a+J1,rank_b+J2)
                       
                       do J3 = J3min,J3max,2 
                          
                          sm1 = sm1 + (-1) ** (J3/2) * d6ji(J3,J1,rank_a,jp1,ja,jp2) &
                               * d6ji(rank_a,rank_b,rank_c,J2,J1,J3) * f_tensor_elem(p1,a,AA,jbas) &
                               * iso_ladder_elem(a,p2,h1,h2,J3,J2,BB,jbas)* sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 
                       end do

                       sm = sm + (-1) ** ((jp1+jp2+J1+J2+rank_a+rank_c)/2) * sm1 

                    end if
                 end if
              end if

              sm1 = 0.d0
              if (jbas%itzp(a)== tp2 ) then
                 if (mod(jbas%ll(a)+par_a+lp2,2)==0) then 
                    if (triangle(ja,jp2,rank_a) )then

                       J3min = max(abs(ja - jp1),abs(rank_a-J1),abs(rank_b-J2))
                       J3max = min(ja + jp1,rank_a+J1,rank_b+J2)

                       do J3 = J3min,J3max,2 

                          sm1 = sm1 - (-1) ** (J3/2) * d6ji(J3,J1,rank_a,jp2,ja,jp1) &
                               * d6ji(rank_a,rank_b,rank_c,J2,J1,J3) * f_tensor_elem(p2,a,AA,jbas) &
                               * iso_ladder_elem(a,p1,h1,h2,J3,J2,BB,jbas)*sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0))  

                       end do
                       sm = sm + (-1) ** ((J2+rank_a+rank_c)/2) * sm1 
                    end if
                 end if
              end if

              sm1 = 0.d0 
              if (jbas%itzp(a) - BB%dTz*2 == th1 ) then
                 if (mod(jbas%ll(a)+par_b+lh1,2)==0) then 
                    if (triangle(ja,jh1,rank_b) )then

                       J3min = max(abs(ja - jh2),abs(rank_a-J1),abs(rank_b-J2))
                       J3max = min(ja + jh2,rank_a+J1,rank_b+J2)

                       do J3 = J3min,J3max,2
                          sm1 = sm1 - (-1)**(J3/2) * d6ji(J3,J2,rank_b,jh1,ja,jh2)&
                               * d6ji(rank_a,rank_b,rank_c,J2,J1,J3) * f_iso_ladder_elem(a,h1,BB,jbas)&
                               * tensor_elem(p1,p2,h2,a,J1,J3,AA,jbas)*sqrt((J2+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 
                       end do
                       sm = sm + sm1 * (-1)**((J1+rank_b+rank_c)/2) 
                    end if
                 end if
              end if

              sm1 = 0.d0 
              if (jbas%itzp(a)  - BB%dTz*2 == th2 ) then
                 if (mod(jbas%ll(a)+par_b+lh2,2)==0) then 
                    if (triangle(ja,jh2,rank_b) )then

                       J3min = max(abs(ja - jh1),abs(rank_a-J1),abs(rank_b-J2))
                       J3max = min(ja + jh1,rank_a+J1,rank_b+J2)

                       do J3 = J3min,J3max,2
                          sm1 = sm1 + (-1)**(J3/2) * d6ji(J3,J2,rank_b,jh2,ja,jh1)&
                               * d6ji(rank_a,rank_b,rank_c,J2,J1,J3) * f_iso_ladder_elem(a,h2,BB,jbas)&
                               * tensor_elem(p1,p2,h1,a,J1,J3,AA,jbas)*sqrt((J2+1.d0)*(J3+1.d0)*(rank_c+1.d0))  
                       end do
                       sm = sm + sm1 * (-1)**((jh1+jh2+J2+J1+rank_b+rank_c)/2) 
                    end if
                 end if
              end if
           end do


           do ix = 1, holes
              i = jbas%holes(ix)
              ji = jbas%jj(i)

              sm1 = 0.d0
              if (jbas%itzp(i)== th1 ) then
                 if (mod(jbas%ll(i)+par_a+lh1,2)==0) then 
                    if (triangle(ji,jh1,rank_a) )then

                       J3min = max(abs(ji - jh2),abs(rank_a-J2),abs(rank_b-J1))
                       J3max = min(ji + jh2,rank_a+J2,rank_b+J1)

                       do J3 = J3min,J3max,2 

                          sm1 = sm1 + (-1) ** (J3/2) * d6ji(J3,J2,rank_a,jh1,ji,jh2) &
                               * d6ji(rank_a,rank_b,rank_c,J1,J2,J3) * f_tensor_elem(i,h1,AA,jbas) &
                               * iso_ladder_elem(p1,p2,h2,i,J1,J3,BB,jbas)*sqrt((J2+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 

                       end do
                       sm = sm + (-1) ** ((J1+rank_b)/2) * sm1 
                    end if
                 end if
              end if

              sm1 = 0.d0
              if (jbas%itzp(i)== th2 ) then
                 if (mod(jbas%ll(i)+par_a+lh2,2)==0) then 
                    if (triangle(ji,jh2,rank_a) )then

                       J3min = max(abs(ji - jh1),abs(rank_a-J2),abs(rank_b-J1))
                       J3max = min(ji + jh1,rank_a+J2,rank_b+J1)

                       do J3 = J3min,J3max,2 

                          sm1 = sm1 - (-1) ** (J3/2) * d6ji(J3,J2,rank_a,jh2,ji,jh1) &
                               * d6ji(rank_a,rank_b,rank_c,J1,J2,J3) * f_tensor_elem(i,h2,AA,jbas) &
                               * iso_ladder_elem(p1,p2,h1,i,J1,J3,BB,jbas)*sqrt((J2+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 

                       end do
                       sm = sm + (-1) ** ((jh1+jh2+J1+J2+rank_b)/2) * sm1 
                    end if
                 end if
              end if

              sm1 = 0.d0 
              if (jbas%itzp(i)== tp1-BB%dTz*2 ) then
                 if (mod(jbas%ll(i)+par_b+lp1,2)==0) then 
                    if (triangle(ji,jp1,rank_b) )then

                       J3min = max(abs(ji - jp2),abs(rank_a-J2),abs(rank_b-J1))
                       J3max = min(ji + jp2,rank_a+J2,rank_b+J1)

                       do J3 = J3min,J3max,2
                          sm1 = sm1 - (-1)**(J3/2) * d6ji(J3,J1,rank_b,jp1,ji,jp2)&
                               * d6ji(rank_a,rank_b,rank_c,J1,J2,J3) * f_iso_ladder_elem(p1,i,BB,jbas)&
                               * tensor_elem(i,p2,h1,h2,J3,J2,AA,jbas)*sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0)) 
                       end do
                       sm = sm + sm1 * (-1)**((jp1+jp2+J2+J1+rank_a)/2) 
                    end if
                 end if
              end if

              sm1 = 0.d0 
              if (jbas%itzp(i)== tp2-BB%dTz*2 ) then
                 if (mod(jbas%ll(i)+par_b+lp2,2)==0) then 
                    if (triangle(ji,jp2,rank_b) )then

                       J3min = max(abs(ji - jp1),abs(rank_a-J2),abs(rank_b-J1))
                       J3max = min(ji + jp1,rank_a+J2,rank_b+J1)

                       do J3 = J3min,J3max,2
                          sm1 = sm1 + (-1)**(J3/2) * d6ji(J3,J1,rank_b,jp2,ji,jp1)&
                               * d6ji(rank_a,rank_b,rank_c,J1,J2,J3) * f_iso_ladder_elem(p2,i,BB,jbas)&
                               * tensor_elem(i,p1,h1,h2,J3,J2,AA,jbas)*sqrt((J1+1.d0)*(J3+1.d0)*(rank_c+1.d0))  
                       end do
                       sm = sm + sm1 * (-1)**((J2+rank_a)/2) 
                    end if
                 end if
              end if

           end do
           CC%tblck(q)%Xpphh(I_BIG,J_BIG) = sm*pre
        end do
     end do
  end do

  call dTZ_tensor_product_222_pp_hh(AA,BB,CC,jbas)
  call dTz_tensor_product_222_ph(AA,BB,CC,jbas) 
  ! do nothing
end subroutine tensor_dTz_tensor_product
!=========================================================================================================
!=========================================================================================================
subroutine EOM_observables( ladder_ops, iso_ops, O1,HS, Hcm, trans, mom, eom_states , jbas)
  implicit none

  type(sq_op),dimension(:) :: ladder_ops
  type(iso_ladder),dimension(:) :: iso_ops
  type(sq_op) :: O1,HS,Hcm
  type(spd) :: jbas
  type(eom_mgr) :: eom_states
  type(obsv_mgr) :: trans,mom
  integer :: q,Jin,Jout,Pin,Pout,in,out,states,instate,DTz
  integer :: i,a,b,j,jx,ix,ax,bx,ja,ji,jj,jb,rank
  logical :: to_ground
  real(8) :: Mfi,strength_down,strength_up,moment,dcgi,dcgi00
  real(8) :: E_in, E_out,sm,d6ji
  real(8),allocatable,dimension(:) :: STRENGTHS,MOMENTS,ENERGIES
  character(2) :: flts,dTzlab 
  character(1) :: statlab

  Mfi = dcgi00()
  ! CALCULATE TRANSITIONS 
  do q = 1, trans%num
     to_ground = .false.
     dTz = trans%dtz(q)
     read(trans%Jpi1(q)(1:1),'(I1)') Jin
     
     if (trans%Jpi1(q)(2:2) == '+' ) then
        Pin = 0
     else
        Pin = 2
     end if
     IF ( trans%Jpi2(q) == 'GS' ) then
        Jout = 0
        to_ground = .true.
     else
        read(trans%Jpi2(q)(1:1),'(I1)') Jout  
        if (trans%Jpi2(q)(2:2) == '+' ) then
           Pout = 0
        else
           Pout = 2
        end if
     end if


     print*
     print*, '============================================================================='
     print*, '        E_in                E_out         B('&
          //trans%oper//';'//trans%Jpi1(q)//' -> '//trans%Jpi2(q)//&
          ')      B('//trans%oper//';'//trans%Jpi2(q)//' -> '//trans%Jpi1(q)//')' 
     print*, '============================================================================='

     Jin = 2* Jin
     Jout = 2* Jout


     if (dTZ ==0 ) then 
        If (to_Ground) then
           ! count number of states
           states = 0.d0 
           do in = 1,size(ladder_ops)
              IF ( ladder_ops(in)%rank .ne. Jin) cycle
              IF ( ladder_ops(in)%dpar .ne. Pin) cycle
              states = states + 1 
           end do

           allocate(Strengths(states),Energies(states))
           states = 2*states + 1 ! energies, observs, and gs
           write(flts,'(I2)') states       
           states = 0

           do in = 1, size(ladder_ops)

              IF ( ladder_ops(in)%rank .ne. Jin) cycle
              IF ( ladder_ops(in)%dpar .ne. Pin) cycle

              Mfi = transition_to_ground_ME(O1,ladder_ops(in),jbas)

              strength_down = Mfi * Mfi /(ladder_ops(in)%rank+1.d0)
              strength_up = Mfi * Mfi

              E_in = ladder_ops(in)%E0
              E_out = 0.d0
              states = states + 1
              strengths(states) = strength_down
              Energies(states) = ladder_ops(in)%E0

              write(*,'(4(f19.12))') E_in,E_out,Strength_down,Strength_up           

           end do

           open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
                '_energies_strengths_'//trans%oper//'_'//trans%Jpi1(q)//'_'//trans%Jpi2(q)//'.dat')
           write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') nint(HS%hospace),HS%eMax,HS%E0,Energies,strengths
           close(31)
           
           deallocate(Energies,strengths)  
           
        else
           instate = 0
           do In = 1, size(ladder_ops)

              IF ( ladder_ops(in)%rank .ne. Jin) cycle
              IF ( ladder_ops(in)%dpar .ne. Pin) cycle

              states = 0.d0 
              do out = 1,size(ladder_ops)
                 IF ( out == in) cycle
                 IF ( ladder_ops(out)%rank .ne. Jout) cycle
                 IF ( ladder_ops(out)%dpar .ne. Pout) cycle
                 states = states + 1 
              end do

              allocate(Strengths(states),Energies(states))
              states = 2*states + 1 ! energies, observs, and gs
              write(flts,'(I2)') states       
              states = 0

              do out = 1, size(ladder_ops)


                 if (out==in) cycle ! not a transition
                 IF ( ladder_ops(out)%rank .ne. Jout) cycle
                 IF ( ladder_ops(out)%dpar .ne. Pout) cycle
                 states = states + 1
                 Mfi = transition_ME(ladder_ops(out),O1,ladder_ops(in),jbas) 

                 if (ladder_ops(out)%E0 > ladder_ops(in)%E0) then 
                    strength_down = Mfi * Mfi / (ladder_ops(out)%rank+1.d0)
                    strength_up = Mfi * Mfi / (ladder_ops(in)%rank+1.d0)
                    E_in = ladder_ops(out)%E0
                    E_out = ladder_ops(in)%E0
                    strengths(states) = strength_up
                    Energies(states) = E_in 
                    write(*,'(4(f19.12))') E_in,E_out,Strength_up,Strength_down           
                 else
                    strength_down = Mfi * Mfi / (ladder_ops(in)%rank+1.d0)
                    strength_up = Mfi * Mfi / (ladder_ops(out)%rank+1.d0)
                    E_in = ladder_ops(in)%E0
                    E_out = ladder_ops(out)%E0
                    strengths(states) = strength_down
                    Energies(states) = E_out 
                    write(*,'(4(f19.12))') E_in,E_out,Strength_down,Strength_up           
                 end if


              end do
              instate = instate+1
              write(statlab,'(I1)') instate
              write(dTzlab,'(I2)') dTz
              dTzlab = adjustl(dTzlab) 
              open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
                   '_energies_strengths_'//trans%oper//'_'//trans%Jpi1(q)//'_'//statlab//'_'&
                   //trans%Jpi2(q)//'_dTz'//trim(dTzlab)//'.dat')
              write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') &
                   nint(HS%hospace),HS%eMax,ladder_ops(in)%E0,Energies,strengths
              close(31)
              deallocate(Energies,strengths)  

           end do
        end if
     else
        
        instate = 0
        do In = 1, size(iso_ops)
           
           IF ( iso_ops(in)%rank .ne. Jin) cycle
           IF ( iso_ops(in)%dpar .ne. Pin) cycle
           IF ( iso_ops(in)%dTz .ne. dTz ) cycle          
           
           states = 0.d0 
           do out = 1,size(iso_ops)
              IF ( out == in) cycle
              IF ( iso_ops(out)%rank .ne. Jout) cycle
              IF ( iso_ops(out)%dpar .ne. Pout) cycle
              IF ( iso_ops(out)%dTz .ne. dTz ) cycle
              states = states + 1 
           end do

           allocate(Strengths(states),Energies(states))
           states = 2*states + 1 ! energies, observs, and gs
           write(flts,'(I2)') states       
           states = 0

           do out = 1, size(iso_ops)


              if (out==in) cycle ! not a transition
              IF ( iso_ops(out)%rank .ne. Jout) cycle
              IF ( iso_ops(out)%dpar .ne. Pout) cycle
              IF ( iso_ops(out)%dTz  .ne. dTz) cycle
              states = states + 1
              Mfi = transition_ME_Tz_EM_Tz(iso_ops(out),O1,iso_ops(in),jbas) 
       
              if (iso_ops(out)%E0 > iso_ops(in)%E0) then 
                 strength_down = Mfi * Mfi / (iso_ops(out)%rank+1.d0)
                 strength_up = Mfi * Mfi / (iso_ops(in)%rank+1.d0)
                 E_in = iso_ops(out)%E0
                 E_out = iso_ops(in)%E0
                 strengths(states) = strength_up
                 Energies(states) = E_in 
                 write(*,'(4(f19.12))') E_in,E_out,Strength_up,Strength_down           
              else
                 strength_down = Mfi * Mfi / (iso_ops(in)%rank+1.d0)
                 strength_up = Mfi * Mfi / (iso_ops(out)%rank+1.d0)
                 E_in = iso_ops(in)%E0
                 E_out = iso_ops(out)%E0
                 strengths(states) = strength_down
                 Energies(states) = E_out 
                 write(*,'(4(f19.12))') E_in,E_out,Strength_down,Strength_up           
              end if


           end do
           instate = instate+1
           write(statlab,'(I1)') instate
           write(dTzlab,'(I2)') dTz
           dTzlab = adjustl(dTzlab)               
           open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
                '_energies_strengths_'//trans%oper//'_'//trans%Jpi1(q)//'_'//statlab//'_'//&
                trans%Jpi2(q)//'_dTz'//trim(dTzlab)//'.dat')
           write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') nint(HS%hospace),&
                HS%eMax,iso_ops(in)%E0,Energies,strengths
           close(31)
           deallocate(Energies,strengths)  

        end do




        
     end if
     print*         
  end do

  ! CALCULATE MOMENTS 
  do q = 1, mom%num

     read(mom%Jpi1(q)(1:1),'(I1)') Jin
     dTz = mom%dtz(q)
     if (mom%Jpi1(q)(2:2) == '+' ) then
        Pin = 0
     else
        Pin = 2
     end if
     Jin = 2*Jin
     print*

     print*, '======================================='
     print*, '           E              <'//mom%oper//'>('//mom%Jpi1(q)//')'  
     print*, '======================================='

     if (dTz == 0 ) then 
        ! count number of states
        states = 0.d0 
        do in = 1,size(ladder_ops)
           IF ( ladder_ops(in)%rank .ne. Jin) cycle
           IF ( ladder_ops(in)%dpar .ne. Pin) cycle
           states = states + 1 
        end do

        allocate(moments(states),Energies(states))
        states = 2*states + 1 ! energies, observs, and gs
        write(flts,'(I2)') states       
        states = 0


        do In = 1, size(ladder_ops) 

           IF ( ladder_ops(in)%rank .ne. Jin) cycle
           IF ( ladder_ops(in)%dpar .ne. Pin) cycle
           
           Mfi = transition_ME(ladder_ops(in),O1,ladder_ops(in),jbas)  
           moment = Mfi / sqrt(Jin+1.d0)*dcgi(Jin,Jin,O1%rank,0,Jin,Jin)
           E_in = ladder_ops(in)%E0
           write(*,'(2(f19.12))') E_in,moment

           states = states + 1
           moments(states) = moment
           energies(states)=E_in
        end do
        dTZlab = '0 '
        open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
             '_energies_moments_'//mom%oper//'_'//mom%Jpi1(q)//'_dTz'//trim(dTzlab)//'.dat')
        write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') nint(HS%hospace),HS%eMax,HS%E0,Energies,moments
        close(31)
        deallocate(Energies,moments)  
        print* 
     else
        ! count number of states
        states = 0.d0 
        do in = 1,size(iso_ops)
           IF ( iso_ops(in)%rank .ne. Jin) cycle
           IF ( iso_ops(in)%dpar .ne. Pin) cycle
           IF ( iso_ops(in)%dTz .ne. dTz) cycle
           states = states + 1 
        end do

        allocate(moments(states),Energies(states))
        states = 2*states + 1 ! energies, observs, and gs
        write(flts,'(I2)') states       

        states = 0

        
        do In = 1, size(iso_ops) 

           IF ( iso_ops(in)%rank .ne. Jin) cycle
           IF ( iso_ops(in)%dpar .ne. Pin) cycle
           IF ( iso_ops(in)%dTz .ne. dTz) cycle
           
           Mfi = transition_ME_Tz_EM_Tz(iso_ops(in),O1,iso_ops(in),jbas)  

           moment = Mfi / sqrt(Jin+1.d0)*dcgi(Jin,Jin,O1%rank,0,Jin,Jin)
           E_in = iso_ops(in)%E0
           write(*,'(2(f19.12))') E_in,moment
           states = states + 1
           moments(states) = moment
           energies(states)=E_in
        end do
        write(dTzlab,'(I2)') dTz
        dTzlab = adjustl(dTzlab)               
           
        open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
             '_energies_moments_'//mom%oper//'_'//mom%Jpi1(q)//'_dTz'//trim(dTzlab)//'.dat')
        write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') nint(HS%hospace),HS%eMax,HS%E0,Energies,moments
        close(31)
        deallocate(Energies,moments)  
        print* 

     end if
  end do


  if (allocated(Hcm%tblck)) then 
     
     ! CALCULATE Hcm 
     print*
     print*, '======================================='
     print*, '           E              <Hcm>'  
     print*, '======================================='

     do In = 1, size(ladder_ops) 
        Mfi = transition_ME(ladder_ops(in),Hcm,ladder_ops(in),jbas)  
        moment = Mfi/sqrt(ladder_ops(in)%rank+1.d0) + Hcm%E0
        E_in = ladder_ops(in)%E0
        write(*,'(2(f19.12))') E_in,moment
        states = states + 1
     end do
  end if

end subroutine EOM_observables
!=========================================================================================================
!=========================================================================================================
subroutine EOM_beta_observables( ladder_ops, iso_ops, O1,HS, Hcm, trans, mom, eom_states , jbas)
  implicit none

  type(sq_op),dimension(:) :: ladder_ops
  type(iso_ladder),dimension(:) :: iso_ops
  type(sq_op) :: HS,Hcm
  type(iso_operator) :: O1
  type(spd) :: jbas
  type(eom_mgr) :: eom_states
  type(obsv_mgr) :: trans,mom
  integer :: q,Jin,Jout,Pin,Pout,in,out,states,instate,DTz
  logical :: to_ground,from_ground
  real(8) :: Mfi,strength_down,strength_up,moment,dcgi,dcgi00
  real(8) :: E_in, E_out
  real(8),allocatable,dimension(:) :: STRENGTHS,MOMENTS,ENERGIES
  character(2) :: flts,dTzlab 
  character(1) :: statlab

  Mfi = dcgi00()
  ! CALCULATE TRANSITIONS 
  do q = 1, trans%num
     to_ground = .false.
     from_ground = .false. 
     dTz = trans%dtz(q) ! dTz is the dTz of the INCOMING STATE     

     IF ( trans%Jpi1(q) == 'GS' ) then
        Jin = 0 
        from_ground = .true.
     else
        read(trans%Jpi1(q)(1:1),'(I1)') Jin
        if (trans%Jpi1(q)(2:2) == '+' ) then
           Pin = 0
        else
           Pin = 2
        end if
     end if
     IF ( trans%Jpi2(q) == 'GS' ) then
        Jout = 0
        to_ground = .true.
     else
        read(trans%Jpi2(q)(1:1),'(I1)') Jout  
        if (trans%Jpi2(q)(2:2) == '+' ) then
           Pout = 0
        else
           Pout = 2
        end if
     end if

     if (to_Ground .and. from_ground) STOP 'No GT transitions between the same state, silly'

     print*
     print*, '============================================================================='
     print*, '        E_in                E_out         B('&
          //trans%oper//';'//trans%Jpi1(q)//' -> '//trans%Jpi2(q)//&
          ')      B('//trans%oper//';'//trans%Jpi2(q)//' -> '//trans%Jpi1(q)//')' 
     print*, '============================================================================='

     Jin = 2* Jin
     Jout = 2* Jout

     If (to_Ground) then
        ! count number of states
        states = 0
    
        do in = 1,size(iso_ops)          
           IF ( iso_ops(in)%dTz .ne. -1*O1%dTz ) cycle ! need to change dTz to zero 
           IF ( iso_ops(in)%rank .ne. Jin) cycle
           IF ( iso_ops(in)%dpar .ne. Pin) cycle
           states = states + 1 
        end do
        
        allocate(Strengths(states),Energies(states))
        states = 2*states + 1 ! energies, observs, and gs
        write(flts,'(I2)') states       
        states = 0
        do in = 1, size(iso_ops)
           IF ( iso_ops(in)%dTz .ne. -1*O1%dTz ) cycle ! need to change dTz to zero 
           IF ( iso_ops(in)%rank .ne. Jin) cycle
           IF ( iso_ops(in)%dpar .ne. Pin) cycle
           
           Mfi = iso_transition_to_ground_ME(O1,iso_ops(in),jbas)
           
           strength_down = Mfi * Mfi /(iso_ops(in)%rank+1.d0)
           strength_up = Mfi * Mfi
           
           E_in = iso_ops(in)%E0
           E_out = 0.d0
           states = states + 1
           strengths(states) = strength_down
           Energies(states) = iso_ops(in)%E0
           
           write(*,'(4(f19.12))') E_in,E_out,Strength_down,Strength_up           
        end do
        
        open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
             '_energies_strengths_'//trans%oper//'_'//trans%Jpi1(q)//'_'//trans%Jpi2(q)//'.dat')
        write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') nint(HS%hospace),HS%eMax,HS%E0,Energies,strengths
        close(31)
        
        deallocate(Energies,strengths)  
     
     else if ( from_ground) then 
        ! count number of states
        states = 0
            
        do out = 1,size(iso_ops)
           IF ( iso_ops(out)%dTz .ne. O1%dTz ) cycle ! need to change dTz from zero 
           IF ( iso_ops(out)%rank .ne. Jout) cycle
           IF ( iso_ops(out)%dpar .ne. Pout) cycle
           states = states + 1 
        end do
        
        allocate(Strengths(states),Energies(states))
        states = 2*states + 1 ! energies, observs, and gs
        write(flts,'(I2)') states       
        states = 0
        do out = 1, size(iso_ops)
           
           IF ( iso_ops(out)%dTz .ne. O1%dTz ) cycle ! need to change dTz from zero 
           IF ( iso_ops(out)%rank .ne. Jout) cycle
           IF ( iso_ops(out)%dpar .ne. Pout) cycle
           
           Mfi = iso_transition_to_ground_ME(O1,iso_ops(out),jbas)
           
           strength_down = Mfi * Mfi  
           strength_up = Mfi * Mfi/(iso_ops(out)%rank+1.d0) 
           
           E_in =  0.d0 
           E_out = iso_ops(out)%E0

           states = states + 1

           strengths(states) = strength_down
           Energies(states) = iso_ops(out)%E0
           
           write(*,'(4(f19.12))') E_in,E_out,Strength_down,Strength_up           

        end do
        
        open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
             '_energies_strengths_'//trans%oper//'_'//trans%Jpi1(q)//'_'//trans%Jpi2(q)//'.dat')
        write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') nint(HS%hospace),HS%eMax,HS%E0,Energies,strengths
        close(31)
        
        deallocate(Energies,strengths)  

     else

        stop 'transition not currently implemented, must include reference state' 
        if (dTZ == 0) then 
           instate = 0        
           do In = 1, size(ladder_ops)
              
              IF ( ladder_ops(in)%rank .ne. Jin) cycle
              IF ( ladder_ops(in)%dpar .ne. Pin) cycle
              
              states = 0.d0 
              do out = 1,size(iso_ops)

                 IF ( iso_ops(out)%dTz .ne. O1%dTz ) cycle ! change dTz from zero
                 IF ( iso_ops(out)%rank .ne. Jout) cycle
                 IF ( iso_ops(out)%dpar .ne. Pout) cycle
                 states = states + 1 
              end do
              
           allocate(Strengths(states),Energies(states))
           states = 2*states + 1 ! energies, observs, and gs
           write(flts,'(I2)') states       
           states = 0

              do out = 1, size(iso_ops)

                 IF ( iso_ops(out)%dTz .ne. O1%dTz ) cycle ! change dTz from zero
                 IF ( iso_ops(out)%rank .ne. Jout) cycle
                 IF ( iso_ops(out)%dpar .ne. Pout) cycle

                 states = states + 1
 !                Mfi = iso_transition_ME(iso_ops(out),O1,ladder_ops(in),jbas) 

                 if (iso_ops(out)%E0 > ladder_ops(in)%E0) then 
                    strength_down = Mfi * Mfi / (iso_ops(out)%rank+1.d0)
                    strength_up = Mfi * Mfi / (ladder_ops(in)%rank+1.d0)
                    E_in = iso_ops(out)%E0
                    E_out = ladder_ops(in)%E0
                    strengths(states) = strength_up
                    Energies(states) = E_in 
                    write(*,'(4(f19.12))') E_in,E_out,Strength_up,Strength_down           
                 else
                    strength_down = Mfi * Mfi / (ladder_ops(in)%rank+1.d0)
                    strength_up = Mfi * Mfi / (iso_ops(out)%rank+1.d0)
                    E_in = ladder_ops(in)%E0
                    E_out = iso_ops(out)%E0
                    strengths(states) = strength_down
                    Energies(states) = E_out 
                    write(*,'(4(f19.12))') E_in,E_out,Strength_down,Strength_up           
                 end if


              end do
              instate = instate+1
              write(statlab,'(I1)') instate
              write(dTzlab,'(I2)') dTz
              dTzlab = adjustl(dTzlab) 
              open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
                   '_energies_strengths_'//trans%oper//'_'//trans%Jpi1(q)//'_'//statlab//'_'&
                   //trans%Jpi2(q)//'_dTz'//trim(dTzlab)//'.dat')
              write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') &
                   nint(HS%hospace),HS%eMax,ladder_ops(in)%E0,Energies,strengths
              close(31)
              deallocate(Energies,strengths)  

           end do


        else

           if ( dTz + O1%dTz  == 0 ) then

              stop 'transition not implemented yet.' 
              ! the outgoing state is ladder_ops
              instate = 0        
              do In = 1, size(iso_ops)

                 IF ( iso_ops(in)%dTz .ne. O1%dTz ) cycle ! change dTz to zero
                 IF ( iso_ops(in)%rank .ne. Jin) cycle
                 IF ( iso_ops(in)%dpar .ne. Pin) cycle

                 states = 0.d0 
                 do out = 1,size(ladder_ops)

                    IF ( ladder_ops(out)%rank .ne. Jout) cycle
                    IF ( ladder_ops(out)%dpar .ne. Pout) cycle
                    states = states + 1 
                 end do

                 allocate(Strengths(states),Energies(states))
                 states = 2*states + 1 ! energies, observs, and gs
                 write(flts,'(I2)') states       
                 states = 0

                 do out = 1, size(ladder_ops)

                    IF ( ladder_ops(out)%rank .ne. Jout) cycle
                    IF ( ladder_ops(out)%dpar .ne. Pout) cycle

                    states = states + 1
!                    Mfi = reverse_iso_transition_ME(ladder_ops(out),O1,iso_ops(in),jbas) 

                    if (ladder_ops(out)%E0 > iso_ops(in)%E0) then 
                       strength_down = Mfi * Mfi / (ladder_ops(out)%rank+1.d0)
                       strength_up = Mfi * Mfi / (iso_ops(in)%rank+1.d0)
                       E_in = ladder_ops(out)%E0
                       E_out = iso_ops(in)%E0
                       strengths(states) = strength_up
                       Energies(states) = E_in 
                       write(*,'(4(f19.12))') E_in,E_out,Strength_up,Strength_down           
                    else
                       strength_down = Mfi * Mfi / (iso_ops(in)%rank+1.d0)
                       strength_up = Mfi * Mfi / (ladder_ops(out)%rank+1.d0)
                       E_in = iso_ops(in)%E0
                       E_out = ladder_ops(out)%E0
                       strengths(states) = strength_down
                       Energies(states) = E_out 
                       write(*,'(4(f19.12))') E_in,E_out,Strength_down,Strength_up           
                    end if


                 end do
                 instate = instate+1
                 write(statlab,'(I1)') instate
                 write(dTzlab,'(I2)') dTz
                 dTzlab = adjustl(dTzlab) 
                 open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
                      '_energies_strengths_'//trans%oper//'_'//trans%Jpi1(q)//'_'//statlab//'_'&
                      //trans%Jpi2(q)//'_dTz'//trim(dTzlab)//'.dat')
                 write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') &
                      nint(HS%hospace),HS%eMax,iso_ops(in)%E0,Energies,strengths
                 close(31)
                 deallocate(Energies,strengths)  

              end do

             
           else ! both states are iso

              stop 'transition not implemented yet.' 
              instate = 0
              do In = 1, size(iso_ops)

                 IF ( iso_ops(in)%rank .ne. Jin) cycle
                 IF ( iso_ops(in)%dpar .ne. Pin) cycle
                 IF ( iso_ops(in)%dTz .ne. dTz ) cycle          

                 states = 0.d0 
                 do out = 1,size(iso_ops)
                    IF ( out == in) cycle
                    IF ( iso_ops(out)%rank .ne. Jout) cycle
                    IF ( iso_ops(out)%dpar .ne. Pout) cycle
                    IF ( iso_ops(out)%dTz .ne. dTz + O1%dTz) cycle
                    states = states + 1 
                 end do

                 allocate(Strengths(states),Energies(states))
                 states = 2*states + 1 ! energies, observs, and gs
                 write(flts,'(I2)') states       
                 states = 0
                 
                 do out = 1, size(iso_ops)
                    
                    if (out==in) cycle ! not a transition
                    IF ( iso_ops(out)%rank .ne. Jout) cycle
                    IF ( iso_ops(out)%dpar .ne. Pout) cycle
                    IF ( iso_ops(out)%dTz  .ne. dTz + O1%dTz) cycle

                    states = states + 1
!                    Mfi = twoiso_transition_ME_Tz_EM_Tz(iso_ops(out),O1,iso_ops(in),jbas) 
                    
                    if (iso_ops(out)%E0 > iso_ops(in)%E0) then 
                       strength_down = Mfi * Mfi / (iso_ops(out)%rank+1.d0)
                       strength_up = Mfi * Mfi / (iso_ops(in)%rank+1.d0)
                       E_in = iso_ops(out)%E0
                       E_out = iso_ops(in)%E0
                       strengths(states) = strength_up
                       Energies(states) = E_in 
                       write(*,'(4(f19.12))') E_in,E_out,Strength_up,Strength_down           
                    else
                       strength_down = Mfi * Mfi / (iso_ops(in)%rank+1.d0)
                       strength_up = Mfi * Mfi / (iso_ops(out)%rank+1.d0)
                       E_in = iso_ops(in)%E0
                       E_out = iso_ops(out)%E0
                       strengths(states) = strength_down
                       Energies(states) = E_out 
                       write(*,'(4(f19.12))') E_in,E_out,Strength_down,Strength_up           
                    end if


                 end do
                 instate = instate+1
                 write(statlab,'(I1)') instate
                 write(dTzlab,'(I2)') dTz
                 dTzlab = adjustl(dTzlab)               
                 open(unit=31,file=trim(OUTPUT_DIR)//trim(adjustl(prefix))//&
                      '_energies_strengths_'//trans%oper//'_'//trans%Jpi1(q)//'_'//statlab//'_'//&
                      trans%Jpi2(q)//'_dTz'//trim(dTzlab)//'.dat')
                 write(31,'(2(I5),'//trim(adjustl(flts))//'(f25.14))') nint(HS%hospace),&
                      HS%eMax,iso_ops(in)%E0,Energies,strengths
                 close(31)
                 deallocate(Energies,strengths)  
                 
              end do
              
           end if
        end if
        print*

     end if
  end do
     
  
end subroutine EOM_beta_observables


!=========================================================
subroutine initialize_rho21_zerorange(rho21,jbas)
  use gaulag
  implicit none 
  
  type(sq_op) :: rho21
  type(spd) :: jbas
  real(8) :: nu,alpha,BB,AA,sm,angular,d6ji,dcgi,pre,dcgi00,t1,t2
  real(8) :: omp_get_wtime,ang1,ang2
  real(8),dimension(20) :: x,w
  real(8),allocatable,dimension(:) :: Lag_a,Lag_b,Lag_c,Lag_d
  integer :: q,i,ord,a,b,c,d,II,JJ,KK,BigL,weirdx,Lmin,Lmax,qxx
  integer :: ja,jb,jc,jd,la,lb,lc,ld,xmin,xmax,antisym
  integer :: npp,nhh,nph,ta,tb,tc,td,na,nb,nc,nd,Jtot
  integer :: jxstart,c1,c2,n1,n2,g_ix,al,bl,cl,dl,LL,s
  logical :: square


  pre = dcgi00()
  nu = rho21%hospace*0.5d0 /hbarc2_over_mc2

  ord = 20
  AA = 0.0
  BB = 2.0
  
  t1 = omp_get_Wtime()
  if (rho21%Aprot+rho21%Aneut==2) then
     q=block_index(2,0,0)

     nhh = rho21%mat(q)%nhh
     npp = rho21%mat(q)%npp
     nph = rho21%mat(q)%nph
     Jtot = rho21%mat(q)%lam(1)

     ! print*, q,rho21%nblocks
     do g_ix = 1,6 
        n1 = size(rho21%mat(q)%gam(g_ix)%X(:,1))
        n2 = size(rho21%mat(q)%gam(g_ix)%X(1,:))
        if ((n1*n2) == 0) cycle 

        ! read in information about which 
        ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        square = sqs(g_ix) 
        jxstart = jst(g_ix) 


        do ii = 1, n1

           a = rho21%mat(q)%qn(c1)%Y(ii,1)
           b = rho21%mat(q)%qn(c1)%Y(ii,2)
           pre = 1.d0       
           if (a==b) pre = 1/sqrt(2.d0)

           na = jbas%nn(a) 
           nb = jbas%nn(b)
           ja = jbas%jj(a)
           jb = jbas%jj(b)
           la = jbas%ll(a)
           lb = jbas%ll(b)
           ta = jbas%itzp(a)
           tb = jbas%itzp(b) 

           allocate(Lag_a(0:na),Lag_b(0:nb) )
           call generate_associated_laguerre_polynomials(na,la,Lag_a)
           call generate_associated_laguerre_polynomials(nb,lb,Lag_b)

           do jj = min(jxstart,ii),n2 

              c = rho21%mat(q)%qn(c2)%Y(jj,1)
              d = rho21%mat(q)%qn(c2)%Y(jj,2)
              if (c==d) pre = pre/sqrt(2.d0) 

              nc = jbas%nn(c) 
              nd = jbas%nn(d)
              jc = jbas%jj(c)
              jd = jbas%jj(d)
              lc = jbas%ll(c)
              ld = jbas%ll(d)
              tc = jbas%itzp(c)
              td = jbas%itzp(d) 

              if (la + lb == 0 ) then
                 if (lc .ne. ld) cycle
              end if

              if (lc + ld == 0 ) then
                 if (la .ne. lb) cycle
              end if

              !the ls should sum to even due to parity conservation           
              alpha = (la+lb+lc+ld)/2+0.5d0 

              ! generate quadrature rule for \int_0^\infty  x^\alpha e^{-2x} f(x)
              call get_rule(ord,alpha,AA,BB,x,w) 

              allocate(Lag_c(0:nc),Lag_d(0:nd) )

              ! get laguerre polynomial coefficients 
              call generate_associated_laguerre_polynomials(nc,lc,Lag_c)
              call generate_associated_laguerre_polynomials(nd,ld,Lag_d)

              sm = 0.d0 
              do kk = 1, ord

                 sm = sm +  asslag(na,la,lag_a,x(kk))*asslag(nb,lb,lag_b,x(kk))*asslag(nc,lc,lag_c,x(kk)) &
                      *asslag(nd,ld,lag_d,x(kk))*w(kk) 

              end do

              deallocate(Lag_c,Lag_d)

              ! sm is now the integral \int_0^\infty  x^\alpha e^{-2x} f(x)

              sm = sm * (nu*0.5)**(1.5) /Pi_const
              ! sm is now the radial part of the matrix element.          

              ! There is also an angular momentum coupled pre-factor
              angular = 0.d0
              ang2=0.d0 
              ang1=0.d0 

              Lmin = max( abs(la-lb) , abs(lc-ld))
              Lmax = min( la+lb , lc+ld)

              xmin = 1!max(  abs(Jtot-1),abs(2*la-jb),abs(2*lc-jd))          
              xmax = 21!min( Jtot+jb , Jtot+jd)           

              s = 1
              al = 2*la
              bl = 2*lb
              cl = 2*lc
              dl = 2*ld


              if (kron_del(ta,tc)*kron_del(tb,td) == 1) then 
                 do BigL = Lmin,Lmax              
                    LL = 2*BigL
                    do weirdx = xmin,xmax,2

                       ang1 = ang1 + (weirdx + 1.d0)*d6ji(ja,jb,Jtot,weirdx,s,al) * &
                            d6ji(al,bl,LL,s,weirdx,jb) * dcgi(al,0,bl,0,LL,0)*&
                            dcgi(cl,0,dl,0,LL,0)*d6ji(jc,jd,Jtot,weirdx,s,cl) *&
                            d6ji(cl,dl,LL,s,weirdx,jd) 

                    end do
                 end do
              end if

              if (kron_del(ta,td)*kron_del(tb,tc) == 1) then 
                 do BigL = Lmin,Lmax              
                    LL = 2*BigL
                    do weirdx = xmin,xmax,2

                       ang2 = ang2 + (weirdx + 1.d0)*d6ji(ja,jb,Jtot,weirdx,s,al) * &
                            d6ji(al,bl,LL,s,weirdx,jb) * dcgi(al,0,bl,0,LL,0)*&
                            dcgi(cl,0,dl,0,LL,0)*d6ji(jd,jc,Jtot,weirdx,s,dl) *&
                            d6ji(dl,cl,LL,s,weirdx,jc)   ! the permutation comes from the spins, not the
                       ! mls, so, the clebsch with c and d does not change. 

                    end do
                 end do
              end if

              angular = angular + ang1 - ang2 * (-1)**((jc+jd+cl+dl+Jtot+LL)/2) !!! antisymmetry condition              
              angular = angular * (-1) ** ( (ja+jb+jc+jd)/2)    !!! this may not be needed if I got the coupling order wrong.


              angular = angular * sqrt( (ja+1.d0)*(jb+1.d0)*(jc+1.d0)*(jd+1.d0) &
                   *(2*la+1.d0)*(2*lb+1.d0)*(2*lc+1.d0)*(2*ld+1.d0)) / 4.d0 / Pi_const                                                
              rho21%mat(q)%gam(g_ix)%X(II,JJ) = angular*sm *pre
              if (square) rho21%mat(q)%gam(g_ix)%X(JJ,II) = angular*sm *pre 

           end do
           deallocate(Lag_a,Lag_b)
        end do
     end do
     
  else
     !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE), SHARED(RHO21,JBAS)  
     do q= 1, rho21%nblocks

        nhh = rho21%mat(q)%nhh
        npp = rho21%mat(q)%npp
        nph = rho21%mat(q)%nph
        Jtot = rho21%mat(q)%lam(1)

        ! print*, q,rho21%nblocks
        do g_ix = 1,6 
           n1 = size(rho21%mat(q)%gam(g_ix)%X(:,1))
           n2 = size(rho21%mat(q)%gam(g_ix)%X(1,:))
           if ((n1*n2) == 0) cycle 

           ! read in information about which 
           ! array we are using from public arrays
           c1 = sea1(g_ix) 
           c2 = sea2(g_ix) 
           square = sqs(g_ix) 
           jxstart = jst(g_ix) 


           do ii = 1, n1

              a = rho21%mat(q)%qn(c1)%Y(ii,1)
              b = rho21%mat(q)%qn(c1)%Y(ii,2)
              pre = 1.d0       
              if (a==b) pre = 1/sqrt(2.d0)

              na = jbas%nn(a) 
              nb = jbas%nn(b)
              ja = jbas%jj(a)
              jb = jbas%jj(b)
              la = jbas%ll(a)
              lb = jbas%ll(b)
              ta = jbas%itzp(a)
              tb = jbas%itzp(b) 

              allocate(Lag_a(0:na),Lag_b(0:nb) )
              call generate_associated_laguerre_polynomials(na,la,Lag_a)
              call generate_associated_laguerre_polynomials(nb,lb,Lag_b)

              do jj = min(jxstart,ii),n2 

                 c = rho21%mat(q)%qn(c2)%Y(jj,1)
                 d = rho21%mat(q)%qn(c2)%Y(jj,2)
                 if (c==d) pre = pre/sqrt(2.d0) 

                 nc = jbas%nn(c) 
                 nd = jbas%nn(d)
                 jc = jbas%jj(c)
                 jd = jbas%jj(d)
                 lc = jbas%ll(c)
                 ld = jbas%ll(d)
                 tc = jbas%itzp(c)
                 td = jbas%itzp(d) 

                 if (la + lb == 0 ) then
                    if (lc .ne. ld) cycle
                 end if

                 if (lc + ld == 0 ) then
                    if (la .ne. lb) cycle
                 end if

                 !the ls should sum to even due to parity conservation           
                 alpha = (la+lb+lc+ld)/2+0.5d0 

                 ! generate quadrature rule for \int_0^\infty  x^\alpha e^{-2x} f(x)
                 call get_rule(ord,alpha,AA,BB,x,w) 

                 allocate(Lag_c(0:nc),Lag_d(0:nd) )

                 ! get laguerre polynomial coefficients 
                 call generate_associated_laguerre_polynomials(nc,lc,Lag_c)
                 call generate_associated_laguerre_polynomials(nd,ld,Lag_d)

                 sm = 0.d0 
                 do kk = 1, ord

                    sm = sm +  asslag(na,la,lag_a,x(kk))*asslag(nb,lb,lag_b,x(kk))*asslag(nc,lc,lag_c,x(kk)) &
                         *asslag(nd,ld,lag_d,x(kk))*w(kk) 

                 end do

                 deallocate(Lag_c,Lag_d)

                 ! sm is now the integral \int_0^\infty  x^\alpha e^{-2x} f(x)

                 sm = sm * (nu*0.5)**(1.5) /Pi_const
                 ! sm is now the radial part of the matrix element.          

                 ! There is also an angular momentum coupled pre-factor
                 angular = 0.d0
                 ang2=0.d0 
                 ang1=0.d0 

                 Lmin = max( abs(la-lb) , abs(lc-ld))
                 Lmax = min( la+lb , lc+ld)

                 xmin = 1!max(  abs(Jtot-1),abs(2*la-jb),abs(2*lc-jd))          
                 xmax = 21!min( Jtot+jb , Jtot+jd)           

                 s = 1
                 al = 2*la
                 bl = 2*lb
                 cl = 2*lc
                 dl = 2*ld


                 if (kron_del(ta,tc)*kron_del(tb,td) == 1) then 
                    do BigL = Lmin,Lmax              
                       LL = 2*BigL
                       do weirdx = xmin,xmax,2

                          ang1 = ang1 + (weirdx + 1.d0)*d6ji(ja,jb,Jtot,weirdx,s,al) * &
                               d6ji(al,bl,LL,s,weirdx,jb) * dcgi(al,0,bl,0,LL,0)*&
                               dcgi(cl,0,dl,0,LL,0)*d6ji(jc,jd,Jtot,weirdx,s,cl) *&
                               d6ji(cl,dl,LL,s,weirdx,jd) 

                       end do
                    end do
                 end if

                 if (kron_del(ta,td)*kron_del(tb,tc) == 1) then 
                    do BigL = Lmin,Lmax              
                       LL = 2*BigL
                       do weirdx = xmin,xmax,2

                          ang2 = ang2 + (weirdx + 1.d0)*d6ji(ja,jb,Jtot,weirdx,s,al) * &
                               d6ji(al,bl,LL,s,weirdx,jb) * dcgi(al,0,bl,0,LL,0)*&
                               dcgi(cl,0,dl,0,LL,0)*d6ji(jd,jc,Jtot,weirdx,s,dl) *&
                               d6ji(dl,cl,LL,s,weirdx,jc)   ! the permutation comes from the spins, not the
                          ! mls, so, the clebsch with c and d does not change. 

                       end do
                    end do
                 end if

                 angular = angular + ang1 - ang2 * (-1)**((jc+jd+cl+dl+Jtot+LL)/2) !!! antisymmetry condition              
                 angular = angular * (-1) ** ( (ja+jb+jc+jd)/2)    !!! this may not be needed if I got the coupling order wrong.


                 angular = angular * sqrt( (ja+1.d0)*(jb+1.d0)*(jc+1.d0)*(jd+1.d0) &
                      *(2*la+1.d0)*(2*lb+1.d0)*(2*lc+1.d0)*(2*ld+1.d0)) / 4.d0 / Pi_const                                                
                 rho21%mat(q)%gam(g_ix)%X(II,JJ) = angular*sm *pre
                 if (square) rho21%mat(q)%gam(g_ix)%X(JJ,II) = angular*sm *pre 

              end do
              deallocate(Lag_a,Lag_b)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  t2 = omp_get_Wtime()
  print*, 'time:', t2-t1
end subroutine initialize_rho21_zerorange

real(8) function asslag(n,l,lag,x)
  implicit none

  integer :: n,l ,ii
  real(8) :: lag(0:n),x,sm

  sm = 0.d0 
  do ii = 0,n
     sm = sm + lag(ii) * x**ii
  end do

  asslag= sm
end function asslag

 
subroutine generate_associated_laguerre_polynomials(nMax,l,LAGnl)
  implicit none

  integer :: nMax,l,n, II 
  real(8) :: LAG(0:nMax,0:nMax) !!1 coefficients
  real(8) :: norm,LAGnl(0:nMax)

  
  LAG = 0.d0

  LAG(0,0) = 1.d0   ! L_0^(l+1/2) == 1 always

  if (nMax > 0) then 
     LAG(1,0) = 1.5d0 + l  !!! L_1^(l+1/2) = 1.5+l-x 
     LAG(1,1) = -1.d0   
     
     LAG(1,:) = LAG(1,:) 
  end if
  
!!! now use recurrence relation
  do n = 2, nMax !!! loop over Laguerre rank    

     LAG(n,0) = ( (2*n-0.5 + l)* LAG(n-1,0)  - (n-0.5+l) * LAG(n-2,0) ) / n
     do II = 1 , n !!! loop over powers of x        
        LAG(n,II) = ( (2*n-0.5 + l)* LAG(n-1,II)  - (n-0.5+l) * LAG(n-2,II) - LAG(n-1,II-1) )  / n
     end do
  end do

  do n = 0,nMax
     norm = sqrt( 2.d0**( n+ l +3) *fac_over_doublefac(n,2*n+2*l+1) )
     LAG(n,:) = LAG(n,:)*norm
  end do

  LAGnl = LAG(nMax,:)
  
end subroutine generate_associated_laguerre_polynomials

end module
  
  
  
  

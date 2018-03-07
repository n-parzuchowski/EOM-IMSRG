module tensor_products
  use isospin_operators
  use cross_coupled
  ! tensor-tensor product functions 
  !! AA is a HERMITIAN OPERATOR
  !! BB is a LADDER OPERATOR
  ! THE TENSOR MUST BE THE SECOND ARGUMENT
  
contains
!===================================================================
!===================================================================
subroutine tensor_product_222_pp_hh(AA,BB,CC,jbas) 
  !VERIFIED 
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  AA,BB,CC
  integer :: q,q1,q2,J1,J2,Tz,Par,phase,rank_A,rank_B,rank_c,i,Cpar
  integer :: np1,nb1,nh1,np2,nb2,nh2,pm,P_a,P_b,J3min,J3max,J3
  real(8) :: bet_off,al_off,d6ji
  real(8),allocatable,dimension(:,:) :: W 
  
  rank_A = AA%rank
  rank_B = BB%rank
  rank_C = CC%rank 
!construct temporary matrices
  do q = 1, CC%nblocks
     
     J1 = CC%tblck(q)%Jpair(1) 
     J2 = CC%tblck(q)%Jpair(2)

     phase = CC%tblck(q)%lam(1)
     par = CC%tblck(q)%lam(2) 
     Tz = CC%tblck(q)%lam(3)
             
     J3Min = max(abs(rank_a-J1),abs(rank_b-J2))
     J3Max = min(rank_a + J1,rank_b + J2)
     np1 = CC%tblck(q)%npp1
     nh2 = CC%tblck(q)%nhh2

     do J3 = J3Min, J3Max, 2
        if ( J3 > 2*jbas%jtotal_max) cycle
        If ( J3 .ge. J1 )  then

           IF (J2 .ge. J3 ) then 
           
              P_a = CC%tblck(q)%lam(2)
              Tz = CC%tblck(q)%lam(3) 
              P_b = mod(P_a + AA%dpar/2,2)
              
              q1 = tensor_block_index(J1,J3,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J3,J2,rank_b,Tz,P_b) 
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------

               
              np2 = AA%tblck(q1)%npp2
                            
              if (np1*nh2 .ne. 0)  then                       
                 if (np2 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0

                   call dgemm('N','N',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np1,&
                        BB%tblck(q2)%tgam(3)%X,np2,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if
                   
                              
           else ! ( J2< J3) 
                                      
              P_a = CC%tblck(q)%lam(2)
              Tz = CC%tblck(q)%lam(3) 
              P_b = mod(P_a + CC%dpar/2,2) 
              
              q1 = tensor_block_index(J1,J3,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J2,J3,rank_b,Tz,P_b) 
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------
              np2 = AA%tblck(q1)%npp2
              nh1 = BB%tblck(q2)%nhh1
              
              if (np1*nh2 .ne. 0)  then                       
                 if (np2*nh1 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('N','T',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np1,&
                         BB%tblck(q2)%tgam(7)%X,nh1,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if
              
              
           end if

        else ! ( J3 < J1 .le. J2 )
           if ( J2 .ge. J1) then

              P_a = mod(CC%tblck(q)%lam(2)+AA%dpar/2,2) 
              Tz = CC%tblck(q)%lam(3) 
              P_b = P_a
              
              q1 = tensor_block_index(J3,J1,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J3,J2,rank_b,Tz,P_b) 
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------

               
              np2 = AA%tblck(q1)%npp1
                            
              if (np1*nh2 .ne. 0)  then                       
                 if (np2 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*AA%tblck(q1)%lam(1)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0

                   call dgemm('T','N',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np2,&
                        BB%tblck(q2)%tgam(3)%X,np2,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if
                   
                   
              
            end if

              
           
        end if

     end do

     J3Min = max(abs(rank_b-J1),abs(rank_a-J2))
     J3Max = min(rank_b + J1,rank_a + J2)

     
     do J3 = J3Min, J3Max, 2
!        if ( J3 > 2*jbas%jtotal_max) cycle        
        If ( J3 .ge. J1 )  then

           IF (J2 .ge. J3 ) then 


              Tz = CC%tblck(q)%lam(3) 

              !----------------------------------------------------------------------------
              !         Zpphh 
              !----------------------------------------------------------------------------

               

              P_b = CC%tblck(q)%lam(2)
              P_a = mod(P_b + BB%dpar/2,2)

              q1 = tensor_block_index(J1,J3,rank_b,Tz,P_b) 
              q2 = tensor_block_index(J3,J2,rank_a,Tz,P_a) 

              nh1 = BB%tblck(q1)%nhh2     
              !w1pphh = -Bpphh.Ahhhh 

              if (np1*nh2 .ne. 0)  then                       
                 if (nh1 .ne. 0) then 
                    al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*&
                         (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('N','N',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(3)%X,np1,&
                         AA%tblck(q2)%tgam(5)%X,nh1,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if

           else ! ( J2< J3) 

              Tz = CC%tblck(q)%lam(3) 

              !----------------------------------------------------------------------------
              !         Zpphh 
              !----------------------------------------------------------------------------

              P_b = CC%tblck(q)%lam(2)
              P_a = mod(P_b + CC%dpar/2,2)

              q1 = tensor_block_index(J1,J3,rank_b,Tz,P_b) 
              q2 = tensor_block_index(J2,J3,rank_a,Tz,P_a) 

              nh1 = BB%tblck(q1)%nhh2     
              !w1pphh = -Bpphh.Ahhhh 

              if (np1*nh2 .ne. 0)  then                       
                 if (nh1 .ne. 0) then 
                    al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*AA%tblck(q2)%lam(1)*&
                         (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('N','T',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(3)%X,np1,&
                         AA%tblck(q2)%tgam(5)%X,nh2,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if


           end if

        else ! ( J3 < J1 .le. J2 )
           if ( J2 .ge. J1) then
              

              Tz = CC%tblck(q)%lam(3) 
              
           !----------------------------------------------------------------------------
           !         Zpphh 
           !----------------------------------------------------------------------------

               
              
              P_b = mod(CC%tblck(q)%lam(2)+BB%dpar/2,2)
              P_a = P_b

              q1 = tensor_block_index(J3,J1,rank_b,Tz,P_b) 
              q2 = tensor_block_index(J3,J2,rank_a,Tz,P_a) 

              nh1 = BB%tblck(q1)%nhh1     
              !w1pphh = -Bpphh.Ahhhh 

              if (np1*nh2 .ne. 0)  then                       
                 if (nh1 .ne. 0) then 
                    al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*&
                         (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('T','N',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(7)%X,nh1,&
                         AA%tblck(q2)%tgam(5)%X,nh1,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if

           end if



        end if

     end do

!     NOW DO THE OPPOSITE SIDE OF THE MATRIX (J1 > J2)

!    just switching the names so that I don't have to change
!    variables in the expression below. 
     J3 = J1
     J1 = J2
     J2 = J3 
     
     if (J1 > jbas%jtotal_max*2) cycle
     
     J3Min = max(abs(rank_a-J1),abs(rank_b-J2))
     J3Max = min(rank_a + J1,rank_b + J2)
     Cpar= mod(CC%tblck(q)%lam(2)+ CC%dpar/2,2)

     np1 = CC%tblck(q)%npp2
     nh2 = CC%tblck(q)%nhh1

     allocate(W(np1,nh2)) 
     W = 0.d0
     do J3 = J3Min, J3Max, 2
        if ( J3 > 2*jbas%jtotal_max) cycle
        If ( J3 .ge. J1 )  then
           IF (J1 .ge. J2 ) then 
                        
              Tz = CC%tblck(q)%lam(3) 
              
              P_a = Cpar
              P_b = mod(P_a + CC%dpar/2,2)

              q1 = tensor_block_index(J1,J3,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J2,J3,rank_b,Tz,P_b) 
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------

               
              np2 = AA%tblck(q1)%npp2
              nh1 = BB%tblck(q2)%nhh1
                            
              if (np1*nh2 .ne. 0)  then                       
                 if (np2*nh1 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    
                    call dgemm('N','T',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np1,&
                         BB%tblck(q2)%tgam(7)%X,nh1,bet_off,W,np1)
                    
                 end if
              end if
              
           end if
        else ! J3 < J1 
           if (J2 .ge. J3)  then 
                                                    
              Tz = CC%tblck(q)%lam(3) 
              
              P_a = mod(Cpar+AA%dpar/2,2)
              P_b = P_a

              q1 = tensor_block_index(J3,J1,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J3,J2,rank_b,Tz,P_b) 
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------
              np2 = AA%tblck(q1)%npp1
              
              if (np1*nh2 .ne. 0)  then                       
                 if (np2 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*AA%tblck(q1)%lam(1)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('T','N',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np2,&
                         BB%tblck(q2)%tgam(3)%X,np2,bet_off,W,np1)
                 end if
              end if
              
              
              
           else ! ( J3 < J1 ,  J3 >= J2  )
                 
              P_a = mod(Cpar+AA%dpar/2,2) 
              P_b = mod(Cpar+CC%dpar/2,2) 
              Tz = CC%tblck(q)%lam(3) 
              
              q1 = tensor_block_index(J3,J1,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J2,J3,rank_b,Tz,P_b) 
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------

               
              np2 = AA%tblck(q1)%npp1
              nh1 = BB%tblck(q2)%nhh1
                            
              if (np1*nh2 .ne. 0)  then                       
                 if (np2*nh1 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*AA%tblck(q1)%lam(1)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0

                   call dgemm('T','T',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np2,&
                        BB%tblck(q2)%tgam(7)%X,nh1,bet_off,W,np1)
                 end if
              end if
                   
                   
              
            end if

              
           
         end if
         
      end do

      J3Min = max(abs(rank_b-J1),abs(rank_a-J2))
      J3Max = min(rank_b + J1,rank_a + J2)
      
      do J3 = J3Min, J3Max, 2
         if ( J3 > 2*jbas%jtotal_max) cycle
         If ( J3 .ge. J1 )  then
            IF (J1 .ge. J2 ) then


               Tz = CC%tblck(q)%lam(3) 

               !----------------------------------------------------------------------------
               !         Zpphh 
               !----------------------------------------------------------------------------

                

               P_b = Cpar
               P_a = mod(P_b + CC%dpar/2,2)

               q1 = tensor_block_index(J1,J3,rank_b,Tz,P_b) 
               q2 = tensor_block_index(J2,J3,rank_a,Tz,P_a) 

               nh1 = BB%tblck(q1)%nhh2                  
               !w1pphh = -Bpphh.Ahhhh 

               if (np1*nh2 .ne. 0)  then                       
                  if (nh1 .ne. 0) then
                     al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*AA%tblck(q2)%lam(1)*&
                          (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                     bet_off = 1.d0
                     call dgemm('N','T',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(3)%X,np1,&
                          AA%tblck(q2)%tgam(5)%X,nh2,bet_off,W,np1)
                  end if
               end if
            end if
         else ! ( J2< J3) 
            if (J2 .ge. J3)  then 
               Tz = CC%tblck(q)%lam(3) 

               q1 = tensor_block_index(J3,J1,rank_a,Tz,P_a) 
               q2 = tensor_block_index(J3,J1,rank_b,Tz,P_b) 

               !----------------------------------------------------------------------------
               !         Zpphh 
               !----------------------------------------------------------------------------

               P_b = mod(Cpar + BB%dpar/2,2)
               P_a = P_b

               q1 = tensor_block_index(J3,J1,rank_b,Tz,P_b) 
               q2 = tensor_block_index(J3,J2,rank_a,Tz,P_a) 

               nh1 = BB%tblck(q1)%nhh1     
               !w1pphh = -Bpphh.Ahhhh 

               if (np1*nh2 .ne. 0)  then                       
                  if (nh1 .ne. 0) then 
                     al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)* &
                          (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                     bet_off = 1.d0
                     call dgemm('T','N',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(7)%X,nh1,&
                          AA%tblck(q2)%tgam(5)%X,nh1,bet_off,W,np1)
                  end if
               end if
             
            else ! (J2 <  J3 < J1)
               


               Tz = CC%tblck(q)%lam(3) 

               !----------------------------------------------------------------------------
               !         Zpphh 
               !----------------------------------------------------------------------------                

               P_b = mod(Cpar+BB%dpar/2,2)
               P_a = mod(Cpar+CC%dpar/2,2)

               q1 = tensor_block_index(J3,J1,rank_b,Tz,P_b) 
               q2 = tensor_block_index(J2,J3,rank_a,Tz,P_a) 

               nh1 = BB%tblck(q1)%nhh1     
               !w1pphh = -Bpphh.Ahhhh 

               if (np1*nh2 .ne. 0)  then                       
                  if (nh1 .ne. 0) then 
                     al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*AA%tblck(q2)%lam(1)*&
                          (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                     bet_off = 1.d0
                     call dgemm('T','T',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(7)%X,nh1,&
                          AA%tblck(q2)%tgam(5)%X,nh2,bet_off,W,np1)
                  end if
               end if

            end if



         end if

      end do

      CC%tblck(q)%tgam(7)%X = CC%tblck(q)%tgam(7)%X + Transpose(W) 

     deallocate(W) 
  end do
  
end subroutine tensor_product_222_pp_hh
!=================================================================
!=================================================================
 subroutine tensor_product_222_ph(AA,BB,CC,jbas) 
   ! VERIFIED ph channel 2body tensor_product. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(sq_op) :: AA,BB,CC
   real(8),allocatable,dimension(:,:) :: PANDYA_A,PANDYA_B,PANDYA_AB
   type(cc_mat) :: LCC
   integer :: nh,np,nb1,nb2,q,IX,JX,i,r1,r2,Tz,PAR_J3,JTM,q1,q2,J3,J4,rank_a,a,p1,p2,h1,h2
   integer :: ji,jh1,jh2,ti,th1,th2,lh1,lh2,li,n1,n2,c1,c2,jxstart,J4min,J4max,ja,jp1,jp2
   integer :: J1,J2, Jtot,Ntot,qx,J3min,J3max,ril,rjk,rli,rkj,g_ix,thread,total_threads
   integer :: phase1,phase2,phase3,rik,rki,rjl,rlj,PAR2,J1min,J2min,J1max,J2max,n_J3,n_J4,n_J5
   integer :: phase_34,phase_abcd,phase_ac,phase_bc,j5,PAR_J4,PAR_J5,rank_b,rank_c
   integer :: phase_bd,phase_ad,nj_perm,full_int_phase,J5min,J5max,tp1,tp2,lp1,lp2
   integer,allocatable,dimension(:,:) :: qn_J3,qn_J4,qn_J5
   real(8) :: sm ,pre,pre2,omp_get_wtime ,t1,t2,coef9,factor,sm_ex,nj1,nj2  
   real(8) :: prefac_12,prefac_1,nj,Xelem,Yelem,V,al_off,d6ji
   logical :: square
   
   Ntot = CC%Nsp
   JTM = jbas%Jtotal_max*2
   total_threads = size(CC%direct_omp) - 1
   rank_a = AA%rank
   rank_b = BB%rank
   rank_c = CC%rank 
   ! construct intermediate matrices

   ! sum over all cross coupled channels for both matrices 
   do J3 = 0, JTM  ,2
      J4min = abs(rank_a -J3)
      J4max = min(rank_a + J3,JTM)

      do Tz = 0,2,2 ! this is the cross coupled TZ (no negative)
         
         do PAR_J3 = 0,1

            n_J3 = count_configs( J3 ,Tz, PAR_J3 , jbas, .true. ) !boolean true: ph, false: hp configs
            if (n_J3 == 0) cycle
            allocate(qn_J3(n_J3,2)) 
            n_J3 = count_configs( J3 ,Tz, PAR_J3 , jbas, .true.,qn_J3) ! this fills the cc- basis descriptor

            PAR_J4 = mod(PAR_J3 + AA%dpar/2,2)
            PAR_J5 = mod(PAR_J4 + BB%dpar/2,2)
            
            do J4 = J4min,J4max ,2

               J5min = max(abs(rank_b-J4),abs(rank_c-J3)) 
               J5max = min(rank_b+J4,rank_c+J3,JTM)

               n_J4 = count_configs( J4 ,Tz, PAR_J4 , jbas, .true. )
               if (n_J4 == 0) cycle
               allocate(qn_J4(n_J4,2)) 
               n_J4 = count_configs( J4 ,Tz, PAR_J4 , jbas, .true.,qn_J4)

               allocate( PANDYA_A( n_J3, n_J4) )
               PANDYA_A = 0.d0 
               call fill_generalized_pandya_matrix(J3,J4,PANDYA_A,qn_J3,qn_J4,AA,jbas)
               
               do J5 = J5min,J5max,2 
                  
                  n_J5 = count_configs( J5 ,Tz, PAR_J5 , jbas, .false. )
                  if (n_J5 == 0) cycle
                  allocate(qn_J5(n_J5,2)) 
                  n_J5 = count_configs( J5 ,Tz, PAR_J5 , jbas, .false.,qn_J5)

                  allocate(PANDYA_B(n_J4, n_J5),PANDYA_AB(n_J3,n_J5)) 
                  PANDYA_B = 0.d0 
                  PANDYA_AB = 0.d0 
                  call fill_generalized_pandya_matrix(J4,J5,PANDYA_B,qn_J4,qn_J5,BB,jbas)                  

                  al_off = d6ji(rank_a,rank_b,rank_c,J5,J3,J4) * (-1)** ((rank_C + J3)/2) * &
                       sqrt((J3+1.d0)*(J5+1.d0)*(rank_c+1.d0)) 
                  call dgemm('N','N',n_J3,n_J5,n_J4,al_off,PANDYA_A,n_J3,PANDYA_B,n_J4,bet,PANDYA_AB,n_J3) 


                  
                  
                  do JX = 1,n_J5

                     ! GET KET 
                     h2 = qn_J5(JX,1)
                     p2 = qn_J5(JX,2)
                     
                     jp2 = jbas%jj(p2)
                     lp2 = jbas%ll(p2)
                     tp2 = jbas%itzp(p2)
                     jh2 = jbas%jj(h2)
                     lh2 = jbas%ll(h2)
                     th2 = jbas%itzp(h2)

                     do IX = 1, n_J3 

                        ! GET BRA
                        p1 = qn_J3(IX,1)
                        h1 = qn_J3(IX,2)

                        jp1 = jbas%jj(p1)
                        jh1 = jbas%jj(h1)
                        lp1 = jbas%ll(p1)
                        tp1 = jbas%itzp(p1)
                        lh1 = jbas%ll(h1)
                        th1 = jbas%itzp(h1)

                        if ( (tp1 +tp2) .ne. (th1+th2) ) cycle
                        if ( mod(lp1 +lp2,2) .ne. mod(lh1+lh2+CC%dpar/2,2) ) cycle
                        
                        ! CALCULATE X CONTRIBUTIONS

                        J1min = abs(jp1-jp2)
                        J1max = jp2+jp1

                        ! these are the results of the Matmuls 
                        Xelem = PANDYA_AB(IX,JX)

                        do J1 = J1min,J1max,2
                           if ((p1==p2).and.(mod(J1/2,2)==1)) cycle                           
                           J2min = max(abs(jh1-jh2),abs(rank_c-J1))
                           J2max = min(jh1+jh2 ,rank_c+J1) 
                           prefac_1 = sqrt(J1+1.d0)

                           do J2 = J2min,J2max,2
                              if ((h1==h2).and.(mod(J2/2,2)==1)) cycle                           
                              nj = ninej(CC%xindx,jp1,jh1,J3,jp2,jh2,J5,J1,J2,rank_c)
                              prefac_12 = prefac_1 *sqrt(J2+1.d0)

                              V = prefac_12*nj*(-1)**((jp2-jh2)/2)*Xelem
                              
                              if ( h2 .ge. h1 ) then
                                 if ( p2 .ge. p1 ) then

                                    
                                    if (J1.le.J2) then
                                       call add_elem_to_tensor(V,p1,p2,h1,h2,J1,J2,CC,jbas)                                         
                                    end if
                                    if (J1.ge.J2) then
                                       call add_elem_to_tensor(V,h1,h2,p1,p2,J2,J1,CC,jbas)                                         
                                    end if

                                 else


                                    V = prefac_12*nj*(-1)**((jp1-jh2)/2)*Xelem*(-1)**(J1/2)
                                    if (J1.le.J2) then
                                       call add_elem_to_tensor(V,p2,p1,h1,h2,J1,J2,CC,jbas)                                         
                                    end if
                                    if (J1.ge.J2) then

                                       call add_elem_to_tensor(V,h1,h2,p2,p1,J2,J1,CC,jbas)                                         
                                    end if

                                 end if
                              else
                                 if (p1 > p2) then

                                    V = prefac_12*nj*(-1)**((jp1-jh1)/2)*Xelem*(-1)**((J1+J2)/2)
                                    if (J1.le.J2) then

                                       call add_elem_to_tensor(V,p2,p1,h2,h1,J1,J2,CC,jbas)                                         
                                    end if
                                    if (J1.ge.J2) then

                                       call add_elem_to_tensor(V,h2,h1,p2,p1,J2,J1,CC,jbas)                                       
                                    end if
                                 else

                                    V = prefac_12*nj*(-1)**((jp2-jh1)/2)*Xelem*(-1)**((J2)/2)
                                    if (J1.le.J2) then
                                       call add_elem_to_tensor(V,p1,p2,h2,h1,J1,J2,CC,jbas)                                         
                                    end if
                                    if (J1.ge.J2) then
                                       call add_elem_to_tensor(V,h2,h1,p1,p2,J2,J1,CC,jbas)                                       
                                    end if
                                 end if
                              end if
                           end do
                        end do
                     end do
                  end do

                  deallocate(PANDYA_B,qn_J5,PANDYA_AB) 
               end do
               deallocate(PANDYA_A,qn_J4)
            end do
            deallocate(qn_J3)
         end do
      end do
   end do
 end subroutine tensor_product_222_ph
 !===================================================================================
 !===================================================================================
 real(8) function count_configs(J1,Tz,PAR,jbas,ph,qn)
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
                  if ( (jbas%con(i) == 1 ).or. (jbas%con(j) == 0)) cycle
               else
                  if ( (jbas%con(i) == 0 ).or. (jbas%con(j) == 1)) cycle
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

   count_configs = r1 
 end function count_configs


 subroutine fill_generalized_pandya_matrix(J1,J2,MAT,qn1,qn2,OP,jbas)
   ! CALCULATES THE CROSS GENERALIZED PANDYA MATRIX ELEMENTS OF
   ! OP FOR A GIVEN CHANNEL
   implicit none

   type(spd) :: jbas
   type(sq_op) :: OP 
   integer,dimension(:,:) :: qn1,qn2
   real(8),dimension(:,:) :: MAT
   integer :: a,b,c,d,J1,J2,N1,N2,II,JJ

   N1  = size(MAT(:,1))
   N2  = size(MAT(1,:))   

   do JJ = 1,N2
      do II = 1,N1

         a=qn1(II,1)
         b=qn1(II,2)
         c=qn2(JJ,1)
         d=qn2(JJ,2)

         MAT(II,JJ) = Vgenpandya(a,b,c,d,J1,J2,Op,jbas)

      end do
   end do

 end subroutine fill_generalized_pandya_matrix
!===================================================================
!===================================================================
subroutine dTz_tensor_product_222_pp_hh(AA,BB,CC,jbas) 
  !VERIFIED 
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  AA
  type(iso_ladder) :: BB,CC
  integer :: q,q1,q2,J1,J2,Tz,Par,phase,rank_A,rank_B,rank_c,i,Cpar
  integer :: np1,nb1,nh1,np2,nb2,nh2,pm,P_a,P_b,J3min,J3max,J3,Tz1,Tz2
  real(8) :: bet_off,al_off,d6ji
  real(8),allocatable,dimension(:,:) :: W 
  
  rank_A = AA%rank
  rank_B = BB%rank
  rank_C = CC%rank 
!construct temporary matrices
  do q = 1, CC%nblocks
     
     J1 = CC%tblck(q)%Jpair(1) 
     J2 = CC%tblck(q)%Jpair(2)

     phase = CC%tblck(q)%lam(1)
     par = CC%tblck(q)%lam(2) 
     Tz1 = CC%tblck(q)%lam(3)
     Tz2 = CC%tblck(q)%lam(4)

     J3Min = max(abs(rank_a-J1),abs(rank_b-J2))
     J3Max = min(rank_a + J1,rank_b + J2)
     np1 = CC%tblck(q)%npp1
     nh2 = CC%tblck(q)%nhh2

     do J3 = J3Min, J3Max, 2

        if ( J3 > 2*jbas%jtotal_max) cycle
        
        If ( J3 .ge. J1 )  then
           
              P_a = CC%tblck(q)%lam(2)
              Tz = CC%tblck(q)%lam(3) 
              P_b = mod(P_a + AA%dpar/2,2)
              
              q1 = tensor_block_index(J1,J3,rank_a,Tz,P_a) 
              q2 = iso_ladder_block_index(J3,J2,rank_b,Tz,P_b)                            
              np2 = AA%tblck(q1)%npp2
                            
              if (np1*nh2 .ne. 0)  then                       
                 if (np2 .ne. 0) then 
                    ! Apppp . Bpphh 

                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)

                    bet_off = 1.d0
                    
                    call dgemm('N','N',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np1,&
                         BB%tblck(q2)%Xpphh,np2,bet_off,CC%tblck(q)%Xpphh,np1)

                 end if
              end if

           else

              P_a = mod(CC%tblck(q)%lam(2)+AA%dpar/2,2) 
              Tz = CC%tblck(q)%lam(3) 
              P_b = P_a
              q2 = iso_ladder_block_index(J3,J2,rank_b,Tz,P_b)                            

              q1 = tensor_block_index(J3,J1,rank_a,Tz,P_a) 
              
              np2 = AA%tblck(q1)%npp1
              
              if (np1*nh2 .ne. 0)  then                       
                 if (np2 .ne. 0) then 
                    ! Apppp . Bpphh 
                    
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*&
                         (-1)**((J2+J3+rank_c)/2)*sqrt(rank_c+1.d0)* AA%herm 
                    
                    bet_off = 1.d0

                    call dgemm('T','N',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np2,&
                         BB%tblck(q2)%Xpphh,np2,bet_off,CC%tblck(q)%Xpphh,np1)

                 end if
              end if
              
           end if

        end do


        J3Min = max(abs(rank_b-J1),abs(rank_a-J2))
        J3Max = min(rank_b + J1,rank_a + J2)
        
        do J3 = J3Min, J3Max, 2


           if ( J3 > 2*jbas%jtotal_max) cycle        


           IF (J2 .ge. J3 ) then 


              Tz = CC%tblck(q)%lam(3) 

              !----------------------------------------------------------------------------
              !         Zpphh 
              !----------------------------------------------------------------------------
                            
              P_b = CC%tblck(q)%lam(2)
              P_a = mod(P_b + BB%dpar/2,2)

              q1 = iso_ladder_block_index(J1,J3,rank_b,Tz1,P_b) 
              q2 = tensor_block_index(J3,J2,rank_a,Tz2,P_a) 

              nh1 = BB%tblck(q1)%nhh2     
              !w1pphh = -Bpphh.Ahhhh 

              if (np1*nh2 .ne. 0)  then                       
                 if (nh1 .ne. 0) then 
                    al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*&
                         (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('N','N',np1,nh2,nh1,al_off,BB%tblck(q1)%Xpphh,np1,&
                         AA%tblck(q2)%tgam(5)%X,nh1,bet_off,CC%tblck(q)%Xpphh,np1)
                 end if
              end if

           else ! ( J2< J3) 

              Tz = CC%tblck(q)%lam(3) 

              !----------------------------------------------------------------------------
              !         Zpphh 
              !----------------------------------------------------------------------------

              P_b = CC%tblck(q)%lam(2)
              P_a = mod(P_b + CC%dpar/2,2)

              q1 = iso_ladder_block_index(J1,J3,rank_b,Tz1,P_b)                            
              q2 = tensor_block_index(J2,J3,rank_a,Tz2,P_a) 

              nh1 = BB%tblck(q1)%nhh2     
              !w1pphh = -Bpphh.Ahhhh 

              if (np1*nh2 .ne. 0)  then                       
                 if (nh1 .ne. 0) then 
                    al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*AA%tblck(q2)%lam(1)*&
                         (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)*AA%herm 
                    bet_off = 1.d0
                    call dgemm('N','T',np1,nh2,nh1,al_off,BB%tblck(q1)%Xpphh,np1,&
                         AA%tblck(q2)%tgam(5)%X,nh2,bet_off,CC%tblck(q)%Xpphh,np1)
                 end if
              end if


           end if
           
           
        end do


     end do
   end subroutine dTz_tensor_product_222_pp_hh
!=================================================================
!=================================================================
 subroutine dTz_tensor_product_222_ph(AA,BB,CC,jbas) 
   ! VERIFIED ph channel 2body tensor_product. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(sq_op) :: AA
   type(iso_ladder) :: BB,CC
   real(8),allocatable,dimension(:,:) :: PANDYA_A,PANDYA_B,PANDYA_AB
   type(cc_mat) :: LCC
   integer :: nh,np,nb1,nb2,q,IX,JX,i,r1,r2,Tz1_cc,Tz2_cc,PAR_J3,JTM,q1,q2,J3,J4,rank_a,a,p1,p2,h1,h2
   integer :: ji,jh1,jh2,ti,th1,th2,lh1,lh2,li,n1,n2,c1,c2,jxstart,J4min,J4max,ja,jp1,jp2
   integer :: J1,J2, Jtot,Ntot,qx,J3min,J3max,ril,rjk,rli,rkj,g_ix,thread,total_threads
   integer :: phase1,phase2,phase3,rik,rki,rjl,rlj,PAR2,J1min,J2min,J1max,J2max,n_J3,n_J4,n_J5
   integer :: phase_34,phase_abcd,phase_ac,phase_bc,j5,PAR_J4,PAR_J5,rank_b,rank_c
   integer :: phase_bd,phase_ad,nj_perm,full_int_phase,J5min,J5max,tp1,tp2,lp1,lp2
   integer,allocatable,dimension(:,:) :: qn_J3,qn_J4,qn_J5
   real(8) :: sm ,pre,pre2,omp_get_wtime ,t1,t2,coef9,factor,sm_ex,nj1,nj2  
   real(8) :: prefac_12,prefac_1,nj,Xelem,Yelem,V,al_off,d6ji
   logical :: square

   Ntot = CC%Nsp
   JTM = jbas%Jtotal_max*2
!   total_threads = size(CC%direct_omp) - 1
   rank_a = AA%rank
   rank_b = BB%rank
   rank_c = CC%rank 
   ! construct intermediate matrices

   ! sum over all cross coupled channels for both matrices 
   do J3 = 0, JTM  ,2
      J4min = abs(rank_a -J3)
      J4max = min(rank_a + J3,JTM)

      do Tz1_cc = 0,2,2 ! this is the cross coupled TZ (no negative)
         Tz2_cc = abs(Tz1_cc - 2*abs(CC%dTz))
         
         do PAR_J3 = 0,1

            n_J3 = count_configs( J3 ,Tz1_cc, PAR_J3 , jbas, .true. ) !boolean true: ph, false: hp configs
            if (n_J3 == 0) cycle
            allocate(qn_J3(n_J3,2)) 
            n_J3 = count_configs( J3 ,Tz1_cc, PAR_J3 , jbas, .true.,qn_J3) ! this fills the cc- basis descriptor

            PAR_J4 = mod(PAR_J3 + AA%dpar/2,2)
            PAR_J5 = mod(PAR_J4 + BB%dpar/2,2)

            do J4 = J4min,J4max ,2

               J5min = max(abs(rank_b-J4),abs(rank_c-J3)) 
               J5max = min(rank_b+J4,rank_c+J3,JTM)

               n_J4 = count_configs( J4 ,Tz1_cc, PAR_J4 , jbas, .true. )
               if (n_J4 == 0) cycle
               allocate(qn_J4(n_J4,2)) 
               n_J4 = count_configs( J4 ,Tz1_cc, PAR_J4 , jbas, .true.,qn_J4)

               allocate( PANDYA_A( n_J3, n_J4) )
               PANDYA_A = 0.d0 
               call fill_generalized_pandya_matrix(J3,J4,PANDYA_A,qn_J3,qn_J4,AA,jbas)
               
               do J5 = J5min,J5max,2 

                  n_J5 = count_configs( J5 ,Tz2_cc, PAR_J5 , jbas, .false. )
                  if (n_J5 == 0) cycle
                  allocate(qn_J5(n_J5,2)) 
                  n_J5 = count_configs( J5 ,Tz2_cc, PAR_J5 , jbas, .false.,qn_J5)

                  allocate(PANDYA_B(n_J4, n_J5),PANDYA_AB(n_J3,n_J5)) 
                  PANDYA_B = 0.d0 
                  PANDYA_AB = 0.d0 
                  call fill_isoladder_pandya_matrix(J4,J5,PANDYA_B,qn_J4,qn_J5,BB,jbas)                  

                  al_off = d6ji(rank_a,rank_b,rank_c,J5,J3,J4) * (-1)** ((rank_C + J3)/2) * &
                       sqrt((J3+1.d0)*(J5+1.d0)*(rank_c+1.d0)) 
                  call dgemm('N','N',n_J3,n_J5,n_J4,al_off,PANDYA_A,n_J3,PANDYA_B,n_J4,bet,PANDYA_AB,n_J3) 

                                  
                  do JX = 1,n_J5

                     ! GET KET 
                     h2 = qn_J5(JX,1)
                     p2 = qn_J5(JX,2)
                     
                     jp2 = jbas%jj(p2)
                     lp2 = jbas%ll(p2)
                     tp2 = jbas%itzp(p2)
                     jh2 = jbas%jj(h2)
                     lh2 = jbas%ll(h2)
                     th2 = jbas%itzp(h2)

                     do IX = 1, n_J3 

                        ! GET BRA
                        p1 = qn_J3(IX,1)
                        h1 = qn_J3(IX,2)

                        jp1 = jbas%jj(p1)
                        jh1 = jbas%jj(h1)
                        lp1 = jbas%ll(p1)
                        tp1 = jbas%itzp(p1)
                        lh1 = jbas%ll(h1)
                        th1 = jbas%itzp(h1)

                        if ( (tp1 +tp2) -2*CC%dTz .ne. (th1+th2) ) cycle
                        if ( mod(lp1 +lp2,2) .ne. mod(lh1+lh2+CC%dpar/2,2) ) cycle
                        
                        ! CALCULATE X CONTRIBUTIONS

                        J1min = abs(jp1-jp2)
                        J1max = jp2+jp1

                        ! these are the results of the Matmuls 
                        Xelem = PANDYA_AB(IX,JX)

                        do J1 = J1min,J1max,2
                           if ((p1==p2).and.(mod(J1/2,2)==1)) cycle                           
                           J2min = max(abs(jh1-jh2),abs(rank_c-J1))
                           J2max = min(jh1+jh2 ,rank_c+J1) 
                           prefac_1 = sqrt(J1+1.d0)

                           do J2 = J2min,J2max,2
                              if ((h1==h2).and.(mod(J2/2,2)==1)) cycle                           
                              nj = ninej(CC%xindx,jp1,jh1,J3,jp2,jh2,J5,J1,J2,rank_c)
                              prefac_12 = prefac_1 *sqrt(J2+1.d0)



                              if ( h2 .ge. h1 ) then
                                 if ( p2 .ge. p1 ) then                                   
                                    V = prefac_12*nj*(-1)**((jp2-jh2)/2)*Xelem

                                    call add_elem_to_ladder(V,p1,p2,h1,h2,J1,J2,CC,jbas)                                         
                                 else

                                    V = prefac_12*nj*(-1)**((jp1-jh2)/2)*Xelem*(-1)**(J1/2)

                                    call add_elem_to_ladder(V,p2,p1,h1,h2,J1,J2,CC,jbas)                                         

                                 end if
                              else
                                 if (p1 > p2) then
                                    
                                    V = prefac_12*nj*(-1)**((jp1-jh1)/2)*Xelem*(-1)**((J1+J2)/2)
                                    
                                    call add_elem_to_ladder(V,p2,p1,h2,h1,J1,J2,CC,jbas)                                         
                                 else
                                    
                                    V = prefac_12*nj*(-1)**((jp2-jh1)/2)*Xelem*(-1)**((J2)/2)
                                    
                                    call add_elem_to_ladder(V,p1,p2,h2,h1,J1,J2,CC,jbas)                                         
                                 end if
                              end if
                           end do
                        end do
                     end do
                  end do

                  deallocate(PANDYA_B,qn_J5,PANDYA_AB) 
               end do
               deallocate(PANDYA_A,qn_J4)
            end do
            deallocate(qn_J3)
         end do
      end do
   end do
 end subroutine dTz_tensor_product_222_ph
 !===================================================================================
 !===============================================================================
 subroutine fill_isoladder_pandya_matrix(J1,J2,MAT,qn1,qn2,OP,jbas)
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
         c=qn2(JJ,1)
         d=qn2(JJ,2)

         MAT(II,JJ) = Visopandya(a,b,c,d,J1,J2,Op,jbas)

      end do
   end do

 end subroutine fill_isoladder_pandya_matrix
 
end module 
  
  
  

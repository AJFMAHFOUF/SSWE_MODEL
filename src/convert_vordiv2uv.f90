subroutine convert_vordiv2uv

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: i1, i2, ms, js, jsm, jsp, j_index2
   
! Physical fields required for non linear terms   
   
 call legt_i(phi_m,phi_mn(:,2),0)
 call fft_i(phi,phi_m) 
  
 call legt_i(vor_m,vor_mn(:,2),0)
 call fft_i(vor,vor_m)   
     
! Stream function and velocity potential (spectral)

 psi_mn(:) = (0.,0.)
 khi_mn(:) = (0.,0.) 
 do i1 = 0,mm
   ms=abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     if (i2 > 0) then 
       psi_mn(js) = -a*a/(i2*(i2+1.0))*vor_mn(js,2)
       khi_mn(js) = -a*a/(i2*(i2+1.0))*div_mn(js,2)   
     endif  
   enddo
 enddo      
!
!  U and V wind components (spectral)
!
 u_mn(:) = (0.,0.)
 v_mn(:) = (0.,0.)
 do i1 = 0,mm
   ms=abs(i1)
   do i2 = ms,mm
     js  = j_index2(mm,ms,i2)
     jsm = j_index2(mm,ms,max0(0,i2-1))
     jsp = j_index2(mm,ms,i2+1)
     u_mn(js) = j*i1*khi_mn(js) + (i2-1)*eps(i1,i2)*psi_mn(jsm) - & 
                   (i2+2)*eps(i1,i2+1)*psi_mn(jsp)   
     v_mn(js) = j*i1*psi_mn(js) - (i2-1)*eps(i1,i2)*khi_mn(jsm) + &
                   (i2+2)*eps(i1,i2+1)*khi_mn(jsp)  
   enddo
   i2 = mm + 1
   js  = j_index2(mm,ms,i2)
   jsm = j_index2(mm,ms,i2-1)
   u_mn(js) =   (i2-1)*eps(i1,i2)*psi_mn(jsm)    
   v_mn(js) = - (i2-1)*eps(i1,i2)*khi_mn(jsm) 
 enddo  

! Back to physical space

  call legt_i(u_m,u_mn,1)
  call fft_i(u,u_m) 
  
  call legt_i(v_m,v_mn,1)
  call fft_i(v,v_m) 
 
 return

end subroutine convert_vordiv2uv

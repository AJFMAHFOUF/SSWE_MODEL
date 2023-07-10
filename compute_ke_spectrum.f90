subroutine compute_ke_spectrum(nstep)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer, intent(in)   :: nstep
 integer :: i1, i2, ms, js, j_index2, ihour
 character(len=1)      :: ichst1
 character(len=2)      :: ichst2
 character(len=3)      :: ichst
 character(len=2)      :: tt
 real, dimension(0:mm) :: ke_spectrum, zcoef
 
 ihour = nstep*dt/3600.0
 
 if (ihour < 10) then 
   write(ichst1,'(i1)') ihour
   ichst = '00'//ichst1
 elseif (ihour < 100) then
   write(ichst2,'(i2)') ihour
   ichst = '0'//ichst2
 else
   write(ichst,'(i3)') ihour
 endif      
 
 write(tt,'(i2)') mm 
   
! Compute kinetic energy spectrum (formulation from Lambert, 1984)   
 
 zcoef = 0.0
 do i2 = 1,mm
   zcoef(i2) = a*a/float(i2*(i2+1))
 enddo   
 
 ke_spectrum(:) = 0.0
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     if (i1 == 0) then
       ke_spectrum(i2) = ke_spectrum(i2) + 0.5*zcoef(i2)*(abs(vor_mn(js,2))**2 + abs(div_mn(js,2))**2)
     else
       ke_spectrum(i2) = ke_spectrum(i2) + zcoef(i2)*(abs(vor_mn(js,2))**2 + abs(div_mn(js,2))**2)
     endif 
   enddo   
 enddo
 
! Store results in file 
 
 open (unit=25,file='../DATA/SSWE_model_spectrum_T'//tt//'_step_'//ichst//'_expid_'//expid//'.dat',status='unknown')
 
 do i1=1,mm
   write(25,*) i1,ke_spectrum(i1)
 enddo  
 
 close(unit=25)
      
 return

end subroutine compute_ke_spectrum

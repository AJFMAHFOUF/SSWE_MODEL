subroutine init
 
 use fft99_mod
 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer            :: i1, i2, j1
 real               :: zlon, zlat, zzz, zweight 
 real, dimension(4) :: zfield
 character(len=3)   :: tt
 character(len=2)   :: tt1

 
! Define character for truncation

 if (mm < 99) then    
   write(tt1,'(i2)') mm 
   tt = '0'//tt1
 else
   write(tt,'(i3)') mm
 endif     
 
! Initialisation prior FFT

 call set99(trigs,ifax,nlon)
 
! Define array for recurrence formula (Legendre polynomials)

 eps(:,:) = 0.0
 do i1 = -mm,mm
   do i2 = abs(i1),mm+1
     if (i1 /= i2) eps(i1,i2) = sqrt(float(i2*i2 - i1*i1)/float(4*i2*i2 - 1))
   enddo
 enddo 

! Read initial physical fields (vorticity, divergence, geopotential, u and v winds)  
     
 open (unit=10,file='../data_in/VOR_'//cdate//'_00_T'//tt//'gg.dat',status='old')
 open (unit=11,file='../data_in/U_'//cdate//'_00_T'//tt//'gg.dat',status='old')
 open (unit=12,file='../data_in/V_'//cdate//'_00_T'//tt//'gg.dat',status='old') 
 open (unit=13,file='../data_in/PHI_'//cdate//'_00_T'//tt//'gg.dat',status='old') 

 phibar = 0.0
 zweight = 0.0
 do j1 = 1,nlat  
   do i1 = 1,nlon+1
     read(10,*) zlon, zlat, zfield(1)
     read(11,*) zlon, zlat, zfield(2)
     read(12,*) zlon, zlat, zfield(3)
     read(13,*) zlon, zlat, zfield(4)
     if (i1 /= nlon+1) then
       zzz=cos(pi*zlat/180.)
       vor(i1,j1) = zfield(1)
       utr(i1,j1) = zfield(2)
       vtr(i1,j1) = zfield(3) 
       phibar = phibar + zfield(4)*cos(pi/180.*zlat)
       zweight = zweight + cos(pi/180.*zlat)
       if (.not.l_real_ic) then
         utr(i1,j1) = 25.0*zzz - 30.0*zzz**3 + 300.0*(1.0-zzz**2)*zzz**6
         vtr(i1,j1) = 0.0
       endif
     endif
   enddo
 enddo
 phibar = phibar/zweight
 
 print *,'Initial fields have been read from input files for vor, u and v ',phibar
 
 close(unit=10)
 close(unit=11)
 close(unit=12)
 close(unit=13)

! Read Gaussian latitudes (sin) a and Gaussian weights w 
 
 open (20,file='../data_in/gaulat_legpol_T'//tt//'a.dat',status='old')
 
 read(20,*) (x(j1),j1=1,nlat)
 read(20,*) (w(j1),j1=1,nlat)
 
 print *,'Gaussian latitudes and weights have been read'
 
 read(20,*) ((pl_legendr(i1,j1),j1=1,nlat),i1=1,mmax)
 
 print *,'Legendre polynomials have been read'    
       
 close(unit=20)
 
! Compute Coriolis parameter and transformed winds
 
 do j1 = 1,nlat
   f(j1) = 2.*omega*x(j1)
   u(:,j1) = utr(:,j1)*a*sqrt(1.0 - x(j1)*x(j1))
   v(:,j1) = vtr(:,j1)*a*sqrt(1.0 - x(j1)*x(j1))
 enddo
 
! Spectral coefficients for vorticity (assumes vorticity field defined)

 call fft_d(vor,vor_m)    
 call legt_d(vor_m,vor_mn(:,2))
                                                                                                                                   
! Fill other arrays
   
 vor_mn(:,1) = vor_mn(:,2)
 vor_mn(:,3) = vor_mn(:,2)   
 
 print *,'Exit from initialisation subroutine'      
      
 return     
end subroutine init

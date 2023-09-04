subroutine init
 
 use fft99_mod
 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer            :: i1, i2, j1
 real               :: zlon, zlat, zweight
 real               :: zfield, zu0, zh0, zr0, zr, zr2, zlonc, zlatc, zzz
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

! Read Gaussian latitudes (sin) a and Gaussian weights w 
 
 open (20,file='../data_in/gaulat_legpol_T'//tt//'a.dat',status='old')
 
 read(20,*) (x(j1),j1=1,nlat)
 read(20,*) (w(j1),j1=1,nlat)
 
 print *,'Gaussian latitudes and weights have been read'
 
 read(20,*) ((pl_legendr(i1,j1),j1=1,nlat),i1=1,mmax)
 
 print *,'Legendre polynomials have been read'    
       
 close(unit=20)

! Define analytical initial fields for Case 5 from Williamson et al. (1992)  

 phi_bar = 5500.0*g
 zu0 = 2.0*pi*a/(12.0*86164.1)
 zu0 = 0.5*zu0 ! necessary to recover published results
 zh0 = 2000.0
 zr0 = pi/9.0
 zlonc = 0.5*pi
 zlatc = pi/6.0
 zzz = a*omega*zu0 + 0.5*zu0*zu0
 
 do j1 = 1,nlat  
   do i1 = 1,nlon
     zlat = asin(real(x(j1)))*180.0/pi
     zlon = 360./float(nlon)*(i1-1)
     vor(i1,j1) = 2.0*zu0/a*sin(zlat*pi/180.)
     div(i1,j1) = 0.0
     phi(i1,j1) = 5960.0*g - zzz*(sin(zlat*pi/180.))**2 
     utr(i1,j1) = zu0*cos(zlat*pi/180.)
     vtr(i1,j1) = 0.0
     zr2 = min(zr0*zr0, (zlon*pi/180. - zlonc)**2 + (zlat*pi/180. - zlatc)**2)
     zr = sqrt(zr2)
     phis(i1,j1) = zh0*g*(1.0 - zr/zr0)
     phi(i1,j1) = phi(i1,j1) - phis(i1,j1) 
   enddo
 enddo
 
 print *,'Initial fields have been defined for vor, div, phi, u and v',phi_bar
 
! Compute Coriolis parameter and transformed winds
 
 do j1 = 1,nlat
   f(j1) = 2.*omega*x(j1)
   u(:,j1) = utr(:,j1)*a*sqrt(1.0 - x(j1)*x(j1))
   v(:,j1) = vtr(:,j1)*a*sqrt(1.0 - x(j1)*x(j1))
 enddo
 
! Spectral coefficients for orography

 call fft_d(phis,phis_m)
 call legt_d(phis_m,phis_mn) 
 
! Spectral coefficients for vorticity (assumes vorticity field defined)

 call fft_d(vor,vor_m)    
 call legt_d(vor_m,vor_mn(:,2))

! Spectral coefficients for divergence (assumes divergence field defined)
 
 call fft_d(div,div_m) 
 call legt_d(div_m,div_mn(:,2))
                                                                                                                                
! Spectral coefficients for geopotential
 
 call fft_d(phi,phi_m) 
 call legt_d(phi_m,phi_mn(:,2))
 
! Spectral coefficients for tracer
 
 call fft_d(qv,qv_m) 
 call legt_d(qv_m,qv_mn(:,2)) 
   
! Fill other arrays
   
 vor_mn(:,1) = vor_mn(:,2)
 vor_mn(:,3) = vor_mn(:,2)   
 div_mn(:,1) = div_mn(:,2)
 div_mn(:,3) = div_mn(:,2)
 phi_mn(:,1) = phi_mn(:,2)
 phi_mn(:,3) = phi_mn(:,2)
 qv_mn(:,1)  = qv_mn(:,2)
 qv_mn(:,3)  = qv_mn(:,2)
 
 print *,'Exit from initialisation subroutine'      
      
 return     
end subroutine init

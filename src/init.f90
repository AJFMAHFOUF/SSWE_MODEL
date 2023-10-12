subroutine init
 
 use fft99_mod
 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer            :: i1, i2, j1
 real               :: zlon, zlat, zweight
 real               :: zomega, zk, zr, zphi0, za, zb, zc 
 real               :: zlatr, zlonr, zclatr, zslatr
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

 phi_bar = 0.0
 zweight = 0.0
!
! Set-up parameters for Haurwitz wave 
! 
 zomega = 7.848E-6
 zk = zomega
 zr = 4.0
 zphi0 = g*8.0E3
! 
! Set-up analytical form for initial conditions 
!
 do j1 = 1,nlat  
   do i1 = 1,nlon
     zlat = asin(real(x(j1)))*180.0/pi
     zlon = 360.0/float(nlon)*(i1-1)
     zlatr = zlat*pi/180.0
     zlonr = zlon*pi/180.0    
     zclatr = cos(zlatr)
     zslatr = sin(zlatr)
     
     vor(i1,j1) = 2.0*zomega*zslatr - zk*zslatr*zclatr**zr*(zr*zr + 3.0*zr + 2.0)*cos(zr*zlonr)
                
     div(i1,j1) = 0.0
     
     za = 0.5*zomega*(2.0*omega + zomega)*zclatr**2 + &
        & 0.25*zk*zk*zclatr**(2.0*zr)*((zr + 1.0)*zclatr**2 + & 
        & (2.0*zr*zr - zr - 2.0) - 2.0*zr*zr/zclatr**2)   
     zb = 2.0*(omega + zomega)*zk/((zr + 1.0)*(zr + 2.0))*zclatr**zr* &
        & ((zr*zr + 2.0*zr + 2.0) - ((zr + 1.0)*zclatr)**2)         
     zc = 0.25*zk*zk*zclatr**(2.0*zr)*((zr + 1.0)*zclatr**2 - (zr + 2.0))     
     phi(i1,j1) = zphi0 + a*a*(za + zb*cos(zr*zlonr) + zc*cos(2.0*zr*zlonr))   
     
     utr(i1,j1) = a*zomega*zclatr + a*zk*zclatr**(zr - 1.0)* &
                & (zr*zslatr**2 - zclatr**2)*cos(zr*zlonr)
                
     vtr(i1,j1) = -a*zk*zr*zclatr**(zr - 1.0)*zslatr*sin(zr*zlonr)  
     
     qv(i1,j1) = 0.0
       
     zweight = zweight + abs(zclatr)    
     phi_bar = phi_bar + phi(i1,j1)*abs(zclatr)               
   enddo
 enddo
 phi_bar = phi_bar/zweight
 
 print *,'Initial fields have been initialized for vor, div, phi, u and v',phi_bar
 
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

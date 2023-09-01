program read_era5_for_gg
 integer, parameter :: nx1=1440, ny1=721 ! regular 0.25°x0.25° lat/lon grid
 integer, parameter :: nx2=128,   ny2=64  ! Gaussian grid for T42 triagular truncation
 !integer, parameter :: nx2=192,   ny2=96  ! Gaussian grid for T63 triagular truncation
 !integer, parameter :: nx2=160,   ny2=80  ! Gaussian grid for T53 triagular truncation
 !integer, parameter :: nx2=96,    ny2=48  ! Gaussian grid for T31 triagular truncation
 !integer, parameter :: nx2=320,   ny2=160  ! Gaussian grid for T106 triagular truncation
 real*8, dimension(ny2) :: xg, wg
 real, dimension(nx1,ny1,5) :: field1
 real, dimension(nx2,ny2,5) :: field2
 real, dimension(nx2,ny2)   :: lat1, lon1
 real :: lat_0, lon_0, dlat, dlon, pi, zdummy
 real :: x1, x2, delta
 integer :: i, j, ii, jj, k, kk
 character (3) :: trunc
 character (8) :: date
 
 date='21121978'
 
 trunc='042'
!
! Open input files - ERA5 regular lat/lon grid 0.25°x0.25° resolution
!  
 open (unit=55,file='../data_in/vorticity_'//date//'_00_025res.dat',status='old')
 open (unit=65,file='../data_in/divergence_'//date//'_00_025res.dat',status='old')
 open (unit=75,file='../data_in/geopotential_'//date//'_00_025res.dat',status='old')
 open (unit=85,file='../data_in/uwind_'//date//'_00_025res.dat',status='old')
 open (unit=95,file='../data_in/vwind_'//date//'_00_025res.dat',status='old')
! 
! Open output files - Gaussian (quadratic) grid on selected truncation
! 
 open (unit=56,file='../data_in/VOR_'//date//'_00_T'//trunc//'gg.dat',status='unknown')
 open (unit=66,file='../data_in/DIV_'//date//'_00_T'//trunc//'gg.dat',status='unknown')
 open (unit=76,file='../data_in/PHI_'//date//'_00_T'//trunc//'gg.dat',status='unknown')
 open (unit=86,file='../data_in/U_'//date//'_00_T'//trunc//'gg.dat',status='unknown')
 open (unit=96,file='../data_in/V_'//date//'_00_T'//trunc//'gg.dat',status='unknown')
!
! Define the elements of the original grid
!  
 lat_0 = -90.0
 lon_0 = 0.0
 dlat = 0.25
 dlon = 0.25
!
! Define constants
! 
 pi = acos(-1.0)
!
! Read initial field
! 
 k = 1
 do kk=0,40,10
   read (55+kk,*) ! skip header
   do j = 1,ny1
     do i = 1,nx1
       read (55+kk,*) zdummy,zdummy,field1(i,ny1-j+1,k)
     enddo  
   enddo 
   k = k + 1
 enddo    
!
! Find the gaussian latitudes and associated weights
!
 call gauleg(-1.0D0,1.0D0,xg,wg,ny2)
!
! Parameters for a regular grid along a longitude circle
! Gaussian grid being close to a regular grid (for plotting purposes)
! 
 x1 = asin(real(xg(1)))*180.0/pi
 x2 = asin(real(xg(ny2)))*180.0/pi
 delta = (x2 - x1)/float(ny2-1) 
!
! Find the closest point for the target grid
! 
 k = 1
 do kk=0,40,10
   do j=1,ny2
     do i = 1,nx2
       lat1(i,j) = asin(real(xg(j)))*180.0/pi
       lon1(i,j) = lon_0 + 360./float(nx2)*(i-1)
       ii = INT((lon1(i,j) - lon_0)/dlon + 1)
       jj = INT((lat1(i,j) - lat_0)/dlat + 1)
       field2(i,j,k) = field1(ii,jj,k)
       write (56+kk,*) lon1(i,j),lat1(i,j),field2(i,j,k)
       if (i == nx2) then ! a longitude is added to insure periodicity (for plotting purposes)
         write(56+kk,*) 360.0,lat1(i,j),field2(1,j,k)
       endif   
     enddo
   enddo
   k = k + 1
 enddo      
 stop  
end program read_era5_for_gg

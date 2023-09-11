module params

 implicit none

 integer, parameter :: mm = 63                ! maximum wave number
 integer, parameter :: nlat = (3*mm+1)/2 +1   ! number of latitudes
 integer, parameter :: nlon = 2*nlat          ! number of longitudes
 integer, parameter :: mmax = (mm+1)*(mm+4)/2 ! number of stored wavenumbers
 integer, parameter :: nfft = 1               ! number of FFT to be done
 
 complex, parameter :: j = (0,1)              ! square root of -1 
 real, parameter    :: a = 6371.22E3          ! Earth radius
 real, parameter    :: pi = acos(-1.0)        ! Pi constant
 real, parameter    :: g = 9.80616            ! Earth gravitational acceleration
 real, parameter    :: omega = 2.0*pi/86164.1 ! Earth angular speed (stellar day)
 real, parameter    :: nu = 0.02, wk = 0.53   ! tunable parameters for 2*dt filter 
 real, parameter    :: kdiff = 1.0E15         ! Coefficient for horizontal diffusion
 real, parameter    :: dt = 900.0             ! model time step
 integer, parameter :: nhtot = 120            ! number of hours of model integration
 integer, parameter :: npdt = nhtot*3600/dt   ! number of model time steps
 integer, parameter :: nfreq = 24*3600/dt     ! hourly output archiving frequency
 character(len=3)   :: expid='t10'            ! experiment identifier
 logical            :: lreaduv=.true.         ! logical to use u v at initial time
 logical            :: lsemimp=.true.         ! semi-implicit scheme

 
end module params 
 
module model_vars  
    
 use params   
 
 implicit none
 
 real, dimension (nlon,nlat) :: vor, div, phi, qv  ! prognostic variables in physical space
 real, dimension (nlon,nlat) :: psi, khi, u, v, ke ! diagnostic variables in physical space
 real, dimension (nlon,nlat) :: utr, vtr           ! u and v winds in geographical coordinates
 real, dimension (nlon,nlat) :: phis               ! surface geopotential
 real                        :: phi_bar            ! mean value of geopotential (for SI scheme)
 
 complex, dimension(nlat,-mm:mm) :: vor_m, div_m, phi_m, qv_m, phis_m
 complex, dimension(nlat,-mm:mm) :: uvor_m, vvor_m
 complex, dimension(nlat,-mm:mm) :: psi_m, khi_m, u_m, v_m, ke_m  
 
 complex, dimension(mmax,3) :: vor_mn, div_mn, phi_mn, qv_mn  ! prognostic variables in spectral space (3 time steps)
 complex, dimension(mmax)   :: psi_mn, khi_mn, u_mn, v_mn, ke_mn, phis_mn
 
end module model_vars

module spectral_vars

 use params
 
 implicit none
  
 real, dimension(mmax,nlat) :: pl_legendr ! Legendre polynomials
 real, dimension(nlat)      :: w, x       ! Gaussian weights and latitudes (sinus value)
 real, dimension(nlat)      :: f          ! Coriolis parameter
 
 real, dimension(-mm:mm,0:mm+1) :: eps    ! array for Legendre polynomials derivatives
 
! arrays for FFT991 

 integer, dimension(13)         :: ifax
 real, dimension(3*nlon/2+1)    :: trigs 
 real, dimension(nfft*(nlon+2)) :: acoef
 real, dimension(nfft*(nlon+1)) :: work
 
end module spectral_vars

program main_sswm
!=================================================================================
!
! Global spectral shallow water model in Fortran 90
! with vorticity, divergence and geopotential as prognostic variables
! Explicit leapfrog scheme - Triangular truncation (defined in module params)
!
! Equations taken from Coiffier (2011) - pages 128-129
!
! Jean-François Mahfouf (22/06/2023)
!
! Use of FFT99 defined by Temperton (ECMWF) 
!
! Option in order to apply a DFI to the initial fields (Lynch and Huang, 1992)
!
!==================================================================================
 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: nstep, npdt_max
 real    :: t1, t2, dt1 
 logical :: ldfi, loutput
 
 type (prog_var) :: xin, xin0
 type (prog_var) :: xout, xout_filtered

! Required initialisations and read fields in physical space

 call init

 if (lreaduv) call convert_uv2vordiv
 
 call cpu_time(time=t1)
 
! Set-up initial fields in spectral space

 xin%vormn = vor_mn(:,1)
 xin%divmn = div_mn(:,1) 
 xin%phimn = phi_mn(:,1)
 xin%qvmn  = qv_mn(:,1)
 
! Save initial conditions 
 
 xin0 = xin
 
 if (linit) then 
 
! DFI Forward integration  
 
   dt1 = dt
   ldfi = .true.
   npdt_max = ndfi_win*1800/dt
   loutput = .false.
 
   call model(xin,xout,dt1,npdt_max,ldfi,loutput)
 
! Save output

   xout_filtered = xout
 
! DFI Backward integration

   dt1 = -dt
   xin = xin0 
 
   call model(xin,xout,dt1,npdt_max,ldfi,loutput)
 
! Save output

   xout_filtered%vormn = xout_filtered%vormn + xout%vormn
   xout_filtered%divmn = xout_filtered%divmn + xout%divmn 
   xout_filtered%phimn = xout_filtered%phimn + xout%phimn  
   xout_filtered%qvmn  = xout_filtered%qvmn  + xout%qvmn
   
   xin = xout_filtered
   
 else ! no DFI
 
   xin = xin0
   
 endif    
 
! Model integration with raw or filtered initial fields

 dt1 = dt
 npdt_max = npdt + 1
 ldfi = .false.
 loutput = .true. 
    
 call model(xin,xout,dt1,npdt_max,ldfi,loutput)
 
 call cpu_time(time=t2)
 
 print *,'Model execution CPU time =',t2-t1
 
 stop

end program main_sswm

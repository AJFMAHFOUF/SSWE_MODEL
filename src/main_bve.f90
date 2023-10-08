program main_bve
!=================================================================================
!
! Global spectral non-divergent vorticity model in Fortran 90
! with relative vorticity as prognostic variable
!
! Explicit leapfrog scheme - Triangular truncation (defined in module params)
!
! Derived from a shallow water model
!
! Jean-François Mahfouf (08/10/2023)
!
! Use of FFT99 defined by Temperton (ECMWF) 
!
!==================================================================================
 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: nstep, npdt_max
 real    :: t1, t2, dt1 
 logical :: loutput
 
 type (prog_var) :: xin
 type (prog_var) :: xout

! Required initialisations and read fields in physical space

 call init

 if (lreaduv) call convert_uv2vor
 
 call cpu_time(time=t1)
 
! Set-up initial fields in spectral space

 xin%vormn = vor_mn(:,1)     
 
! Model integration with initial fields

 dt1 = dt
 npdt_max = npdt + 1
 loutput = .true. 
    
 call model(xin,xout,dt1,npdt_max,loutput)
 
 call cpu_time(time=t2)
 
 print *,'Model execution CPU time =',t2-t1
 
 stop

end program main_bve

program test
  real :: ke, ket
  read (200,*)
  ket = 0.0
  read (*,*) nn
  do i=1,945
    read (200,*) m,n,ke
    if (n == nn) then
       ket = ket + ke
    endif
  enddo    
  print *,'ket for n=',nn,'  ',ket
end program test

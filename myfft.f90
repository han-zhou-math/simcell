module myfft
  use, intrinsic :: iso_c_binding
  !include '/apps/fftw3/intel/18.0/mvapich2/2.3.1/3.3.8/include/fftw3.f03' !  pitzer
  !include '/usr/local/fftw3/intel/18.0/mvapich2/2.3/3.3.8/include/fftw3.f03'!ovens
  !include '/Users/yaol/local/fftw/include/fftw3.f03'
  !include '/home/lyao/local/fftw/gfc/include/fftw3.f03'
  !include '/Users/yaol/local/fftw/ifc/include/fftw3.f03'
  !include '/Users/yaol/local/fftw/gfc/include/fftw3.f03'
  include 'fftw3.f03'
  !!include '/Users/yaol/local/fftw/include/fftw3.f'
end module 

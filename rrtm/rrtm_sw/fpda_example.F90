program main

      use parkind, only: im => kind_im, rb => kind_rb

      use m_fpda_rrtm_sw
      implicit none

      integer(im),parameter :: nlay=20
      real(rb),dimension(1,nlay+1) :: plev
      real(rb),dimension(1,nlay)   :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

      integer(im) :: nbands
      real(rb),allocatable,target,dimension(:)   :: band_lbound,band_ubound,weights
      real(rb),allocatable,target,dimension(:,:) :: tau_gas,tau_ray

      integer(im) :: k

      do k=1,nlay+1
        plev(1,k)=1e3 - (k-1)*1e3/nlay
      enddo
      do k=1,nlay
        tlay(1,k)=50 - (k-1)*50/nlay
      enddo

      h2ovmr=9e-6
      o3vmr =5e-9
      co2vmr=400e-6
      ch4vmr=10e-6
      n2ovmr=320e-9
      o2vmr =.209

      call fpda_rrtm_sw(nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                        nbands, band_lbound,band_ubound, weights, tau_gas, tau_ray)

end program

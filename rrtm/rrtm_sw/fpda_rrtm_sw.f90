module m_fpda_rrtm_sw
      use rrtmg_sw_init, only: rrtmg_sw_ini
      use parkind, only: im => kind_im, rb => kind_rb
      use rrsw_wvn
      use parrrsw, only: nbndsw,naerec,jpb1, jpb2
      use rrtmg_sw_rad, only: rrtmg_sw
      use rrtmg_sw_spcvrt, only: fpda_taugas,fpda_tauray,fpda_solsrc
      
      use iso_c_binding
      implicit none

      logical :: linit=.False.


    contains
      subroutine f2c_fpda_rrtm_sw(nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, c_nbands,c_band_lbound,c_band_ubound, c_weights, c_tau_gas, c_tau_ray) bind(C,name="fpda_rrtm_sw")

          integer(c_int),intent(in) :: nlay

          real(c_double),dimension(1,nlay+1),intent(in) :: plev
          real(c_double),dimension(1,nlay),intent(in)   :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

          integer(c_int),intent(out) :: c_nbands
          type(c_ptr),intent(out) :: c_band_lbound,c_band_ubound
          type(c_ptr),intent(out) :: c_weights,c_tau_gas,c_tau_ray

          real(rb),dimension(1,nlay+1) :: plev_r
          real(rb),dimension(1,nlay)   :: tlay_r, h2ovmr_r, o3vmr_r, co2vmr_r, ch4vmr_r, n2ovmr_r, o2vmr_r

          integer(im) :: nbands,ib
          real(rb),allocatable,save,target,dimension(:)   :: band_lbound,band_ubound,weights  ! [      nbands ]
          real(rb),allocatable,save,target,dimension(:,:) :: tau_gas,tau_ray                  ! [ nlay,nbands ]

          plev_r         (1,:)= rev( plev         (1,:) )
          tlay_r         (1,:)= rev( tlay         (1,:) )
          h2ovmr_r       (1,:)= rev( h2ovmr       (1,:) )
          o3vmr_r        (1,:)= rev( o3vmr        (1,:) )
          co2vmr_r       (1,:)= rev( co2vmr       (1,:) )
          ch4vmr_r       (1,:)= rev( ch4vmr       (1,:) )
          n2ovmr_r       (1,:)= rev( n2ovmr       (1,:) )
          o2vmr_r        (1,:)= rev( o2vmr        (1,:) )

          call fpda_rrtm_sw(nlay,plev_r,tlay_r, h2ovmr_r, o3vmr_r, co2vmr_r, ch4vmr_r, n2ovmr_r, o2vmr_r, nbands, band_lbound, band_ubound, weights, tau_gas, tau_ray)

          do ib = 1,ngptsw
            tau_gas        (:,ib)= rev( tau_gas      (:,ib) )
            tau_ray        (:,ib)= rev( tau_ray      (:,ib) )
          enddo

          c_nbands = ngptsw
          c_band_lbound = c_loc(band_lbound)
          c_band_ubound = c_loc(band_ubound)
          c_weights     = c_loc(weights)
          c_tau_gas     = c_loc(tau_gas)
          c_tau_ray     = c_loc(tau_ray)
      end subroutine

      subroutine fpda_rrtm_sw(nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, nbands, band_lbound,band_ubound, weights, tau_gas, tau_ray)
        integer(im),parameter :: ncol=1,dyofyr=0,inflgsw=0,iceflgsw=0,liqflgsw=0
        real(rb)   ,parameter :: adjes=1, scon=1.36822e+03

        integer(im),intent(in) :: nlay

        real(rb),dimension(:,:),intent(in) :: plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

        integer(im),intent(out) :: nbands
        real(rb),allocatable,target,dimension(:),intent(out)   :: band_lbound,band_ubound,weights
        real(rb),allocatable,target,dimension(:,:),intent(out) :: tau_gas,tau_ray

        real(rb),dimension(ncol,nlay+1) :: tlev,hhl

        real(rb),dimension(ncol,nlay  ) :: play, cldfr, cicewp, cliqwp, reice, reliq 

        real(rb),dimension(nbndsw, ncol,nlay  ) :: taucld, ssacld, asmcld, fsfcld
        real(rb),dimension(ncol, nlay, nbndsw ) :: tauaer, ssaaer, asmaer
        real(rb),dimension(ncol, nlay, naerec ) :: ecaer

        real(rb),dimension(ncol     ) :: tsfc, asdir, aldir, asdif, aldif, coszen

        real(rb),dimension(ncol,nlay+1) :: swuflx,swdflx,swuflxc,swdflxc
        real(rb),dimension(ncol,nlay  ) :: swhr,swhrc

        integer(im) :: iv,ig,k,icol,iw

        integer(kind=im) :: icld=0         ! Cloud overlap method
        integer(kind=im) :: iaer=0         ! Aerosol option flag

        do icol=1,ubound(plev,1)
          call hydrostat_lev(plev(icol,:),tlay(icol,:), hhl(icol,:))
          play(icol,:)      = .5_rb*(plev(icol,1:nlay)+plev(icol,2:nlay+1))

          tlev(icol,nlay+1) = tlay(icol,nlay)
          tsfc(icol)        = tlev(icol,nlay+1)
          tlev(icol,1)      = tlay(icol,1)

          do k = 1,nlay-1
            tlev(icol,k+1)  = .5_rb * (tlay(icol,k+1)+tlay(icol,k) )
          enddo
        enddo

        cldfr  = 0; taucld = 0; ssacld  = 0; asmcld  = 0;
        fsfcld = 0; cicewp = 0; cliqwp  = 0; reice   = 0;
        reliq  = 0; tauaer = 0; ssaaer  = 0; asmaer  = 0;
        ecaer  = 0; coszen = 1; asdir   = 0; aldir   = 0;
        asdif  = 0; aldif  = 0; swdflxc = 0; swuflxc = 0;

        if(.not.linit) then
          call rrtmg_sw_ini(1006._rb)
          linit = .True.

!          print *,'hhl',hhl(1,:)
!          print *,'tlev',tlev(1,:)
!          print *,'tlay',tlay(1,:)
!          print *,'plev',plev(1,:)
!          print *,'play',play(1,:)
        endif

        call rrtmg_sw &
            (ncol    ,nlay    ,icld    ,iaer    , &
            play    ,plev    ,tlay    ,tlev    ,tsfc    , &
            h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
            asdir   ,asdif   ,aldir   ,aldif   , &
            coszen  ,adjes   ,dyofyr  ,scon    , &
            inflgsw ,iceflgsw,liqflgsw,cldfr   , &
            taucld  ,ssacld  ,asmcld  ,fsfcld  , &
            cicewp  ,cliqwp  ,reice   ,reliq   , &
            tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
            swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc)

        nbands = ngptsw
        if(.not.allocated(band_lbound)) allocate(band_lbound(ngptsw) )
        if(.not.allocated(band_ubound)) allocate(band_ubound(ngptsw) )
        if(.not.allocated(weights    )) allocate(weights    (ngptsw) )
        if(.not.allocated(tau_gas    )) allocate(tau_gas    (nlay,ngptsw) )
        if(.not.allocated(tau_ray    )) allocate(tau_ray    (nlay,ngptsw) )

        iw=0
        do iv=1,nbndsw
          do ig=1,ngc(iv)
            iw = iw+1
            band_lbound(iw) = wavenum1(iv+jpb1-1)
            band_ubound(iw) = wavenum2(iv+jpb1-1)
            weights(iw)     = fpda_solsrc(iw)
            tau_gas(:,iw)   = fpda_taugas(:,iw)
            tau_ray(:,iw)   = fpda_tauray(:,iw)      
          enddo
        enddo
      end subroutine

      subroutine hydrostat_lev(plev,tlay, hhl)
        real(rb),intent(in) :: plev(:),tlay(:)
        real(rb),intent(out) :: hhl(size(plev))
        real(rb) :: dp,dz,rho
        integer(im) :: k
        hhl(1) = 0._rb
        do k=1,size(tlay)
          dp  = abs( plev(k)-plev(k+1) ) 
          rho = ( plev(k)+dp/2._rb  ) / 287.058_rb / tlay(k)
          dz  = dp / rho / 9.8065_rb
          hhl(k+1) = hhl(k) + dz
        enddo
      end subroutine

      function rev(inp)
          real(rb),intent(in) :: inp(:)
          real(rb) :: rev(size(inp))
          rev = inp(ubound(inp,1):lbound(inp,1):-1)
      end function

end module


module aero_sectional
   use shr_kind_mod, only: r8 => shr_kind_r8
   !use ppgrid,          only:  pcols, pver, pverp
   use chem_mods,    only: gas_pcnst
   use aerosoldef,   only: l_soa_a1, l_so4_a1

  public :: aerosect_register



   save
        integer, public, parameter                                  ::  secNrBins=5 ! nr of bins
        integer, public, parameter                                  ::  secNrSpec=2  ! number of condensing species
        character(len=20), public, parameter,dimension(secNrSpec)   ::  secspecNames=(/ 'SO4_SEC','SOA_SEC'/)  ! names of condensing species
        character(len=20), public, parameter,dimension(secNrSpec)   ::  specNames=(/ 'SO4','SOA'/)  !condensing species

        real(r8), parameter                                         :: max_diameter=39.6e-9_r8! [m] volume median? diameter 1st mode23.6e-9_r8 ! m
        real(r8), parameter                                         :: min_diameter=5.0e-9_r8 ! m
        integer, public, dimension(secNrSpec)                       ::  secCoagulate_receiver!=(/ l_so4_a1,l_soa_a1/)  ! coagulate receivers

        real(r8),public, dimension(secNrBins)                       :: secMeanD ![m] mean diameter
        real(r8),public,parameter, dimension(secNrSpec)             :: rhopart_sec = (/ 1769.0_r8,1500.0_r8 /)
        !holds the chemistry indices of the sectional scheme:
        integer, public, dimension(secNrSpec, secNrBins)            :: secConstIndex
        real(r8), public, dimension(secNrBins)                      :: sec_rhopart_vols

        integer, public, dimension(secNrBins, secNrBins)            :: autocoag_receiver_index




contains
subroutine aerosect_init()     
        use cam_history,        only: addfld, add_default, outfld,fieldname_len!, phys_decomp
        use ppgrid,             only : pcols, pver

        implicit none
        real(r8)  :: d_rat
        integer   :: i,j
        character(len=20)  :: field_name

        secCoagulate_receiver=(/ l_so4_a1,l_soa_a1/)  ! number of coagulate receivers
        
        ! Use discrite geometric distribution/volume-ratio size distrib:
        d_rat=(max_diameter/min_diameter)**(1._r8/(secNrBins))
        secMeanD(1)=min_diameter
        do i=2,secNrBins
                secMeanD(i) = secMeanD(i-1)*d_rat 
        end do




end subroutine aerosect_init



subroutine aerosect_register()

   use constituents, only: cnst_get_ind
        implicit none
        integer ::          secInd,specInd
        character(len=20)   :: cnst_name
        call aerosect_init()  !
        do secInd=1,secNrBins
                do specInd=1,secNrSpec !names constituents as 'SOA_SEC01'/'SO4_SEC01'
                        WRITE(cnst_name,'(A,I2.2)') trim(secspecNames(specInd)),secInd
                        call cnst_get_ind(trim(cnst_name), secConstIndex(specInd,secInd),abort=.true.)
                end do
        end do




        

end subroutine aerosect_register


! Calculate properties:
subroutine aerosect_write2file(q, lchnk,ncol, pmid, temperature)
        use cam_history,    only: addfld, add_default, outfld,fieldname_len!, phys_decomp
        use ppgrid,         only:  pcols, pver, pverp
         use aerosoldef,    only : chemistryIndex
        use physconst,      only: rair
        implicit none
        integer,  intent(in)    :: lchnk                      ! chunk identifier
        integer,  intent(in)    :: ncol                       ! number of columns

        real(r8), intent(in)    :: q(pcols,pver,gas_pcnst) ! tmr [kg/kg]including moisture
        real(r8), intent(in)    :: pmid(pcols,pver) ! tmr [kg/kg]including moisture
        real(r8), intent(in)    :: temperature(pcols,pver) ! tmr [kg/kg]including moisture
        character(len=20)       :: field_name
        real(r8)                :: rhoAir

        integer :: indBin, indVol, ind_sec,i,k
        real(r8)   :: dummy(pcols,pver)
        do indBin = 1, secNrBins
        !Go through all core species in that mode
                do indVol = 1,secNrSpec
                        ind_sec=chemistryIndex(secConstIndex(indVol,indBin))
                         do k=1,pver
                               do i=1,ncol
                                        rhoAir = pmid(i,k)/rair/temperature(i,k)
                                        call sec_numberConc(q(i,k,ind_sec),indVol, indBin,rhoAir,  dummy(i,k))

                                end do 
                        end do

                        !WRITE(field_name,'(A,I2.2)') 'nr',trim(secspecNames(indVol)),indBin
                        WRITE(field_name,'(A,A,I2.2)') 'nr',trim(secspecNames(indVol)),indBin
                        call outfld(trim(field_name),dummy, pcols, lchnk) !#
                end do 
        end do


end subroutine aerosect_write2file



subroutine sec_numberConc(mass, volNr, binNr , rhoAir, numberConc)
        implicit none

        real(r8), intent(in)    :: mass !kg/kg
        integer, intent(in)     :: binNr ! bin_index
        real(r8), intent(in)    :: rhoAir
        real(r8), intent(out)   :: numberConc !#/m3_air
        integer, intent(in)     :: volNr
        !local:
        integer                 :: specInd
        real(r8), parameter     :: pi=3.141592654_r8

        numberConc = mass/rhopart_sec(volNr)*rhoAir/ & ![kg_aer/kg_air]/[kg_aer/m3_aer]*[kg_air/m3_air]
                              (secMeanD(binNr)**3*pi/6._r8) ! /[m3_aer/#]--> #/m3_air
        if (mass .lt. 1.e-35) then
                numberConc=0.0_r8
        end if






end subroutine sec_numberConc



subroutine sec_moveMass(massDistrib, numberConc_old, leave_sec, rhoAir, modeDiam, decrease_dt)!,rhopart)
        use aerosoldef, only : chemistryIndex
        implicit none

        real(r8),dimension(gas_pcnst),intent(inout)         :: massDistrib
        real(r8), dimension(secNrBins),intent(in)           :: numberConc_old
        real(r8), dimension(secNrSpec), intent(out)         :: leave_sec
        logical, intent(out)                                :: decrease_dt
        real(r8)                                            :: rhoAir
        real(r8)                                            :: modeDiam


        real(r8), dimension(secNrSpec, secNrBins)           :: numberConc_new
        real(r8), dimension( secNrBins)                     :: volume
        real(r8), dimension( secNrBins)                     :: volume_new
        integer                                             :: indBin,indVol
        real(r8), parameter                                 :: pi = 3.141592654_r8
        real(r8)                                            :: xfrac
        real(r8),dimension(secNrSpec,secNrBins)             :: volfrac  !
        real(r8)                                            ::sumUp
        real(r8),dimension(secNrSpec)                       ::moveMass
        real(r8), dimension(secNrSpec)                      :: volume_end_new


        decrease_dt=.FALSE.
        !compute volume in each bin with condensation (by mass) and by
        !numberconcentration 
        do indBin = 1, secNrBins
                volume_new(indBin) = 0.0_r8
                volfrac(:,indBin)=0.0_r8
                do indVol = 1, secNrSpec! calculate volume in each bin by mass/density! m3
                        if (numberConc_old(indBin)<1.e-30_r8) then
                                volume_new(indBin)=0.0_r8
                        else
                            volume_new(indBin) = volume_new(indBin) + massDistrib(chemistryIndex(secConstIndex(indVol,indBin)))/&
                                    rhopart_sec(indVol) * rhoAir/&
                                    (numberConc_old(indBin))
                        end if

                        !if (numberConc_old(indBin)<1.e-30_r8) then
                        !        volume_new(indBin)=0.0_r8
                        !end if
                        volfrac(indVol,indBin)=massDistrib(chemistryIndex(secConstIndex(indVol,indBin)))/&
                                rhopart_sec(indVol)*rhoAir
                                !kg/kg(air)*[kg(air)/m3(air)][kg/m3]--> m3/m3(air)
                        !write(*,*) 'SMB: volume_new',indBin, volume_new(indBin)
                        !write(*,*) 'SMB: rhopart', rhopart_sec(indVol)!volume_new(indBin)
                end do ! calculate volume in each bin by numberconcentration (volume from before condenstion) 
                volfrac(:,indBin)=volfrac(:,indBin)/(sum(volfrac(:,indBin))+1.E-50_r8)
                if (sum(volfrac(:,indBin))<1.e-50_r8) then
                        volfrac(:,indBin)=0.0_r8
                end if
                !write(*,*) 'volfrac', volfrac
                !do indVol = 1, secNrSpec! calculate volume in each bin by mass/density! m3
                volume(indBin) =  pi * secMeanD(indBin)**3/6._r8 
                !write(*,*) 'SMB: volume_old',indBin, volume(indBin)
                !volume(indBin) = numberConc_old(indBin) * pi * secMeanD(indBin)**3/6._r8 
                !end do ! calculate volume in each bin by numberconcentration (volume from before condenstion) 
        end do
        numberConc_new(:,:) = 0._r8
        !numberConc_new(:,secNrBins) = 0.0_r8!volfrac(:,secNrBins)*numberConc_old(secNrBins) 
        do indBin =  1, secNrBins-1
                ! fraction to stay in bin
                xfrac=(volume(indBin+1)-volume_new(indBin)) &
                                /(volume(indBin+1)-volume(indBin))
                if (numberConc_old(indBin)<1.e-30) then
                        xfrac=1.0_r8
                end if
                if (xfrac .le. 0._r8) then      ! if the fraction to stay is equal to
                                                ! less than zero, then the
                                                ! aerosols have grown too large
                                                ! for the next bin and we will
                                                ! want to decrease the time step
                                                ! to avoid this. 
                        decrease_dt=.TRUE.
                end if

                if (xfrac .le. 0._r8) then
                        decrease_dt=.TRUE.
                end if
                xfrac=max(0._r8, min(1._r8,xfrac))
                !nwrite(*,*) 'SMB: xfrac',xfrac
                do indVol= 1, secNrSpec
                        numberConc_new(indVol, indBin) = numberConc_new(indVol, indBin) + &
                                        xfrac*numberConc_old(indBin) &
                                        *volfrac(indVol,indBin)
                        numberConc_new(indVol, indBin+1) = numberConc_new(indVol, indBin+1) + &
                                        (1-xfrac)*numberConc_old(indBin)   &
                                        *volfrac(indVol,indBin)


                end do
                
        end do

        xfrac = (max_diameter**3 * pi/6.0_r8 - volume_new(secNrBins)) &
                                /(max_diameter**3*pi/6.0_r8-volume(secNrBins))

        if (xfrac .le. 0._r8) then
                decrease_dt=.TRUE.
        end if
        
        xfrac=max(0._r8, min(1._r8,xfrac))
        
        do indVol=1, secNrSpec
                numberConc_new(indVol,secNrBins) = numberConc_new(indVol,secNrBins) + &
                                                xfrac * numberConc_old(secNrBins) &
                                                * volfrac(indVol,secNrBins)
                leave_sec(indVol) = & !massDistrib(chemistryIndex(secConstIndex(indVol, secNrBins)))*(1-xfrac) 
                        pi * max_diameter**3 / 6.0_r8 * rhopart_sec(indVol)/rhoAir &   ! [m3_aer/#]*[kg_aer/m3_aer]/[kg_air/m3_air]--> [kg_aer/kg_air/#][m3_air]
                                                * (1-xfrac) * numberConc_old(secNrBins) &          ! *[#/m3_air] --> kg_aer/kg_air
                                                * volfrac(indVol,secNrBins)
        end do 
        do indBin=1,secNrBins
                do indVol=1,secNrSpec !Assume
                        massDistrib(chemistryIndex(secConstIndex(indVol,indBin))) = &! &!massDistrib(secConstIndex(indVol,indBin))+&
                                rhopart_sec(indVol)/rhoAir &!* massfrac(indVol,indBin)* numberConc_new(indBin)! &
                                * numberConc_new(indVol, indBin) * pi * secMeanD(indBin)**3/6.0_r8 !&
                end do
        end do

        !bulk++


end subroutine sec_moveMass





end module aero_sectional

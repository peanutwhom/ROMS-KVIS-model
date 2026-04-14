#include "cppdefs.h"

#if defined NONLINEAR && (defined LMD_SKPP || defined SOLAR_SOURCE) && \
    defined SOLVE3D

     SUBROUTINE lmd_swfrac_tile (ng, tile,                              &
     &                            LBi, UBi, LBj, UBj,k, KVIS_3d,          &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Zscale, Z,swdk)
!
! SVN $Id: lmd_swfrac.F 995 2020-01-10 04:01:28Z arango $
!------------------------------------------------- Hernan G. Arango ---
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                     !
!  This subroutine computes the fraction of solar shortwave flux      !
!  penetrating to specified depth due to exponential decay in         !
!  Jerlov water type. Also calculates solar elevation angle.          !
!                                                                     !
!=======================================================================

      USE mod_param
      USE mod_mixing
      USE mod_scalars
      USE mod_iounits
      USE mod_grid
     ! USE mod_forces
      ! USE mod_ocean
     ! USE mod_stepping
     ! USE mod_ncparam
      implicit none

!  Imported variable declarations.
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj,k
      real(r8), intent(in) :: KVIS_3d(IminS:ImaxS,JminS:JmaxS)
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      real(r8), intent(in) :: Zscale

      real(r8), intent(in) ::  Z(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: swdk(IminS:ImaxS,JminS:JmaxS)


!  Local variable declarations.
      integer :: i, j

      real(r8) :: F_IR = 0.53
      real(r8) :: F_VIS 
      
      real(r8) :: Tr_IR, Tr_VIS,K_IR
      real(r8) :: fac1,fac2,fac3,fac4
      real(r8) :: theta_a,theta_water
  !    real(r8), pointer :: lat
  !    real(r8), parameter :: attenuation_length = 10.0_r8 

# include "set_bounds.h"
!-----------------------------------------------------------------------
!  Loop over each grid
!-----------------------------------------------------------------------
    F_VIS = 1-F_IR
    DO j = Jstr, Jend
        DO i = Istr, Iend
            theta_a = calculate_solar_elevation(ng,GRID(ng)%latr(i,j))
            !fac4 = 0.27+0.0518*(KVIS_3d(i,j)**0.428)
            Tr_VIS = EXP(KVIS_3d(i,j)*Z(i,j)*Zscale)

            fac3 = (0.001_r8 + ABS(Z(i,j) * Zscale))**0.65_r8
            fac2 = 0.560_r8 + (2.304_r8 / fac3)
            K_IR = fac2 * (1.0_r8 + 0.002_r8 * theta_a)
            Tr_IR = EXP(K_IR * Z(i,j) * Zscale)

            swdk(i,j) = F_VIS * Tr_VIS + F_IR * Tr_IR
           !   swdk(i,j) = 1
              IF (i == 3 .AND. j == 3 .AND. k == 10) THEN
              PRINT *, 'Z(i,j): ', Z(i,j)
              PRINT *, 'KVIS_3d', KVIS_3d(i,j)
              PRINT *, 'lat', GRID(ng)%latr(i,j)
              PRINT *, 'swdk',swdk(i,j)
              END IF
         END DO
      END DO
    
 
      RETURN
      
      CONTAINS

      FUNCTION calculate_solar_elevation(ng,lat) RESULT(zenith)
     
        USE mod_scalars   
        IMPLICIT NONE

        integer, intent(in) :: ng
        real(r8), intent(in) :: lat 
        real(r8) :: theta_a,zenith,theta_water                      

        real(r8) :: julian_day               
        real(r8) :: n_days                  
        real(r8) :: zenith_angle             
        real(r8) :: hour_angle,current_hour
        real(r8), parameter :: dtor = pi / 180.0_r8
        real(r8), parameter :: rtod = 180.0_r8 / pi
        real(r8), parameter :: n_air = 1.0003_r8     
        real(r8), parameter :: n_water = 1.333_r8    
        
        julian_day = time_ref + tdays(ng)
       
        n_days = MOD(julian_day -1.0_r8, 365.0_r8)

        zenith_angle = 23.45_r8 * SIN(2.0 * pi * (n_days - 81) / 365)

        current_hour = MOD(tdays(ng),1.0_r8) * 24.0_r8

        hour_angle = 15.0_r8 * (current_hour - 12.0_r8)

        theta_a = ASIN( SIN(dtor * lat) * SIN(dtor * zenith_angle) +    &
      &      COS(dtor * lat) * COS(dtor * zenith_angle) *               &
       &     COS(dtor * hour_angle) )

       theta_a = theta_a * rtod
    
       theta_a = MAX(theta_a, 0.0_r8)
       zenith = 90-theta_a  

      END FUNCTION calculate_solar_elevation

      END SUBROUTINE lmd_swfrac_tile
#else
      SUBROUTINE lmd_swfrac
      RETURN
      END SUBROUTINE lmd_swfrac
#endif



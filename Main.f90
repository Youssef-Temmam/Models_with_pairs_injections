PROGRAM InjectionCalculation

!----------------------------------------------------------------------------------------------------------------------------
!///////////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!------------------------------------------------------ DESCRIPTION GENERAL -------------------------------------------------
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////////////////////////
!----------------------------------------------------------------------------------------------------------------------------

! This program is built in the aim of calculate the pairs injection term produce by photon-photon <--> e_+ + e_- interraction.
! The energetics photon are emmited by an accretion disk located in the equatorial plane. The ambiant medium is approximated 
! as transparent and pairs production negligible in comparison of the photons numbers. Following these assumptions the photons 
! distribution function is transported along the geodesics (Eq.279a).  

! The code could be split in three par : 
!    - One part located in GeodesicsRoutine.f90 and BesselFunction.f90 are dedicated to calculation related to geodesics 
!    trajectory which crosses the axis
!    - a second part located in different routine modules are dedicated to the disk emission.
!    - the last part is located in Injection.f90 calculate the collision integrals of pairs production (Use open MPI).

! Other part of the programm could be used for other aims (In and Out data in EntreSortie.f90 / gestion and calculation of 
! public variables in PublicVariable.f90 / RoutineOthers.f90 for all the rest).

! export OMP_NUM_THREADS=5
!EXECUTION : ./Main 



!EXECUTION With terminal log: ./Main < Trm_CommonReading.in
  
!----------------------------------------------------------------------------------------------------------------------------

!TOUTES LES GRANDEURS SONT ADIMENSIONNES

!----------------------------------------------------------------------------------------------------------------------------
!///////////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!--------------------------------------------------- ROUTINES ET PROGRAMMES -------------------------------------------------
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////////////////////////
!----------------------------------------------------------------------------------------------------------------------------

USE EntreSortie             !Entree et sortie
USE PublicVariable          !Variable public
USE BesselFunction          !Bessels function
USE Routine_other           !Module containing useful routines
USE GeodesicsRoutine        !Modules link to geodesics calculations
USE IntegralsCalculations   !Integration modules
USE NarayanYi               !Quantities link to accretion disk of Narayan-Yi 1995
! USE NovikovThorne         !Quantities link to accretion disk of Novikov Thorne 1963
USE KinoKaburaki          !Quantities link to accretion disk of Kino and Kaburacki 2000
USE Injection               !Modules link to injection collision integrals calculation

!----------------------------------------------------------------------------------------------------------------------------
!///////////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!------------------------------------------------------ DECLARATIONS --------------------------------------------------------
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////////////////////////
!----------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
REAL(KIND=wp)                                               :: total_time          ! Temps de calculs
INTEGER                                                     :: count_time_av

REAL(KIND=wp)                                               :: time_compare
INTEGER                                                     :: count_compare
! REAL(KIND=wp),DIMENSION(:),ALLOCATABLE                      :: ra_tab

REAL(KIND=wp)                                               :: r_ax,psi,R_eq
REAL(KIND=wp),DIMENSION(2)                                  :: PsiB,Inj_r,psi_bounds_disk

REAL(KIND=wp)                                               :: Integ,ra,psim,psip
INTEGER                                                     :: j,n_try
REAL(KIND=wp)                                               :: tau_num,tau_den,tau_0
!----------------------------------------------------------------------------------------------------------------------------
!///////////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////////////////////////
!----------------------------------------------------------------------------------------------------------------------------

CALL SYSTEM_CLOCK(count_rate=rate)
CALL SYSTEM_CLOCK(count_time_av)

! ------ Terminal output or no
WRITE(*,'(a)',advance='no')' Terminal output (T/F) : '
READ(*,*) terminal_log
WRITE(n_log,*)
IF (terminal_log) THEN
    n_log = 0
ELSE 
    n_log = 10
    OPEN(n_log, file = "LogFile/basics_output.log", status ='replace') 
END IF

! ----- Print the header
CALL HEADER()

! ------ Read parameter
CALL READ_BH_PARAM()
CALL READ_INJECTIONPROFILE()
CALL READ_DISK_PARAM()

! ------ Create output folder of quantiest related to geodesics
CALL OUTPUTFOLDER_GEOD(path_geod,sol_folder_name_geod)

! ------ Create output folder of quantiest related to disk type
CALL OUTPUTFOLDER_DISK(path_disk,sol_folder_name_disk)

! ------ Test geodesics routine R_eq/r_a for different value of r_a
ra_tab = Stretched_tab(rh+drh_prof,r_ax_max_prof,drh_prof,n_ax_prof)
CALL TEST_Reqsra(n_ax_prof,ra_tab,n_psi)

! Comparison between ordinary calculation and interpolation
! ra  = 3.496_wp
! psi = 2.021_wp
! n_try = 10000


! ------- Test spectrum profil on the axis.
IF (kind_emmit_prof .EQ. 1) THEN
    ! ----- Test NY N spectrum profil on the axis
    CALL ANALYSIS_NY_N()
ELSE IF (kind_emmit_prof .EQ. 3) THEN
    ! ----- Test KK uncomptonized spectrum profil on the axis
    CALL ANALYSIS_KK_U()
ENDIF

! print*, "ee_brem_param",ee_brem_param

! tau_0 = F_adim_disk_kk_u(10.4_wp,10._wp**(-9))
! print*, "value F_adim = ",tau_0
! tau_0 = F_adim_disk_kk_u(10.4_wp,0.85_wp)
! print*, "value F_adim = ",tau_0

! ! ------- InjectionOnePoint(r,kind_emmit,Rm,Rp,Theta_0,rel_err)

! Inj_r = InjectionOnePoint(10._wp*rh,kind_emmit_prof,R_isco,R_ext_ny_n,1._wp,0.05_wp)


! ----- Final time
CALL PRINT_TIME(count_time_av,.TRUE.,'Total time             :',total_time)
! -----
CLOSE(n_log)

END PROGRAM InjectionCalculation



! ! ---- Test integrals

! Integ = GAUSS_LEGENDRE_1D_8(cube,0._wp,1._wp)
! PRINT*,"Integ",Integ
! Integ = GAUSS_LEGENDRE_1D_16(cube,0._wp,1._wp)
! PRINT*,"Integ",Integ
! Integ = GAUSS_LEGENDRE_1D(cube,0._wp,1._wp,0.00000001_wp)
! PRINT*,"Integ",Integ

! Integ = GAUSS_LEGENDRE_api_8(gauss,0._wp)
! PRINT*,"Integ_gauss",Integ
! Integ = GAUSS_LEGENDRE_api_16(gauss,0._wp)
! PRINT*,"Integ_gauss",Integ
! Integ = GAUSS_LEGENDRE_api(gauss,0._wp,)
! PRINT*,"Integ_gauss",Integ,abs(Integ-0.5_wp*sqrt(pi_circ))
! Integ = GAUSS_LEGENDRE_ra_X_psi_psim_psip_d_psi_8(test_int,1._wp,0.5_wp,1._wp,2._wp)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*exp(-0.25_wp)*31._wp/5._wp
! Integ = GAUSS_LEGENDRE_ra_X_psi_psim_psip_d_psi_16(test_int,1._wp,0.5_wp,1._wp,2._wp)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*exp(-0.25_wp)*31._wp/5._wp
! Integ = GAUSS_LEGENDRE_ra_X_psi_psim_psip_d_psi(test_int,1._wp,0.5_wp,1._wp,2._wp,0.000000001_wp)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*exp(-0.25_wp)*31._wp/5._wp

! Integ = GAUSS_LEGENDRE_ra_X_api_dX_8(test_int1,1._wp,2._wp,1.0_wp,0._wp)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*sqrt(pi_circ)/2._wp
! Integ = GAUSS_LEGENDRE_ra_X_api_dX_16(test_int1,1._wp,2._wp,1.0_wp,0._wp)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*sqrt(pi_circ)/2._wp
! Integ = GAUSS_LEGENDRE_ra_X_api_dX_32(test_int1,1._wp,2._wp,1.0_wp,0._wp)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*sqrt(pi_circ)/2._wp
! Integ = GAUSS_LEGENDRE_ra_X_api_dX_64(test_int1,1._wp,2._wp,1.0_wp,0._wp)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*sqrt(pi_circ)/2._wp
! Integ = GAUSS_LEGENDRE_ra_X_api_dX_127(test_int1,1._wp,2._wp,1.0_wp,3._wp,0._wp)
! ! GAUSS_LEGENDRE_ra_X_api_dX_127(f,psim,psip,ra,k,Xm)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*sqrt(pi_circ)/2._wp
! Integ = GAUSS_LEGENDRE_ra_X_api_dX(test_int1,1._wp,2._wp,1.0_wp,0._wp,0.00000001_wp)
! PRINT*,"Integ_gauss",Integ,cos(1._wp)*sqrt(pi_circ)/2._wp
! Integ = Integrand_Flux_axis_X_Psi_ny_n(3.44928769219581947392_wp,0.292984500595344155623_wp&
!     &,1.57079632228393202792_wp,1.57079634016717765803_wp,1.57079632669304705228_wp)
! PRINT*,"Integ",Integ

! ra = 3.44928769219581947392_wp
! !PsiB    = PsiBounds(ra)
! psim    = 0.993396787406716099355_wp !PsiFromRad(ra,R_ext_ny_n,PsiB(1),PsiB(2),PsiB)    ! In ny_n models the exterior bounds of the disk is an entered parameter
! psip    = 1.59345893684167663447_wp  ! PsiFromRad(ra,R_isco,PsiB(1),PsiB(2),PsiB)        ! In ny_n models the interior bounds of the disk is ISCO radius
! psi_bounds_disk     = [psim,psip]
! R_eq = RadiusEqu(3.44928769219581947392_wp,1.57079632669304705228_wp,psi_bounds_disk)

! !print*,psiB
! print*,psi_bounds_disk

! PsiB    = PsiBounds(ra)
! CALL SYSTEM_CLOCK(count_compare)
! DO j=1,n_try
!     R_eq = RadiusEqu(ra,psi,PsiB)
! ENDDO
! CALL PRINT_TIME(count_compare,.TRUE.,'Ordinary routine       :',time_compare)
! PRINT*,'R_eq = ',R_eq
! PRINT*,'Full time =', time_compare

! CALL ACTUALISE_REQ_TAB_PUBLIC(ra)
! CALL SYSTEM_CLOCK(count_compare)
! DO j=1,n_try
!     R_eq = RadiusEqu_interp(ra,psi)
! ENDDO
! CALL PRINT_TIME(count_compare,.TRUE.,'Interpol routine       :',time_compare)
! PRINT*,'R_eq = ',R_eq
! PRINT*,'Full time =', time_compare

! ! PRINT*,Req_tab_public
! ! PRINT*,psi_tab_public


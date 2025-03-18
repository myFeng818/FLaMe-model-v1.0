!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------
      PROGRAM RUNLAKE
C
C=======================================================================
C     * THE CANADIAN SMALL LAKE MODEL STANDALONE DRIVER
C     *
C     * AUTHOR: MURRAY MACKAY 		JAN 2014
C     *
C=======================================================================
C     * FLaMe v1.0 model built on Canadian small lake model
C     * 
C     * Main program of FLaMe         
C     *
C     * AUTHOR: Manon Maisonnier, Maoyuan Feng, ... ,Pierre Regnier.         
C     * Final revision: March 2025
C     * 
C     * Reference: A new biogeochemical modelling framework (FLaMe v1.0) 
C     *            for lake methane emissions on the regional scale: 
C     *            Development and application to the European domain
C     *         
C     * Correspondance: maoyuan.feng@ulb.be; maoyuanfeng93@gmail.com         
C=======================================================================
C     * Primary subroutines of FLaMe biogeochemical modules
C     *         Clabile_module.f: Organic carbon subroutine
C     *         CH4_module.f:     Methane production and emission subroutine
C     *         O2_module.f:      Oyxgen dynamics subroutine
C     *         oxidation.f:      Oxidation of methane in water column
C     *         transport.f:      Diffusive transport of gases within water         
C     *         geometry_bgc.f:   Setting of geometry of biogeochemical modules
C     *         search_zebmin.f:  Solving the threshold sediment depth for the
C     *                           split between diffusion and ebullition  
C     *         Writing.f:        Writing the key fluxes and stocks into output
C                                 files
C     *              
C     *         CLASSL.f: a key module from CSLM that calls the biogeochemical
C     *                   subroutines
C     *         
C     *         The other subroutines are from CSLM              
C=======================================================================              

      IMPLICIT NONE
C

      CHARACTER(41) ,PARAMETER :: GEN_DIR = 
     1 "/Users/fengmaoyuan/Documents/Maoyuan_ULB/"
        ! "D:/"
                                                          
      CHARACTER(33),PARAMETER :: FILE_DIR =
     1 "current_project/FLaMe/Europe_test"

     
      CHARACTER(len=*),PARAMETER ::SET="/variations_site_1/"   ! variations/        set_Rinta/Boreal/
 
      CHARACTER(len=*),PARAMETER :: OP=
     1  "output/"      !base/



      INTEGER, PARAMETER :: x_per_m = 4    ! 1 for one output by month, 4 for 4 output by month
 
      INTEGER :: DELTAT_READING
      INTEGER :: N_READING,N_step_met,NSTEP_write
      
      INTEGER,PARAMETER :: NLAT=4,NMOS=1,ILG=NLAT*NMOS    !MM NMOS =3 originaly   ! europe : 28 (29)  NLAT=28 (2/5/23)
      INTEGER,PARAMETER :: NLAKMAX=410
      INTEGER,PARAMETER :: NBS=4,IGL=1,IRSTRT=0
      INTEGER :: MIDROT (NLAT,NMOS)
      INTEGER :: N,NLTEST,NMTEST,I,J,K,L,M,JL1,JL2,
     1        NCOUNT,NDAY,IHOUR,IMIN,IDAY,IYEAR,NML_,NMW,
     2        IDISP,IZREF,ISLFD,IPCP,ITG,IALS,ISNOALB
      INTEGER, DIMENSION(ILG) :: NMTEST_I     ! NMTEST -> DIMENSION    NMTEST only is max NMTEST_I
C
C-- GATHER-SCATTER INDEX ARRAYS
      INTEGER,DIMENSION(ILG) :: ILMOS,JLMOS,IWMOS,JWMOS
C
C-- LAKE TILE PARAMETERS                                            
      INTEGER NLAKGAT(ILG), NLAKROT(NLAT,NMOS)
      REAL,   DIMENSION(NLAT,NMOS) :: HLAKROT, LLAKROT, BLAKROT,    
     1           HFSLROT, HEVLROT, FSGLROT, FLGLROT, HMFLROT,
     2           ASVDROT, ASIDROT,REFROT,BCSNROT,limP_ROT
      REAL,   DIMENSION(ILG) :: HLAKGAT, LLAKGAT, BLAKGAT   
      REAL,   DIMENSION(ILG) :: ASVLGAT, ASILGAT,BCSNLGAT,REFLGAT,
     1           ZDMLGAT,ZDHLGAT
C-- LAKE PROGNOSTIC VARIABLES                                       
      REAL,   DIMENSION(NLAT,NMOS,NLAKMAX) :: TLAKROT               
      REAL,   DIMENSION(ILG,NLAKMAX) :: TLAKGAT                     
      REAL,   DIMENSION(ILG) :: EXPW,DTEMP,HDPTH,DELU,GRED,         
     1                    TKELAK,T0LAK,LKICEH,RHOMIX,TSED,SNICEH,ROFICEH
C-- LAKE DIAGNOSTIC VARIABLES                                       
      REAL,   DIMENSION(ILG) :: HFSLGAT, HEVLGAT, FSGLGAT, FLGLGAT, 
     1                       HMFLGAT,HTCLGAT,FICEGAT,FLSGAT,G0SLGAT,
     2                       FSGSLGAT,FLGSLGAT,HFSSLGAT,HEVSLGAT,
     3                       HMFNLGAT,HTCSLGAT
      REAL,   DIMENSION(ILG) :: PCPLGAT,PCPNLGAT,QFLGAT,QFNLGAT,
     1                       ROFNLGAT,SFTLGAT,SFULGAT,SFVLGAT,SFQLGAT,
     2                       SFHLGAT,QLWOLGAT,ALVLGAT,ALILGAT,EFLGAT,
     3                       GTLGAT,QGLGAT,DRLGAT,PETLGAT,QSENLGAT,
     4                       TFXLGAT,QEVPLGAT,QFSLGAT,QFXLGAT,SNOLGAT,
     5                       RHOSLGAT,TSNOLGAT,ALBSLGAT,WSNOLGAT
      REAL,DIMENSION(NLAT,NMOS) :: PCPLROT,PCPNLROT,QFLROT,QFNLROT,
     1                       ROFNLROT,FSGSLROT,FLGSLROT,HFSSLROT,
     2                       HEVSLROT,HMFNLROT,HTCSLROT,HTCLROT,
     3                       SFTLROT,SFULROT,SFVLROT,SFQLROT,SFHLROT,
     4                       QLWOLROT,ALVLROT,ALILROT,EFLROT,GTLROT,
     5                       QGLROT,DRLROT,PETLROT,QSENLROT,TFXLROT,
     6                       QEVPLROT,QFSLROT,QFXLROT,SNOLROT,RHOSLROT,
     7                       TSNOLROT,ALBSLROT,WSNOLROT
CC    INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(15,307)
      INTEGER, PARAMETER :: SP=kind(1.0)
      INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(2*precision(1.0_sp))
      REAL(dp),   DIMENSION(ILG) :: CTLSTP 
C
      REAL,DIMENSION(ILG,NBS) :: 
     1        FSDBLGAT,   FSFBLGAT,   FSSBLGAT
      REAL,DIMENSION(NLAT,NBS) ::
     1        FSDBROL,   FSFBROL,   FSSBROL
C
C--   ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES
C
      REAL,DIMENSION(NLAT,NMOS) :: FAREROT
      REAL,DIMENSION(NLAT) ::
     1      ZRFMROW,   ZRFHROW,   ZDMROW ,   ZDHROW ,  
     2      ZBLDROW,   FSVHROW,   FSIHROW,   RADJROW,
     3      CSZROW ,   FDLROW ,   ULROW  ,   VLROW  ,   
     4      TAROW  ,   QAROW  ,   PRESROW,   PREROW ,  
     5      PADRROW,   VPDROW ,   TADPROW,   RHOAROW,  
     6      UVROW  ,   GCROW  ,   RPCPROW,   TRPCROW,
     7      SPCPROW,   TSPCROW,   RHSIROW,   RPREROW,
     8      SPREROW  

      REAL,DIMENSION(ILG) :: ATM_CH4
 
C
      REAL,DIMENSION(ILG) ::
     1      ZRFMGAT,   ZRFHGAT,   ZDMGAT ,   ZDHGAT ,  
     2      ZBLDGAT,   FSVHGAT,   FSIHGAT,   RADJGAT,
     3      CSZGAT ,   FDLGAT ,   ULGAT  ,   VLGAT  ,   
     4      TAGAT  ,   QAGAT  ,   PRESGAT,   PREGAT ,  
     5      PADRGAT,   VPDGAT ,   TADPGAT,   RHOAGAT,
     6      RPCPLGAT,  TRPCLGAT,  SPCPLGAT,  TSPCLGAT,
     7      RHSILGAT,  RADJLGAT,  PADRLGAT
C
C-- LAND SURFACE DIAGNOSTIC VARIABLES
C
      REAL,DIMENSION(NLAT,NMOS) :: CDHROT ,   CDMROT 
      REAL,DIMENSION(ILG) ::       CDHGAT ,   CDMGAT 
C
C-- ARRAYS USED FOR OUTPUT AND DISPLAY PURPOSES.
C
      CHARACTER     TITLE1*4,     TITLE2*4,     TITLE3*4,
     1              TITLE4*4,     TITLE5*4,     TITLE6*4
      CHARACTER     NAME1*4,      NAME2*4,      NAME3*4,
     1              NAME4*4,      NAME5*4,      NAME6*4
      CHARACTER     PLACE1*4,     PLACE2*4,     PLACE3*4,
     1              PLACE4*4,     PLACE5*4,     PLACE6*4
C 
C-- LOCAL CONSTANTS AND VARIABLES.
C
      REAL FSDOWN,DAY,DECL,HOUR,COSZ,          ! M.M. errase DEGLAT and DEGLON (see below)
     1     EA,CA,CB,EASAT,CONST,HCAP,ZTOP,ZBOT,Z,QSUML,
     2     RHOIW,ICETOP,ICEBOT

      REAL,DIMENSION(ILG) :: DEGLAT,DEGLON    ! M.Maisonnier  (add ILG dimension to DEGLAT and DEGLON compare to original)
      REAL,DIMENSION(NLAT) :: LAT_GRID,LON_GRID
      integer :: index_start, index_end

C
C-- COMMON BLOCK PARAMETERS.
C
      REAL DELT,TFREZ, CLHVAP,PI,
     1     RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,TCSAND,TCCLAY,
     2     TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,
     3     HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,
     4     CKARM,CGRAV,CPD,DELTA,AS,ASX,ANGMAX,BS,BETA,CI,FACTN,HMIN
C
C    * LAKE TILE COMMON BLOCK PARAMETERS                            
C    *                                                              
      REAL TKECN,TKECF,TKECE,TKECS,HDPTHMIN, DUMAX, 
     1     TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,TKECL


C=======================================================================
C     * BIOCHEMICAL VARIABLES (Manon Maisoniier)

      INTEGER,PARAMETER :: NLAKMAX_bgc = 820
	  
      INTEGER, DIMENSION(ILG) :: NLAK_bgc,HLAK_bgc  !!!!!!!!!!!!!!!!!!!!!!!!!! DEPTH BGC !!!!!!
      
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: CH4,O2
      REAL, DIMENSION (ILG,0:NLAKMAX_bgc+1) :: Area
      REAL, DIMENSION (ILG) :: Clabile,FPP,FRESP,FBUR
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Clabs_1,Clabs_2,Clabs_3
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FsCH4_C_1,FsCH4_C_2,FsCH4_C_3
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Fsed_1,Fsed_2,Fsed_3
      REAL :: Pchl
      REAL, DIMENSION (ILG) :: FatmCH4,Hs,limP,ext_C
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FsCH4,FsDiff,FsEbul,zebmin
      REAL, DIMENSION (ILG) :: FatmO2,Fw
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FPPWC,FSAER,FswDiff

      REAL, DIMENSION (ILG,NLAKMAX_bgc+1) :: K_diff
      REAL, DIMENSION (ILG) :: ZPhotic,ZRESP,ZPhotic_m,ZRESP_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FRESPWC,FRESPWC_m
      
	  REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Fibs
      INTEGER, DIMENSION (ILG) :: previous_ice,ZPROD
	  
C=======================================================================
      REAL, DIMENSION (ILG) ::  G0_m,F0_m,FSGL_m,
     1     Q0_m,FLGL_m,HFSL_m,HEVL_m,T0_m,LKICEH_m,SNO_m,
     2     RHOSNI_m,S_m,R_m,OSW_m,RHOSNO_m,SNICEH_m
     
      REAL, DIMENSION (ILG,NLAKMAX) :: TLAK_m
      REAL, DIMENSION (ILG) :: FQU_m,BFLX_m,DISS_m,TKE_m,DELU_m,HDPTH_m
      INTEGER :: JMIX_m
      REAL :: ZMIX_m,TMIX_m,DELTHRM_m
      REAL, DIMENSION (ILG) :: DTEMP_m,
     2     FSHEAR_m,FENTRA_m,USTAR_m,EXPW_m,WEDB_m,
     1     ZREFM_m,VA_m,CDML_m,TSNOW_m,
     1     WSNOW_m,ZSNOW_m,FLS_m,Hs_m,
     1     QSWNS_m,QTRANSL_m,QLWIN_m,QLWOS_m,QSENSS_m,QEVAPS_m,
     2     QMELTS_m,GZEROSL_m,CFLUXS_m,QSWIN_m,TSED_m
      REAL, DIMENSION (ILG,NLAKMAX) :: QFLX_m, FFLX_m
	  
      REAL, DIMENSION (ILG,NLAKMAX_bgc) ::
     1     O2_m,CH4_m,
     2     FPPWC_m,FSAER_m,
     6     FsCH4_C_m,FsDiff_m,
     1     FsEbul_m,zebmin_m
      REAL, DIMENSION (ILG) :: Clabile_m,FPP_m,FRESP_m,FatmO2_m,
     1                         FatmCH4_m,Fw_m,FBUR_m

      REAL :: DELT_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: 
     1    Clabs_1_m,Clabs_2_m,Clabs_3_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FsCH4_C_1_m,FsCH4_C_2_m,
     1                                 FsCH4_C_3_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Fsed_1_m,Fsed_2_m,Fsed_3_m      
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FswDiff_m,OxCH4_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FatmEbul , FatmEbul_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc+1) :: K_diff_m

      REAL,DIMENSION (ILG,NLAKMAX_bgc) :: TLAK_bgc_m
      REAL,DIMENSION (ILG) :: HDPTH_bgc_m
      
      REAL :: year_old

      
      REAL, DIMENSION(ILG) :: Fstor_m,DELT_stor_m

      REAL,DIMENSION (ILG) :: DELT_ice_m

      REAL, PARAMETER :: KsP = 1.*30.97d-9*1.0d6*3
C=======================================================================
C     * PHYSICAL CONSTANTS.
C     * THE FOLLOWING COMMON BLOCKS ARE DEFINED SPECIFICALLY FOR USE 
C     * IN CLASS, AND ARE RETAINED HERE TO MAINTAIN COMPATIBILITY
C
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
      COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
C
C    * LAKE TILE COMMON BLOCK                                       
C    *                                                              
      COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,            
     2             TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,  
     3             TKECL,DUMAX                                      
C
C     * ADDITIONAL VALUES FOR RPN AND GCM COMMON BLOCKS.
C
      COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX              
C
      DATA      DELTA,      AS,         ASX,        ANGMAX
     1       /  0.608,      12.0,       4.7,        0.85   /

      DATA      CI,         BS,         BETA,       FACTN,      HMIN 
     1       /  40.0,       1.0,        1.0,        1.2,        40.   /
C
C    * PARAMETERS ORIGINALLY FROM CLASS COMMON BLOCK DATA           
C    *                                                              
      DATA      VKC,        CT,         VMIN,     EMSW
     1       /  0.40,       1.15E-3,    0.1,      0.97     /
      DATA      TCW,        TCICE,      TCSAND,    TCCLAY,  TCOM
     1       /  0.57,       2.24,       2.5,       2.5,     0.25  /
      DATA      TCDRYS,     RHOSOL,     RHOOM
     1       /  0.275,      2.65E3,     1.30E3  /    
      DATA      HCPW,       HCPICE,     HCPSOL,    HCPOM
     1       /  4.187E6,    1.9257E6,   2.25E6,    2.50E6 /
      DATA      HCPSND,     HCPCLY,     SPHW,      SPHICE,  SPHVEG
     1       /  2.13E6,     2.38E6,     4.186E3,   2.10E3,  2.70E3 /
      DATA      RHOW,       RHOICE,     TCGLAC,    CLHMLT,  CLHVAP
     1       /  1.0E3,      0.917E3,    2.24,      0.334E6, 2.501E6/
C
C--Lake TKE process efficiencies                                    
      DATA  TKECN,      TKECF,      TKECE,      TKECS,      TKECL   
     1/     1.33,       0.25,       1.15,       0.20,      0.2350/ 
CToo?1/     1.33,       0.25,       1.15,       0.20,      0.2250/
CELA 1/     1.33,       0.25,       1.15,       0.20,      0.2350/
C
C--Lake process parameter limits, and grid spacing                  
      DATA  HDPTHMIN, TKEMIN,  DELMAX,  DELMIN,  DELZLK,  DELSKIN   
     1/     0.5,      1.0E-12, 5.0,     0.5,     0.5,     0.050 /   
      DATA  DHMAX,    DUMAX                                         
     1/     2.0,      0.1   /                                       
C
C     * ASSIGN VALUES NORMALLY SPECIFIED WITHIN THE GCM.
C
      DATA      PI   /  3.1415926535898    /
      DATA      GRAV,       RGAS,       SPHAIR,    RGASV
     1/         9.80616,    287.04,     1.00464E3,  461.50   /
      DATA      TFREZ,      SBC
     1/         273.16,     5.66796E-8  /
      DATA      DELT
     1/         600.0  /
C     1/         600.0   /       ! 10 minutes (eg. Raksjon, ELA)
C    1/         300.0   /       ! 5 minutes (eg. Toolik)
C    1/         1800.0  /       ! half hrly
      CGRAV=GRAV
      CKARM=VKC
      CPD=SPHAIR
C

C      
C=======================================================================
C     * CLASS SWITCHES select different runtime options

C     * IF IDISP=0, VEGETATION DISPLACEMENT HEIGHTS ARE IGNORED,
C     * BECAUSE THE ATMOSPHERIC MODEL CONSIDERS THESE TO BE PART
C     * OF THE "TERRAIN".
C     * IF IDISP=1, VEGETATION DISPLACEMENT HEIGHTS ARE CALCULATED.

C     * IF IZREF=1, THE BOTTOM OF THE ATMOSPHERIC MODEL IS TAKEN
C     * TO LIE AT THE GROUND SURFACE.
C     * IF IZREF=2, THE BOTTOM OF THE ATMOSPHERIC MODEL IS TAKEN
C     * TO LIE AT THE LOCAL ROUGHNESS HEIGHT.

C     * IF ISLFD=0, DRCOEF IS CALLED FOR SURFACE STABILITY CORRECTIONS
C     * AND THE ORIGINAL GCM SET OF SCREEN-LEVEL DIAGNOSTIC CALCULATIONS 
C     * IS DONE.
C     * IF ISLFD=1, DRCOEF IS CALLED FOR SURFACE STABILITY CORRECTIONS
C     * AND SLDIAG IS CALLED FOR SCREEN-LEVEL DIAGNOSTIC CALCULATIONS. 
C     * IF ISLFD=2, FLXSURFZ IS CALLED FOR SURFACE STABILITY CORRECTIONS
C     * AND DIASURF IS CALLED FOR SCREEN-LEVEL DIAGNOSTIC CALCULATIONS. 

C     * IF IPCP=1, THE RAINFALL-SNOWFALL CUTOFF IS TAKEN TO LIE AT 0 C.
C     * IF IPCP=2, A LINEAR PARTITIONING OF PRECIPITATION BETWEEEN 
C     * RAINFALL AND SNOWFALL IS DONE BETWEEN 0 C AND 2 C.
C     * IF IPCP=3, RAINFALL AND SNOWFALL ARE PARTITIONED ACCORDING TO
C     * A POLYNOMIAL CURVE BETWEEN 0 C AND 6 C.
C     * IF IPCP=4, THE RAINFALL, SNOWFALL AND TOTAL PRECIPITATION RATES
C     * ARE READ IN DIRECTLY.

C     * ITC, ITCG AND ITG ARE SWITCHES TO CHOOSE THE ITERATION SCHEME TO
C     * BE USED IN CALCULATING THE CANOPY OR GROUND SURFACE TEMPERATURE
C     * RESPECTIVELY.  IF THE SWITCH IS SET TO 1, A BISECTION METHOD IS
C     * USED; IF TO 2, THE NEWTON-RAPHSON METHOD IS USED.
    
C     * IF IPAI, IHGT, IALC, IALS AND IALG ARE ZERO, THE VALUES OF 
C     * PLANT AREA INDEX, VEGETATION HEIGHT, CANOPY ALBEDO, SNOW ALBEDO
C     * AND SOIL ALBEDO RESPECTIVELY CALCULATED BY CLASS ARE USED.
C     * IF ANY OF THESE SWITCHES IS SET TO 1, THE VALUE OF THE
C     * CORRESPONDING PARAMETER CALCULATED BY CLASS IS OVERRIDDEN BY
C     * A USER-SUPPLIED INPUT VALUE.



      IDISP=1             !1 for stanalone with field obs.  Not used for lake model
      IZREF=1             !1 for standalone with field obs.  Used in TLSPREP 
      ISLFD=0
      IPCP=1
      ITG=1               !std
      IALS=0
      ISNOALB=0
C=======================================================================
C     * OPEN FILES FOR READING AND WRITING.

      OPEN(UNIT=50,FILE=GEN_DIR//FILE_DIR//SET//
     1     "input/lake_characteristic.ini",  !  'inp/PB_lat_4_avril_NM_4_limP.ini',
     1     STATUS='OLD')
      OPEN(UNIT=51, FILE=GEN_DIR//FILE_DIR//SET//
     1     "input/meteorological_data.met",  !  'inp/PB_lat_4_avril_yr2002to2018.met',
     1     STATUS='OLD')

      OPEN(UNIT=70,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of0")!,STATUS='OLD')
      OPEN(UNIT=71,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of1")!,STATUS='OLD')     
      OPEN(UNIT=72,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of2")!,STATUS='OLD')                 
      OPEN(UNIT=73,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of3")!,STATUS='OLD')                  
      OPEN(UNIT=74,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of4")!,STATUS='OLD')                    
      OPEN(UNIT=75,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of5")!,STATUS='OLD')
      OPEN(UNIT=76,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of6")!,STATUS='OLD')          
      OPEN(UNIT=77,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of7")!,STATUS='OLD')

      OPEN(UNIT=78,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of8")!,STATUS='OLD')
      OPEN(UNIT=79,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of9")!,STATUS='OLD')

      OPEN(UNIT=80,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of10")
      OPEN(UNIT=81,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of11")

      OPEN(UNIT=82,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of12")
      OPEN(UNIT=83,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of13")
C      OPEN(UNIT=84,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of14")
C      OPEN(UNIT=90,FILE=GEN_DIR//FILE_DIR//SET//OP//"LAKE.of20")!,STATUS='OLD')




C     * READ AND PROCESS INITIALIZATION AND BACKGROUND INFORMATION.
C     * FIRST, MODEL RUN SPECIFICATIONS.

      READ (50,5010) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
      READ (50,5010) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      READ (50,5010) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
C                                                                   
C   * LAKE DIAGNOSTIC OUTPUT FILES                                  
C                                                                   

      WRITE(70,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(70,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
      WRITE(71,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(71,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
      WRITE(72,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(72,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
      WRITE(73,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(73,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
      WRITE(74,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(74,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
      WRITE(75,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(75,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
      WRITE(76,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(76,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
      WRITE(77,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(77,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6          
      WRITE(78,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(78,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      WRITE(79,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(79,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      WRITE(80,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(80,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      WRITE(81,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
      WRITE(81,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

      
      
C

      WRITE(70,7010)                                                
7010  FORMAT('I_max (-) dt (s) DELZLK (m) NLAKMAX_bgc (-) ')        ! NLAKMAX (-)   
      WRITE(71,7011)                                                
7011  FORMAT('NLAK_bgc (-)  Area (m²)')
      WRITE(72,7012)                                                
7012  FORMAT('YEAR (yr)   DAY (DOY)    HR (hr)    MIN (min)   I (-) ')             
      WRITE(73,7013)                                                
7013  FORMAT('HDTPH (m)   LKICEH (m)    TSED (°C)')                                        
      WRITE(74,7014)                                                
7014  FORMAT('[Clabile] (g/m³)    FPP (g/m³/Dt)    FRESP (g/m³/Dt) ')                     
      WRITE(75,7015)                                                
7015  FORMAT('[Clabile1] (g/m³)    [Clabile2] (g/m³)    ', 
     1       '[Clabile3] (g/m³)    Fsed1 (g/m³/Dt)  ',
     2       'Fsed2 (g/m³/Dt)    Fsed3 (g/m³/Dt) ',
     3       'FsCH41 (g/m³/Dt)   FsCH42 (g/m³/Dt)  FsCH43 (g/m³/Dt)')              
      WRITE(76,7016)                                                
7016  FORMAT('[O2] (g/m³)    FatmO2 (g/m²/Dt)    FPPWC (g/m³/Dt)  ', 
     1       '  Fw (g/m³/Dt)    Fsaer (g/m³/Dt)   ')
      WRITE(77,7017)                                                
7017  FORMAT('FsCH4 (g/m³/Dt)   FsDiff (g/m³/Dt)  ',
     1       '   FsEbul (g/m³/Dt) ')
      WRITE(78,7018)                                                
7018  FORMAT('[CH4] (g/m³)  FatmCH4 (g/m²/Dt)  FswDiff (g/m³/Dt) ' ,
     1       '   OxCH4 (g/m³/Dt) ')

      WRITE(79,7019)                                                
7019  FORMAT('T_bgc (°C)   ')     ! NLAK  (-)

      WRITE(80,7020)                                                
7020  FORMAT('k_diff (m2/s)')

      WRITE(81,7021)                                                
7021  FORMAT('ZPhotic (m)   ZRESP (m)   FRESPWC (g/m³/Dt)')

      WRITE(82,7022)                                                
7022  FORMAT('FatmEbu (g/m2/Dt)')

      WRITE(83,7023)                                                
7023  FORMAT(' DELT_stor (s)   Fstor (g/m2/Dt)')

C      WRITE(84,7024)                                                
C7024  FORMAT(' DELT_m (s)   DELT_ice (s)')


      
C                                                                   
C=======================================================================
C   * READ GEOPHYSICAL DATA
C
        
      READ(50,5020) ZRFMROW(1),ZRFHROW(1),ZBLDROW(1),    ! New version M.Maisonnier (replace the two previous lines)
     1     GCROW(1),NLTEST,NMTEST

      ZRFMROW(1:NLAT) = ZRFMROW(1)
      ZRFHROW(1:NLAT) = ZRFHROW(1)
      ZBLDROW(1:NLAT) = ZBLDROW(1)             ! New version M.Maisonnier (replace the two previous lines)
      GCROW(1:NLAT) = GCROW(1)
      
      DO 50 I=1,NLTEST  
         READ(50,5021) LAT_GRID(I),LON_GRID(I),NMTEST_I(I)   ! M.Maisonnier               DO 51 M=1,NMTEST_I(I)
         if (I==1) THen
             DEGLAT(1:NMTEST_I(I)) = LAT_GRID(I)
             DEGLON(1:NMTEST_I(I)) = LON_GRID(I)

         ELSE
             index_start = int(SUM(NMTEST_I(1:I-1)) + 1)
             index_end   = int(SUM(NMTEST_I(1:I)) )
                 
             DEGLAT(index_start:index_end) = LAT_GRID(I)
             DEGLON(index_start:index_end) = LON_GRID(I)

         ENDif

         DO 51 M=1,NMTEST_I(I)
          READ(50,5011) TITLE1
          READ(50,5045) MIDROT(I,M),FAREROT(I,M)                    
        IF (MIDROT(I,M) .GT. 0.0 ) THEN                           
           CALL XIT('RUNLAKE',-1)
        ELSE
C-- READ IN WATER TILE INFORMATION                              
           READ(50,5095) HLAKROT(I,M),LLAKROT(I,M),
     1                   limP_ROT(I,M)
           NLAKROT(I,M)=NINT(HLAKROT(I,M)/DELZLK)
           limP_ROT(I,M) = limP_ROT(I,M)/
     1     ( KsP + limP_ROT(I,M) )
           BLAKROT(I,M) = limP_ROT(I,M)*4.+0.1
           READ(50,5096) (TLAKROT(I,M,J),J=1,NLAKROT(I,M))
        ENDIF
 51	  CONTINUE
 50   CONTINUE

 

      DO 100 I=1,NLTEST
      DO 101 M=1,NMTEST_I(I)   !!!!!!!!!!!!MM add _I(I)
          IF (NLAKROT(I,M) .GE. 1) THEN                             
          DO 73 J=1,NLAKROT(I,M)                                    
            IF (TLAKROT(I,M,J).LE.200.0) THEN                       
              TLAKROT(I,M,J)=TLAKROT(I,M,J)+TFREZ                   
            ENDIF                                                   
73        CONTINUE                                                  
          ENDIF                                                     
C-- INITIALIZE LAKE SNOW CONDITIONS
          SNOLROT(I,M)=0.0
          ALBSLROT(I,M)=0.0
          RHOSLROT(I,M)=0.0
          TSNOLROT(I,M)=0.0
          WSNOLROT(I,M)=0.0
101	  CONTINUE
100   CONTINUE
C
C=======================================================================
5010  FORMAT(2X,6A8)
5011  FORMAT(8X,A8)
5020  FORMAT(3F10.2,F7.1,3I5)   ! M.Maisonnier -> 2 first F10.2 are remove below -> next line) (replace 3F10.2 by 5F10.2 if you want to go back to the 
5021  FORMAT(2F10.2,I3)     ! M.Maisonnier
5045  FORMAT(I8,F8.3)                                               
5095  FORMAT(2F10.2,F12.4)                                      
5096  FORMAT(410F6.2)                                               
5300  FORMAT(1X,I2,I3,I5,I6,2F9.2,E14.4,F9.2,E12.3,F8.2,F12.2,3F9.2,
     1       F9.4)
6001  FORMAT('CLASS TEST RUN:     ',6A4)
6002  FORMAT('RESEARCHER:         ',6A4)
6003  FORMAT('INSTITUTION:        ',6A4)


C
C------------------------------------------------------------------
C THE CLASS GATHER-SCATTER MACHINERY IS RETAINED HERE FOR
C COMPATIBILITY
C------------------------------------------------------------------
      NML_=0
      NMW=0
      DO 110 I=1,NLTEST
c      PRINT*, NMTEST_I(I)
              DO 120 J=1,NMTEST_I(I)
                  IF(FAREROT(I,J).GT.0.0)      THEN
                      IF(MIDROT(I,J).GT.0)    THEN
                          NML_=NML_+1
                          ILMOS(NML_)=I
                          JLMOS(NML_)=J
c                      PRINT*, '1'
                      ELSE
                          NMW=NMW+1
                          IWMOS(NMW)=I
                          JWMOS(NMW)=J
c                      PRINT*, '2'
                      ENDIF
                  ENDIF
  120         CONTINUE
  110 CONTINUE
      print*, "NML=",NML_, NMW 

C=======================================================================
C     * LAUNCH RUN.

      N=0
      NCOUNT=1
      NDAY=86400/NINT(DELT)
      N_READING=0               ! M. Maisonnier
      DELTAT_READING = 60*60*24
      N_step_met = 0
      NSTEP_write = 1
      IYEAR = 0
      
200   CONTINUE
C
C========================================================================
C     * READ IN METEOROLOGICAL FORCING DATA FOR CURRENT TIME STEP;
C     * CALCULATE SOLAR ZENITH ANGLE AND COMPONENTS OF INCOMING SHORT-
C     * WAVE RADIATION FLUX; ESTIMATE FLUX PARTITIONS IF NECESSARY.

      N=N+1
      N_READING=N_READING+1
      N_step_met = N_step_met +1
      
C      print*, 'OK here runlake.f'
      
      IF ( N_READING .gt. (DELTAT_READING/DELT) ) THEN      ! M.Maisonnier :  allow to have a model time step lower than the time resolution of meteorological data
          N_READING = 1
      ENDIF
      
      year_old = IYEAR
      
      IF (N_READING == 1) THEN
      DO 250 I=1,NLTEST 
          READ(51,5300,END=999) IHOUR,IMIN,IDAY,IYEAR,FSDOWN,FDLROW(I),
     1        PREROW(I),TAROW(I),QAROW(I),UVROW(I),PRESROW(I) 
          
     	if ((IYEAR .ne. year_old) .and. (I==1)) then
            print*, IYEAR
        endif
C
          FSVHROW(I)=0.5*FSDOWN
          FSIHROW(I)=0.5*FSDOWN
          TAROW(I)=TAROW(I)+TFREZ
          ULROW(I)=UVROW(I)
          VLROW(I)=0.0
          UVROW(I)=MAX(VMIN,UVROW(I))
          ZDMROW(I)=10.0
          ZDHROW(I)=2.0
          FSSBROL(I,1)=FSVHROW(I)
          FSSBROL(I,2)=FSIHROW(I)
          RPREROW(I)=0.0           !only needed for IPC=4
          SPREROW(I)=0.0           !only needed for IPC=4
 250   CONTINUE
       ELSE
          IMIN = IMIN + DELT/60
          IHOUR = IHOUR + FLOOR(REAL(IMIN)/60)
          IMIN = IMIN - FLOOR(REAL(IMIN)/60)*60
       ENDIF

      DAY=REAL(IDAY)+(REAL(IHOUR)+REAL(IMIN)/60.)/24.
      DECL=SIN(2.*PI*(284.+DAY)/365.)*23.45*PI/180.
      HOUR=(REAL(IHOUR)+REAL(IMIN)/60.)*PI/12.-PI
      DO 300 I=1,NLTEST
          RADJROW(I)=LAT_GRID(I)*PI/180.
          COSZ=SIN(RADJROW(I))*SIN(DECL)+COS(RADJROW(I))*
     >                     COS(DECL)*COS(HOUR)
          CSZROW(I)=SIGN(MAX(ABS(COSZ),1.0E-3),COSZ)
300   CONTINUE
C================================================================
C
C     * CALCULATION OF ATMOSPHERIC INPUT VARIABLES.

      CALL CLASSI(VPDROW,TADPROW,PADRROW,RHOAROW,RHSIROW,
     1            RPCPROW,TRPCROW,SPCPROW,TSPCROW,TAROW,QAROW,
     2            PREROW,RPREROW,SPREROW,PRESROW,
     3            IPCP,NLAT,1,NLTEST)

C======================================================================
C GATHER CODE FROM CLASSG
C--   GATHER LAKE RELEVANT VARIABLES INTO WATER-TILE GAT ARRAYS
C
      DO 700 K=1,NMW                                                
          HLAKGAT(K+NML_)=HLAKROT(IWMOS(K),JWMOS(K))                 
          LLAKGAT(K+NML_)=LLAKROT(IWMOS(K),JWMOS(K))                 
          BLAKGAT(K+NML_)=BLAKROT(IWMOS(K),JWMOS(K))                 
          NLAKGAT(K+NML_)=NLAKROT(IWMOS(K),JWMOS(K))
          limP(K+NML_) =limP_ROT(IWMOS(K),JWMOS(K))
          DO 705 L=1,NLAKGAT(K+NML_)                                 
            TLAKGAT(K+NML_,L)=TLAKROT(IWMOS(K),JWMOS(K),L)           
705       CONTINUE                                                  
          ASVLGAT(K+NML_)=ASVDROT(IWMOS(K),JWMOS(K))
          ASILGAT(K+NML_)=ASIDROT(IWMOS(K),JWMOS(K))
          BCSNLGAT(K+NML_)=BCSNROT(IWMOS(K),JWMOS(K))
          REFLGAT(K+NML_)=REFROT(IWMOS(K),JWMOS(K))
          SNOLGAT(K+NML_)=SNOLROT(IWMOS(K),JWMOS(K))
          RHOSLGAT(K+NML_)=RHOSLROT(IWMOS(K),JWMOS(K))
          TSNOLGAT(K+NML_)=TSNOLROT(IWMOS(K),JWMOS(K))
          ALBSLGAT(K+NML_)=ALBSLROT(IWMOS(K),JWMOS(K))
          WSNOLGAT(K+NML_)=WSNOLROT(IWMOS(K),JWMOS(K))
700   CONTINUE                                                      
C
C-- ATMOSPHERIC FORCING VARIABLES NEEDED FOR LAKE TILES         
C   GATHERED ON TOP OF LAND TILES                               
C
      DO 800 K=1,NMW                                                
          FSVHGAT(K+NML_)=FSVHROW(IWMOS(K))                          
          FSIHGAT(K+NML_)=FSIHROW(IWMOS(K))                          
          CSZGAT (K+NML_)=CSZROW (IWMOS(K))                          
          FDLGAT (K+NML_)=FDLROW (IWMOS(K))                          
          ULGAT  (K+NML_)=ULROW  (IWMOS(K))                          
          VLGAT  (K+NML_)=VLROW  (IWMOS(K))                           
          TAGAT  (K+NML_)=TAROW  (IWMOS(K))                          
          QAGAT  (K+NML_)=QAROW  (IWMOS(K))                          
          PRESGAT(K+NML_)=PRESROW(IWMOS(K))                          
          RHOAGAT(K+NML_)=RHOAROW(IWMOS(K))                          
          ZRFMGAT(K+NML_)=ZRFMROW(IWMOS(K))                          
          ZRFHGAT(K+NML_)=ZRFHROW(IWMOS(K))                          
          DO L=1,NBS
            FSDBLGAT(K+NML_,L)=FSDBROL(IWMOS(K),L) 
            FSFBLGAT(K+NML_,L)=FSFBROL(IWMOS(K),L) 
            FSSBLGAT(K+NML_,L)=FSSBROL(IWMOS(K),L) 
          ENDDO
          ZDMLGAT(K+NML_)=ZDMROW(IWMOS(K))
          ZDHLGAT(K+NML_)=ZDHROW(IWMOS(K))
          RPCPLGAT(K+NML_)=RPCPROW(IWMOS(K))
          TRPCLGAT(K+NML_)=TRPCROW(IWMOS(K))
          SPCPLGAT(K+NML_)=SPCPROW(IWMOS(K))
          TSPCLGAT(K+NML_)=TSPCROW(IWMOS(K))
          RHSILGAT(K+NML_)=RHSIROW(IWMOS(K))
          RADJLGAT(K+NML_)=RADJROW(IWMOS(K))
          PADRLGAT(K+NML_)=PADRROW(IWMOS(K))
 800   CONTINUE

       
C========================================================================
C CHECK ENERGY BALANCE - start of timestep
C ICEBOT doesn't include weight of snow
C
    
      RHOIW=RHOICE/RHOW
      JL1=1+NML_                                                     
      JL2=NMW+NML_                                                  
      DO 320 I=JL1,JL2
        ICEBOT=RHOIW*LKICEH(I)
        ICETOP=LKICEH(I)-ICEBOT
        IF (ICEBOT .GE. DELSKIN) THEN 
          HCAP=HCPICE
        ELSE IF (LKICEH(I) .LE. 0.0) THEN
          HCAP=HCPW
        ELSE 
          HCAP=(LKICEH(I)*HCPICE + (DELSKIN-ICEBOT)*HCPW)/DELSKIN
        ENDIF
        IF (N .EQ. 1) THEN 
          CTLSTP(I)= -HCAP*TLAKGAT(I,1)*DELSKIN
        ELSE
          CTLSTP(I)= -HCAP*T0LAK(I)*DELSKIN
        ENDIF
C
        DO 330, J=1,NLAKGAT(I)
          ZTOP=DELSKIN + DELZLK*(J -1)
          ZBOT=DELSKIN + DELZLK*J
          IF (ICEBOT .GE. ZBOT) THEN 
            HCAP=HCPICE
          ELSE IF (ICEBOT .LE. ZTOP) THEN
            HCAP=HCPW
          ELSE 
            Z=ICEBOT-ZTOP
            HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK
          ENDIF
          CTLSTP(I)=CTLSTP(I) - HCAP*TLAKGAT(I,J)*DELZLK
330     CONTINUE
320   CONTINUE
C
C========================================================================
C          * LAKE MODEL
C
C       print*, 'before Class L'
     
      ATM_CH4 = 1800.0*1.0e-9 
      CALL CLASSL (HLAKGAT, LLAKGAT, BLAKGAT, NLAKGAT, TLAKGAT,   
     1 T0LAK, HDPTH, LKICEH, SNICEH, ROFICEH,
     2 SNOLGAT,RHOSLGAT,TSNOLGAT,ALBSLGAT,WSNOLGAT,
     3 CDHGAT,CDMGAT,QSENLGAT,TFXLGAT,QEVPLGAT,QFSLGAT,QFXLGAT,
     4 PETLGAT, EFLGAT, GTLGAT, QGLGAT, DRLGAT,
     5 SFTLGAT,SFULGAT,SFVLGAT,SFQLGAT,SFHLGAT,QLWOLGAT,ALVLGAT,ALILGAT,
     6 FSGLGAT, FLGLGAT, HFSLGAT, HEVLGAT, HMFLGAT, HTCLGAT,
     7 FSGSLGAT, FLGSLGAT, HFSSLGAT, HEVSLGAT, HMFNLGAT, HTCSLGAT,
     8 PCPLGAT,PCPNLGAT,QFLGAT,QFNLGAT,ROFNLGAT,FICEGAT,FLSGAT,G0SLGAT,
     9 EXPW, DTEMP, TKELAK, DELU, GRED, RHOMIX,   
     1 FSVHGAT, FSIHGAT, FDLGAT, ULGAT, VLGAT, TAGAT, QAGAT,
     2 RHOAGAT, PADRLGAT, PRESGAT, CSZGAT, ZRFMGAT, ZRFHGAT, 
     3 ZDMLGAT,ZDHLGAT,RPCPLGAT,TRPCLGAT,SPCPLGAT,TSPCLGAT,RHSILGAT,
     4   RADJLGAT,
     5 ASVLGAT,ASILGAT,FSDBLGAT, FSFBLGAT, FSSBLGAT, REFLGAT, BCSNLGAT,
     6 ILG, JL1, JL2, NLAT, NLAKMAX, ISLFD, IZREF, ITG,
     7 IALS, NBS, ISNOALB, IGL, IRSTRT,
     8     N, IYEAR, IDAY, IHOUR, IMIN, TSED,
     9 Clabile,FPP,FRESP,FBUR,CH4,O2,Area,limP ,N_step_met,
     1                        FatmO2,Fw,FPPWC,FSAER,
     6        FswDiff,FatmCH4,Hs,Pchl,FsCH4,FsDiff,FsEbul,zebmin,
     7   Clabs_1,Clabs_2,Clabs_3,FsCH4_C_1,FsCH4_C_2,FsCH4_C_3,
     8   Fsed_1,Fsed_2,Fsed_3,K_diff,	 FatmEbul ,
     1     G0_m,F0_m,FSGL_m,
     1     Q0_m,FLGL_m,HFSL_m,HEVL_m,T0_m,LKICEH_m,SNO_m,
     2     RHOSNI_m,S_m,R_m,OSW_m,RHOSNO_m,SNICEH_m,
     3     TLAK_m,FQU_m,BFLX_m,DISS_m,TKE_m,DELU_m,HDPTH_m,
     3     JMIX_m,ZMIX_m,TMIX_m,DTEMP_m,
     2     FSHEAR_m,FENTRA_m,USTAR_m,EXPW_m,WEDB_m,
     1     DELTHRM_m,ZREFM_m,VA_m,CDML_m,TSNOW_m,
     1     WSNOW_m,ZSNOW_m,FLS_m,
     1     QSWNS_m,QTRANSL_m,QLWIN_m,QLWOS_m,QSENSS_m,QEVAPS_m,
     2     QMELTS_m,GZEROSL_m,CFLUXS_m,QSWIN_m,QFLX_m,
     1     FFLX_m,TSED_m,
     1     Clabile_m,FPP_m,FRESP_m,FBUR_m,O2_m,CH4_m,
     2     FatmO2_m,Fw_m,FPPWC_m,FSAER_m,
     6     FatmCH4_m,Hs_m,FsCH4_C_m,FsDiff_m,
     1     FsEbul_m,zebmin_m,NSTEP_write,
     2     Clabs_1_m,Clabs_2_m,Clabs_3_m,FsCH4_C_1_m,FsCH4_C_2_m,
     3     FsCH4_C_3_m,Fsed_1_m,Fsed_2_m,Fsed_3_m,
     4     FswDiff_m,OxCH4_m,K_diff_m, FatmEbul_m ,
     5     ZPhotic,ZRESP,FRESPWC,ZPhotic_m,ZRESP_m,FRESPWC_m,
     6	   Fibs,previous_ice,
     7     NLAKMAX_bgc,NLAK_bgc,TLAK_bgc_m,HDPTH_bgc_m, 
     8     Fstor_m, DELT_stor_m, HLAK_bgc,
     8     DELT_ice_m,  x_per_m , DEGLAT, ZPROD,ext_C,ATM_CH4)
	 
C        print*, 'after class L'
      
C========================================================================
C CHECK ENERGY BALANCE - end of timestep
C
     
      DO 325 I=JL1,JL2
        ICEBOT=RHOIW*LKICEH(I)
        ICETOP=LKICEH(I)-ICEBOT
        IF (ICEBOT .GE. DELSKIN) THEN 
          HCAP=HCPICE
        ELSE IF (LKICEH(I) .LE. 0.0) THEN
          HCAP=HCPW
        ELSE 
          HCAP=(LKICEH(I)*HCPICE + (DELSKIN-ICEBOT)*HCPW)/DELSKIN
        ENDIF
        CTLSTP(I)= CTLSTP(I) + HCAP*T0LAK(I)*DELSKIN
        DO 335, J=1,NLAKGAT(I)
          ZTOP=DELSKIN + DELZLK*(J -1)
          ZBOT=DELSKIN + DELZLK*J
          IF (ICEBOT .GE. ZBOT) THEN 
            HCAP=HCPICE
          ELSE IF (ICEBOT .LE. ZTOP) THEN
            HCAP=HCPW
          ELSE 
            Z=ICEBOT-ZTOP
            HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK
          ENDIF
          CTLSTP(I)=CTLSTP(I) + HCAP*TLAKGAT(I,J)*DELZLK
335     CONTINUE
        CTLSTP(I)=CTLSTP(I)/DELT
        QSUML=FSGLGAT(I)+FLGLGAT(I)-HFSLGAT(I)-HEVLGAT(I)-HMFLGAT(I)


c        WRITE(79,6019) IYEAR,IDAY,IHOUR,IMIN,QSUML-CTLSTP(I),QSUML
c6019    FORMAT(I4,1X,3(I3,1X),2F8.2)
CCC     IF(ABS(CTLSTP(I)-QSUML).GE. 10.0) THEN
CCC           WRITE(6,6440) N,DAY,CTLSTP(I),QSUML,QSUML-CTLSTP(I)
CCC6440          FORMAT(2X,'LAKE ENERGY BALANCE  ',I8,F8.3,2F20.8,F8.2)
Cmdm          WRITE(6,6450) FSGLGAT(I),FLGLGAT(I),HFSLGAT(I),
Cmdm 1             HEVLGAT(I),HMFLGAT(I),LKICEH(I),T0LAK(I)
Cmdm6450          FORMAT(2X,7F15.6)
Ctest         STOP
CCC     ENDIF
325   CONTINUE
C
C========================================================================
C SCATTER CODE FROM CLASSS 
C--   LAKE DIAGNOSTIC VARIABLES SPLIT OUT

      
      DO 701 K=1,NMW
          HLAKROT(IWMOS(K),JWMOS(K))=HLAKGAT(K+NML_)                 
          LLAKROT(IWMOS(K),JWMOS(K))=LLAKGAT(K+NML_)                 
          BLAKROT(IWMOS(K),JWMOS(K))=BLAKGAT(K+NML_)                 
          NLAKROT(IWMOS(K),JWMOS(K))=NLAKGAT(K+NML_)                 
          DO 706 L=1,NLAKGAT(K+NML_)                                 
            TLAKROT(IWMOS(K),JWMOS(K),L)=TLAKGAT(K+NML_,L)           
706       CONTINUE                                                  
          ASVDROT(IWMOS(K),JWMOS(K))=ASVLGAT(K+NML_)
          ASIDROT(IWMOS(K),JWMOS(K))=ASILGAT(K+NML_)
          SNOLROT(IWMOS(K),JWMOS(K))=SNOLGAT(K+NML_)
          RHOSLROT(IWMOS(K),JWMOS(K))=RHOSLGAT(K+NML_)
          TSNOLROT(IWMOS(K),JWMOS(K))=TSNOLGAT(K+NML_)
          ALBSLROT(IWMOS(K),JWMOS(K))=ALBSLGAT(K+NML_)
          WSNOLROT(IWMOS(K),JWMOS(K))=WSNOLGAT(K+NML_)
701   CONTINUE                                                      
C
C-- FLUX DIAGNOSTIC VARIABLES SPLIT OUT                     
      DO 385 K=1,NMW                                                
          CDHROT (IWMOS(K),JWMOS(K))=CDHGAT (K+NML_)                 
          CDMROT (IWMOS(K),JWMOS(K))=CDMGAT (K+NML_)                 
          HFSLROT (IWMOS(K),JWMOS(K))=HFSLGAT (K+NML_)               
          HEVLROT (IWMOS(K),JWMOS(K))=HEVLGAT(K+NML_)                
          FSGLROT (IWMOS(K),JWMOS(K))=FSGLGAT(K+NML_)                
          FLGLROT (IWMOS(K),JWMOS(K))=FLGLGAT(K+NML_)                
          HMFLROT (IWMOS(K),JWMOS(K))=HMFLGAT(K+NML_)                
          PCPLROT (IWMOS(K),JWMOS(K))=PCPLGAT(K+NML_)                
          PCPNLROT (IWMOS(K),JWMOS(K))=PCPNLGAT(K+NML_)                
          QFLROT (IWMOS(K),JWMOS(K))=QFLGAT(K+NML_)                
          QFNLROT (IWMOS(K),JWMOS(K))=QFNLGAT(K+NML_)                
          ROFNLROT (IWMOS(K),JWMOS(K))=ROFNLGAT(K+NML_)                
          HFSSLROT (IWMOS(K),JWMOS(K))=HFSSLGAT (K+NML_)               
          HEVSLROT (IWMOS(K),JWMOS(K))=HEVSLGAT(K+NML_)                
          FSGSLROT (IWMOS(K),JWMOS(K))=FSGSLGAT(K+NML_)                
          FLGSLROT (IWMOS(K),JWMOS(K))=FLGSLGAT(K+NML_)                
          HMFNLROT (IWMOS(K),JWMOS(K))=HMFNLGAT(K+NML_)                
          HTCSLROT (IWMOS(K),JWMOS(K))=HTCSLGAT(K+NML_)                
          HTCLROT (IWMOS(K),JWMOS(K))=HTCLGAT(K+NML_)                
          SFTLROT (IWMOS(K),JWMOS(K))=SFTLGAT(K+NML_)                
          SFULROT (IWMOS(K),JWMOS(K))=SFULGAT(K+NML_)                
          SFVLROT (IWMOS(K),JWMOS(K))=SFVLGAT(K+NML_)                
          SFQLROT (IWMOS(K),JWMOS(K))=SFQLGAT(K+NML_)                
          SFHLROT (IWMOS(K),JWMOS(K))=SFHLGAT(K+NML_)                
          QLWOLROT (IWMOS(K),JWMOS(K))=QLWOLGAT(K+NML_)                
          ALVLROT (IWMOS(K),JWMOS(K))=ALVLGAT(K+NML_)                
          ALILROT (IWMOS(K),JWMOS(K))=ALILGAT(K+NML_)                
          PETLROT (IWMOS(K),JWMOS(K))=PETLGAT(K+NML_)                
          EFLROT (IWMOS(K),JWMOS(K))=EFLGAT(K+NML_)                
          GTLROT (IWMOS(K),JWMOS(K))=GTLGAT(K+NML_)                
          QGLROT (IWMOS(K),JWMOS(K))=QGLGAT(K+NML_)                
          DRLROT (IWMOS(K),JWMOS(K))=DRLGAT(K+NML_)                
          QSENLROT (IWMOS(K),JWMOS(K))=QSENLGAT(K+NML_)                
          TFXLROT (IWMOS(K),JWMOS(K))=TFXLGAT(K+NML_)                
          QEVPLROT (IWMOS(K),JWMOS(K))=QEVPLGAT(K+NML_)                
          QFSLROT (IWMOS(K),JWMOS(K))=QFSLGAT(K+NML_)                
          QFXLROT (IWMOS(K),JWMOS(K))=QFXLGAT(K+NML_)                
385   CONTINUE                                                      
C
C-- DAY COUNTER
      NCOUNT=NCOUNT+1
      IF(NCOUNT.GT.NDAY) THEN
          NCOUNT=1
      ENDIF
C=======================================================================
C
      GO TO 200
 999   CONTINUE
C C

      CLOSE(UNIT=50)
      CLOSE(UNIT=51)

      CLOSE(UNIT=70)
      CLOSE(UNIT=71)
      CLOSE(UNIT=72)
      CLOSE(UNIT=73)
      CLOSE(UNIT=74)
      CLOSE(UNIT=75)
      CLOSE(UNIT=76)
      CLOSE(UNIT=77)
      CLOSE(UNIT=78)
      CLOSE(UNIT=79)
      CLOSE(UNIT=80)
      CLOSE(UNIT=81)
      CLOSE(UNIT=82)
      CLOSE(UNIT=83)
C      CLOSE(UNIT=84)
C      CLOSE(UNIT=90)


      STOP
      END

      SUBROUTINE Clabile_module(NSTEP,ILG,I,NLAKMAX,NLAK,
     1                      DELZLK,TLAK,Clabile,HDPTH, Area,          
     2                      FPP,FMIN,FBUR,
     3                      IDAY,IHOUR,IMIN,
     5                      LKICEH,limP,shortwave,DEGLAT,
     6                      ext_C,ZPROD, VOL_phot, VOL_tot)
      
C============================================================================
C    ORGANIC MATTER SUBROUTINE
C     * AUTHOR: Manon Maisonnier, Maoyuan Feng, ... ,Pierre Regnier.         
C     * Final revision: March 2025
C     * 
C     * Reference: A new biogeochemical modelling framework (FLaMe v1.0) 
C     *            for lake methane emissions on the regional scale: 
C     *            Development and application to the European domain
C     *         
C     * Correspondance: maoyuan.feng@ulb.be; maoyuanfeng93@gmail.com         
C 
C     NSTEP:      ith time step
C     ILG:        Number of total lakes simulated
C     I:          Index of lakes
C     NLAKMAX:    Max number of layers across all lakes 
C     NLAK:       Number of layers of all lakes 
C     DELZLK:     Water layer depth
C     Clabile:    Organic matter concentration
C     HDPTH:      Mixing depth
C     Area:       Area of different layers for different lakes
C     FPP:        Autochthonous primary production
C     FMIN:       Mineralization rate of organic carbon
C     FBUR:       Burial rate of organic carbon
C     IDAY,IHOUR,IMIN: ith day in the year, ith hour in the day, and ith minutes
C                      in the hour 
C     LKICEH:     Ice cover thickness
C     limP:       Limitation factor of Phorsphorus
C     shortwave:  shortwave radiation
C     DEGLAT:     Latitude of the lake
C     ext_C:      light attenuation
C     ZPROD:      photic depth
C     VOL_phot:   Photic volume
C     VOL_tot:    Total water volume      
C============================================================================    
C DLL=The hours of daylight following Forsythe et al (Ecological Modelling, 1997)
C BC= Average depth-integrated chlorophyll concentration (mgChla/m3), based on Mavaara et al (2017)
C PPRP= The ratio of maximum gross photosynthesis to algal respiration per unit chlorophyll, fixed at 15
C DOB= Bottom oxygen at the first cell (mg/l)
C DF= Clabile based on calibration (HURTADO CAICEDO; mg C/m3)
C P=  Diffusive CH4 production based on calibrated Clabile (mg CH4/m3/d)
C Clab= Clabile based on primary production (BC)
C PP= Diffusive CH4 production based on calculated Clabile (mg CH4/m3/d)
C Pchl= maximum, annually integrated, chlorophyll specific carbon fixation rate (2.5 gC/gChla/hr)
C HOD=  hypolimnetic oxygen demand based on the depth of hypolimnion (0.003*HDPTH(I)+0.021; g/m3/day)
C SOD= sediment oxygen demand (g/m2/day)
C IR= the site-specific PAR hour by hour
C P0Ol = 2 x 10^-8 mol m-3 s-1 and alpha0Ol = 3  m-1 for the oligotrophic state
C P0Eu = 1 x 10^-6 mol m-3 s-1 and alpha0Eu = 30 m-1 for the eutrophic state.
C hs is the depth of the entire sediment column= 5 m for the oligotrophic state and 2 m for the eutrophic state
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ::  NLAKMAX,ILG,IL1,IL2,JMIX,IDAY,IHOUR,IMIN
      INTEGER :: IRSTRT,NSTEP,io    
      INTEGER,DIMENSION(ILG) :: NLAK
      REAL,DIMENSION(ILG) :: LLAK,LKICEH,VA,T0,CH0,USTAR,USTARB,EB
      REAL,DIMENSION(ILG) :: O0,CH4DELTA,IR,PMAX,Clabile,DClab,DOB 
      REAL,DIMENSION(ILG) :: KDGDyn,FPP,FMIN,FBUR,TSUM,TAVE,HDPTHSUM            
      REAL,DIMENSION(ILG) :: Theta,Phi,DLL,ACS,HOD,SOD
      REAL, DIMENSION(ILG) :: HDPTH,limP,ext_C
      REAL,DIMENSION(ILG) :: QSWIN,QLWIN,LI0,DELTHRM,UMLGR,TF,P,BC,M
      REAL,DIMENSION(ILG) :: Theta2,Phi2,Dur,ZEBCAL,kge,FatmCH4
      REAL,DIMENSION(ILG,NLAKMAX) :: FsCH4,FsDiff,FswDiff,FsEbul,
     1                               DCH4s,CH4s         
      REAL,DIMENSION(ILG,NLAKMAX) :: TLAK,CH4,DCH4LAK,
     1								 CH4FLUX,LI !,QFLX
      REAL,DIMENSION(ILG,NLAKMAX) :: O2,DOY,AF,ProdEu,ProdOl
      REAL,DIMENSION(ILG,NLAKMAX) :: OMAX,OxCH4
      REAL,DIMENSION(ILG,NLAKMAX) :: KZ_DIFF,DCH4LAKback,DCH4LAKfor
      REAL :: DELSKIN,DELZLK,TIMSTP,RHOW,TCW
c     REAL,DIMENSION(100) :: HDPTHR
      REAL,DIMENSION(ILG,NLAKMAX) :: Prod,Zebmin, P0
      REAL,DIMENSION(ILG) :: alpha0, Hs, VOL_phot, VOL_tot,DAYT

      REAL, DIMENSION (ILG,NLAKMAX) :: Clabs_1,Clabs_2,Clabs_3
      REAL, DIMENSION (ILG,NLAKMAX) :: Fsed_1,Fsed_2,Fsed_3


      REAL,DIMENSION(ILG,0:NLAKMAX+1) :: Area
      REAL, PARAMETER :: P_Chl_max = 1.*0.5
      REAL,DIMENSION(ILG) :: shortwave
      REAL,DIMENSION(ILG) :: DEGLAT,ZPROD
      
C ----* COMMON BLOCKS *------------------------------------------------
 
      REAL DELT
      COMMON /CLASS1/ DELT
 
C ----* LOCAL VARIABLES *---------------------------------------------
 
      INTEGER :: I,J,K,N,ok
      REAL :: FW,C1,C2,UC,UT,D,nu,ze,gamma,lambda,Sc,
     1     Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,
     2     Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19,Q20,Q21,
     3     fg, A1,A2,A3,A4,S,B1,B2,B3,KCH4,KO2,Vmax,
     4     RC,PQ10,TPR,k20,kbur,Nk,ScCH4,
     5     CH4CSTAR,CH4USTAR,CH4EX,U6X,KZ,
     6     cp, ror,vk,
     7     KC,PPRP,DL,IK,KDW,KDP,KDG,Pchl,CF,Pi,Vol,
     8     P0Ol,P0Eu,alpha0Ol,alpha0Eu,K10,HsOl,HsEu,
     9     HO2,PO2,HCH4,PCH4,Patm,alphaO2,Ck1,Ck2,
     1     Tmp,rho0, g, a, b,PchlEu,PchlOl,
     2     z1,z2,z3,zerror

C
C---- zebmin
      
      DATA Tmp ,     rho0 ,   g
     1/    273.15  , 1000 ,  9.81  /    ! K , kg/m3 , m/s2
   
      
C	  
C-1--Wind speed and shear velocity coefficients
      DATA  Q1,      Q2,      Q3       
     1/     3.35E-2, 9.6E-4,  4.1E-5/
C
C-2--Ebullition coefficients
      DATA  Q4,       Q5,    Q6,     Q7       
     1/     -3.6E-4,  0.16,  -0.94,  1.48/
C    1/     4.6771,   0.1703, 0.00,  0.00/
C
C-3--Methane flux coefficients
C-FW
      DATA  Q8,       Q9
     1/     5.8E-2,   4.7E-2/
C-Methane exchange coefficients
      DATA  C1,          C2,          D,        nu
     1/     2.55237E-2,  5.47723E-2,  1.62E-9,  1E-6/
      DATA  ror,       	 vk
     1/     0.0012,		 0.4/
      DATA  UC,      UT,    ze,      gamma,  lambda
     1/     1.1E-1,  1E-1,  2.5E-3,  6.5,    3.0/                   
C    1/     1.1E-1,  1E-1,  3.5E-3,  5.0,    3.0/                    !Large scale Toolik?
C    1/     6.2E-2,  9E-2,  1.5E-3,  2.5,    2.0/                    !Large scale Toolik?
      DATA  Q10,        Q11,        Q12,        Q13,  Q14   
     1/     1.60626E+2, 9.0909E+0,  5.7735E-4,  10,   1.33333E-4/
C-Methane solubility coefficients  
      DATA  fg,      A1,           A2,          A3         
     1/     1.9E-6, -4.152807E+2, 5.968104E+2, 3.792599E+2/
      DATA  A4,            S,  B1,         B2,         B3          
     1/     -6.20757E+1,   0.0,  -5.916E-2,  3.2174E-2,  -4.8198E-3/ !No salinity Toolik   
      DATA  Q15        
     1/     1.604246E-5/                                           
C
C-4--Productivity coefficients
      DATA  Q16,       Q17,      Q18    
     1/     8.4737E0,  -2.83E-1, 1E-6/
      DATA  RC,    PQ10,  TPR,    Q19,  k20  , kbur  
c     1/     0.02, 2.2, 3.5, 10., 0.0082, 0.0192 /	
     1/     0.02, 2.2, 3.5, 10., 0.008, 0.004 /						         !Maavara et al. (2017) up  boost
c     1/     0.02,  2.2,   3.5,	 10.,	 0.05/						         !Maavara et al. (2017) up
c     1/     0.02,  2.2,   3.5,	 10.,	 0.03/						         !Maavara et al. (2017)
c     1/     0.02,  2.2,   3.5,	 10.,	 0.09/						         !Kessler et al. (2012)
c     1/     0.02,  3.5,   -3.0,  10./                        		 !Toolik-D Zhuang et al. (2004)
C    1/     0.02,  4.0,   -5.5,  10./        						 !Toolik-W Zhuang et al. (2004)
C
C-5--Methane concentration coefficients
C-Methane oxidation coefficients
      DATA  Q20,        Q21    												
     1/     1.604246E1, 1.59994E1/
C-Diffusion coefficient
      DATA  KZ     
     1/     1.0E-7/
C-Clabile coefficient (Lewis (2011) & Maavara et al (2017))
      DATA  KC,    PPRP, DL,  IK,  KDW,  KDP,  KDG, CF    
     1/     0.014, 15,   12,  120, 0.13, 0.06, 0.24, 0.5/
      DATA  Pi     
     1/     3.14159265359/   
C-Lake Volume (surface area:1.49 km3, mean depth: 7.4 m)
      DATA  Vol     
     1/     11026000/
C-Methane Production Parameters (The oligotrophic values from Stepanenko et al. 2016 and the eutrophic values from Langenegger et al. 2019)
      DATA  P0Ol,  P0Eu,    alpha0Ol,alpha0Eu, K10, HsOl,HsEu
     1/     2.0E-8,1.0E-6,  3,		 30,	   2.3, 5,	 2/
          DATA PchlEu, PchlOl
     1/         0.15 , 0.01/     
      DATA  HO2,	PO2,   HCH4,   PCH4,       Patm
     1/     1.3E-5, 0.2095,1.3E-5, 1800.0E-9  ,101325/	! Atmospheric pressure=101325 (Pa), HO2: Henry's constant for O2 at reference temperature (mol m-3 Pa-1), PO2: Atmospheric partial pressure of O2 (Atmospheric partial pressure)
      DATA  Ck1,     Ck2,     Nk,   ScCH4									! Ck1 = 5.75d-6 (m s-1) Ck2 = 5.97d-6 (m s-1)^{-Nk}
     1/     5.75d-6, 5.97d-7, 1.7  ,570/


          
          
C-----set----------------------------------------------------------------
       
       Pchl = P_Chl_max*limP(I)   ! ????   Calibtated ... should be 2.5
          
C-----------------------------------------------------------------------

      	P(I)=0.0
      	EB(I)=0.0

        IF ( (-tan((2.*Pi/360.)*DEGLAT(I))*tan((2.*Pi/360.)*
     2          23.45*SIN((2.*Pi/360.)*360.*(IDAY-81.0)/365.)
     3          )) .GE. 1.0) THEN
        DAYT(I) = 0.0

        IR(I) = 0.0
 
        ELSEIF ( (-tan((2.*Pi/360.)*DEGLAT(I))*tan((2.*Pi/360.)*
     2          23.45*SIN((2.*Pi/360.)*360.*(IDAY-81.0)/365.)
     3          )) .LE. -1.0) THEN
        DAYT(I) = 24.0

        IR(I) = shortwave(I)*0.58*4.6/(DAYT(I)/24)
 
        ELSE
        DAYT(I) = REAL((24./Pi)*ACOS(
     1          -tan((2.*Pi/360.)*DEGLAT(I))*tan((2.*Pi/360.)*
     2          23.45*SIN((2.*Pi/360.)*360.*(IDAY-81.0)/365.)
     3          )))

        IR(I) = shortwave(I)*0.58*4.6/(DAYT(I)/24)

      ENDIF     

  

        IF (IR(I) < 0.0) then
            PRINT*, IR(I)
            STOP 
        ENDIF

        Theta(I)=0.2163108+2*ATAN(0.9671396*TAN(0.00860*(IDAY-186)))
        Phi(I)=ASIN(0.39795*COS(Theta(I)))
            ACS(I)=( SIN(0.833*Pi/180)+SIN(68.63*Pi/180)*SIN(Phi(I)))/
     1         ( COS(68.63*Pi/180)*COS(Phi(I)) ) 
          IF (ACS(I) .GE. 1) THEN  
            ACS(I)=1
          ELSE IF (ACS(I) .LE. -1) THEN
            ACS(I)=-1
          ENDIF
        DLL(I)=24-24/3.14159*ACOS(ACS(I))	

        if (Clabile(I)/1000
     1       .lt. 1.d-30) then 
            KDGDyn(I)=0.0
        elseif ( Clabile(I)/1000 .le. 400. ) THEN
        KDGDyn(I)= EXP(-4.44+1.8*LOG(Clabile(I)/1000.)-
     1                 0.149*(LOG(Clabile(I)/1000.))**2) 
        else
            KDGDyn(I)= EXP(-4.44+1.8*LOG(400.)-
     1         0.149*(LOG(400.))**2) 
        endif
 
        IF  (( LKICEH(I).GT. 1.d-30 ) .OR. (IR(I) < 1.d-4)) THEN 
                BC(I)=0.0
        ELSE    
                 BC(I)=(1/KC)*
     1    ( 0.75*PPRP*LOG(0.7*IR(I)/(0.5*IK))*(1/ZPROD(I))*DAYT(I)/24
     2            -(KDW+KDP+KDGDyn(I)) )

               BC(I)=min(BC(I),600.0)

        END IF
        IF(BC(I).LT. 0.0) THEN  
          BC(I)=0.0
        END IF

        ext_C(I) = KDW+KDP+KDGDyn(I) 


        ZPROD(I) = 4.6/(KC*BC(I) + ext_C(I))

        ZPROD(I) = MAX(ZPROD(I),HDPTH(I))

        ZPROD(I) = MIN(ZPROD(I),DELZLK*NLAK(I))
        VOL_tot(I) = SUM(Area(I,1:NLAK(I))*DELZLK)
        VOL_phot(I) = 0.
        DO 10500 J=1,max(int(ZPROD(I)/DELZLK),1)
           VOL_phot(I) = VOL_phot(I) + Area(I,J)*DELZLK !!!!
10500     CONTINUE
        VOL_phot(I) = min(VOL_tot(I), VOL_phot(I))

               
C     C-------Depth Average Temp------------------------------------


        TSUM(I)=0.
        DO 410 J=1,max(int(ZPROD(I)/DELZLK),2)
            TSUM(I)=TSUM(I)+TLAK(I,J)-273.15
 410    CONTINUE
        TAVE(I)=TSUM(I)/max(int(ZPROD(I)/DELZLK),2) !smax(int(ZPROD(I)/DELZLK),2)
        IF(TAVE(I).LT. 28.) THEN
         M(I)=2.**((TAVE(I)-28.)/10.)
        ELSE
         M(I)=1.
        ENDIF 

C     Calculate DClabile based on primary production

        FPP(I)=BC(I)*Pchl*M(I)*(DELT/3600.) 
     1         *VOL_phot(I)/VOL_tot(I) 

                               !      *min(ZPROD(I)/(NLAK(I)*DELZLK),1.) ! mgC/m3/dt
        FPP(I)= MAX(FPP(I),0.0)

        FMIN(I)=k20*Clabile(I)*(DELT/(3600.*24.))*			     ! mgC/m3/dt
     1        1.02**(TAVE(I)-20)
        FMIN(I) = MAX(FMIN(I),0.0)

        FBUR(I)=kbur*Clabile(I)*(DELT/(3600.*24.))			     ! mgC/m3/dt
        FBUR(I) = MAX(FBUR(I),0.0)

            IF((Clabile(I)+DClab(I)).LT.0.0) THEN
               DClab(I) = - Clabile(I) 
               Clabile(I)= 0.0
            ELSE
               DClab(I)=FPP(I) - FMIN(I) - FBUR(I)
               Clabile(I)=Clabile(I) + DClab(I)
	        END IF

        
  
      END

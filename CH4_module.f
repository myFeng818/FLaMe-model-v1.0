      SUBROUTINE CH4_module(NSTEP,ILG,I,NLAKMAX_bgc,NLAK_bgc,DELZLK,
     1                      TLAK_bgc,Clabile,CH4,O2,Area,HDPTH_bgc,           
     2                      FPP,FRESP,
     3                      IDAY,IHOUR,
     4                      JMIX,VA,
     5                      LKICEH,limP,
     4     FswDiff,
     6               FatmCH4,Hs,Pchl,FsCH4,FsDiff,FsEbul,zebmin,
     7    Clabs_1,Clabs_2,Clabs_3,FsCH4_C_1,FsCH4_C_2,FsCH4_C_3,
     3     Fsed_1,Fsed_2,Fsed_3, K_diff , DTHERMO , FatmEbul,
     7      Fibs,previous_ice, f_sed_ , 
     6                  VOL_phot, VOL_tot,ATM_CH4,zox)

C============================================================================
C    METHANE SUBROUTINE
C              
C     * AUTHOR: Manon Maisonnier, Maoyuan Feng, ... ,Pierre Regnier.         
C     * Final revision: March 2025
C     * 
C     * Reference: A new biogeochemical modelling framework (FLaMe v1.0) 
C     *            for lake methane emissions on the regional scale: 
C     *            Development and application to the European domain
C     *         
C     * Correspondance: maoyuan.feng@ulb.be; maoyuanfeng93@gmail.com   

C============================================================================
C
C DLL=The hours of daylight following Forsythe et al (Ecological Modelling, 1997)
C BC= Average depth-integrated chlorophyll concentration (mgChla/m3), based on Mavaara et al (2017)
C PPRP= The ratio of maximum gross photosynthesis to algal respiration per unit chlorophyll, fixed at 15
C DOB= Bottom oxygen at the first cell (mg/l)
C DF= Clabile based on calibration (HURTADO CAICEDO; mg C/m3)
C P=  Diffusive CH4 production based on calibrated Clabile (mg CH4/m3/d)
C Clab= Clabile based on primary production (BC)
C PP= Diffusive CH4 production based on calculated Clabile (mg CH4/m3/d)
C Pchl= maximum, annually integrated, chlorophyll specific carbon fixation rate (2.5 gC/gChla/hr)
C HOD=  hypolimnetic oxygen demand based on the depth of hypolimnion (0.003*HDPTH_bgc(I)+0.021; g/m3/day)
C SOD= sediment oxygen demand (g/m2/day)
C IMAX= the maximum site-specific PAR based on empirical relationship (Lewis 2011)
C P0Ol = 2 x 10^-8 mol m-3 s-1 and alpha0Ol = 3  m-1 for the oligotrophic state  
C P0Eu = 1 x 10^-6 mol m-3 s-1 and alpha0Eu = 30 m-1 for the eutrophic state.
C hs is the depth of the entire sediment column= 5 m for the oligotrophic state and 2 m for the eutrophic state
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ::  NLAKMAX_bgc,ILG,IL1,IL2,JMIX,IDAY,IHOUR
      INTEGER :: IRSTRT,NSTEP,io    , test
      INTEGER,DIMENSION(ILG) :: NLAK_bgc
      REAL,DIMENSION(ILG) :: LLAK,LKICEH,VA,T0,CH0,USTAR,USTARB,EB
      REAL,DIMENSION(ILG) :: O0,CH4DELTA,IMAX,PMAX,Clabile,DClab,DOB 
      REAL,DIMENSION(ILG) :: KDGDyn,FPP,FRESP,TSUM,TAVE,HDPTHSUM            
      REAL,DIMENSION(ILG) :: Theta,Phi,DLL,ACS,HOD,SOD
      REAL, DIMENSION(ILG) :: HDPTH_bgc,limP
      REAL,DIMENSION(ILG) :: QSWIN,QLWIN,LI0,DELTHRM,UMLGR,TF,P,BC,M
      REAL,DIMENSION(ILG) :: Theta2,Phi2,Dur,ZEBCAL,kge,FatmCH4
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: FsCH4,FsDiff,FswDiff,FsEbul,
     1                           DCH4s,CH4s, FatmEbul, FED, FsDiff_         
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: TLAK_bgc,CH4,DCH4LAK,CH4FLUX,
     1									 LI !,QFLX
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: O2,DOY,AF,ProdEu,ProdOl
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: OMAX,OxCH4
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: KZ_DIFF,DCH4LAKback,DCH4LAKfor
      REAL :: DELSKIN,DELZLK,TIMSTP,RHOW,TCW,DEGLAT
c     REAL,DIMENSION(100) :: HDPTHR
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: Prod,Zebmin, P0
      REAL,DIMENSION(ILG) :: Hs, alpha0, VOL_phot, VOL_tot
      !REAL,DIMENSION(ILG,NLAKMAX_bgc) :: alpha_0

      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Clabs_1,Clabs_2,Clabs_3
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FsCH4_C_1,FsCH4_C_2,
     1		 FsCH4_C_3,
     3     Fsed_1,Fsed_2,Fsed_3,FsO2
	 
	     REAL, DIMENSION (ILG) :: DTHERMO
             REAL, DIMENSION (ILG) :: ATM_CH4
      
      REAL, DIMENSION(ILG,NLAKMAX_bgc+1):: K_diff


      REAL,DIMENSION(ILG,0:NLAKMAX_bgc+1) :: Area

      REAL :: k_C2,x_C2,t_C2  
      REAL :: k_C3,x_C3,t_C3
	  
	  REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Fibs
      INTEGER, DIMENSION (ILG) :: previous_ice

      REAL,DIMENSION (ILG,NLAKMAX_bgc) :: zox
      REAL,DIMENSION (ILG,NLAKMAX_bgc) :: Ox_red_sed
	

      REAL, DIMENSION(ILG,NLAKMAX_bgc) :: Dox_sed
      
      REAL, PARAMETER :: K_ox = 0.67 !gO2/m3
      REAL, PARAMETER :: f_sed =1.*0.25
      REAL, PARAMETER :: f_ox = 1.*5.0
      REAL, PARAMETER :: f_bubble = 1.*1.0
      REAL, PARAMETER :: f_ox_ice = 1.*0.5
      REAL, PARAMETER :: f_ice_diss = 1.*0.5
      
      REAL, PARAMETER :: alpha = 1.0*1.0

      REAL :: min_alpha != 3.0/param_alpha
      REAL :: var_alpha != param_alpha*1.0

      REAL :: f_sed_ , f_ox_, f_ox_ice_ , f_ice_diss_

C ----* COMMON BLOCKS *------------------------------------------------
 
      REAL DELT
      COMMON /CLASS1/ DELT
 
C ----* LOCAL VARIABLES *---------------------------------------------
 
      INTEGER :: I,J,K,N,ok
      REAL :: FW,C1,C2,UC,UT,D,nu,ze,gamma,lambda,Sc,
     1     Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,
     2     Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19,Q20,Q21,
     3     fg, A1,A2,A3,A4,S,B1,B2,B3,KCH4,KO2,Vmax,
     4     Nk,ScCH4,
     5     CH4CSTAR,CH4USTAR,CH4EX,U6X,KZ,
     6     cp, ror,vk,
     7     KC,PPRP,DL,IK,KDW,KDP,KDG,Pchl,CF,Pi,Vol,
     8     P0Ol,P0Eu,alpha0Ol,alpha0Eu,K10,HsOl,HsEu,
     9     HO2,PO2,HCH4,PCH4,Patm,alphaO2,Ck1,Ck2,
     1     Tmp,rho0, g, a, b,PchlEu,PchlOl,
     2     z1,z2,z3,zerror , fracbubble

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
     1/     0.15 , 0.01/     
 	    DATA  HO2,	PO2,   HCH4,   PCH4,       Patm
     1/     1.3E-5, 0.2095,1.3E-5, 1800.0E-9  ,101325/	! Atmospheric pressure=101325 (Pa), HO2: Henry's constant for O2 at reference temperature (mol m-3 Pa-1), PO2: Atmospheric partial pressure of O2 (Atmospheric partial pressure)
 	    DATA  Ck1,     Ck2,     Nk,   ScCH4									! Ck1 = 5.75d-6 (m s-1) Ck2 = 5.97d-6 (m s-1)^{-Nk}
     1/     5.75d-6, 5.97d-7, 1.7  ,677 /    ! Sc from   "Advances in Quantifying Air-Sea Gas Exchange and Environmental Forcing* " of RikWanninkhof et at 2009,  ! sc other 570/

        

  
c-----Sensibility set--------------------------

      if (f_sed .gt. 1.0) then
          f_sed_ = 1.0
      elseif (f_sed .lt. 1.d-10) then
          f_sed_ = 0.0
      else
          f_sed_ = f_sed
      endif
      !PRINT*, f_sed

      if (f_ox .gt. 95) then
          f_ox_ = 95
      elseif (f_ox .lt. 1.d-10) then
          f_ox_ = 0.0
      else
          f_ox_ = f_ox
      endif

      if (f_ox_ice .gt. 1.0) then
          f_ox_ice_ = 1.0
      elseif (f_ox_ice .lt. 1.d-10) then
          f_ox_ice_ = 0.0
      else
          f_ox_ice_ = f_ox_ice
      endif

      if (f_ice_diss .gt. 1.0) then
          f_ice_diss_ = 1.0
      elseif (f_ice_diss .lt. 1.d-10) then
          f_ice_diss_ = 0.0
      else
          f_ice_diss_ = f_ice_diss
      endif     


      ! var_alpha = alpha*1.0
      ! if (alpha .lt. 1.e-10) then
          ! min_alpha = 1000000.0
      ! else
          ! min_alpha = 3.0/alpha
      ! endif
     
C-----set----------------------------------------------------------------
          !f_sed_ = f_sed

          Hs(I) = 5. !% HsOl

          t_C2 = 3.*30.*24.*60.*60.   
          x_C2 = 50./100.
          k_C2 = LOG(x_C2)/t_C2
       
          t_C3 = 6.*30.*24.*60.*60.   
          x_C3 = 50./100.
          k_C3 = LOG(x_C3)/t_C3
        
C-----Ebullition based on Total Production - Diffusion

        alpha0(I)= alpha*(10+FPP(I)*VOL_tot(I)/VOL_phot(I))

        DO 103 J=1,NLAK_bgc(I)

       alphaO2  =  0. ! 10. !  316.8   ! m³/mol  from LAKE 2.0  
	   

        IF ( Area(I,J)-Area(I,J+1) .lt. 1.e-30 ) THEN
              FsCH4(I,J) = 0.d0
              Fsed_1(I,J) = 0.d0
              Fsed_2(I,J) = 0.d0
              Fsed_3(I,J) = 0.d0
        ELSE   
            Fsed_1(I,J) = (100./100.)*f_sed
     1             *FRESP(I)*sum(Area(I,1:NLAK_bgc(I)))*DELZLK
     2             /(Area(I,0)*Hs(I))   

            Fsed_2(I,J) =  (0./100.)*f_sed
     1             *FRESP(I)*sum(Area(I,1:NLAK_bgc(I)))*DELZLK
     2              /(Area(I,0)*Hs(I))  

            Fsed_3(I,J) =  (0./100.)*f_sed
     1                *FRESP(I)*sum(Area(I,1:NLAK_bgc(I)))*DELZLK	
     2              /(Area(I,0)*Hs(I))  
 
           ENDIF


           Clabs_1(I,J) = Clabs_1(I,J) + Fsed_1(I,J)     ! mgC/m3
           Clabs_2(I,J) = Clabs_2(I,J) + Fsed_2(I,J)
           Clabs_3(I,J) = Clabs_3(I,J) + Fsed_3(I,J)

           FsCH4_C_1(I,J) = Clabs_1(I,J)         ! mgC/m3
           FsCH4_C_2(I,J) = Clabs_2(I,J)*(1-exp(k_C2*DELT))
           FsCH4_C_3(I,J) = Clabs_3(I,J)*(1-exp(k_C3*DELT))

           Clabs_1(I,J) = Clabs_1(I,J) - FsCH4_C_1(I,J)
           Clabs_2(I,J) = Clabs_2(I,J) - FsCH4_C_2(I,J)
           Clabs_3(I,J) = Clabs_3(I,J) - FsCH4_C_3(I,J)
           
           FsCH4(I,J) = ( FsCH4_C_1(I,J) + FsCH4_C_2(I,J) +
     1                    FsCH4_C_3(I,J) ) *16./12. ! mgCH4/m3/dt  
           


         IF (FsCH4(I,J) .eq. 0.0) THEN
            P0(I,J) = 0.
         ELSE   
            P0(I,J) = FsCH4(I,J)
     1           *(Hs(I)*alpha0(I))/(1-exp(-alpha0(I)*Hs(I)))        !mgCH4/m3/dt  

         ENDIF

         
c     ---- ZEBMIN ---



         IF (P0(I,J) .lt. 1.d-20) THEN
              Zebmin(I,J) = Hs(I)
         ELSE
         CALL search_zebmin(J,Zebmin(I,J),P0(I,J),alpha0(I),Hs(I),
     1                         CH4(I,J),TLAK_bgc(I,J),DELT,DELZLK)

         ENDIF



      IF (FsCH4(I,J) .eq. 0.0) THEN
         FsDiff(I,J) = 0.
         FsEbul(I,J) = 0.
      ELSE
       FsDiff(I,J)=P0(I,J)					! mgCH4/m3/dt; for P0:    gC/m3/s   x 0.5molCH4/molC x molC/12grC x 1000mgCH4/molCH4 x 300s/dt
     2        	 *(1-exp(-alpha0(I)*Zebmin(I,J)))/
     3        	 (alpha0(I)*Hs(I))
       FsEbul(I,J)= FsCH4(I,J) - FsDiff(I,J)      


      ENDIF

           FsDiff_(I,J) = FsDiff(I,J)*(1.0-O2(I,J)/(O2(I,J)+K_ox))

             
          Ox_red_sed(I,J) = 2.* (FsDiff(I,J)-FsDiff_(I,J)) 
     1	    * (32./12.)/1000 + FRESP(I)*(1.-f_sed)*(32./12.)/1000.    !gO2/m3/dt
          ! oxygen consumption rate

          Ox_red_sed(I,J) = max(Ox_red_sed(I,J),0.000001)

       Dox_sed(I,J) = 1E-9*10**(3.672-984.26/(273.+TLAK_bgc(I,J)))
     1             *DELT  ! m2/dt, oxygen diffusion coefficient in sediment

          zox(I,J) = (2.*Dox_sed(I,J)*O2(I,J)/Ox_red_sed(I,J))**0.5 !
                      ! (m2/dt * gO2/m3/ (gO2/m3/dt))**0.5 = m oxygen
                      ! penetration depth
                              
          zox(I,J) = min(zox(I,J),50.)

     	  O2(I,J) = O2(I,J) - 
     1	   2.* (FsDiff(I,J)-FsDiff_(I,J)) 
     1	      * (32./12.) * (Hs(I)/DELZLK)
     1           *((Area(I,J)-Area(I,J+1))/Area(I,J))/1000.     ! gO2/m3
           
 
         FsDiff(I,J) = FsDiff_(I,J)


      IF ( Area(I,J)-Area(I,J+1) .gt. 1.d-30 ) THEN
	      DCH4s(I,J)=(FsCH4(I,J)-FsDiff(I,J)-FsEbul(I,J))	! mgCH4/m3/dt 
      ELSE
            DCH4s(I,J) = 0.d0
      ENDIF   
            CH4s(I,J)=CH4s(I,J)+DCH4s(I,J)
  
 

C------------ Proportion of methane bubbling from sediment reaching the surface as bubble --------


	
        fracbubble = f_bubble*(  
     1     - (0.0000003*(J*DELZLK)**4)  
     1     - (0.00043*(J*DELZLK)**3)
     2     + (0.075*(J*DELZLK)**2) 
     3     - (4.55*(J*DELZLK))
     4	    + 100. ) /100.
        IF (fracbubble .gt. 100.0/100) THEN
            fracbubble = 100.0/100
        ELSEIF (fracbubble .lt. 1.d-10/100) THEN
            fracbubble = 0.0 
        ENDIF
	 
        FatmEbul(I,J) = fracbubble*FsEbul(I,J)
     	FatmEbul(I,J) = FatmEbul(I,J)*Hs(I)*(Area(I,J)-Area(I,J+1))  
        FatmEbul(I,J) = FatmEbul(I,J)/(Area(I,0)*1000.)                   	 ! gCH4 / m² /dt
			
			
		IF ((LKICEH(I).LE.1.e-30) .AND. (previous_ice(I)==1)) THEN
           FatmEbul(I,J) = FatmEbul(I,J) + Fibs(I,J)*f_ox_ice
		   Fibs(I,J) = 0.0
        ENDIF
			
        IF(LKICEH(I).GT. 1.e-30) THEN
           CH4(I,1) =  CH4(I,1) + f_ice_diss*FatmEbul(I,J)/DELZLK
			Fibs(I,J) = Fibs(I,J) + (1-f_ice_diss)*FatmEbul(I,J)
           FatmEbul(I,J) = 0.
           IF (J == NLAK_bgc(I) ) THEN
            previous_ice(I) = 1
           ENDIF
		ELSEIF (J == NLAK_bgc(I) ) THEN 
		   previous_ice(I) = 0
		ENDIF
			
        FED(I,J)=(1.-fracbubble)*FsEbul(I,J)
        FED(I,J)=FED(I,J)*(Hs(I)*(Area(I,J)-Area(I,J+1)))/(J*1000.)   ! gCH4/dt (number of g of CH4 that are going to be diluted in each layer


            FswDiff(I,J)= (FsDiff(I,J)*(Hs(I)/DELZLK)
     1           *(Area(I,J)-Area(I,J+1))/Area(I,J))/1000. ! gCH4/m3/dt
  
C
CC--------Bottom CH4--------------------	    
C

       CH4(I,J)=CH4(I,J)+FswDiff(I,J) ! mgCH4/L   = gCH4/m3

       CH4(I,1:J)  = CH4(I,1:J) 
     1   		  +FED(I,J)/(Area(I,1:J)*DELZLK)


 103     CONTINUE
 

C     C--------Atmospheric Fluxes (Flame; Maisonnier)------------
            
	      kge(I) = (Ck1 + Ck2*VA(I)**Nk)*(sqrt(600/ScCH4)) ! m/s
            
        IF(LKICEH(I).GT. 1.e-30) THEN
        	FatmCH4(I)=0.0
        ELSEIF (TLAK_bgc(I,1) == 0.0 ) THEN
          FatmCH4(I)= MIN( CH4(I,1)*DELZLK ,        !                   mg/l * m/s *s/dt = g/m³ * m/s *s/dt = g/(m²dt)
     1        kge(I)*(CH4(I,1)-4.08d-5)*DELT)
        ELSE
	      FatmCH4(I) = MIN( CH4(I,1)*DELZLK ,               ! mg/l * m/s * s/dt = g/(m²dt)
     1        	 kge(I)*(CH4(I,1)-(HCH4
     2          *exp( 1800.*((1/TLAK_bgc(I,1))-(1/298.15)) )
     3           *ATM_CH4(I)*1.0e-9*Patm*16))*DELT )
	      CH4(I,1)=  CH4(I,1)-FatmCH4(I)/DELZLK
        END IF


C
C     C-------Transport---------------------

            CALL transport(I,ILG,CH4,DELT,TLAK_bgc,
     1                  NLAK_bgc,NLAKMAX_bgc,HDPTH_bgc,JMIX,VA,Area,  
     2                    K_diff , DTHERMO
     3                    )

        


      RETURN
                
    
      END

      SUBROUTINE O2_module(NSTEP,ILG,I,NLAKMAX_bgc,NLAK_bgc,DELZLK,
     1                     TLAK_bgc,O2,Area,HDPTH_bgc,FPP,FRESP,
     4                     JMIX,VA,
     5                     LKICEH,
     6                     FatmO2,Rw,FPPWC,RSAER,
     5     ZPhotic,ZRESP,FRESPWC,QFLX_bgc, K_diff,DTHERMO,f_sed,
     6     ZPROD, VOL_phot, VOL_tot,zox)


C=======================================================================
C     * MAR  31/16 - M.MACKAY.  PLACED IN SEPARATE SUBROUTINE
C     * JAN  22/16 - L.BOEGMAN,  LAKE DISSOLVED OXYGEN CALCULATIONS
C     *             A.JABBARI

C    OXYGEN SUBROUTINE
C     * AUTHOR: Manon Maisonnier, Maoyuan Feng, ... ,Pierre Regnier.         
C     * Final revision: March 2025
C     * 
C     * Reference: A new biogeochemical modelling framework (FLaMe v1.0) 
C     *            for lake methane emissions on the regional scale: 
C     *            Development and application to the European domain
C     *         
C     * Correspondance: maoyuan.feng@ulb.be; maoyuanfeng93@gmail.com    

C=======================================================================
C
C------------------------
C DISSOLVED OXYGEN BUDGET  
C  - FOLLOWING FROM HAMILTON AND SCHLADOW ENVIRON. MOD. (1997)
C  - HOD FROM BOUFFARD ET AL WATER RESOUR RES (2013)
C  - JUNE 2015, L Boegman 
C
C CALCULATE THE EFFECT OF SURFACE AERATION ON OXYGEN. USE THE METHOD OF 
C PATTERSON ET AL. (1985). INCORRECT COEFFFICIENTS FROM THE 1985 PAPER HAVE BEEN ADJUSTED.
C CALCULATE THE WIND-AIR SHEAR AT THE SURFACE (FROM IVEY AND PATTERSON 1984)
C
C OCSTAR=OXYGEN CONCENTRATION AT SATURATION
C DELTA=CHANGE IN OXYGEN CONCENTRATION OVER THE TIME STEP
C CHANGE IN OXYGEN FOR THE SURFACE LAYER OVER THE TIME STEP
C OUSTAR=WIND SHEAR AT THE SURFACE
C OWG=OXYGEN TRANSFER VELOCITY AT THE SURFACE
C U6X=WIND SPEED ADJUSTED BY MULTIPLICATIVE FACTOR
C WINDX=WIND MULTIPLICATION FACTOR (WIND SHELTERING): 0 < XWIND < 1.0 
C
C HSCCDO=HALF SATURATION CONSTANT FOR BIOLOGICAL LIMITATION BY DO IN THE SEDIMENT (MG L**-1)
C Half saturation constant for DO effect on SOD decay (mean=1.0,range=0.1-2.0): D.H. USED HSCCDO=1.5
C
C BIOLC=COEFFICIENT FOR BIOCHEMICAL SEDIMENT OXYGEN UPTAKE (G M**-2 D**-1)
C Oxygen demand of sub-euphotic sediments/g/m2/day (0.02-50.0 range, 5.0 mean): D.H. USED BIOLC=0.2
C USING BIOLC=0.46 FROM BOUFFARD ET AL WAT RESOUR RES (2013).  SHOULD TRY NADER'S VALUE FOR EAGLE LAKE
C OLAK=DISSOLVED OXYGEN profile of lake tile [mg/l] (DO)
C O0=Lake surface DO [mg/l]
C SOD=SEDIMENT OXYGEN DEMAND (G M**2 D**-1), *LAYER DEPTH*TIMESTEP TO GET G M**3
C SEDTEM=sediment temperature multiplier for release rate: D.H. USED SEDTEM=1.08
C OFLUX=OXYGEN FLUX from the bed computed from Fick’s Law J=-KZ dDO/dZ 𝐽 = −𝐾 𝑑𝐷𝑂 ⁄𝑑𝑧.
C KZ=the turbulent diffusivity
C OBOD=CHANGE IN OXYGEN DUE TO BIOCHEMICAL OXYGEN DEMAND (MG L**-1)
C HSCDO=HALF SATURATION CONSTANT FOR BIOLOGICAL LIMITATION BY DO IN THE SEDIMENT (MG L**-1)
C Half saturation constant for DO effect on BOD decay (mean=1.0,range=0.1-2.DELZLK0): D.H. USED HSCDO=1.5
C PB and RSP= Production and Respiration (mg/L) based on Patterson et al, Freshwater Biology (1985) 
C PBS= Maximum potential specific productivity (41.5 mg O2/mg Chla/h based on Patterson et al, 1985)
C THETABO=temperature adjustment for decay of organic matter  
C CONSTBO=temperature-dependent organic decay rate (day**-1) 
C
      IMPLICIT NONE
      INTEGER ::  NLAKMAX_bgc,ILG,IL1,IL2,JMIX,NSTEP
      INTEGER,DIMENSION(ILG) :: NLAK_bgc, NLAK
      REAL,DIMENSION(ILG) :: LKICEH,VA,T0,O0,USTAR,USTARB,Chla
      REAL,DIMENSION(ILG) :: QSWIN,QLWIN,LI0,DELTHRM,UMLGR,SOD,HOD
      REAL,DIMENSION(ILG) :: ZPhotic,FPP,FRESP,HDPTH_bgc,ZPROD,
     1   				HDPTH_eff,RW
      REAL,DIMENSION(ILG) :: FPPW,ODELTA,FatmO2,kge,ZRESP,RW_1,RW_2
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: TLAK_bgc,O2,DOLAK,OFLUX,LI
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: DOLAKback,DOLAKfor,KZ_diff
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: PB,RSP,OBOD,FPPWC,OxCH4,
     1 FRESPWC
      REAL :: DELSKIN,DELZLK,TIMSTP

      REAL, DIMENSION(ILG,NLAKMAX_bgc) :: RSAER,QFLX_bgc

      REAL, DIMENSION(ILG,NLAKMAX_bgc+1):: K_diff

      REAL, DIMENSION(ILG,0:NLAKMAX_bgc+1) :: Area
      REAL,DIMENSION(ILG) :: VOL_eff, VOL_tot, VOL_phot
	  REAL, DIMENSION (ILG) :: DTHERMO

      REAL :: f_sed

      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: zox,Fox_sed
      REAL,DIMENSION(ILG) :: FRESP_ratio,VOL_ox,FRESP_sed
      
C ----* COMMON BLOCKS *------------------------------------------------
C
      REAL DELT
      COMMON /CLASS1/ DELT
C
C ----* LOCAL VARIABLES *---------------------------------------------
C
      INTEGER I,J,K
      REAL OMIX,SODGMD,OBAR,WINDX,HSCCDO,HSCDO,BIOLC,SEDTEM,THETABO,
     1     CONSTBO,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,
     2     OCSTAR,OUSTAR,OWG,U6X,KZ,KZTT,Z,
     3	   PBS,alpha,beta,alphaO2,	
     4	   HO2,PO2,HCH4,PCH4,Patm,Ck1,Ck2,Nk,ScCH4,ScO2
C
C--OXYGEN MODEL COEFFICIENTS
      DATA  Q1,      Q2,      Q3,      Q4,      Q5
     1/     3.35E-2, 9.6E-4,  4.1E-5,  1.1E-1,  6.4022E-5/
      DATA  Q6,        Q7,         Q8,       Q9,        Q10
     1/     1.4101E+5, 5.24838E+3, 7.7117E0, 1.31403E0, 4.593E+1/
      DATA  WINDX, HSCCDO, HSCDO, BIOLC, SEDTEM, THETABO, CONSTBO
     1/     1.0,   1.5,    1.5,   0.462, 1.08,   1.08,    0.005/ 
      DATA  KZ,			PBS,		alpha,		beta
     1/     1.0E-7,		41.5,		0.20,		4.2/
 	    DATA  HO2,	PO2,   HCH4,   PCH4,       Patm
     1/     1.3E-5, 0.2095,1.3E-5, 1800.0E-9  ,101325/	! Atmospheric pressure=101325 (Pa), HO2: Henry's constant for O2 at reference temperature (mol m-3 Pa-1), PO2: Atmospheric partial pressure of O2 (Atmospheric partial pressure)
 	    DATA  Ck1,     Ck2,     Nk,   ScCH4,   ScO2									! Ck1 = 5.75d-6 (m s-1) Ck2 = 5.97d-6 (m s-1)^{-Nk}
     1         /     5.75d-6, 5.97d-7, 1.7  , 677 , 589 /    ! Sc from   "Advances in Quantifying Air-Sea Gas Exchange and Environmental Forcing* " of RikWanninkhof et at 2009,  ! sc other 570,     441/

C--------Atmospheric Fluxes (M. Maisonnier)------------	 
   
	     kge(I) = (Ck1 + Ck2*VA(I)**Nk)*
     1       (sqrt(600/
     2       (13750*(0.10656*exp(-0.0627*TLAK_bgc(I,1))+0.00495))             !schmidt come from "dissolved oxygen model" from stefan and fang 1994
     3       ))		 ! m/s
             
         IF(LKICEH(I).GT. 0.05) THEN
        	 FatmO2(I)=0.0
         ELSE
	     FatmO2(I) = MIN( O2(I,1)*DELZLK ,        !               mg/l * m/s *s/dt = g/m³ * m/s *s/dt = g/(m²dt)
     1        kge(I)*(O2(I,1)-(HO2
     2        *exp(1700.*((1./TLAK_bgc(I,1))-(1./298.15)))
     3        *PO2*Patm*32.))*DELT)
	     O2(I,1)=  O2(I,1)-FatmO2(I)/DELZLK
        ENDIF 

c----------------Bubble---------------------------

        DO J = 1,NLAK_bgc(I)    
           
			if ( O2(I,J) .gt. 
     2        3.*(HO2*exp(1700.*((1./TLAK_bgc(I,J))-(1./298.15)))
     3        *PO2*(Patm+1000.*9.81*(J*DELZLK))*32.) ) THEN

            FatmO2(I) = FatmO2(I) +  O2(I,J) - 
     2        3.*(HO2*exp(1700.*((1./TLAK_bgc(I,J))-(1./298.15)))
     3        *PO2*(Patm+1000.*9.81*(J*DELZLK))*32.)
            O2(I,J) =  
     2        3.*(HO2*exp(1700.*((1./TLAK_bgc(I,J))-(1./298.15)))
     3        *PO2*(Patm+1000.*9.81*(J*DELZLK))*32.)

            endif

        ENDDO    


C-----------PRIMARY PRODUCTION--------------------------

        DO 425 J=1,NLAK_bgc(I)
        	HOD(I)=0.0
        	SOD(I)=0.0
        	Chla(I)=0.0
        	OBOD(I,J)=0.0
            FPPWC(I,J)=0.0
            RSP(I,J)=0.0
            PB(I,J)=0.0 
           
 425    CONTINUE  

        
        ZPhotic(I) = ZPROD(I)
        
        FPPW(I)=FPP(I)*32./12./1000.* ! FPPW=FPP=BC(I)*Pchl*M(I) applied to ZPhotic
     1          VOL_tot(I)/VOL_phot(I)           

        DO 415 J=1,NLAK_bgc(I)                 
   
            Z=DELZLK*J
                   
            IF ( Z .LE. ZPhotic(I)) THEN	! Production applied to ZPhotic        !!!! 
                   FPPWC(I,J)=FPPW(I) !mgO2 L-1
    		ELSE
                   FPPWC(I,J)=0.0
            ENDIF
 
            O2(I,J)=O2(I,J)+FPPWC(I,J)               
														
 415     CONTINUE
         
C-----------RESPIRATION---------------------------------

          HDPTH_eff(I) = int(HDPTH_bgc(I)/DELZLK)*DELZLK
            ZRESP(I) = max(ZPhotic(I),HDPTH_eff(I))
            VOL_eff(I) = Area(I,1)*DELZLK
        DO 106 J=2,int(ZRESP(I)/DELZLK)
           VOL_eff(I) = VOL_eff(I) + Area(I,J)*DELZLK
 106    CONTINUE

C-----------In the sediment------------------------------
        DO J=1,NLAK_bgc(I)
            Fox_sed(I,J) = FRESP(I)*(1.-f_sed)*
     1                   (Area(I,J)-Area(I,J+1))*zox(I,J)
     1                   /(Area(I,J)*DELZLK)   
     1                   *(32./12./1000) ! g O2/m3

            Fox_sed(I,J) = min(O2(I,J),Fox_sed(I,J))
                
            O2(I,J) = O2(I,J) - Fox_sed(I,J)
            
        enddo
C----------O2 consumption in water column------------------
        RW_1(I)= (20./100)*FRESP(I)*(1.-f_sed)*(32./12.)/1000.
        RW_2(I)= (80./100.)*FRESP(I)*(1.-f_sed)*(32./12.)/1000.
   
         
        RW_1(I) = RW_1(I)!*VOL_tot(I)/VOL_phot(I)  
        RW_2(I) = RW_2(I)*VOL_tot(I)/VOL_eff(I)
        !    PRINT*, RW(I)			
        DO 105 J = 1,NLAK_bgc(I)
            IF (J*DELZLK .le. ZRESP(I) ) THEN
            FRESPWC(I,J) = RW_2(I)
            ELSE
            FRESPWC(I,J) = 0.
            ENDIF
            If  (J .LE. NLAK_bgc(I)) THEN   !(J*DELZLK .le. ZPhotic(I)) THEN
                FRESPWC(I,J) = FRESPWC(I,J) + RW_1(I)
            endif
        
        FRESPWC(I,J) = FRESPWC(I,J)-  FRESP(I)*(1.-f_sed)*
     1                   (Area(I,J)-Area(I,J+1))*zox(I,J)
     1                   /(Area(I,J)*DELZLK)   
     1                   *(32./12./1000) ! g O2/m3

        FRESPWC(I,J) = min(O2(I,J),FRESPWC(I,J))
        
 105    CONTINUE
        RW(I) = sum(FRESPWC(I,1:NLAK_bgc(I))
     1              *AREA(I,1:NLAK_bgc(I))*DELZLK)/VOL_tot(I)
            
            O2(I,1:NLAK_bgc(I))=O2(I,1:NLAK_bgc(I))
     1	                         -FRESPWC(I,1:NLAK_bgc(I))


C-----------CONTROLE---------------------------------

         DO 416 J=1,NLAK_bgc(I)
            	IF (O2(I,J) .LT. 0) THEN
                   O2(I,J)=0.d0
                   PRINT*, 'O2 < 0  ->  O2 = 0'
            	ENDIF	
 416     CONTINUE            

C----------------Transport------------------------


      CALL transport(I,ILG,O2,DELT,TLAK_bgc,
     1               NLAK_bgc,NLAKMAX_bgc,HDPTH_bgc,JMIX,VA,Area,
     2               K_diff,DTHERMO)

C------------------------

        DO J = 1,NLAK_bgc(I)

        IF (ABS(O2(I,J)) .gt. HUGE(1.0)) THEN
            PRINT*, 'O2' , I, J , O2(I,J)
        ENDIF

        IF (ABS(FRESPWC(I,J)) .gt. HUGE(1.0)) THEN
            PRINT*, 'FRESPWC',I,J,FRESPWC(I,J),Rw(I),FRESP(I),f_sed
        ENDIF

        IF (ABS(FPPWC(I,J)) .gt. HUGE(1.0)) THEN
            PRINT*, 'FPPWC' , I, J , FPPWC(I,J), FPPW(I),FPP(I), 
     1           VOL_tot(I),VOL_phot(I),VOL_phot(I)/VOL_tot(I)
        ENDIF


        ENDDO

      
      RETURN
      END

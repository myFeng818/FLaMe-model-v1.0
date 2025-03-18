       SUBROUTINE oxidation(NSTEP,ILG,I,NLAKMAX_bgc,NLAK_bgc,O2,CH4,
     4                     TLAK_bgc,OxCH4)

C============================================================================
C    OXIDATION SUBROUTINE
C               
C     * AUTHOR: Manon Maisonnier, Maoyuan Feng, ... ,Pierre Regnier.         
C     * Final revision: March 2025
C     * 
C     * Reference: A new biogeochemical modelling framework (FLaMe v1.0) 
C     *            for lake methane emissions on the regional scale: 
C     *            Development and application to the European domain
C     *         
C     * Correspondance: maoyuan.feng@ulb.be; maoyuanfeng93@gmail.com 
C               
C----------------------------------------------------------------------

      IMPLICIT NONE
C
C ----* LAKE MODEL VARIABLES *----------------------------------------
C
      INTEGER :: NLAKMAX_bgc,NSTEP
      INTEGER,DIMENSION(ILG) :: NLAK_bgc
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: TLAK_bgc
      
C
C ----* INTEGER CONSTANTS *-------------------------------------------
C
      INTEGER :: ILG,IL1,IL2      
      
C ----* COMMON BLOCKS *------------------------------------------------
 
      REAL DELT,TFREZ,GRAV,DELZLK
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /LAKECON/ DELZLK
      
C ----* LOCAL VARIABLES *---------------------------------------------

      INTEGER :: I,J

C=======================================================================
C     * BIOCHEMICAL VARIABLES (Manon Maisoniier)

      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: CH4,O2
      
C ----* CH4 Oxydation VARIABLES *----------------------------------

      REAL,DIMENSION (ILG,NLAKMAX_bgc) :: OxCH4
      
	  REAL,PARAMETER :: V_max = 1.*0.5E-6  ! mol CH4 m3 s-1   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL :: KCH4, KO2
      REAL :: Q10,Tr
      
          DATA   KCH4,       KO2
     1     /    3.75E-2 , 2.1E-2 / ! molCH4 m-3 , molO2 m-3 , molCH4 m-3 s-1
          DATA   Q10 , Tr
     1     /     1.2  , 283.15 / ! Q10 in [1.4,3.5]  q10 : Thottathil et al. 2019 (deduce)
          
      
C=======================================================================

c --------------------------* Main *-------------------------------          

c      DO 600 I=IL1,IL2    
        DO 610 J=1,NLAK_bgc(I)      
        IF (O2(I,J) .GT. 0. .AND. CH4(I,J) .GT. 0.) THEN                      
            OxCH4(I,J) =V_max*(Q10**(TLAK_bgc(I,J)-Tr)/10)*
     1          ((O2(I,J)/32)/((O2(I,J)/32)+(KO2)))*
     2          ((CH4(I,J)/16)/((CH4(I,J)/16)+(KCH4)))
     3          *16.*DELT                                   ! g/m³/s * s/dt
     
     	ELSE
     		OXCH4(I,J)=0.d0
        END IF  

        OxCH4(I,J) = min(CH4(I,J),(1./2.)*O2(I,J)*16./32.,OxCH4(I,J))

           CH4(I,J) = CH4(I,J) - OxCH4(I,J)   !  MAX(0.,CH4(I,J) - OxCH4(I,J)*16*DELT)
           O2(I,J)  = O2(I,J) - 2.*OxCH4(I,J)*32./16.    !  MAX(0., O2(I,J) - 2*OxCH4(I,J)*32*DELT)

          ! PRINT*, OxCH4(I,J)
           
610     CONTINUE  
c     600   CONTINUE
        
      END SUBROUTINE

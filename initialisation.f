      SUBROUTINE initialisation(ILG,IL1,IL2,NLAK_bgc,NLAKMAX_bgc,
     1                          Clabile,CH4,O2,Clabs_1,Clabs_2,Clabs_3,
     2                          FsCH4_C_1,FsCH4_C_2,FsCH4_C_3,
     3                          Fsed_1,Fsed_2,Fsed_3,Fibs)	     
C============================================================================
C            Initialisation of Species  (Manon Maisonnier) 
C===========================================================================
C----------------------------------------------------------------------

      IMPLICIT NONE
C
C ----* LAKE MODEL VARIABLES *----------------------------------------
C
      INTEGER NLAKMAX_bgc
      INTEGER,DIMENSION(ILG) :: NLAK_bgc

C
C ----* INTEGER CONSTANTS *-------------------------------------------
C
      INTEGER :: ILG,IL1,IL2            
      
C ----* LOCAL VARIABLES *---------------------------------------------

      INTEGER :: I

C=======================================================================
C     * BIOCHEMICAL VARIABLES (Manon Maisoniier)

      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: CH4,O2,Fibs
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Clabs_1,Clabs_2,Clabs_3
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FsCH4_C_1,FsCH4_C_2,FsCH4_C_3
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Fsed_1,Fsed_2,Fsed_3

      REAL, DIMENSION (ILG) :: Clabile
      
C ===========================================================================

      DO I=IL1,IL2
         CH4(I,1:NLAK_bgc(I)) = 0.d0
         O2(I,1:NLAK_bgc(I)) = 0.d0
         Clabile(I) = 0.d0

         Clabs_1(I,1:NLAK_bgc(I)) = 0.d0
         Clabs_2(I,1:NLAK_bgc(I)) = 0.d0
         Clabs_3(I,1:NLAK_bgc(I)) = 0.d0

         FsCH4_C_1(I,1:NLAK_bgc(I)) = 0.d0
         FsCH4_C_2(I,1:NLAK_bgc(I)) = 0.d0
         FsCH4_C_3(I,1:NLAK_bgc(I)) = 0.d0

         Fsed_1(I,1:NLAK_bgc(I)) = 0.d0
         Fsed_2(I,1:NLAK_bgc(I)) = 0.d0
         Fsed_3(I,1:NLAK_bgc(I)) = 0.d0

         Fibs(I,1:NLAK_bgc(I)) = 0.d0


      ENDDO   
      
      
      RETURN
                
    
      END

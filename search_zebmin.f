      SUBROUTINE search_zebmin(J,zebmin,P0,alpha0,Hs,
     1                         CH4,TLAK,TIMSTP,DELZLK)

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
C=============================================================================

      INTEGER :: ok,J
      REAL :: DELZLK,TIMSTP
      REAL :: alpha0, Hs
      REAL :: zebmin, P0, CH4, TLAK

      REAL :: z1, z2, z3,z_err, PH2O, Pres, K_diff_CH4

      REAL :: HCH4, Patm

    

      Patm = 101325.   ! Pa
      
      
      HCH4 = 1.4d-5*exp( (1600.) * ( (1./TLAK) - (1./298.15) ) )
      

      K_diff_CH4 = (9.798e-10 + 2.986e-11*(TLAK-273.15)    ! m²s⁻¹ 
     1               + 4.381e-13*(TLAK-273.15)**2)
     2            *(1.-log(0.5**2))
      

        IF (TLAK .gt. 273.15 ) THEN
           PH2O = 610.94
     1          *exp((17.625*(TLAK-273.15))
     2               /(TLAK-273.15+243.04))
        ELSE
           PH2O = 611.21
     1          *exp((22.587*(TLAK-273.15))
     2               /(TLAK+0.71))
        ENDIF
        
           Pres = 1000.*9.81*(J*DELZLK) + (1.-0.78)*Patm - PH2O
     1            -(CH4/16.)/HCH4     

           
           z1 = 0.d0
           z2 = Hs
           z_err = 0.005
           ok = 0
           Zebmin = Hs
           z3 = Hs

         IF (Fzebmin(Hs)*Fzebmin(0.) > 0.) THEN
            IF (Fzebmin(Hs) < 0.) THEN 
                z3 = Hs
            ELSE 
                z3 = 0.
            ENDIF
            ok = 1
         ENDIF

11       IF (ok .eq. 0) THEN
           IF (Fzebmin(z1)*Fzebmin(z2) .lt. 0.) THEN
              z3 = (z1+z2)/2.
           ELSEIF (Fzebmin(z1)*Fzebmin(z2) == 0.) THEN
               IF ( Fzebmin(z1) == 0. ) THEN
                    z3 = z1
               ELSE
                    z3 = z2
               ENDIF
               ok = 1
           ELSE 
             PRINT*, 'bug zebmin research'
             STOP
           ENDIF
         ENDIF

         IF (ok .eq. 0) THEN
            IF ((ABS(z3-z1) < z_err) .OR. (ABS(z3-z1) < z_err)) THEN
                ok = 1
            ELSE
				IF (Fzebmin(z1)*Fzebmin(z3) .lt. 0.) THEN
                 z2 = z3
				ELSE
                 z1 = z3
				ENDIF
                GOTO 11
            ENDIF
         ENDIF
          
       
        zebmin = z3
         
      RETURN
      
      CONTAINS

      REAL FUNCTION Fzebmin(z_)
      IMPLICIT NONE
      REAL :: z_
      Fzebmin = 1. - z_*alpha0*exp(-z_*alpha0) - exp(-z_*alpha0)
     1     -(alpha0**2)*K_diff_CH4*Pres*HCH4/
     2       (P0/(1000.*16.*TIMSTP))
      END FUNCTION Fzebmin

      END

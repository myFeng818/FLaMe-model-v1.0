      SUBROUTINE transport(I,ILG,X,TIMSTP,TLAK_bgc,NLAK_bgc,
     1      NLAKMAX_bgc,HDPTH_bgc,JMIX,VA,Area,K_diff,DTHERMO)

	 

C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     Transport routine   (MANON MAISONNIER)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------

      
C     ----------------------------------------------------------------
C     Variable declaration
C     ----------------------------------------------------------------

      IMPLICIT NONE
      

C     Basic 
      
      INTEGER :: NLAKMAX_bgc, JMIX
      INTEGER, DIMENSION(ILG) :: NLAK_bgc
      REAL, DIMENSION(ILG) :: HDPTH_bgc
      REAL, DIMENSION(ILG,NLAKMAX_bgc) :: TLAK_bgc
      REAL,DIMENSION(ILG) :: VA
      REAL :: TIMSTP
      INTEGER :: I,J,ILG,IL1,IL2  !,X_number
  
      REAL, DIMENSION(ILG,0:NLAKMAX_bgc+1) :: Area
      
C     Transport
      
      REAL, DIMENSION(ILG,NLAKMAX_bgc) :: XDELTAback,XDELTAfor,XFLUX,X
      REAL, DIMENSION(ILG) :: X_mix

C     K_diff
      REAL, DIMENSION (ILG) :: DTHERMO
      INTEGER :: DTH

      REAL, DIMENSION(ILG,NLAKMAX_bgc+1) :: K_diff
      REAL :: K_diff_0, N_0

c      REAL :: rho_JMIX_plus, rho_JMIX_minus, rho_JMIX
      REAL :: mixing = 0.15             ! mixing efficiency    include in [0.05 ; 0.25]
      REAL :: C10 = 0.001               ! wind stress coefficient set for w10 <= 7 m s-1
      REAL :: k_karman = 0.41           ! Karman constant
      REAL :: rhoa = 1.225              ! kg/m3  rho aire      
      REAL :: Tmp = 273.15   ! K    melting point temperature

      REAL, DIMENSION(ILG,NLAKMAX_bgc) :: rho
	
      REAL, DIMENSION(ILG) ::   X_before, X_after
	  
      REAL :: Z2,ZBOT,ZMIX,Z

      REAL :: TKECN,TKECF,TKECE,TKECS,TKECL,HDPTHMIN,
     1   TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,DUMAX
	 
	 
	     COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,        
     2                 TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,
     3                 DHMAX,TKECL,DUMAX

C     ----------------------------------------------------------------
C     Some files
C     ----------------------------------------------------------------


      
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     Main
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------

       X_before(I)= sum(X(I,1:NLAK_bgc(I))*Area(I,1:NLAK_bgc(I))*DELZLK)

       Z2=DELZLK+DELSKIN
       ZBOT=DELSKIN+DELZLK*NLAK_bgc(I)
       IF (HDPTH_bgc(I) .LT. Z2) THEN	!fully stratified: no mixing
         JMIX=1
         ZMIX=DELZLK+DELSKIN
        ELSE IF (HDPTH_bgc(I) .GE. ZBOT) THEN	!fully mixed
         JMIX=NLAK_bgc(I)
         ZMIX=ZBOT
        ELSE
         DO 510, J=2,NLAK_bgc(I)
           ZMIX=DELSKIN + DELZLK*J
           JMIX=J
           IF (HDPTH_bgc(I) .LT. ZMIX) EXIT
510      CONTINUE
       ENDIF
			

C     ----------------------------------------------------------------
C     K_diff
C     ----------------------------------------------------------------

c	CALL density(I,ILG,NLAKMAX_bgc,TLAK_bgc,NLAK_bgc,rho)


         
        DTH = max(NINT(DTHERMO(I)/DELZLK),1)


        K_diff(I,1:JMIX) = 1.E-4
       IF (JMIX .lt. NLAK_bgc(I)) THEN         
         K_diff(I,JMIX+1:NLAK_bgc(I)+1) = 1.E-7
       ENDIF   
       IF ( (JMIX+ DTH ).lt. (NLAK_bgc(I)+1) ) THEN
         K_diff(I,JMIX+1+DTH :NLAK_bgc(I)+1) = 1.E-5
       ENDIF
         
        
      
C     ----------------------------------------------------------------
C     Transport Scheme
C     ----------------------------------------------------------------

        ! K_diff(I,1:NLAK_bgc(I)+1) = 1.e-6


       DO 420 J=1,NLAK_bgc(I)

        IF (J==1) THEN
         XDELTAback(I,J)=0.d0
         ELSE
          XDELTAback(I,J) = X(I,J)-X(I,J-1)
        ENDIF

        IF (J==NLAK_bgc(I)) THEN
         XDELTAfor(I,J)=0.
        ELSE
          XDELTAfor(I,J) = X(I,J+1)-X(I,J)
        ENDIF
            ! XFLUX(I,J)=( K_diff(I,J+1)*XDELTAfor(I,J)*Area(I,J+1)*DELZLK
       ! 1                 - K_diff(I,J)*XDELTAback(I,J)*Area(I,J)*DELZLK )
       ! 2           * TIMSTP/(2*DELZLK*DELZLK)
        XFLUX(I,J)=( K_diff(I,J+1)*XDELTAfor(I,J) *Area(I,J+1)*DELZLK
     1              - K_diff(I,J)*XDELTAback(I,J) *Area(I,J)*DELZLK )
     2              *( 1/(DELZLK*DELZLK))/(Area(I,J)*DELZLK)
 
420    CONTINUE

       DO 103 J=1,NLAK_bgc(I)
          X(I,J) = max(0., X(I,J)
     1               + XFLUX(I,J)*TIMSTP )
103    CONTINUE

c       PRINT*,  'back : ',(XDELTAback(1,J),J=1,NLAK_bgc(I))

c       PRINT*,  'for : ',(XDELTAfor(1,J),J=1,NLAK_bgc(I))
       
c       PRINT*,  'stabilité : ',
c     1  (K_diff(1,J)*TIMSTP/(DELZLK*DELZLK),J=1,NLAK_bgc(I)) 

c       PRINT*,  'flux : ',(XFLUX(1,J),J=1,NLAK_bgc(I))
       
c       PRINT*, 'stuff : ', (K_diff(1,J),J=1,NLAK_bgc(I)+1) ,
c     1          DELZLK, 1/(DELZLK*DELZLK)  , TIMSTP ,
c     2        TIMSTP/(DELZLK*DELZLK) 

c --------- MELANGE MIXING  ----------

         X_mix(I) = 0.
        ! IF (((TLAK_bgc(I,1) .gt. 3.7+Tmp).AND.(TLAK_bgc(I,1) 
       ! 1         .lt. 4.2+Tmp)) ) THEN
        X_mix(I) = SUM(X(I,1:JMIX)*AREA(I,1:JMIX))/SUM(AREA(I,1:JMIX))
        X(I,1:JMIX) = X_mix(I)
         ! ENDIF



C     ----------------------------------------------------------------


c        X_after(I) = sum (X(I,1:NLAK_bgc(I))*Area(I,1:NLAK_bgc(I))*DELZLK )

c        IF  (X_before(I) .ne. X_after(I)) THEN
c         If (X_before(I) == 0.0 ) THEN
c            X(I,1:NLAK_bgc(I)) = 0.0
c         ELSE
c         X(I,1:NLAK_bgc(I)) = X(I,1:NLAK_bgc(I))  * X_before(I)/X_after(I)
c         ENDIF
c        ENDIF

      RETURN
      
      END SUBROUTINE 



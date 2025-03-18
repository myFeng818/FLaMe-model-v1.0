      SUBROUTINE geometry_bgc(ILG,I,DELZLK,HLAK,HLAK_bgc,
     1						HDPTH, HDPTH_bgc, TLAK, TLAK_bgc, 
     1 				      NLAKMAX,NLAK,NLAKMAX_bgc,NLAK_bgc,
     2                    LLAK,Area,NSTEP)
	 
	 
      INTEGER :: ILG, I, IL1, IL2, NLAKMAX, NLAKMAX_bgc
      INTEGER, DIMENSION(ILG) :: NLAK, NLAK_bgc
      REAL :: DELZLK
      REAL,DIMENSION(ILG) :: LLAK
      REAL,DIMENSION(ILG,0:NLAKMAX_bgc+1) :: Area


      REAL,DIMENSION(ILG) :: HLAK,HLAK_bgc

!      INTEGER :: J , DN
      REAL :: Pi = 3.1416
	  
      REAL,DIMENSION(ILG,NLAKMAX) :: TLAK
      REAL,DIMENSION(ILG,NLAKMAX_bgc) :: TLAK_bgc
      REAL,DIMENSION(ILG) :: HDPTH , HDPTH_bgc
	 
      REAL :: V1, Vv, Vvv
      INTEGER :: n , l, j, indice

	  
      IF ( NSTEP == 1 ) THEN
	  
	        HLAK_bgc(I) = 2.*HLAK(I) ! 1.7302*HLAK(I)  !Pi/2)*HLAK(I) 
            NLAK_bgc(I)=NINT(HLAK_bgc(I)/DELZLK)
		
	  
	       CALL area_layer(ILG,I,NLAKMAX_bgc,NLAK_bgc,DELZLK,
     1                 LLAK,Area)
        
      ENDIF
	  
	  
      TLAK_bgc(I,1:NLAK(I))=TLAK(I,1:NLAK(I))
      TLAK_bgc(I,NLAK(I):NLAK_bgc(I))=TLAK(I,NLAK(I))


	
	  
	  
      HDPTH_bgc(I) = max(HDPTH(I),0.5)
	  
      IF (HDPTH(I) .gt. HLAK(I)- DELZLK ) THEN
        HDPTH_bgc(I) = HLAK_bgc(I) 
      ENDIF

       IF ((TLAK_bgc(I,1) .gt. 3.8+Tmp).AND.(TLAK_bgc(I,1)               ! 3.98275
     1        .lt. 4.1+Tmp)) THEN
            HDPTH_bgc(I) = JMIX*DELZLK
       ENDIF
	
  
  
      END
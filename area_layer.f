      SUBROUTINE area_layer(ILG,I,NLAKMAX_bgc,NLAK_bgc,DELZLK,
     1                      LLAK,Area)

      INTEGER :: ILG, I,IL1,IL2,NLAKMAX_bgc
      INTEGER, DIMENSION(ILG) :: NLAK_bgc
      REAL :: DELZLK
      REAL,DIMENSION(ILG) :: LLAK
      REAL,DIMENSION(ILG,0:NLAKMAX_bgc+1) :: Area

      INTEGER :: J 
      REAL, PARAMETER :: Pi = 3.1416
      REAL,PARAMETER :: g = 2.0

      
c      OPEN(UNIT=10,FILE='area.txt')
      
      
c      DO 100 I=IL1,IL2

      Area(I,:) = 0.d0
      Area(I,0) = 1.d6*(LLAK(I)/2)*(LLAK(I)/2)/Pi      ! 0 index is for the skin
      Area(I,1) = Area(I,0)
      
         DO J=2,NLAK_bgc(I)
            Area(I,J)=Area(I,0)
     1              *(1./(2.*g))*(g+(ATANH(
     2              ((1./2.)-(REAL(J-1)/NLAK_bgc(I)))
     3               *2.*tanh(g)
     4               ))) 

c     1        !*COS(((J-1)*((DELZLK)/(DELZLK*NLAK_bgc(I)))**2)*(Pi/2))

         ENDDO
            Area(I,NLAK_bgc(I)+1)= 0.0

          

c         WRITE(10,*)  I, (Area(I,J),J=0,NLAK_bgc(I))     
	
c       PRINT*, I, NLAKMAX_bgc, NLAK_bgc(I), 'Area', 
c     1             (Area(I,j), j = 0,NLAK_bgc(I)+1)

c 100  ENDDO

      END

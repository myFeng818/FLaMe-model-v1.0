      SUBROUTINE density(I,ILG,NLAKMAX,TLAK,NLAK,rho)
      
      INTEGER :: I,ILG,NLAKMAX
      INTEGER, DIMENSION(ILG) :: NLAK
      REAL :: Tmp
      REAL, DIMENSION (ILG,NLAKMAX) :: rho, TLAK
      
      Tmp = 273.15  ! K

	DO J=1,NLAK(I)
         rho(I,J) = (999.83952 + 16.945176*(TLAK(I,J)-Tmp)
     1       -7.9870401d-3*(TLAK(I,J)-Tmp)**2
     2        - 46.170461d-6*(TLAK(I,J)-Tmp)**3
     3        +105.56302d-9*(TLAK(I,J)-Tmp)**4
     4        -280.54253d-12*(TLAK(I,J)-Tmp)**5)
     5        /(1+16.897850d-3*(TLAK(I,J)-Tmp))
	ENDDO
      
      END

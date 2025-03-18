      SUBROUTINE Writing(NSTEP,NSTEP_write,ILG,IL1,IL2,JL, 
     1     NLAKMAX,NLAK,IYEAR,IDAY,IHOUR,IMIN,LKICEH,
     3     TLAK,HDPTH,TSED,
     1     Clabile,FPP,FRESP,FBUR,Area,O2,CH4,
     2     FatmO2,Fw,FPPWC,FSAER,
     6     FatmCH4,FsCH4,FsDiff,FsEbul,
     7     Clabs_1,Clabs_2,Clabs_3,FsCH4_C_1,FsCH4_C_2,FsCH4_C_3,
     3     Fsed_1,Fsed_2,Fsed_3,
     4     FswDiff,OxCH4, K_diff,  FatmEbul ,
     1     LKICEH_m,
     3     TLAK_m,HDPTH_m,TSED_m,
     1     Clabile_m,FPP_m,FRESP_m,FBUR_m,O2_m,CH4_m,
     2     FatmO2_m,Fw_m,FPPWC_m,FSAER_m,
     6     FatmCH4_m,FsCH4_C_m,FsDiff_m,
     1     FsEbul_m,
     7     Clabs_1_m,Clabs_2_m,Clabs_3_m,FsCH4_C_1_m,FsCH4_C_2_m,
     3     FsCH4_C_3_m,Fsed_1_m,Fsed_2_m,Fsed_3_m,
     4     FswDiff_m,OxCH4_m,K_diff_m, FatmEbul_m ,
     5     ZPhotic,ZRESP,FRESPWC,ZPhotic_m,ZRESP_m,FRESPWC_m,
     6	   NLAKMAX_bgc,NLAK_bgc,TLAK_bgc,HDPTH_bgc,
     7	   TLAK_bgc_m,HDPTH_bgc_m, Fstor_m , DELT_stor_m , 
     8     DELT_ice_m, x_per_m )
C============================================================================
C            Writing  (Manon Maisonnier) 
C============================================================================
C
C  Write the value by week of key variables
C
C----------------------------------------------------------------------


      IMPLICIT NONE

C
C ----* LAKE MODEL VARIABLES *----------------------------------------
C
      INTEGER :: NLAKMAX,N_step_met,NSTEP_write,NLAKMAX_bgc
      INTEGER,DIMENSION(ILG) :: NLAK,NLAK_bgc
      REAL,DIMENSION(ILG) :: HDPTH,LKICEH,TSED
      REAL,DIMENSION(ILG,NLAKMAX) :: TLAK
 
 
      REAL,DIMENSION (ILG,NLAKMAX_bgc) :: TLAK_bgc
      REAL,DIMENSION (ILG) :: HDPTH_bgc
	

C                                                                       
C ----* COMMON BLOCK PARAMETERS *--------------------------------------
C
      REAL DELT,TFREZ                     
      REAL TKECN,TKECF,TKECE,TKECS,TKECL,HDPTHMIN,
     1     TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,DUMAX
C
C ----* CLASS COMMON BLOCKS *------------------------------------------
C
      COMMON /CLASS1/ DELT,TFREZ       
      COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,        
     2                 TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,
     3                 DHMAX,TKECL,DUMAX


C
C ----* LOCAL VARIABLES *---------------------------------------------
C
      INTEGER I,J

C
C ----* INTEGER CONSTANTS *-------------------------------------------
C
      INTEGER :: ILG,IL1,IL2,JL,IALS 
      INTEGER :: NSTEP,IYEAR,IDAY,IHOUR,IMIN

C=======================================================================
C     * BIOCHEMICAL VARIABLES (Manon Maisoniier)

      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: CH4,O2
      REAL, DIMENSION (ILG,0:NLAKMAX_bgc+1) :: Area
      REAL, DIMENSION (ILG) :: Clabile,FPP,FRESP,FBUR
      REAL :: Pchl
      REAL, DIMENSION (ILG) :: FatmCH4,H_s,limP
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FsCH4,FsDiff,FsEbul,zebmin
      REAL, DIMENSION (ILG) :: FatmO2,Fw
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FPPWC,FSAER      
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Clabs_1,Clabs_2,Clabs_3
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FsCH4_C_1,FsCH4_C_2,
     1   FsCH4_C_3
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Fsed_1,Fsed_2,Fsed_3
      REAL, DIMENSION (ILG,NLAKMAX_bgc) ::  FswDiff,OxCH4


      REAL, DIMENSION(ILG,NLAKMAX_bgc+1) :: K_diff 

      REAL, DIMENSION(ILG) :: ZPhotic,ZRESP,ZPhotic_m,ZRESP_m
      REAL, DIMENSION(ILG,NLAKMAX_bgc) :: FRESPWC,FRESPWC_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FatmEbul , FatmEbul_m
     
C
C ----* Variable printed *--------------------------------------------
C

      REAL,DIMENSION(ILG) :: HDPTH_m,LKICEH_m,TSED_m!,HDPTH_p
      REAL,DIMENSION(ILG,NLAKMAX) :: TLAK_m

      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: O2_m,CH4_m,
     2     FPPWC_m,FSAER_m,
     6     FsCH4_C_m,FsDiff_m,
     1     FsEbul_m,zebmin_m
      REAL, DIMENSION (ILG) :: Clabile_m,FPP_m,FRESP_m,FatmO2_m,
     1                         FatmCH4_m,Fw_m,FBUR_m
      
      REAL, DIMENSION(ILG) :: Fstor_m,DELT_stor_m

      REAL :: DELT_m      
      
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Clabs_1_m,Clabs_2_m,
     1  Clabs_3_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: FsCH4_C_1_m,FsCH4_C_2_m,
     1  FsCH4_C_3_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) :: Fsed_1_m,Fsed_2_m,Fsed_3_m
      REAL, DIMENSION (ILG,NLAKMAX_bgc) ::  FswDiff_m,OxCH4_m


      REAL, DIMENSION(ILG,NLAKMAX_bgc+1):: K_diff_m

      INTEGER :: bis
      INTEGER :: MONTH , X , DAY_MONTH_X , DAY_MONTH_X_MAX , 
     1           N_DELT_MONTH_X , DAY_MONTH_X_FIRST

	 
      REAL,DIMENSION (ILG,NLAKMAX_bgc) :: TLAK_bgc_m
      REAL,DIMENSION (ILG) :: HDPTH_bgc_m
      


      REAL :: Tmp = 273.15   ! K    melting point temperature
	 
      INTEGER :: x_per_m
      REAL,DIMENSION (ILG) :: DELT_ice_m

C -----------------------------------------------------------------------------
C ------------------------  WRITING (Manon Maisonnier) ------------------------
C -----------------------------------------------------------------------------

       !PRINT*, 'writing', ZPhotic(6), ZRESP(6)

C
C     -------------------- DAY BY MONTH -------------------------------------------
C

C --- Bisextile year ? ---

      SELECT CASE (MODULO(IYEAR,4))
         CASE (0)
         
            SELECT CASE (MODULO(IYEAR,100))
               CASE(0)
               
                  SELECT CASE (MODULO(IYEAR,400))
                     CASE(0)
                        bis = 1
                     CASE DEFAULT
                        bis = 0
                  END SELECT
         
               CASE DEFAULT
                  bis = 1
            END SELECT
               
         CASE DEFAULT
            bis = 0
      END SELECT


C --- number of day in this month/x_per_m ---

      IF (x_per_m == 1) THEN 

      SELECT CASE (bis)
         CASE (0)

            SELECT CASE (IDAY)
               CASE(1:31)
                  MONTH = 1
                  X = 1
                  DAY_MONTH_X = IDAY
                  DAY_MONTH_X_FIRST = 1
                  DAY_MONTH_X_MAX = 31
               CASE(32:59)
                  MONTH = 2
                  X = 1
                  DAY_MONTH_X = IDAY-31
                  DAY_MONTH_X_FIRST = 32
                  DAY_MONTH_X_MAX = 28
               CASE(60:90)
                  MONTH = 3
                  X = 1
                  DAY_MONTH_X = IDAY-59
                  DAY_MONTH_X_FIRST = 60
                  DAY_MONTH_X_MAX = 31
               CASE(91:120)
                  MONTH = 4
                  X = 1
                  DAY_MONTH_X = IDAY-90
                  DAY_MONTH_X_FIRST = 91
                  DAY_MONTH_X_MAX = 30
               CASE(121:151)
                  MONTH = 5
                  X = 1
                  DAY_MONTH_X = IDAY-120
                  DAY_MONTH_X_FIRST = 121
                  DAY_MONTH_X_MAX = 31
               CASE(152:181)   
                  MONTH = 6
                  X = 1
                  DAY_MONTH_X = IDAY-151
                  DAY_MONTH_X_FIRST = 152
                  DAY_MONTH_X_MAX = 30
               CASE(182:212)
                  MONTH = 7
                  X = 1
                  DAY_MONTH_X = IDAY-181
                  DAY_MONTH_X_FIRST = 182
                  DAY_MONTH_X_MAX = 31
               CASE(213:243)
                  MONTH = 8
                  X = 1
                  DAY_MONTH_X = IDAY-212
                  DAY_MONTH_X_FIRST = 213
                  DAY_MONTH_X_MAX = 31
               CASE(244:273)
                  MONTH = 9
                  X = 1
                  DAY_MONTH_X = IDAY-243
                  DAY_MONTH_X_FIRST = 244
                  DAY_MONTH_X_MAX = 30
               CASE(274:304)
                  MONTH = 10
                  X = 1
                  DAY_MONTH_X = IDAY-273
                  DAY_MONTH_X_FIRST = 274
                  DAY_MONTH_X_MAX = 31
               CASE(305:334)
                  MONTH = 11
                  X = 1
                  DAY_MONTH_X = IDAY-304
                  DAY_MONTH_X_FIRST = 305
                  DAY_MONTH_X_MAX = 30
               CASE(335:365)
                  MONTH = 12
                  X = 1
                  DAY_MONTH_X = IDAY-334
                  DAY_MONTH_X_FIRST = 335
                  DAY_MONTH_X_MAX = 31
               CASE DEFAULT
                  PRINT*, 'error day'
                  STOP
            END SELECT
                  
         CASE (1)

            SELECT CASE (IDAY)
               CASE(1:31)
                  MONTH = 1
                  X = 1
                  DAY_MONTH_X = IDAY
                  DAY_MONTH_X_FIRST = 1
                  DAY_MONTH_X_MAX = 31
               CASE(32:60)
                  MONTH = 2
                  X = 1
                  DAY_MONTH_X = IDAY-31
                  DAY_MONTH_X_FIRST = 32
                  DAY_MONTH_X_MAX = 29
               CASE(61:91)
                  MONTH = 3
                  X = 1
                  DAY_MONTH_X = IDAY-60
                  DAY_MONTH_X_FIRST = 61
                  DAY_MONTH_X_MAX = 31
               CASE(92:121)
                  MONTH = 4
                  X = 1
                  DAY_MONTH_X = IDAY-91
	              DAY_MONTH_X_FIRST = 92
                  DAY_MONTH_X_MAX = 30
               CASE(122:152)
                  MONTH = 5
                  X = 1
                  DAY_MONTH_X = IDAY-121
                  DAY_MONTH_X_FIRST = 122
                  DAY_MONTH_X_MAX = 31
               CASE(153:182)   
                  MONTH = 6
                  X = 1
                  DAY_MONTH_X = IDAY-152
		          DAY_MONTH_X_FIRST = 153
                  DAY_MONTH_X_MAX = 30
               CASE(183:213)
                  MONTH = 7
                  X = 1
                  DAY_MONTH_X = IDAY-182
		          DAY_MONTH_X_FIRST = 183
                  DAY_MONTH_X_MAX = 31
               CASE(214:244)
                  MONTH = 8
                  X = 1
                  DAY_MONTH_X = IDAY-213
		          DAY_MONTH_X_FIRST = 214
                  DAY_MONTH_X_MAX = 31
               CASE(245:274)
                  MONTH = 9
                  X = 1
                  DAY_MONTH_X = IDAY-244
                  DAY_MONTH_X_FIRST = 245
                  DAY_MONTH_X_MAX = 30
               CASE(275:305)
                  MONTH = 10
                  X = 1
                  DAY_MONTH_X = IDAY-274
                  DAY_MONTH_X_FIRST = 275
                  DAY_MONTH_X_MAX = 31
               CASE(306:335)
                  MONTH = 11
                  X = 1
                  DAY_MONTH_X = IDAY-305
                  DAY_MONTH_X_FIRST = 306
                  DAY_MONTH_X_MAX = 30
               CASE(336:366)
                  MONTH = 12
                  X = 1
                  DAY_MONTH_X = IDAY-335
                  DAY_MONTH_X_FIRST = 336
                  DAY_MONTH_X_MAX = 31
               CASE DEFAULT
                  PRINT*, 'error day'
                  STOP
          END SELECT
         
       END SELECT

      ELSEIF (x_per_m == 4) THEN

      SELECT CASE (bis)
         CASE (0)

            SELECT CASE (IDAY)

               CASE(1:8)
                  MONTH = 1
                  X = 1
                  DAY_MONTH_X = IDAY
                  DAY_MONTH_X_FIRST = 1
                  DAY_MONTH_X_MAX = 8
			   CASE(9:16)
                  MONTH = 1
                  X = 2
                  DAY_MONTH_X = IDAY-8
                  DAY_MONTH_X_FIRST = 9
                  DAY_MONTH_X_MAX = 8
               CASE(17:24)
                  MONTH = 1
                  X = 3
                  DAY_MONTH_X = IDAY-16
                  DAY_MONTH_X_FIRST = 17
                  DAY_MONTH_X_MAX = 8
               CASE(25:31)
                  MONTH = 1
                  X = 4
                  DAY_MONTH_X = IDAY-24
                  DAY_MONTH_X_FIRST = 25
                  DAY_MONTH_X_MAX = 7

               CASE(32:38)
                  MONTH = 2
                  X = 1
                  DAY_MONTH_X = IDAY-31
                  DAY_MONTH_X_FIRST = 32
                  DAY_MONTH_X_MAX = 7
               CASE(39:45)
                  MONTH = 2
                  X = 2
                  DAY_MONTH_X = IDAY-38
                  DAY_MONTH_X_FIRST = 39
                  DAY_MONTH_X_MAX = 7
               CASE(46:52)
                  MONTH = 2
                  X = 3
                  DAY_MONTH_X = IDAY-45
                  DAY_MONTH_X_FIRST = 46
                  DAY_MONTH_X_MAX = 7
               CASE(53:59)
                  MONTH = 2
                  X = 4
                  DAY_MONTH_X = IDAY-52
                  DAY_MONTH_X_FIRST = 53
                  DAY_MONTH_X_MAX = 7

               CASE(60:67)
                  MONTH = 3
                  X = 1
                  DAY_MONTH_X = IDAY-59
                  DAY_MONTH_X_FIRST = 60
                  DAY_MONTH_X_MAX = 8
               CASE(68:75)
                  MONTH = 3
                  X = 2
                  DAY_MONTH_X = IDAY-67
                  DAY_MONTH_X_FIRST = 68
                  DAY_MONTH_X_MAX = 8
               CASE(76:83)
                  MONTH = 3
                  X = 3
                  DAY_MONTH_X = IDAY-75
                  DAY_MONTH_X_FIRST = 76
                  DAY_MONTH_X_MAX = 8
               CASE(84:90)
                  MONTH = 3
                  X = 4
                  DAY_MONTH_X = IDAY-83
                  DAY_MONTH_X_FIRST = 84
                  DAY_MONTH_X_MAX = 7

               CASE(91:98)
                  MONTH = 4
                  X = 1
                  DAY_MONTH_X = IDAY-90
                  DAY_MONTH_X_FIRST = 91
                  DAY_MONTH_X_MAX = 8
               CASE(99:106)
                  MONTH = 4
                  X = 2
                  DAY_MONTH_X = IDAY-98
                  DAY_MONTH_X_FIRST = 99
                  DAY_MONTH_X_MAX = 8
               CASE(107:113)
                  MONTH = 4
                  X = 3
                  DAY_MONTH_X = IDAY-106
                  DAY_MONTH_X_FIRST = 107
                  DAY_MONTH_X_MAX = 7
               CASE(114:120)
                  MONTH = 4
                  X = 4
                  DAY_MONTH_X = IDAY-113
                  DAY_MONTH_X_FIRST = 114
                  DAY_MONTH_X_MAX = 6

               CASE(121:128)
                  MONTH = 5
                  X = 1
                  DAY_MONTH_X = IDAY-120
                  DAY_MONTH_X_FIRST = 121
                  DAY_MONTH_X_MAX = 8
               CASE(129:136)
                  MONTH = 5
                  X = 2
                  DAY_MONTH_X = IDAY-128
                  DAY_MONTH_X_FIRST = 129
                  DAY_MONTH_X_MAX = 8
               CASE(137:144)
                  MONTH = 5
                  X = 3
                  DAY_MONTH_X = IDAY-136
                  DAY_MONTH_X_FIRST = 137
                  DAY_MONTH_X_MAX = 8
               CASE(145:151)
                  MONTH = 5
                  X = 4
                  DAY_MONTH_X = IDAY-144
                  DAY_MONTH_X_FIRST = 145
                  DAY_MONTH_X_MAX = 7

               CASE(152:159)   
                  MONTH = 6
                  X = 1
                  DAY_MONTH_X = IDAY-151
                  DAY_MONTH_X_FIRST = 152
                  DAY_MONTH_X_MAX = 8
               CASE(160:167)   
                  MONTH = 6
                  X = 2
                  DAY_MONTH_X = IDAY-159
                  DAY_MONTH_X_FIRST = 160
                  DAY_MONTH_X_MAX = 8
               CASE(168:174)   
                  MONTH = 6
                  X = 3
                  DAY_MONTH_X = IDAY-167
                  DAY_MONTH_X_FIRST = 168
                  DAY_MONTH_X_MAX = 7
               CASE(175:181)   
                  MONTH = 6
                  X = 4
                  DAY_MONTH_X = IDAY-174
                  DAY_MONTH_X_FIRST = 175
                  DAY_MONTH_X_MAX = 7

               CASE(182:189)
                  MONTH = 7
                  X = 1
                  DAY_MONTH_X = IDAY-181
                  DAY_MONTH_X_FIRST = 182
                  DAY_MONTH_X_MAX = 8
               CASE(190:197)
                  MONTH = 7
                  X = 2
                  DAY_MONTH_X = IDAY-189
                  DAY_MONTH_X_FIRST = 191
                  DAY_MONTH_X_MAX = 8
               CASE(198:205)
                  MONTH = 7
                  X = 3
                  DAY_MONTH_X = IDAY-197
                  DAY_MONTH_X_FIRST = 198
                  DAY_MONTH_X_MAX = 8
               CASE(206:212)
                  MONTH = 7
                  X = 4
                  DAY_MONTH_X = IDAY-205
                  DAY_MONTH_X_FIRST = 206
                  DAY_MONTH_X_MAX = 7

               CASE(213:220)
                  MONTH = 8
				  X = 1
                  DAY_MONTH_X = IDAY-212
                  DAY_MONTH_X_FIRST = 213
                  DAY_MONTH_X_MAX = 8
               CASE(221:228)
                  MONTH = 8
				  X = 2
                  DAY_MONTH_X = IDAY-220
                  DAY_MONTH_X_FIRST = 221
                  DAY_MONTH_X_MAX = 8
               CASE(229:236)
                  MONTH = 8
				  X = 3
                  DAY_MONTH_X = IDAY-228
                  DAY_MONTH_X_FIRST = 229
                  DAY_MONTH_X_MAX = 8
               CASE(237:243)
                  MONTH = 8
				  X = 4
                  DAY_MONTH_X = IDAY-236
                  DAY_MONTH_X_FIRST = 237
                  DAY_MONTH_X_MAX = 7

               CASE(244:251)
                  MONTH = 9
				  X = 1
                  DAY_MONTH_X = IDAY-243
                  DAY_MONTH_X_FIRST = 244
                  DAY_MONTH_X_MAX = 8
               CASE(252:259)
                  MONTH = 9
				  X = 2
                  DAY_MONTH_X = IDAY-251
                  DAY_MONTH_X_FIRST = 252
                  DAY_MONTH_X_MAX = 8
               CASE(260:266)
                  MONTH = 9
				  X = 3
                  DAY_MONTH_X = IDAY-259
                  DAY_MONTH_X_FIRST = 260
                  DAY_MONTH_X_MAX = 7
               CASE(267:273)
                  MONTH = 9
				  X = 4
                  DAY_MONTH_X = IDAY-266
                  DAY_MONTH_X_FIRST = 267
                  DAY_MONTH_X_MAX = 7

               CASE(274:281)
                  MONTH = 10
				  X = 1
                  DAY_MONTH_X = IDAY-273
                  DAY_MONTH_X_FIRST = 274
                  DAY_MONTH_X_MAX = 8
               CASE(282:289)
                  MONTH = 10
				  X = 2
                  DAY_MONTH_X = IDAY-281
                  DAY_MONTH_X_FIRST = 282
                  DAY_MONTH_X_MAX = 8
               CASE(290:297)
                  MONTH = 10
				  X = 3
                  DAY_MONTH_X = IDAY-289
                  DAY_MONTH_X_FIRST = 290
                  DAY_MONTH_X_MAX = 8
               CASE(298:304)
                  MONTH = 10
				  X = 4
                  DAY_MONTH_X = IDAY-297
                  DAY_MONTH_X_FIRST = 298
                  DAY_MONTH_X_MAX = 7

               CASE(305:312)
                  MONTH = 11
				  X = 1
                  DAY_MONTH_X = IDAY-304
                  DAY_MONTH_X_FIRST = 305
                  DAY_MONTH_X_MAX = 8
               CASE(313:320)
                  MONTH = 11
				  X = 2
                  DAY_MONTH_X = IDAY-312
                  DAY_MONTH_X_FIRST = 313
                  DAY_MONTH_X_MAX = 8
               CASE(321:327)
                  MONTH = 11
				  X = 3
                  DAY_MONTH_X = IDAY-320
                  DAY_MONTH_X_FIRST = 321
                  DAY_MONTH_X_MAX = 7
               CASE(328:334)
                  MONTH = 11
				  X = 4
                  DAY_MONTH_X = IDAY-327
                  DAY_MONTH_X_FIRST = 328
                  DAY_MONTH_X_MAX = 7

               CASE(335:342)
                  MONTH = 12
				  X = 1
                  DAY_MONTH_X = IDAY-334
                  DAY_MONTH_X_FIRST = 335
                  DAY_MONTH_X_MAX = 8
               CASE(343:350)
                  MONTH = 12
				  X = 2
                  DAY_MONTH_X = IDAY-342
                  DAY_MONTH_X_FIRST = 343
                  DAY_MONTH_X_MAX = 8
               CASE(351:358)
                  MONTH = 12
				  X = 3
                  DAY_MONTH_X = IDAY-350
                  DAY_MONTH_X_FIRST = 351
                  DAY_MONTH_X_MAX = 8
               CASE(359:365)
                  MONTH = 12
				  X = 4
                  DAY_MONTH_X = IDAY-358
                  DAY_MONTH_X_FIRST = 359
                  DAY_MONTH_X_MAX = 7

               CASE DEFAULT
                  PRINT*, 'error day'
                  STOP
            END SELECT
                  
         CASE (1)

            SELECT CASE (IDAY)

               CASE(1:8)
                  MONTH = 1
				  X = 1
                  DAY_MONTH_X = IDAY
                  DAY_MONTH_X_FIRST = 1
                  DAY_MONTH_X_MAX = 8
			   CASE(9:16)
                  MONTH = 1
				  X = 2
                  DAY_MONTH_X = IDAY-8
                  DAY_MONTH_X_FIRST = 9
                  DAY_MONTH_X_MAX = 8
               CASE(17:24)
                  MONTH = 1
				  X = 3
                  DAY_MONTH_X = IDAY-16
                  DAY_MONTH_X_FIRST = 17
                  DAY_MONTH_X_MAX = 8
               CASE(25:31)
                  MONTH = 1
				  X = 4
                  DAY_MONTH_X = IDAY-24
                  DAY_MONTH_X_FIRST = 25
                  DAY_MONTH_X_MAX = 7

               CASE(32:39)
                  MONTH = 2
				  X = 1
                  DAY_MONTH_X = IDAY-31
                  DAY_MONTH_X_FIRST = 32
                  DAY_MONTH_X_MAX = 8
               CASE(40:46)
                  MONTH = 2
				  X = 2
                  DAY_MONTH_X = IDAY-39
                  DAY_MONTH_X_FIRST = 40
                  DAY_MONTH_X_MAX = 7
               CASE(47:53)
                  MONTH = 2
				  X = 3
                  DAY_MONTH_X = IDAY-46
                  DAY_MONTH_X_FIRST = 47
                  DAY_MONTH_X_MAX = 7
               CASE(54:60)
                  MONTH = 2
				  X = 4
                  DAY_MONTH_X = IDAY-53
                  DAY_MONTH_X_FIRST = 54
                  DAY_MONTH_X_MAX = 7

               CASE(61:68)
                  MONTH = 3
				  X = 1
                  DAY_MONTH_X = IDAY-60
                  DAY_MONTH_X_FIRST = 60
                  DAY_MONTH_X_MAX = 8
               CASE(69:76)
                  MONTH = 3
				  X = 2
                  DAY_MONTH_X = IDAY-68
                  DAY_MONTH_X_FIRST = 69
                  DAY_MONTH_X_MAX = 8
               CASE(77:84)
                  MONTH = 3
				  X = 3
                  DAY_MONTH_X = IDAY-76
                  DAY_MONTH_X_FIRST = 77
                  DAY_MONTH_X_MAX = 8
               CASE(85:91)
                  MONTH = 3
				  X = 4
                  DAY_MONTH_X = IDAY-84
                  DAY_MONTH_X_FIRST = 85
                  DAY_MONTH_X_MAX = 7

               CASE(92:99)
                  MONTH = 4
				  X = 1
                  DAY_MONTH_X = IDAY-91
	              DAY_MONTH_X_FIRST = 92
                  DAY_MONTH_X_MAX = 8
               CASE(100:107)
                  MONTH = 4
				  X = 2
                  DAY_MONTH_X = IDAY-99
	              DAY_MONTH_X_FIRST = 100
                  DAY_MONTH_X_MAX = 8
               CASE(108:114)
                  MONTH = 4
				  X = 3
                  DAY_MONTH_X = IDAY-107
	              DAY_MONTH_X_FIRST = 108
                  DAY_MONTH_X_MAX = 7
               CASE(115:121)
                  MONTH = 4
				  X = 4
                  DAY_MONTH_X = IDAY-114
	              DAY_MONTH_X_FIRST = 115
                  DAY_MONTH_X_MAX = 7

               CASE(122:129)
                  MONTH = 5
				  X = 1
                  DAY_MONTH_X = IDAY-121
                  DAY_MONTH_X_FIRST = 122
                  DAY_MONTH_X_MAX = 8
               CASE(130:137)
                  MONTH = 5
				  X = 2
                  DAY_MONTH_X = IDAY-129
                  DAY_MONTH_X_FIRST = 130
                  DAY_MONTH_X_MAX = 8
               CASE(138:145)
                  MONTH = 5
				  X = 3
                  DAY_MONTH_X = IDAY-137
                  DAY_MONTH_X_FIRST = 138
                  DAY_MONTH_X_MAX = 8
               CASE(146:152)
                  MONTH = 5
				  X = 4
                  DAY_MONTH_X = IDAY-145
                  DAY_MONTH_X_FIRST = 146
                  DAY_MONTH_X_MAX = 7

               CASE(153:160)   
                  MONTH = 6
				  X = 1
                  DAY_MONTH_X = IDAY-152
		          DAY_MONTH_X_FIRST = 153
                  DAY_MONTH_X_MAX = 8
               CASE(161:168)   
                  MONTH = 6
				  X = 2
                  DAY_MONTH_X = IDAY-160
		          DAY_MONTH_X_FIRST = 161
                  DAY_MONTH_X_MAX = 8
               CASE(169:175)   
                  MONTH = 6
				  X = 3
                  DAY_MONTH_X = IDAY-168
		          DAY_MONTH_X_FIRST = 169
                  DAY_MONTH_X_MAX = 7
               CASE(176:182)   
                  MONTH = 6
				  X = 4
                  DAY_MONTH_X = IDAY-175
		          DAY_MONTH_X_FIRST = 176
                  DAY_MONTH_X_MAX = 7

               CASE(183:190)
                  MONTH = 7
				  X = 1
                  DAY_MONTH_X = IDAY-182
		          DAY_MONTH_X_FIRST = 183
                  DAY_MONTH_X_MAX = 8
               CASE(191:198)
                  MONTH = 7
				  X = 2
                  DAY_MONTH_X = IDAY-190
		          DAY_MONTH_X_FIRST = 191
                  DAY_MONTH_X_MAX = 8
               CASE(199:206)
                  MONTH = 7
				  X = 3
                  DAY_MONTH_X = IDAY-198
		          DAY_MONTH_X_FIRST = 199
                  DAY_MONTH_X_MAX = 8
               CASE(207:213)
                  MONTH = 7
				  X = 4
                  DAY_MONTH_X = IDAY-206
		          DAY_MONTH_X_FIRST = 207
                  DAY_MONTH_X_MAX = 7

               CASE(214:221)
                  MONTH = 8
				  X = 1
                  DAY_MONTH_X = IDAY-213
		          DAY_MONTH_X_FIRST = 214
                  DAY_MONTH_X_MAX = 8
               CASE(222:229)
                  MONTH = 8
				  X = 2
                  DAY_MONTH_X = IDAY-221
		          DAY_MONTH_X_FIRST = 222
                  DAY_MONTH_X_MAX = 8
               CASE(230:237)
                  MONTH = 8
				  X = 3
                  DAY_MONTH_X = IDAY-229
		          DAY_MONTH_X_FIRST = 230
                  DAY_MONTH_X_MAX = 8
               CASE(238:244)
                  MONTH = 8
				  X = 4
                  DAY_MONTH_X = IDAY-237
		          DAY_MONTH_X_FIRST = 238
                  DAY_MONTH_X_MAX = 7

               CASE(245:252)
                  MONTH = 9
				  X = 1
                  DAY_MONTH_X = IDAY-244
                  DAY_MONTH_X_FIRST = 245
                  DAY_MONTH_X_MAX = 8
               CASE(253:260)
                  MONTH = 9
				  X = 2
                  DAY_MONTH_X = IDAY-252
                  DAY_MONTH_X_FIRST = 253
                  DAY_MONTH_X_MAX = 8
               CASE(261:267)
                  MONTH = 9
				  X = 3
                  DAY_MONTH_X = IDAY-260
                  DAY_MONTH_X_FIRST = 261
                  DAY_MONTH_X_MAX = 7
               CASE(268:274)
                  MONTH = 9
				  X = 4
                  DAY_MONTH_X = IDAY-267
                  DAY_MONTH_X_FIRST = 268
                  DAY_MONTH_X_MAX = 7

               CASE(275:282)
                  MONTH = 10
				  X = 1
                  DAY_MONTH_X = IDAY-274
                  DAY_MONTH_X_FIRST = 275
                  DAY_MONTH_X_MAX = 8
               CASE(283:290)
                  MONTH = 10
				  X = 2
                  DAY_MONTH_X = IDAY-282
                  DAY_MONTH_X_FIRST = 283
                  DAY_MONTH_X_MAX = 8
               CASE(291:298)
                  MONTH = 10
				  X = 3
                  DAY_MONTH_X = IDAY-290
                  DAY_MONTH_X_FIRST = 291
                  DAY_MONTH_X_MAX = 8
               CASE(299:305)
                  MONTH = 10
				  X = 4
                  DAY_MONTH_X = IDAY-298
                  DAY_MONTH_X_FIRST = 299
                  DAY_MONTH_X_MAX = 7

               CASE(306:313)
                  MONTH = 11
				  X = 1
                  DAY_MONTH_X = IDAY-305
                  DAY_MONTH_X_FIRST = 306
                  DAY_MONTH_X_MAX = 8
               CASE(314:321)
                  MONTH = 11
				  X = 2
                  DAY_MONTH_X = IDAY-313
                  DAY_MONTH_X_FIRST = 314
                  DAY_MONTH_X_MAX = 8
               CASE(322:328)
                  MONTH = 11
				  X = 3
                  DAY_MONTH_X = IDAY-321
                  DAY_MONTH_X_FIRST = 322
                  DAY_MONTH_X_MAX = 7
               CASE(329:335)
                  MONTH = 11
				  X = 4
                  DAY_MONTH_X = IDAY-328
                  DAY_MONTH_X_FIRST = 329
                  DAY_MONTH_X_MAX = 7

               CASE(336:343)
                  MONTH = 12
				  X = 1
                  DAY_MONTH_X = IDAY-335
                  DAY_MONTH_X_FIRST = 336
                  DAY_MONTH_X_MAX = 8
               CASE(344:351)
                  MONTH = 12
				  X = 2
                  DAY_MONTH_X = IDAY-343
                  DAY_MONTH_X_FIRST = 344
                  DAY_MONTH_X_MAX = 8
               CASE(352:359)
                  MONTH = 12
				  X = 3
                  DAY_MONTH_X = IDAY-351
                  DAY_MONTH_X_FIRST = 352
                  DAY_MONTH_X_MAX = 8
               CASE(360:366)
                  MONTH = 12
				  X = 4
                  DAY_MONTH_X = IDAY-359
                  DAY_MONTH_X_FIRST = 360
                  DAY_MONTH_X_MAX = 7

               CASE DEFAULT
                  PRINT*, 'error day'
                  STOP
          END SELECT
         
       END SELECT
      ELSE
        PRINT*, 'error timing printing : not 1 or 4 per month'
      ENDIF
       
c --- Number of step in this month ---
      
      DELT_m = DAY_MONTH_X_MAX*24.*60.*60. ! * number of second in one month/x
      N_DELT_MONTH_X =   DELT_m/DELT   ! the number of step in one month/x


c --- Write commun basic values ---

      IF (NSTEP==1) THEN
        NSTEP_write = 1
        WRITE(70,7000) IL2-IL1+1,DELT,DELZLK,NLAKMAX_bgc
        !stor_p = "no"
        HDPTH_m = NLAK_bgc*DELZLK
      ENDIF

C
C ------------------ Monthly/X value ------------------------------------------
C

C --- INITIALISATION at the beginning of the month's part ---
	
      IF (
     1     (DAY_MONTH_X==1) .AND. (IHOUR==0) .AND. 
     2       (IMIN==0) )  THEN
        !HDPTH_p = HDPTH_bgc
        LKICEH_m = 0.
        TLAK_bgc_m = 0.
        HDPTH_bgc_m = 0.
        TSED_m = 0.
        Clabile_m = 0.
        FPP_m = 0.
       	FRESP_m = 0.
        O2_m = 0. 
        CH4_m = 0.
        FatmO2_m = 0.
        Fw_m = 0.
        FSAER_m = 0.
        FPPWC_m = 0.
        FatmCH4_m = 0.
        FsCH4_C_m = 0.
        FsDiff_m = 0.
        FsEbul_m = 0.
        FatmEbul_m = 0.
        Clabs_1_m = 0.
        Clabs_2_m = 0.
        Clabs_3_m = 0.
        FsCH4_C_1_m = 0.
        FsCH4_C_2_m = 0.
        FsCH4_C_3_m = 0.
        Fsed_1_m = 0.
        Fsed_2_m = 0.
        Fsed_3_m = 0.
        FswDiff_m = 0.
        OxCH4_m = 0.
        K_diff_m = 0.
        ZPhotic_m = 0. 
        ZRESP_m = 0. 
        FRESPWC_m = 0.
        Fstor_m = 0.
        DELT_stor_m = 0.
        DELT_ice_m = 0.
      ENDIF


C --- COMPUTE VALUE FOR THE MONTH's PART ---

        LKICEH_m = LKICEH_m + LKICEH/N_DELT_MONTH_X
        TLAK_bgc_m = TLAK_bgc_m + TLAK_bgc/N_DELT_MONTH_X
        HDPTH_bgc_m = HDPTH_bgc_m + HDPTH_bgc/N_DELT_MONTH_X
        TSED_m = TSED_m + TSED/N_DELT_MONTH_X
        Clabile_m = Clabile_m + Clabile/N_DELT_MONTH_X

        FPP_m = FPP_m + FPP
        FRESP_m = FRESP_m + FRESP

        O2_m = O2_m + O2/N_DELT_MONTH_X
        CH4_m = CH4_m + CH4/N_DELT_MONTH_X

        FatmO2_m = FatmO2_m + FatmO2
        Fw_m = Fw_m + Fw
        FPPWC_m = FPPWC_m + FPPWC
        FSAER_m = FSAER_m + FSAER
        FatmCH4_m = FatmCH4_m + FatmCH4
        FsCH4_C_m = FsCH4_C_m + FsCH4
        FsDiff_m = FsDiff_m + FsDiff
        FsEbul_m = FsEbul_m + FsEbul
        FatmEbul_m = FatmEbul_m + FatmEbul

        Clabs_1_m = Clabs_1_m + Clabs_1/N_DELT_MONTH_X
        Clabs_2_m = Clabs_2_m + Clabs_2/N_DELT_MONTH_X
        Clabs_3_m = Clabs_3_m + Clabs_3/N_DELT_MONTH_X

        FsCH4_C_1_m = FsCH4_C_1_m + FsCH4_C_1
        FsCH4_C_2_m = FsCH4_C_2_m + FsCH4_C_2
        FsCH4_C_3_m = FsCH4_C_3_m + FsCH4_C_3
        Fsed_1_m = Fsed_1_m + Fsed_1
        Fsed_2_m = Fsed_2_m + Fsed_2
        Fsed_3_m = Fsed_3_m + Fsed_3
        FswDiff_m = FswDiff_m + FswDiff
        OxCH4_m = OxCH4_m + OxCH4	

        K_diff_m = K_diff_m + K_diff/N_DELT_MONTH_X

        ZPhotic_m = ZPhotic_m + ZPhotic/N_DELT_MONTH_X
        ZRESP_m = ZRESP_m + ZRESP/N_DELT_MONTH_X
        FRESPWC_m = FRESPWC_m + FRESPWC

		
        DO I=IL1,IL2
            IF ( ( (TLAK_bgc(I,1) .gt. 3.5+Tmp ).AND. 
     1                ( TLAK_bgc(I,1) .lt. 4.5+Tmp) ) )  THEN
              IF ( FatmCH4(I) .gt. 0.0 ) THEN
              Fstor_m(I) = Fstor_m(I) + FatmCH4(I)
              DELT_stor_m(I) = DELT_stor_m(I) + DELT	
		      ENDIF
            ENDIF
            IF ( LKICEH(I) > 0. ) THEN
              DELT_ice_m(I) = DELT_ice_m(I) + DELT
            ENDIF
        ENDDO


C --- WRITE AT THE END OF THE MONTH ---

      IF ( 
     1     (DAY_MONTH_X==DAY_MONTH_X_MAX) .AND. (IHOUR==23) .AND. 
     2     (IMIN==60-INT(DELT/60)) ) THEN


       DO 666 I = IL1,IL2
            
        IF (NSTEP_write==1) THEN
	     IF (I==IL1) THEN
	        WRITE(71,7100) NLAK_bgc(I),(Area(I,J),J=0,NLAK_bgc(I)+1),  ! !!!!!!!!!! WE SHOULD WRITE THE DIFFERENCE BETWEEN TO LAYE
     1               (6666.6666,J=NLAK_bgc(I)+1,NLAKMAX_bgc)
C                WRITE(90,9000) (Area(I,J)-Area(I,J+1),J=1,NLAK_bgc(I)),  ! !!!!!!!!!! WE SHOULD WRITE THE DIFFERENCE BETWEEN TO LAYER TO ENSURE THAT WE HAVE THATNO MATHER THE PRESICION AND THE SIZE OF THE LAKE !!!
C     1               (6666.6666,J=NLAK_bgc(I)+1,NLAKMAX_bgc)
          ELSE
            WRITE(71,7100) NLAK_bgc(I), (Area(I,J),J=0,NLAK_bgc(I)+1)
C            WRITE(90,9000) (Area(I,J)-Area(I,J+1),J=1,NLAK_bgc(I))
         ENDIF
        ENDIF

         WRITE(72,7200) IYEAR,IDAY,IHOUR,IMIN,I
         WRITE(73,7300) HDPTH_bgc_m(I),LKICEH_m(I),TSED_m(I)-TFREZ
         WRITE(74,7400) Clabile_m(I)/1000.,FPP_m(I)/1000.,
     1                   FRESP_m(I)/1000.

       IF ( (NSTEP_write==1) .AND. (I==IL1) ) THEN
       WRITE(75,7500) (Clabs_1_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     1                   (Clabs_2_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     2                   (Clabs_3_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     3                   (Fsed_1_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     4                   (Fsed_2_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     5                   (Fsed_3_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     6                   (FsCH4_C_1_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     7                   (FsCH4_C_2_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     8                   (FsCH4_C_3_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     9			 (6666.6666,J=9*NLAK_bgc(I)+1,9*NLAKMAX_bgc)
       WRITE(76,7600) (O2_m(I,J),J=1,NLAK_bgc(I)),FatmO2_m(I),
     1                   (FPPWC_m(I,J),J=1,NLAK_bgc(I)),Fw_m(I),
     2                   (Fsaer_m(I,J),J=1,NLAK_bgc(I)),
     9			 (6666.6666,J=3*NLAK_bgc(I)+3,3*NLAKMAX_bgc+2 )
       WRITE(77,7700) (FsCH4_C_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     1                   (FsDiff_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     2                   (FsEbul_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     3                   (6666.6666,J=3*NLAK_bgc(I)+1,3*NLAKMAX_bgc)
       WRITE(78,7800) (CH4_m(I,J),J=1,NLAK_bgc(I)),FatmCH4_m(I),
     1                   (FswDiff_m(I,J),J=1,NLAK_bgc(I)),
     2                   (OxCH4_m(I,J),J=1,NLAK_bgc(I)),
     2                   (6666.6666,J=3*NLAK_bgc(I)+2,3*NLAKMAX_bgc+1)
       WRITE(79,7900) (TLAK_bgc_m(I,J)-TFREZ,J=1,NLAK_bgc(I)),
     1         (6666.6666,J=NLAK_bgc(I)+1,NLAKMAX_bgc)
       WRITE(80,8000) (K_diff_m(I,J),J=1,NLAK_bgc(I)+1),
     2          (6666.6666,J=NLAK_bgc(I)+2,NLAKMAX_bgc+1)
       WRITE(81,8010) ZPhotic_m(I),ZRESP_m(I),
     1                   (FRESPWC_m(I,J), J=1,NLAK_bgc(I)),
     2          (6666.6666,J=NLAK_bgc(I)+1,NLAKMAX_bgc)
       WRITE(82,8020) (FatmEbul_m(I,J), J=1,NLAK_bgc(I)),
     2          (6666.6666,J=NLAK_bgc(I)+1,NLAKMAX_bgc)
       WRITE(83,8030) DELT_stor_m(I),Fstor_m(I),DELT_m, DELT_ice_m(I)

       ELSE
       WRITE(75,7500) (Clabs_1_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     1                   (Clabs_2_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     2                   (Clabs_3_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     3                   (Fsed_1_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     4                   (Fsed_2_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     5                   (Fsed_3_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     6                   (FsCH4_C_1_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     7                   (FsCH4_C_2_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     8                   (FsCH4_C_3_m(I,J)/1000.,J=1,NLAK_bgc(I))
       WRITE(76,7600) (O2_m(I,J),J=1,NLAK_bgc(I)),FatmO2_m(I),
     1                   (FPPWC_m(I,J),J=1,NLAK_bgc(I)),Fw_m(I),
     2                   (Fsaer_m(I,J),J=1,NLAK_bgc(I))
       WRITE(77,7700) (FsCH4_C_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     1                   (FsDiff_m(I,J)/1000.,J=1,NLAK_bgc(I)),
     2                   (FsEbul_m(I,J)/1000.,J=1,NLAK_bgc(I))
       WRITE(78,7800) (CH4_m(I,J),J=1,NLAK_bgc(I)),FatmCH4_m(I),
     1                   (FswDiff_m(I,J),J=1,NLAK_bgc(I)),
     2                   (OxCH4_m(I,J),J=1,NLAK_bgc(I))
       WRITE(79,7900) (TLAK_bgc_m(I,J)-TFREZ,J=1,NLAK_bgc(I))
       WRITE(80,8000) (K_diff_m(I,J),J=1,NLAK_bgc(I)+1)
       WRITE(81,8010) ZPhotic_m(I),ZRESP_m(I),
     1                   (FRESPWC_m(I,J), J=1,NLAK_bgc(I))
        WRITE(82,8020) (FatmEbul_m(I,J), J=1,NLAK_bgc(I))
        WRITE(83,8030) DELT_stor_m(I),Fstor_m(I),DELT_m, DELT_ice_m(I)

      ENDIF

 666  CONTINUE

        NSTEP_write = NSTEP_write + 1
      
      ENDIF

C --- DECLARE FORMATS ---

7000    FORMAT( I5,1X,F5.1,1X,F4.2,1X,I3 )
7100    FORMAT( I5,1X,823(G20.3E3,1X) )    ! 820 + 3
7200    FORMAT( I4,1X,I3,1X,I2,1X,I2,1X,I5 )
7300    FORMAT( F10.2,1X,G10.3E2,1X, G12.4E3 )
7400    FORMAT( 3(G12.4E3,1X) )
7500    FORMAT( 7380(G12.4E3,1X) )      ! 9*NLAKMAX = 9*205     9*NLAKMAX_bgc = 9*322  9*360
7600    FORMAT( 2462(G12.4E3,1X) )       ! 3*NLAKMAX+2 = 617     3*NLAKMAX_bgc+2 = 968
7700    FORMAT( 2464(G12.4E3,1X) )
7800    FORMAT( 2461(G12.4E3,1X) )       ! 3*NLAKMAX+1 = 411     3*NLAKMAX_bgc+1 = 967
7900    FORMAT( 820(G12.4E3,1X) )     ! 820
8000    FORMAT( 822(G12.4E3,1X) )     ! 820 + 2    
8010    FORMAT( 823(G12.4E3,1X) )     ! 820 + 3
8020    FORMAT( 820(G12.4E3,1X) )     ! 820
8030    FORMAT( F15.1,1X,G14.4E3, F15.1 ,1X, F15.1  )
C8040    FORMAT( F12.1 ,1X, F12.1 )

      
C 7000    FORMAT( I5,1X,F5.1,1X,F4.2,1X,I3 )
C 7100    FORMAT( I5,1X,363(G26.20E3,1X) )
C 7200    FORMAT( I4,1X,I3,1X,I2,1X,I2,1X,I5 )
C 7300    FORMAT( F7.2,1X,G10.3E2,1X, G12.4E3 )
C 7400    FORMAT( 3(G12.4E3,1X) )
C 7500    FORMAT( 3240(G12.4E3,1X) )      ! 9*NLAKMAX = 9*205     9*NLAKMAX_bgc = 9*322  9*360
C 7600    FORMAT( 1082(G12.4E3,1X) )       ! 3*NLAKMAX+2 = 617     3*NLAKMAX_bgc+2 = 968
C 7700    FORMAT( 1084(G12.4E3,1X) )
C 7800    FORMAT( 1081(G12.4E3,1X) )       ! 3*NLAKMAX+1 = 411     3*NLAKMAX_bgc+1 = 967
C 7900    FORMAT( 360(G12.4E3,1X) )
C 8000    FORMAT( 362(G12.4E3,1X) )
C 8010    FORMAT( 363(G12.4E3,1X) )
C 8020    FORMAT( 360(G12.4E3,1X) )
C 8030    FORMAT( F12.1,1X,G14.4E3 )
C 8040    FORMAT( F12.1 ,1X, F12.1 )

C 9000    FORMAT( 360(G26.20E3,1X) )


        RETURN
        END



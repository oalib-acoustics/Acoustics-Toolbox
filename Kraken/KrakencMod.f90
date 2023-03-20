MODULE KrakencMod

  USE MathConstants
  SAVE

  INTEGER,          PARAMETER   :: ENVFile = 5, PRTFile = 6, MODFile = 20, EVMFile = 22, MaxM = 20000, MaxMedium = 500, NSets = 5
  
  INTEGER                       :: FirstAcoustic, LastAcoustic, NV( NSets ), ISet, M, &
                                   LRecordLength, IRecProfile = 1, ModeCount, Mode, IProf, ifreq
  REAL      (KIND=8)            :: ET( NSets ), hV( NSets ), VG( MaxM ), cMin, cLow, cHigh, freq, omega, omega2, RMax
  COMPLEX   (KIND=8)            :: k( MaxM ), EVMat( NSets, MaxM ), Extrap( NSets, MaxM )
  CHARACTER (LEN= 8)            :: TopOpt, BotOpt
  CHARACTER (LEN=80)            :: Title

  ! finite-difference grid
  INTEGER                       :: Loc( MaxMedium ), NG( MaxMedium ), N( MaxMedium )
  REAL      (KIND=8)            :: h(   MaxMedium )

  ! storage for finite-difference equations
  REAL    (KIND=8), ALLOCATABLE :: rho( : )
  COMPLEX (KIND=8), ALLOCATABLE :: B1( : ), B2( : ), B3( : ), B4( : )

END MODULE KrakencMod

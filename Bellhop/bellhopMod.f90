MODULE bellhopMod

  USE MathConstants
  INTEGER, PARAMETER :: ENVFile = 5, PRTFile = 6, RAYFile = 21, SHDFile = 25, ARRFile = 36, SSPFile = 40, MaxN = 100000

  ! Reduce MaxN (= max # of steps along a ray) to reduce storage
  ! Note space is wasted in NumTopBnc, NumBotBnc ...

  LOGICAL            :: ThreeD   ! flag to indicate BELLHOP vs BELLHOP3D run
  INTEGER            :: Nrz_per_range, istep
  REAL    ( KIND= 8) :: freq, omega, SrcDeclAngle, SrcAzimAngle, xs_3D( 3 )
  CHARACTER (LEN=80) :: Title

  ! *** Beam structure ***

  TYPE rxyz
     REAL (KIND=8) :: r, x, y, z
  END TYPE rxyz

  TYPE BeamStructure
     INTEGER           :: NBeams, Nimage = 1, Nsteps, iBeamWindow = 5
     REAL     (KIND=8) :: deltas, epsMultiplier = 1, rLoop
     CHARACTER (LEN=1) :: Component = 'P'              ! Pressure or displacement
     CHARACTER (LEN=4) :: Type = 'G S '
     CHARACTER (LEN=7) :: RunType
     TYPE( rxyz )      :: Box
  END TYPE BeamStructure

  TYPE( BeamStructure ) :: Beam

  ! *** ray structure ***

  TYPE ray2DPt
     INTEGER          :: NumTopBnc, NumBotBnc
     REAL   (KIND=8 ) :: x( 2 ), t( 2 ), p( 2 ), q( 2 ), c, Amp, Phase
     COMPLEX (KIND=8) :: tau
  END TYPE ray2DPt
  TYPE( ray2DPt )     :: ray2D( MaxN )

  !!! uncomment COMPLEX below in place of REAL if using paraxial beams
  TYPE ray3DPt
     REAL    (KIND=8) :: p_tilde( 2 ), q_tilde( 2 ), p_hat( 2 ), q_hat( 2 ), DetQ
     REAL    (KIND=8) :: x( 3 ), t( 3 ), phi, c, Amp, Phase
     INTEGER          :: NumTopBnc, NumBotBnc
     ! COMPLEX (KIND=8) :: p_tilde( 2 ), q_tilde( 2 ), p_hat( 2 ), q_hat( 2 ), f, g, h, DetP, DetQ
     COMPLEX (KIND=8) :: tau

  END TYPE ray3DPt
  TYPE( ray3DPt )     :: ray3D( MaxN )

END MODULE bellhopMod

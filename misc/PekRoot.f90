MODULE PekRoot

  ! evaluates a particular branch of the square root that exposes leaky modes

  IMPLICIT NONE

  CONTAINS
    FUNCTION PekerisRoot( z )

      ! At one time, this was used to return the 'Pekeris branch cut'     
      ! which is just a particular branch of the square root that         
      ! exposes many 'leaky' or 'virtual' modes.                          

      ! The current version implements a particular branch that was convenient                                                   

      COMPLEX (KIND=8), INTENT(  IN ) :: z
      COMPLEX (KIND=8)                :: PekerisRoot

      IF ( REAL( z ) >= 0.0D0 ) THEN 
         PekerisRoot = SQRT( z ) 
      ELSE 
         PekerisRoot = ( 0.0D0, 1.0D0 ) * SQRT( -z ) 
      ENDIF

    END FUNCTION PekerisRoot
END MODULE PekRoot

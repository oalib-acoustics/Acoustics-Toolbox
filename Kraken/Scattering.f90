MODULE Scattering

  IMPLICIT NONE

  CONTAINS
    FUNCTION KupIng( sigma, eta1Sq, rho1, eta2Sq, rho2, P, U ) 

      ! Evalutates the imaginary perturbation due to interfacial roughness using the Kuperman-Ingenito formulation                           
      ! P is the pressure at the interface                                
      ! U is P'/rho       "   "      "                                    

      REAL    (KIND=8), INTENT( IN ) :: sigma, rho1, rho2
      COMPLEX (KIND=8), INTENT( IN ) :: P, U, eta1Sq, eta2Sq
      COMPLEX (KIND=8)               :: KupIng, Del, eta1, eta2, A11, A12, A21, A22
      COMPLEX (KIND=8), PARAMETER    :: i = ( 0.0D0, 1.0D0 )

      KupIng = 0.0D0 
      IF ( sigma == 0.0D0 ) RETURN 
      eta1 = ScatterRoot( eta1Sq ) 
      eta2 = ScatterRoot( eta2Sq ) 
      Del  = rho1 * eta2 + rho2 * eta1 

      IF ( Del /= 0.0D0 ) THEN
         A11 = 0.5D0 * ( eta1Sq - eta2Sq ) - ( rho2 * eta1Sq - rho1 * eta2Sq ) * ( eta1 + eta2 ) / Del
         A12 =   i * ( rho2 - rho1 ) ** 2 * eta1 * eta2 / Del
         A21 =  -i * ( rho2 * eta1Sq - rho1 * eta2Sq ) ** 2 / ( rho1 * rho2 * Del )
         A22 = 0.5D0 * ( eta1Sq - eta2Sq ) + ( rho2 - rho1 ) * eta1 * eta2 * ( eta1 + eta2 ) / Del

         KupIng = -sigma**2 * ( -A21 * P ** 2 + ( A11 - A22 ) * P * U + A12 * U ** 2 )
      ENDIF

    END FUNCTION KupIng

    !**********************************************************************

    FUNCTION ScatterRoot( z ) 

      ! Root for interfacial scatter                                      

      COMPLEX ( KIND=8 ) :: ScatterRoot, z

      IF ( REAL( z ) >= 0.0D0 ) THEN 
         ScatterRoot = SQRT( z ) 
      ELSE 
         ScatterRoot = -( 0.0D0, 1.0D0 ) * SQRT( -z ) 
      ENDIF

    END FUNCTION ScatterRoot

  END MODULE Scattering

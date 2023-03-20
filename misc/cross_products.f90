MODULE cross_products

! mbp 1/2012

  IMPLICIT NONE
  
  INTERFACE cross_product
     MODULE PROCEDURE cross_product_sngl, cross_product_dble
  END INTERFACE cross_product

CONTAINS
  FUNCTION cross_product_sngl( a, b )
    REAL (KIND=4), DIMENSION( 3 ) :: cross_product_sngl
    REAL (KIND=4), DIMENSION( 3 ), INTENT(  IN ) :: a, b

    cross_product_sngl( 1 ) = a( 2 ) * b( 3 ) - a( 3 ) * b( 2 )
    cross_product_sngl( 2 ) = a( 3 ) * b( 1 ) - a( 1 ) * b( 3 )
    cross_product_sngl( 3 ) = a( 1 ) * b( 2 ) - a( 2 ) * b( 1 )

  END FUNCTION cross_product_sngl

  FUNCTION cross_product_dble( a, b )
    REAL (KIND=8), DIMENSION( 3 ) :: cross_product_dble
    REAL (KIND=8), DIMENSION( 3 ), INTENT(  IN ) :: a, b

    cross_product_dble( 1 ) = a( 2 ) * b( 3 ) - a( 3 ) * b( 2 )
    cross_product_dble( 2 ) = a( 3 ) * b( 1 ) - a( 1 ) * b( 3 )
    cross_product_dble( 3 ) = a( 1 ) * b( 2 ) - a( 2 ) * b( 1 )

  END FUNCTION cross_product_dble
END MODULE cross_products

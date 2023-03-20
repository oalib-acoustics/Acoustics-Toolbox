!     Last change:  MBP  20 Jun 1999    7:24 am
      SUBROUTINE RFT(A, B, N, ISN)       
!    C, 10 / 11 / 80.V1 RTD <PN> :FFT LIBRARY      
!    IF ISN = 1, THIS SUBROUTINE COMPLETES THE FOURIER TRANSFORM
!     OF 2 * N REAL DATA VALUES, WHERE THE ORIGINAL DATA VALUES ARE
!    STORED ALTERNATELY IN ARRAYS A AND B, AND ARE FIRST
!    TRANSFORMED BY A COMPLEX FOURIER TRANSFORM OF DIMENSION N.
!    THE COSINE COEFFICIENTS ARE IN A(1), A(2), ...A(N + 1) AND
!     THE SINE COEFFICIENTS ARE IN B(1), B(2), ...B(N + 1).
!    A TYPICAL CALLING SEQUENCE IS
!         CALL FFT(A, B, N, N, N, 1)
!         CALL RFT(A, B, N, 1)
!    THE RESULTS SHOULD BE MULTIPLIED BY 0.5  / N TO GIVE THE
!    USUAL SCALING OF COEFFICIENTS.
!    IF ISN = -1, THE INVERSE TRANSFORMATION IS DONE, THE FIRST
!    STEP IN EVALUATING A REAL FOURIER SERIES.
!    A TYPICAL CALLING SEQUENCE IS
!         CALL RFT(A, B, N,  - 1)
!         CALL FFT(A, B, N, N, N,  - 1)
!    THE RESULTS SHOULD BE MULTIPLIED BY 0.5 TO GIVE THE USUAL
!    SCALING , AND THE TIME DOMAIN RESULTS ALTERNATE IN ARRAYS
!     A AND B, I.E., A(1), B(1), A(2), B(2), ...A(N), B(N).
!    THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX
!     ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO
!    GIVE THE CORRECT INDEXING INCREMENT AND A(2) USED TO PASS
!    THE INITIAL ADDRESS FOR THE SEQUENCE OF IMAGINARY
!    VALUES, E.G.
!          CALL FFT(A, A(2), N, N, N, 2)
!          CALL RFT(A, A(2), N, 2)
!    IN THIS CASE, THE COSINE AND SINE COEFFICIENTS ALTERNATE IN
!    A,  BY R.C. SINGLETON, STANFORD RESEARCH INSTITUTE, OCT. 1968
!
! NOTE FROM THE USERS:
!
!     BUGS WERE FOUND IN THIS VERSION OF RFT, WHICH HAS BEEN USED IN
!    BASIC I, II, AND III.   THESE ARE:
!
!      1)  TWO VALUES OUTSIDE THE USER'S BUFFER, X(N + 1) AND X(N + 2) FOR
!     A TRANSFORM OF N DATA POINTS, ARE CLOBBERED IN AN ATTEMPT AT
!     SCRATCH STORAGE.  THIS HAS BEEN CORRECTED IN THE LATER VERSION, AN
!     THE OLD PROGRAM, IDENTICAL IN EXECUTION, BUT WHICH DID NOT CONTAIN
!     CORREKT THINKING IN ITS NONUSE OF COMMENTS, HAS BEEN DESTROYED.
!
!     2)  THE INVERSE TRANSFORM STATEMENTS SHOULD BE CHANGED, SO THAT
!     A(K), A(K + 1) (BETWEEN STATEMENTS 5 & 10) IS DEFINED.
!
!     3)   ROUNDED ARITHMETIC SHOULD BE USED.
!
      DIMENSION A(1), B(1)
      INC = IABS(ISN)        
      NK = N * INC + 2 
      NH = NK / 2    
      SD = -2.0 * ATAN(1.0) / FLOAT(N) 
      CD = 2.0 * SIN(SD) **2    
      SD = SIN(SD + SD)        
      IF ( ISN .LT. 0 ) GOTO 30   
      CN = 1.       
!     DO FIRST LOOP  
!     STORING NYQUIST REAL IN IMAGINARY 0 SAMPLE & VICEVERSA
      AA   = 2. * (A(1) + B(1))
      B(1) = 2. * (A(1) - B(1))   
! 
12    A(1) = AA
! INCREMENT SIN & COSINE    
      SN = SD * CN
      CN = CN * (1. - CD)         
      L = 1 + INC

!     *** MAIN LOOP ***       

      DO J = L, NH, INC      
!        OK, FROM HERE ON IT'S ALMOST THE SAME.    
         K = NK - J
         AA = A(J) + A(K)         
         AB = A(J) - A(K)         
         BA = B(J) + B(K)         
         BB = B(J) - B(K)         
         RE = CN * BA + SN * AB      
         FIM = SN * BA - CN * AB      
         B(K) = FIM - BB
         B(J) = FIM + BB         
         A(K) = AA - RE 
         A(J) = AA + RE 
         AA = CN - (CD * CN + SD * SN)  
         SN = (SD * CN - CD * SN) + SN  
!        THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION         
!        ERROR.  IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE
         CN = AA
       END DO     
!      CN = 0.5 / (AA**2 + SN**2) + 0.5  
!      SN = CN * SN
!      CN = CN * AA

       RETURN

!      *** INVERSE TRANSFORM ***         
 30    CN = -1.0
       SD = -SD     
       AA   = A(1) + B(1)       
       B(1) = A(1) - B(1)       
       GO TO 12
   
       END        
!**********************************************************************C
      SUBROUTINE RFFT(A, NX, IFORW)
      DIMENSION A(NX)
      IF (IFORW.GT.0) THEN
         CALL CFFT(A, NX / 2, 1)
         CALL RFT(A(1), A(2), NX / 2, 2)
      ELSE
         CALL RFT(A(1), A(2), NX / 2,  - 2)
         CALL CFFT(A, NX / 2,  - 1)
      END IF
      RETURN
      END
!**********************************************************************C
      SUBROUTINE RFFTSC(A, NX, IPACK, IFAC)
      DIMENSION A(*)

      IF (IFAC == 1) THEN
         FAC = .5 / NX
      ELSE IF (IFAC == - 1) THEN
         FAC = 0.25 / NX
      ELSE
         FAC = 1E0
      END IF

      IF (IPACK == 2) THEN
         A(2) = 0E0
      ELSE IF (IPACK == 3) THEN
         A(NX + 1) = A(2)
         A(NX + 2) = 0E0
         A(2) = 0E0
      ELSE IF (IPACK == - 3) THEN
         A(2) = A(NX + 1)
      ELSE IF (IPACK == - 2) THEN
         A(2) = 0E0
      END IF

      IF (IFAC /= 0) THEN
        DO I = 1, NX
           A(I) = A(I) * FAC
        END DO
      END IF

      RETURN
      END

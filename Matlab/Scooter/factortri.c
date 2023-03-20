#include <math.h>
#include "mex.h"

void factor(int n,
	    double *di, double *dq,
	    double *ei, double *eq,
	    double *mi, double *mq)
     /********************************************************
      % [ dt, et, mults ] = factor( N, d, e )
      % Gaussian elimination to factor a symmetric tridiagonal linear system
      % 
      % N is the order of the matrix
      % d contains the diagonal elements of the input matrix
      % e              subdiagonal in its last N-1 positions
      % dt contains the diagonal after reduction
      % et contains the upper diagonal
      % mults contains the multipliers used during elimination
      ********************************************************/

/*******************************************************************
 Paul's notes: M-functions factor.m and backsub.m (and their MEX
 counterparts) provide a solution to A*x=b, where A is a COMPLEX,
 SYMMETRIC, and TRI-DIAGONAL square matrix, and b is also COMPLEX.

 In both of these M-functions, d is the main diagonal, e is the
 sub (and super) diagonal with the first element being ignored
 (i.e. elements 2 to N are the subdiagonal elements), and m is the
 set of multipliers (for forward elimination).

 M-function factor.m produces the multipliers and modified d values
 used during forward elimination, which is the first loop in backsub. 
 The second loop in backsub performs the back-substitution.
********************************************************************/

{
  double *pmi, *pmq, *pei, *peq, *pdi, *pdq;
  double lei, leq, ldi, ldq, lmi, lmq;
  double dmag;
  int i;

  /* 
     % LU decomposition without interchanges
     dt( 1 ) = d( 1 );
     if N >= 2
       for I = 1 : N-1
         mults( I ) = e( I+1 ) / dt( I ); % multiplier
         dt( I+1  ) = d( I+1 ) - mults( I ) * e( I+1 );  % new diagonal
       end
     end
     et = e;
  */

  /* dt( 1 ) = d( 1 ); */
  /* do nothing for d(1) or d[0] */

  if (n>=2) {
    pmi=mi; pmq=mq;
    pdi=di;  pdq=dq;
    pei=ei+1; peq=eq+1;
    for (i=0; i<n-1; i++) {

      lei=*(pei++); leq=*(peq++);
      ldi=*(pdi++); ldq=*(pdq++);

      /* mults( I ) = e( I+1 ) / dt( I ); % multiplier */
      dmag=ldi*ldi+ldq*ldq;
      lmi=( lei*ldi + leq*ldq)/dmag;
      lmq=(-lei*ldq + leq*ldi)/dmag;
      *(pmi++)=lmi; *(pmq++)=lmq;

      /* dt( I+1  ) = d( I+1 ) - mults( I ) * e( I+1 );  % new diagonal */
      *(pdi) += - lmi*lei + lmq*leq;
      *(pdq) += - lmi*leq - lmq*lei;

    }
  }
  
}

/* [ mults, dt, et ] = factor( N, d, e ) */

/* Input Arguments */

#define N prhs[0]
#define D prhs[1]
#define E prhs[2]

/* Output Arguments */

#define Mt plhs[0]
#define Dt plhs[1]
#define Et plhs[2]

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int n;
  double *di, *dq, *ei, *eq;
  double *dti, *dtq, *eti, *etq, *mti, *mtq;
  double *pdi, *pdq, *pei, *peq;
  int i;

  /* start of bulletproofing junk from extern/mex/mexfunction.c */

  /* Examine input (right-hand-side) arguments. */
  /*******
  mexPrintf("\nThere are %d right-hand-side argument(s).", nrhs);
  for (i=0; i<nrhs; i++)  {
    mexPrintf("\n\tInput Arg %i is of type:\t%s ",i,mxGetClassName(prhs[i]));
  }
  *******/

  /* Examine output (left-hand-side) arguments. */
  /*******
  mexPrintf("\n\nThere are %d left-hand-side argument(s).\n", nlhs);
  if (nlhs > nrhs)
  mexErrMsgTxt("Cannot specify more outputs than inputs.\n");
  for (i=0; i<nlhs; i++)  {
    plhs[i]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[i])=mxGetNumberOfElements(prhs[i]);
  }
  ********/

  /* end of bulletproofing junk from extern/mex/mexfunction.c */

  /* Check for proper number of arguments */
  if (nrhs != 3) { 
    mexErrMsgTxt("Three input arguments required."); 
  } else if (nlhs > 3) {
    mexErrMsgTxt("Too many output arguments."); 
  } 

  /* could forego passing N and just read size of d */
  n=mxGetScalar(N);

  /* Create a matrix for the return argument */ 
  Mt = mxCreateDoubleMatrix(n, 1, mxCOMPLEX); 
  Dt = mxCreateDoubleMatrix(n, 1, mxCOMPLEX); 
  Et = mxCreateDoubleMatrix(n, 1, mxCOMPLEX); 
    
  /* Assign pointers to the various parameters */ 
  di=mxGetPr(D);
  dq=mxGetPi(D);
  ei=mxGetPr(E);
  eq=mxGetPi(E);
    
  dti=mxGetPr(Dt);
  dtq=mxGetPi(Dt);
  eti=mxGetPr(Et);
  etq=mxGetPi(Et);
  mti=mxGetPr(Mt);
  mtq=mxGetPi(Mt);
        
  /* Do the actual computations in a subroutine */
  factor(n,
	 di, dq,
	 ei, eq,
	 mti, mtq);
  /* now, copy modified d and e into output vectors */
  pdi=dti; pdq=dtq;
  pei=eti; peq=etq;
  /* can increment di, dq, ei, eq, since these will never be used
     after this loop */
  for (i=0; i<n; i++) {
    *(pdi++)=*(di++); *(pdq++)=*(dq++);
    *(pei++)=*(ei++); *(peq++)=*(eq++);
  }

  return;

}

#include <math.h>
#include "mex.h"

void backsub(int n, 
	     double *mi, double *mq, 
	     double *di, double *dq, 
	     double *ei, double *eq, 
	     double *bi, double *bq,
	     double *xi, double *xq)
/**************************
 % Based on TINVIT in EISPACK
 % performs back-substitution for a symmetric tridiagonal linear system
 % N is the order of the matrix
 % d contains the diagonal of the upper triangular matrix
 % e contains the upper diagonal
 % m contains the multipliers used during forward elimination
 % b contains the right hand side (Ax=b)
 % the answer is returned in vector b
 %
 % Michael B. Porter 7/1/85
 (adapted to C by Paul Hursky 7/6/2005)
 **************************/

/*******************************************************************
 Paul's notes: M-functions factor.m and backsub.m (and their MEX
 counterparts) provide a solution to A*x=b, where A is a COMPLEX,
 SYMMETRIC, and TRI-DIAGONAL square matrix, and b is also COMPLEX.

 In both of these M-functions, d is the main diagonal, e is the
 sub (and super) diagonal with the first element being ignored
 (i.e. elements 2 to N are the subdiagonal elements), and m is the
 set of multipliers (for forward elimination).

 M-function factor.m produces the multipliers and modified d values
 used during forward elimination, which is the first loop below. The
 second loop performs the back-substitution.
********************************************************************/

{
  int i;
  double *pbi, *pbq, *pmi, *pmq, *pdi, *pdq, *pei, *peq;
  double bi_last, bq_last, mi_last, mq_last;
  double *pxi, *pxq;
  double dmag, ni, nq;
  double lbi, lbq, ldi, ldq, lei, leq, lxi, lxq;

  /************************
  % Forward elimination
  for I = 2 : N
      b( I ) = b( I ) - mults( I-1 ) * b( I - 1 );
  end
  ************************/
  pbi=bi; pbq=bq;
  pmi=mi; pmq=mq;
  /* must increment pmi and pmq so that first element in loop is 2nd;
     note that this is used to set mi_last and mq_last, so first time
     in the loop we are using 2nd element to set mi_last and mq_last
     for use in the NEXt pass when I=3 */
  mi_last=*(pmi++); mq_last=*(pmq++);
  /* must increment pbi and pbq so that first element in loop is 2nd */
  bi_last=*(pbi++); bq_last=*(pbq++);
  for (i=2; i<=n; i++) {
    *(pbi) += - mi_last*bi_last + mq_last*bq_last;
    *(pbq) += - mi_last*bq_last - mq_last*bi_last;
    mi_last=*(pmi++); mq_last=*(pmq++);
    bi_last=*(pbi++); bq_last=*(pbq++);
  }

  /******************************************************************
  % Back-substitution (result in b)
  x( N ) = b( N ) / d( N );
  if ( N >= 2 )
     for I = N - 1 : -1 : 1
         x( I ) = ( b( I ) - e( I + 1 ) * x( I + 1 ) ) / d( I );
     end
  end
  ******************************************************************/
  pxi=xi+n-1; pxq=xq+n-1;
  pbi=bi+n-1; pbq=bq+n-1;
  pdi=di+n-1; pdq=dq+n-1;

  /*   x( N ) = b( N ) / d( N ); */
  lbi=*(pbi--); lbq=*(pbq--);
  ldi=*(pdi--); ldq=*(pdq--);
  dmag=ldi*ldi+ldq*ldq;
  /* cannot decrement x yet, because we need x(I+1) below */
  *(pxi) = ( lbi*ldi + lbq*ldq)/dmag;
  *(pxq) = (-lbi*ldq + lbq*ldi)/dmag;
  if (n>=2) {
    pei=ei+n-1; peq=eq+n-1;
    /* note that at this point have already decremented b and d once,
       but NOT e */
    for (i=n-1; i>0; i--) {
      /* x( I ) = ( b( I ) - e( I + 1 ) * x( I + 1 ) ) / d( I ); */
      lbi=*(pbi--); lbq=*(pbq--);
      lei=*(pei--); leq=*(peq--);
      ldi=*(pdi--); ldq=*(pdq--);
      /* these x values are from x(I+1), so affter grabbing these x
	 value, can decrement pxi and pxq */
      lxi=*(pxi--); lxq=*(pxq--);
      dmag=ldi*ldi+ldq*ldq;
      ni = lbi - lei*lxi + leq*lxq;
      nq = lbq - lei*lxq - leq*lxi;
      *(pxi) = ( ni*ldi + nq*ldq)/dmag;
      *(pxq) = (-ni*ldq + nq*ldi)/dmag;
    }
  }

}

/* function [x,bt] = backsub( N, mults, d, e, b ) */

/* Input Arguments */

#define N prhs[0]
#define MULTS prhs[1]
#define D prhs[2]
#define E prhs[3]
#define B prhs[4]

/* Output Arguments */

#define X plhs[0]
#define Bt plhs[1]

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int n;
  double *mi, *mq, *di, *dq, *ei, *eq, *bi, *bq, *xi, *xq;
  double *bti, *btq, *pbi, *pbq, *pbia, *pbqa;
  int i;

  /* start of bulletproofing junk from extern/mex/mexfunction.c */

  /* Examine input (right-hand-side) arguments. */
  /*********************
  mexPrintf("\nThere are %d right-hand-side argument(s).", nrhs);
  for (i=0; i<nrhs; i++)  {
    mexPrintf("\n\tInput Arg %i is of type:\t%s ",i,mxGetClassName(prhs[i]));
  }
  **********************/
  
  /* Examine output (left-hand-side) arguments. */
  /*********************
  mexPrintf("\n\nThere are %d left-hand-side argument(s).\n", nlhs);
  if (nlhs > nrhs)
  mexErrMsgTxt("Cannot specify more outputs than inputs.\n");
  for (i=0; i<nlhs; i++)  {
    plhs[i]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[i])=mxGetNumberOfElements(prhs[i]);
  }
  *********************/

  /* end of bulletproofing junk from extern/mex/mexfunction.c */

  /* Check for proper number of arguments */
  if (nrhs != 5) { 
    mexErrMsgTxt("Five input arguments required."); 
  } else if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments."); 
  } 

  /* could have read n from size of d or e */
  n=mxGetScalar(N);
    
  /* Create a matrix for the return argument */ 
  X = mxCreateDoubleMatrix(n, 1, mxCOMPLEX); 
  Bt = mxCreateDoubleMatrix(n, 1, mxCOMPLEX); 
    
  /* Assign pointers to the various parameters */ 
  mi=mxGetPr(MULTS);
  mq=mxGetPi(MULTS);
  di=mxGetPr(D);
  dq=mxGetPi(D);
  ei=mxGetPr(E);
  eq=mxGetPi(E);
  bi=mxGetPr(B);
  bq=mxGetPi(B);
    
  xi=mxGetPr(X);
  xq=mxGetPi(X);
  bti=mxGetPr(Bt);
  btq=mxGetPi(Bt);

  /* Do the actual computations in a subroutine */
  /* copy original b values into output array bt */
  pbi=bti; pbq=btq;
  pbia=bi; pbqa=bq;
  for (i=0; i<n; i++) {
    *(pbi++)=*(pbia++);
    *(pbq++)=*(pbqa++);
  }
  backsub(n, 
	  mi, mq, 
	  di, dq, 
	  ei, eq, 
	  bti, btq,
	  xi, xq);
  /* Hmmm. That's interesting - if I modify the data values of
     a Matlab variable in C, it is a SIDE-EFFECT and the variable
     retains the modified values in the calling routine. */
  /* the following version of the code above demonstrates this,
     or something: the original b values are modified... */
  /**************************
  backsub(n, 
	  mi, mq, 
	  di, dq, 
	  ei, eq, 
	  bi, bq,
	  xi, xq);
  pbi=bti; pbq=btq;
  pbia=bi; pbqa=bq;
  for (i=0; i<n; i++) {
    *(pbi++)=*(pbia++);
    *(pbq++)=*(pbqa++);
  }
  *****************************/
  return;

}

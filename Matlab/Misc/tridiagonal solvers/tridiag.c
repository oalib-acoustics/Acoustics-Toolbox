#include <stdlib.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int am, an, bm, bn, cm, cn, ym, yn, N;  /* dimensions of input arrays */
   double *a, *b, *c, *y, *x, *g;          /* locations of data */
   double beta;

   /********************************************************************/
   /*  Various checks on the input arguments  */

   /*  Check that there are the proper number of arguments on each side */
   if ( nrhs!=4 || nlhs!=1 )
      { char c[200];
        sprintf(c,"Usage:  x = %s( a, b, c, y )",mexFunctionName());
        mexErrMsgTxt(c); }

   /* Check that the input arrays have acceptable sizes, and find N */
   am=mxGetM(prhs[0]); an=mxGetN(prhs[0]);
   bm=mxGetM(prhs[1]); bn=mxGetN(prhs[1]);
   cm=mxGetM(prhs[2]); cn=mxGetN(prhs[2]);
   ym=mxGetM(prhs[3]); yn=mxGetN(prhs[3]);
   if ( ( am<=1 && an<=1 ) || ( bm<=1 && bn<=1 )
             || ( cm<=1 && cn<=1 ) || ( ym<=1 && yn<=1 ) )
      mexErrMsgTxt("a, b, c, y must be vectors");
   N = (ym>1) ? ym : yn;
   if ( ( am<N-1 && an<N-1 ) || ( bm<N && bn<N ) || ( cm<N-1 && cn<N-1 ) )
      mexErrMsgTxt("a, b, c must be vectors of length at least N-1, N, N-1");

   /*  Check that they are double-precision real numbers */
   if (     !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
         || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
         || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
         || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
      mexErrMsgTxt("Does not (yet) work for complex or other input");

   /********************************************************************/
   /*  Allocate space for output and set up pointers */
   plhs[0] = mxCreateDoubleMatrix(ym,yn,mxREAL);
   a=mxGetPr(prhs[0]);  b=mxGetPr(prhs[1]); c=mxGetPr(prhs[2]);
   y=mxGetPr(prhs[3]);  x=mxGetPr(plhs[0]);
   g = (double *) mxMalloc( N * sizeof(double) );

   /********************************************************************/
   /* Solve the problem by LU decomposition and back-substitution */
   beta=b[0];
   if ( beta==0. )
      { mxFree(g);
        mxErrMsgTxt("beta = 0  at j=0"); }
   x[0] = y[0]/beta;
   { int j;
     for ( j=1 ; j<N ; j++ )
      {
      g[j] = c[j-1]/beta;
      beta = b[j] - a[j-1]*g[j];
      if ( beta==0. )
         { char c[200];
           sprintf(c,"beta=0  at  j=%d",j);
           mxFree(g);
           mxErrMsgTxt(c); }
      x[j] = ( y[j] - a[j-1]*x[j-1] )/beta;
      }}
   { int j;
     for ( j=N-2 ; j>=0 ; j-- )   x[j] -= g[j+1]*x[j+1]; }
   mxFree(g);
}

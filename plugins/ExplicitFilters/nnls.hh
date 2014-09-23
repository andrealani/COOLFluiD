#ifndef COOLFluiD_MathTools_NonNegativeLeastSquares_hh
#define COOLFluiD_MathTools_NonNegativeLeastSquares_hh

#include "MathTools/MathChecks.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

class NonNegativeLeastSquares
{
public:
  NonNegativeLeastSquares () {
    c__1 = 1;
    c__0 = 0;
    c__2 = 2;
  }
  // virtual ~NonNegativeLeastSquares ();
  
  void solve(RealMatrix& A, RealVector& b, RealVector& x);
  
private:
  
  /*     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) */

  /*  Algorithm NNLS: NONNEGATIVE LEAST SQUARES */

  /*  The original version of this code was developed by */
  /*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory */
  /*  1973 JUN 15, and published in the book */
  /*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
  /*  Revised FEB 1995 to accompany reprCFuinting of the book by SIAM. */

  /*     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN */
  /*     N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM */

  /*                      A * X = B  SUBJECT TO X .GE. 0 */
  /*     ------------------------------------------------------------------ */
  /*                     Subroutine Arguments */

  /*     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE */
  /*                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N */
  /*                     MATRIX, A.           ON EXIT A() CONTAINS */
  /*                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN */
  /*                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY */
  /*                     THIS SUBROUTINE. */
  /*     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- */
  /*             TAINS Q*B. */
  /*     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL */
  /*             CONTAIN THE SOLUTION VECTOR. */
  /*     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE */
  /*             RESIDUAL VECTOR. */
  /*     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN */
  /*             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0. */
  /*             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z */
  /*     ZZ()     AN M-ARRAY OF WORKING SPACE. */
  /*     INDEX()     AN INT WORKING ARRAY OF LENGTH AT LEAST N. */
  /*                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS */
  /*                 P AND Z AS FOLLOWS.. */

  /*                 INDEX(1)   THRU INDEX(NSETP) = SET P. */
  /*                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z. */
  /*                 IZ1 = NSETP + 1 = NPP1 */
  /*                 IZ2 = N */
  /*     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING */
  /*             MEANINGS. */
  /*             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY. */
  /*             2     THE DIMENSIONS OF THE PROBLEM ARE BAD. */
  /*                   EITHER M .LE. 0 OR N .LE. 0. */
  /*             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. */

  /*     ------------------------------------------------------------------ */
  /* Subroutine */ 
  CFuint nnls(CFreal* a, CFuint mda, CFuint m, CFuint n, 
  	  CFreal* b, CFreal* x, CFreal* rnorm, 
  	  CFreal* w, CFreal* zz, CFuint* index, CFuint* mode);



  /*     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV) */

  /*  CONSTRUCTION AND/OR APPLICATION OF A SINGLE */
  /*  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B */

  /*  The original version of this code was developed by */
  /*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory */
  /*  1973 JUN 12, and published in the book */
  /*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
  /*  Revised FEB 1995 to accompany reprCFuinting of the book by SIAM. */
  /*     ------------------------------------------------------------------ */
  /*                     Subroutine Arguments */

  /*     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a */
  /*            Householder transformation, or Algorithm H2 to apply a */
  /*            previously constructed transformation. */
  /*     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. */
  /*     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO */
  /*            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M */
  /*            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. */
  /*     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot */
  /*            vector.  IUE is the storage increment between elements. */
  /*            On exit when MODE = 1, U() and UP contain quantities */
  /*            defining the vector U of the Householder transformation. */
  /*            on entry with MODE = 2, U() and UP should contain */
  /*            quantities previously computed with MODE = 1.  These will */
  /*            not be modified during the entry with MODE = 2. */
  /*     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH */
  /*            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE */
  /*            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED. */
  /*            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS. */
  /*     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C(). */
  /*     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C(). */
  /*     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0 */
  /*            NO OPERATIONS WILL BE DONE ON C(). */
  /*     ------------------------------------------------------------------ */
  /* Subroutine */ 
  CFuint h12(CFuint mode, CFuint* lpivot, CFuint* l1, 
  	 CFuint m, CFreal* u, CFuint* iue, CFreal* up, CFreal* c__, 
  	 CFuint* ice, CFuint* icv, CFuint* ncv);


      /*     COMPUTE ORTHOGONAL ROTATION MATRIX.. */

      /*  The original version of this code was developed by */
      /*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory 
       */
      /*  1973 JUN 12, and published in the book */
      /*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
      /*  Revised FEB 1995 to accompany reprCFuinting of the book by SIAM. */

      /*     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2)) */
      /*                        (-S,C)         (-S,C)(B)   (   0          ) */
      /*     COMPUTE SIG = SQRT(A**2+B**2) */
      /*        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT */
      /*        SIG MAY BE IN THE SAME LOCATION AS A OR B . */
  CFuint g1(CFreal* a, CFreal* b, CFreal* cterm, CFreal* sterm, CFreal* sig);

private:
  /* Table of constant values */
  CFuint c__1;
  CFuint c__0;
  CFuint c__2;
};


//////////////////////////////////////////////////////////////////////////////

  } // end namespace MathTools

} // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_NonNegativeLeastSquares_hh

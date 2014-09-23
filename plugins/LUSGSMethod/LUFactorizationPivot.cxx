#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/LUFactorizationPivot.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LUFactorizationPivot, LUSGSIteratorData, LUSGSMethodModule> LUFactorizationPivotProvider("LUFactPivot");


//////////////////////////////////////////////////////////////////////////////

LUFactorizationPivot::LUFactorizationPivot(std::string name) :
  LUFactorization(name),
  socket_pivotLUFactorization("pivotLUFactorization")
{
}

//////////////////////////////////////////////////////////////////////////////

void LUFactorizationPivot::execute()
{
  CFAUTOTRACE;

  // Gets the diagonal matrices
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();
  cf_assert(diagBlockJacobMatr.size() == getMethodData().getNbrStatesSets());

  // Gets the number of states sets
  const CFuint nbrStatesSets = diagBlockJacobMatr.size();

  // Resizes pivotLUFactorization
  DataHandle< vector< CFuint > > pivotLUFactorization = socket_pivotLUFactorization.getDataHandle();
  pivotLUFactorization.resize(nbrStatesSets);

  // get isStatesSetParUpdatable data handle
  DataHandle< bool > isStatesSetParUpdatable = socket_isStatesSetParUpdatable.getDataHandle();

  // Loops over the states sets
  for (CFuint iSet = 0; iSet < nbrStatesSets; ++iSet)
  {
    if (isStatesSetParUpdatable[iSet])
    {
      // Dereferences the current matrix
      RealMatrix& currDiagMatrix = diagBlockJacobMatr[iSet];

      //Dereferences the pivot elements
      vector< CFuint >& currPivot = pivotLUFactorization[iSet];

      // number of states in the current set
      const CFuint resSize = currDiagMatrix.nbRows();
      
      // initialize currPivot
      currPivot.resize(resSize);
      for (CFuint iRow = 0; iRow < resSize; ++iRow)
      {
        currPivot[iRow] = iRow;
      }

      // actual LU factorization
      // loop over the diagonal elements
      const CFuint resSizeM1 = resSize - 1;
      for (CFuint iDiag = 0; iDiag < resSizeM1; ++iDiag)
      {
        // PIVOTING
        // find largest element in absolute value in below current diagonal element
        CFreal max = std::abs(currDiagMatrix(iDiag,iDiag));
        CFuint maxValRow = iDiag;
        for (CFuint iRow = iDiag+1; iRow < resSize; ++iRow)
        {
          const CFreal absVal = std::abs(currDiagMatrix(iRow,iDiag));
          if (absVal > max)
          {
            max = absVal;
            maxValRow = iRow;
          }
        }

        // if necessary, update currPivot and swap rows
        if (iDiag != maxValRow)
        {
          const CFuint swap = currPivot[iDiag];
          currPivot[iDiag] = currPivot[maxValRow];
          currPivot[maxValRow] = swap;

          swapRows(iDiag,maxValRow,currDiagMatrix);
        }

        // LU FACTORIZATION
        factorizeMatrix(iDiag,currDiagMatrix);
      }
    }
  }

  // set second element of socket_statesSetIdx to 0 --> updateCoefs are not recomputed
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  statesSetIdx[1] = 0;
}

//////////////////////////////////////////////////////////////////////////////

void LUFactorizationPivot::swapRows(const CFuint row1, const CFuint row2, RealMatrix& matrix)
{
  cf_assert(row1 < matrix.nbRows());
  cf_assert(row2 < matrix.nbRows());
  const CFuint nbrCols = matrix.nbCols();
  for (CFuint iCol = 0; iCol < nbrCols; ++iCol)
  {
    const CFreal swap = matrix(row1,iCol);
    matrix(row1,iCol) = matrix(row2,iCol);
    matrix(row2,iCol) = swap;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LUFactorizationPivot::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  LUFactorization::setup();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > LUFactorizationPivot::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = LUFactorization::needsSockets();

  result.push_back(&socket_pivotLUFactorization);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

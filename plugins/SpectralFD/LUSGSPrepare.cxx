#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/LUSGSPrepare.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LUSGSPrepare, SpectralFDMethodData, SpectralFDModule> lusgsPrepareProvider("LUSGSPrepare");

//////////////////////////////////////////////////////////////////////////////

LUSGSPrepare::LUSGSPrepare(const std::string& name) : SpectralFDMethodCom(name),
    socket_diagBlockJacobMatr("diagBlockJacobMatr"),
    socket_updateCoeff("updateCoeff"),
    socket_gradients("gradients")
{
}

//////////////////////////////////////////////////////////////////////////////

LUSGSPrepare::~LUSGSPrepare()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSPrepare::execute()
{
  CFAUTOTRACE;

  // reset the diagonal block Jacobian matrices
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();
  const CFuint nbrCells = diagBlockJacobMatr.size();
  for (CFuint iCell = 0; iCell < nbrCells; ++iCell)
  {
    diagBlockJacobMatr[iCell] = 0.;
  }

  // reset the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff = 0.0;

  // reset the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  const CFuint nGrads = gradients.size();
  for (CFuint iGrad = 0; iGrad < nGrads; ++iGrad)
  {
    const CFuint nVar = gradients[iGrad].size();
    for (CFuint iVar = 0; iVar < nVar; ++iVar)
    {
      gradients[iGrad][iVar] = 0.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LUSGSPrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_diagBlockJacobMatr);
  result.push_back(&socket_updateCoeff       );
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

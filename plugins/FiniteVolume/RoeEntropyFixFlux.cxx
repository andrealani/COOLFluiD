#include "RoeEntropyFixFlux.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RoeEntropyFixFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
roeEntropyFixFluxProvider("RoeEntropyFix");
      
//////////////////////////////////////////////////////////////////////////////

RoeEntropyFixFlux::RoeEntropyFixFlux(const std::string& name) :
  RoeFlux(name)
{
  addConfigOptionsTo(this);
  _entropyFixID = 1;
  setParameter("entropyFixID",&_entropyFixID);
}

//////////////////////////////////////////////////////////////////////////////

RoeEntropyFixFlux::~RoeEntropyFixFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void RoeEntropyFixFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >
    ("entropyFixID","ID of the entropy correction type.");
}

//////////////////////////////////////////////////////////////////////////////

void RoeEntropyFixFlux::setAbsEigenValues()
{
  //compute eigen values of the left state
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  updateVarSet->computeEigenValues(pdata[0],
				   unitNormal,
				   _leftEvalues);
  
  //compute eigen values of the right state
  updateVarSet->computeEigenValues(pdata[1],
				   unitNormal,
				   _rightEvalues);
  
  CFreal lambdaCorr = 0.0;
  computeLambdaCorr(lambdaCorr);
  lambdaCorr *= 0.5;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  if (_entropyFixID == 1) {
    for (CFuint i = 0; i < nbEqs; ++i) {
      _absEvalues[i] = max(std::abs(_eValues[i]), lambdaCorr);
    }
  }
  else if (_entropyFixID == 2) {
    const CFreal twoLambdaCorr = 2.0*lambdaCorr;
    for (CFuint i = 0; i < nbEqs; ++i) {
      const CFreal absLambda = std::abs(_eValues[i]);
      _absEvalues[i] = (!(std::abs(_eValues[i]) < twoLambdaCorr)) ?
	absLambda : (absLambda*absLambda/(4.0*lambdaCorr) + lambdaCorr);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
    
      } // namespace FiniteVolume

    } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

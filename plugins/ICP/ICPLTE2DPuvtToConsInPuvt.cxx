#include "ICP/ICP.hh"
#include "ICPLTE2DPuvtToConsInPuvt.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ICPLTE2DPuvtToConsInPuvt, 
			    VarSetMatrixTransformer, 
			    ICPModule, 1> 
ICPLTE2DPuvtToConsInPuvtProvider("ICPLTE2DPuvtToConsInPuvt");

//////////////////////////////////////////////////////////////////////////////

ICPLTE2DPuvtToConsInPuvt::ICPLTE2DPuvtToConsInPuvt
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  LTE::Euler2DPuvtLTEToConsInPuvtLTE(model)
{
}

//////////////////////////////////////////////////////////////////////////////

ICPLTE2DPuvtToConsInPuvt::~ICPLTE2DPuvtToConsInPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void ICPLTE2DPuvtToConsInPuvt::setMatrix(const RealVector& state)
{
  LTE::Euler2DPuvtLTEToConsInPuvtLTE::setMatrix(state);
  
  const CFuint totalNbEqsMin2 = state.size()-2;
  const CFuint totalNbEqsMin1 = state.size()-1;
  
  _transMatrix(totalNbEqsMin2,totalNbEqsMin2) = 1.0;
  _transMatrix(totalNbEqsMin1,totalNbEqsMin1) = 1.0;
  
  // CFLog(INFO, "matrix = \n" << _transMatrix << "\n");

  // for (CFuint i = 0; i < state.size(); ++i) {
  //   if (_transMatrix(i,i) < 0.) {
  //     CFLog(WARN, "negative diagonal entry (" << i << "," << i << ") => " << _transMatrix(i,i) << "\n");
  //     abort();
  //   }
  // }
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICP/ICP.hh"
#include "ICPLTE2DdPuvtToPuvt.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ICPLTE2DdPuvtToPuvt, VarSetTransformer, ICPModule,1>
icpLTE2DdPuvtToPuvtProvider("ICPLTE2DdPuvtToPuvt");

//////////////////////////////////////////////////////////////////////

ICPLTE2DdPuvtToPuvt::ICPLTE2DdPuvtToPuvt
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _icpModel(model->getConvectiveTerm().d_castTo<PTERM>()),
  _dhe(3)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
ICPLTE2DdPuvtToPuvt::~ICPLTE2DdPuvtToPuvt()
{
}

//////////////////////////////////////////////////////////////////////

void ICPLTE2DdPuvtToPuvt::transform(const State& state, State& result)
{
  result.copyData(state);
  result[0] += _icpModel->getPressInfComp();
}

//////////////////////////////////////////////////////////////////////

void ICPLTE2DdPuvtToPuvt::transformFromRef(const RealVector& data, State& result)
{
  // here we assume that E's are the last two components
  const CFuint firstScalarVar = _icpModel->getDataSize() - 2;
  result[0] = data[EulerTerm::P] + _icpModel->getPressInfComp();
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = data[EulerTerm::T]; 
  result[4] = data[firstScalarVar];
  result[5] = data[firstScalarVar+1];
}

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

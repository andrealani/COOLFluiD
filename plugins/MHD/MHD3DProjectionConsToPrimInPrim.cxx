#include "MHD/MHD.hh"
#include "MHD3DProjectionConsToPrimInPrim.hh"
#include "MHDProjectionTerm.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionConsToPrimInPrim, VarSetMatrixTransformer, MHDModule, 1> 
mhd3DProjectionConsToPrimInPrimProvider("MHD3DProjectionConsToPrimInPrim");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionConsToPrimInPrim::MHD3DProjectionConsToPrimInPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDProjectionTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionConsToPrimInPrim::~MHD3DProjectionConsToPrimInPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionConsToPrimInPrim::setMatrix(const RealVector& state)
{
  cf_assert(_model.isNotNull());

  const CFreal ovrbar = 1./state[0];
  const CFreal gm1 = _model->getGamma() - 1.;
  const CFreal ubar = state[1];
  const CFreal vbar = state[2];
  const CFreal wbar = state[3];
  const CFreal Bxbar = state[4];
  const CFreal Bybar = state[5];
  const CFreal Bzbar = state[6];
  
  _transMatrix(0,0) = 1.;
  _transMatrix(1,0) = -ubar*ovrbar;
  _transMatrix(1,1) = ovrbar;
  _transMatrix(2,0) = -vbar*ovrbar;
  _transMatrix(2,2) = ovrbar;
  _transMatrix(3,0) = -wbar*ovrbar;
  _transMatrix(3,3) = ovrbar;
  _transMatrix(4,4) = 1.;
  _transMatrix(5,5) = 1.;
  _transMatrix(6,6) = 1.;
  _transMatrix(7,0) = gm1*0.5*(ubar*ubar + vbar*vbar + wbar*wbar);
  _transMatrix(7,1) = -gm1*ubar;
  _transMatrix(7,2) = -gm1*vbar;
  _transMatrix(7,3) = -gm1*wbar;
  _transMatrix(7,4) = -gm1*Bxbar;
  _transMatrix(7,5) = -gm1*Bybar;
  _transMatrix(7,6) = -gm1*Bzbar;
  _transMatrix(7,7) = gm1;
  _transMatrix(8,8) = 1.;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHD.hh"
#include "MHD3DConsToPrimEInRef.hh"
#include "MHDPhysicalModel.hh"
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

Environment::ObjectProvider<MHD3DConsToPrimEInRef, VarSetMatrixTransformer, MHDModule, 1> mhd3DConsToPrimEInRefProvider("MHD3DConsToPrimEInRef");

//////////////////////////////////////////////////////////////////////////////

MHD3DConsToPrimEInRef::MHD3DConsToPrimEInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DConsToPrimEInRef::~MHD3DConsToPrimEInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DConsToPrimEInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());

  const RealVector& linearData = _model->getPhysicalData();

  const CFreal ovrbar = 1./linearData[MHDTerm::RHO];
  const CFreal gm1 = _model->getGamma() - 1.;
  const CFreal ubar = linearData[MHDTerm::VX];
  const CFreal vbar = linearData[MHDTerm::VY];
  const CFreal wbar = linearData[MHDTerm::VZ];
  const CFreal Bxbar = linearData[MHDTerm::BX];
  const CFreal Bybar = linearData[MHDTerm::BY];
  const CFreal Bzbar = linearData[MHDTerm::BZ];

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
  _transMatrix(7,4) = 0.0;
  _transMatrix(7,5) = 0.0;
  _transMatrix(7,6) = 0.0;
  _transMatrix(7,7) = gm1;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

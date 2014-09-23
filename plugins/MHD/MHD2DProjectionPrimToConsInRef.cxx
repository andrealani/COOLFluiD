#include "MHD/MHD.hh"
#include "MHD2DProjectionPrimToConsInRef.hh"
#include "MHDProjection.hh"
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

Environment::ObjectProvider<MHD2DProjectionPrimToConsInRef, VarSetMatrixTransformer, MHDModule, 1> mhd2DProjectionPrimToConsInRefProvider("MHD2DProjectionPrimToConsInRef");

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionPrimToConsInRef::MHD2DProjectionPrimToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDProjectionTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionPrimToConsInRef::~MHD2DProjectionPrimToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrimToConsInRef::setMatrixFromRef()
{
  const RealVector& linearData = _model->getPhysicalData();

  const CFreal rhat = linearData[MHDProjectionTerm::RHO];
  const CFreal gm1 = _model->getGamma() - 1.;
  const CFreal uhat = linearData[MHDProjectionTerm::VX];
  const CFreal vhat = linearData[MHDProjectionTerm::VY];
  const CFreal what = linearData[MHDProjectionTerm::VZ];
  const CFreal Bxhat = linearData[MHDProjectionTerm::BX];
  const CFreal Byhat = linearData[MHDProjectionTerm::BY];
  const CFreal Bzhat = linearData[MHDProjectionTerm::BZ];

  _transMatrix(0,0) = 1.0;
  _transMatrix(1,0) = uhat;
  _transMatrix(1,1) = rhat;
  _transMatrix(2,0) = vhat;
  _transMatrix(2,2) = rhat;
  _transMatrix(3,0) = what;
  _transMatrix(3,3) = rhat;
  _transMatrix(4,4) = 1.0;
  _transMatrix(5,5) = 1.0;
  _transMatrix(6,6) = 1.0;
  _transMatrix(7,0) = 0.5*(uhat*uhat + vhat*vhat + what*what);
  _transMatrix(7,1) = rhat*uhat;
  _transMatrix(7,2) = rhat*vhat;
  _transMatrix(7,3) = rhat*what;
  _transMatrix(7,4) = Bxhat;
  _transMatrix(7,5) = Byhat;
  _transMatrix(7,6) = Bzhat;
  _transMatrix(7,7) = 1./gm1;
  _transMatrix(8,8) = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

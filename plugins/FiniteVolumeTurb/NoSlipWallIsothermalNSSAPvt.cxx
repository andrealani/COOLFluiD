#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolumeTurb/FiniteVolumeSA.hh"
#include "FiniteVolumeTurb/NoSlipWallIsothermalNSSAPvt.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

/////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallIsothermalNSSAPvt,
                      CellCenterFVMData,
                      FiniteVolumeSAModule>
noSlipWallIsothermalNSSAPvtFVMCCProvider("NoSlipWallIsothermalNSSAPvtFVMCC");

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNSSAPvt::NoSlipWallIsothermalNSSAPvt(const std::string& name) :
  NoSlipWallIsothermalNSPvt(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNSSAPvt::~NoSlipWallIsothermalNSSAPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSSAPvt::setup()
{
CFAUTOTRACE;

 NoSlipWallIsothermalNSPvt::setup();

}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSSAPvt::setGhostState(GeometricEntity *const face)
{
  CFAUTOTRACE;
  NoSlipWallIsothermalNSPvt::setGhostState(face);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  cf_assert(m_factor > 0. ||  m_factor < 0.);

  if    (!(m_factor > 0. || m_factor < 0.)) {
    std::cout<<"cf_assertion:: m_factor.in << NoSlipWallIsothermalNSSAPvt.cxx >> m_factor :: " << m_factor << "\n";
    exit(0);
  }

  if (dim == DIM_3D) {
    (*ghostState)[5] = -(*innerState)[5]/m_factor;
  }
  else {

    //  cout << "minus --- |" << (*innerState)[4] << "\n";

    (*ghostState)[4] = -(*innerState)[4]/m_factor;  //to impose zero at the wall
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

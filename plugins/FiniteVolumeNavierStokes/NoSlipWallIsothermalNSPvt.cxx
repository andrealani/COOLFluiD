#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSPvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallIsothermalNSPvt,
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
noSlipWallIsothermalNSPvtFVMCCProvider("NoSlipWallIsothermalNSPvtFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSPvt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("TWall","Wall temperature");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNSPvt::NoSlipWallIsothermalNSPvt(const std::string& name) :
  FVMCC_BC(name),
  m_refTemp(1.)
{
  addConfigOptionsTo(this); 
  
  m_wallTemp = 0.0;
  setParameter("TWall",&m_wallTemp);
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNSPvt::~NoSlipWallIsothermalNSPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSPvt::setup()
{
  FVMCC_BC::setup();
  
  cf_assert(m_wallTemp > 0.0);

  m_tempNode.resize(PhysicalModelStack::getActive()->getDim());
  m_midNode.resize(PhysicalModelStack::getActive()->getDim());
  m_tempGhostNode.resize(PhysicalModelStack::getActive()->getDim());
  m_faceNormal.resize(PhysicalModelStack::getActive()->getDim());

  // adimensionalize the temperature
  SafePtr<EulerTerm> eulerTerm = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();
  cf_assert(eulerTerm.isNotNull());
  m_refTemp = eulerTerm->getTempRef();
  m_wallTemp /= m_refTemp;
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSPvt::setGhostState(GeometricEntity *const face)
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  if (!getPutGhostsOnFace()) {
    // set the current normal
    const CFuint faceID = face->getID();
    const CFuint startID = faceID*dim;
    DataHandle<CFreal> normals = socket_normals.getDataHandle();
    for (CFuint i = 0; i < dim; ++i) {
      m_faceNormal[i] = normals[startID + i];
    }
    
    Node& innerNode = innerState->getCoordinates();
    MathFunctions::computeProjectedPoint(*face->getNode(0), m_faceNormal, innerNode, m_tempGhostNode, 2.);
    
    // this middle node is by construction on the boundary face
    m_midNode = 0.5*(innerNode + m_tempGhostNode);
    
    // here a fix is needed in order to have always ghostT > 0
    // dynamic relocation of the ghost state: the position of the
    // ghost state is locally changed, and the BC is imposed
    // using a weighted average of ghost state (in the new location)
    // and inner state
    const CFreal innerT = (dim == DIM_3D) ? (*innerState)[4] : (*innerState)[3];
    CFreal ghostT = 2.*m_wallTemp - innerT;
    m_factor = 1.0;
    while (ghostT < 0.) {
      m_factor *= m_coeff;
      
      // new position of the ghost node
      m_tempNode = ((m_factor - 1.)*m_midNode  + m_tempGhostNode)/m_factor;
      
      // new temperature in the ghost state
      ghostT = ((m_factor + 1.)*m_wallTemp - innerT)/m_factor;
      
      // move the ghost to the new position
      m_tempGhostNode = m_tempNode;
    }
    
    cf_assert(ghostT > 0.);
    
    // reset the ghost node with the new position
    ghostState->getCoordinates() = m_tempGhostNode;
    
    (*ghostState)[0] = (*innerState)[0];
    (*ghostState)[1] = -(*innerState)[1]/m_factor;
    (*ghostState)[2] = -(*innerState)[2]/m_factor;
    
    if (dim == DIM_3D) {
      (*ghostState)[3] = -(*innerState)[3]/m_factor;
      (*ghostState)[4] = ghostT;
    }
    else {
      (*ghostState)[3] = ghostT;
    }
  }
  else {
    // (*ghostState)[0] = (*innerState)[0];
//     (*ghostState)[1] = 0. ;
//     (*ghostState)[2] = 0. ;
    
//     if (dim == DIM_3D) {
//       (*ghostState)[3] = 0.;
//       (*ghostState)[4] = m_wallTemp;
//     }
//     else {
//       (*ghostState)[3] = m_wallTemp;
//     }
    (*ghostState)[0] = (*innerState)[0];
    (*ghostState)[1] = -(*innerState)[1];
    (*ghostState)[2] = -(*innerState)[2];
    
    
    if (dim == DIM_3D) {
      (*ghostState)[3] = -(*innerState)[3];
      const CFreal tGhost = 2.*m_wallTemp - (*innerState)[4];
      (*ghostState)[4] = (tGhost > 0.) ? tGhost : m_wallTemp;
    }
    else {
      const CFreal tGhost = 2.*m_wallTemp - (*innerState)[3];
      (*ghostState)[3] = (tGhost > 0.) ? tGhost : m_wallTemp;
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

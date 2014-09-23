#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/MirrorVelocity.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MirrorVelocity, CellCenterFVMData,
		      FiniteVolumeModule>
mirrorVelocity("MirrorVelocityFVMCC");

      //////////////////////////////////////////////////////////////////////

MirrorVelocity::MirrorVelocity(const std::string& name) :
  FVMCC_BC(name),
  m_isVelocityComp()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

void MirrorVelocity::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

MirrorVelocity::~MirrorVelocity()
{
}

//////////////////////////////////////////////////////////////////////////////

void MirrorVelocity::setup()
{
  FVMCC_BC::setup();

  if(m_velocityIDs.size() == 0) {
    CFLog(NOTICE, "MirrorVelocity::setup() => choosing default\n");
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    m_velocityIDs.resize(dim);
    for (CFuint i = 0 ; i < dim; ++i) {
      m_velocityIDs[i] = 1 + i;
    }
  }

  m_isVelocityComp.resize(PhysicalModelStack::getActive()->getNbEq());
  m_isVelocityComp = false;
  for (CFuint i = 0 ; i < m_velocityIDs.size(); ++i) {
    m_isVelocityComp[m_velocityIDs[i]] = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void MirrorVelocity::setGhostState
(Framework::GeometricEntity *const face)
{
  using namespace COOLFluiD::Framework;

  // unused // const CFuint dim = PhysicalModelStack::getActive()->getDim();
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  const CFuint faceID = face->getID();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint startID = faceID*dim;
  
  DataHandle< CFreal> normals = socket_normals.getDataHandle();
  // unused // const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused // const CFuint iEqSS = eqSS.getEqSS();

   CFreal vn = 0.;
   CFreal area2 = 0.;
   cf_assert(m_velocityIDs.size() == dim);

   for (CFuint i = 0; i < m_velocityIDs.size(); ++i) {
     // cout << m_velocityIDs[i] << " ";
     const CFreal nxComp = normals[startID + i];
     vn += nxComp*(*innerState)[m_velocityIDs[i]];
     area2 += nxComp*nxComp;
   }

   CFuint jxx = 0;
   for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
     if (!m_isVelocityComp[i]) {
       (*ghostState)[i] = (*innerState)[i];
     }
     else {
       const CFuint nxID = startID + jxx;
       (*ghostState)[i] = (*innerState)[i] - 2.0*vn*normals[nxID]/area2;
       jxx++;
     }
   }

 CFLog(DEBUG_MAX, "MirrorVelocity::setGhostState() => ghostState = " << *ghostState << "\n"); 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

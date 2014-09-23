#include "ArcJetST.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalConsts.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

#include "NavierStokes/EulerTerm.hh"

#include "ArcJet/ArcJetInductionTerm.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ArcJetST, 
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumeArcJetModule> 
arcJetSTFVMCCProvider("ArcJetST");

//////////////////////////////////////////////////////////////////////////////

ArcJetST::ArcJetST(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_LorentzForce("LorentzForce"),
  m_library(CFNULL),
  m_normal(),
  m_curlB(),
  m_LorentzForce(),
  m_B(),
  m_u(),
  m_J(),
  m_gradB()
{
}

//////////////////////////////////////////////////////////////////////////////

ArcJetST::~ArcJetST()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
ArcJetST::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result = 
    ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_LorentzForce);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ArcJetST::setup()
{
  using namespace COOLFluiD::Framework;

  ComputeSourceTermFVMCC::setup();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert (m_library.isNotNull());
  
  const CFuint nbCells = socket_volumes.getDataHandle().size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  socket_LorentzForce.getDataHandle().resize(nbCells*dim);
  
  m_normal.resize(DIM_3D, 0.);
  m_curlB.resize(DIM_3D, 0.);
  m_LorentzForce.resize(DIM_3D, 0.);
  m_B.resize(DIM_3D, 0.);
  m_u.resize(DIM_3D, 0.);
  m_J.resize(DIM_3D, 0.);
  m_gradB.resize(DIM_3D, DIM_3D, 0.);
}

//////////////////////////////////////////////////////////////////////////////
      
void ArcJetST::computeSource
(Framework::GeometricEntity *const element, RealVector& source, RealMatrix& jacobian)
{  
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Physics::NavierStokes;
  
  const EquationSubSysDescriptor& eqSSD = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint nbEqs = eqSSD.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // this is needed for coupling
  if (eqSSD.getEqSS() == 1 || nbEqs == totalNbEqs) {
    CFLogDebugMin( "ArcJetST::computeSource()" << "\n");
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint elemID = element->getID();
    DataHandle<CFreal> LorentzForce = socket_LorentzForce.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
    DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
    DataHandle<RealVector> nstates = this->_sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
    DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();
    
    SafePtr<EulerTerm> eulerTerm = PhysicalModelStack::getActive()->getImplementor()->
      getConvectiveTerm().d_castTo<EulerTerm>();
      
    // please note the reference &
    const State& currState = *element->getState(0);
    this->getMethodData().getUpdateVar()->computePhysicalData(currState, this->m_pdataArray);
    
    // semplification for now, no support for multi-species or turbulence
    typedef Physics::ArcJet::ArcJetInductionTerm<Physics::NavierStokes::EulerTerm> PTERM;
    
    for (CFuint d = 0; d < DIM_3D; ++d) {
      m_u[d] = this->m_pdataArray[PTERM::VX+d];
      m_B[d] = this->m_pdataArray[PTERM::BX+d];
    }
    
    // make sure that uz = 0 in 2D
    if (dim == DIM_2D) {m_u[ZZ] = 0.;}
    
    // compute the gradients of B by applying Green Gauss in the cell volume
    m_gradB = 0.;
    
    // heere it is assumed that Bx,By,Bz,Phi are the last 4 equations
    const CFuint BxID = totalNbEqs-4;
    const CFuint ByID = totalNbEqs-3;
    const CFuint BzID = totalNbEqs-2;
    
    const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
    const CFuint nbFaces = faces.size();
    for (CFuint i = 0; i < nbFaces; ++i) {
      const GeometricEntity *const face = element->getNeighborGeo(i);
      const CFuint faceID = face->getID();
      const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
      const CFuint nbFaceNodes = face->nbNodes();
      const CFreal ovNbFaceNodes = 1./(CFreal)nbFaceNodes;
      
      // store the outward face normal (scaled with the corresponding face area)
      // fill in only the components that are available (==dim)
      // in 2D, m_normal[ZZ] = 0
      const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;
      for (CFuint d = 0; d < dim; ++d) {
	m_normal[d] = normals[startID+d]*factor;
      }
      
      // compute contribution of this face to the Green Gauss integral:
      // \bar{\vec{B}}_f \vec{n}_f = \sum_{i \in f} B_i/N_f \vec{n}_f
      for (CFuint n = 0; n < nbFaceNodes; ++n) {
	const CFuint nodeID = face->getNode(n)->getLocalID();
	const RealVector& nodalState = nstates[nodeID];
	// consider all 3D components here even in 2D
	for (CFuint d = 0; d < DIM_3D; ++d) {
	  const CFreal nd = m_normal[d];
	  m_gradB(XX, d) += nd*nodalState[BxID]*ovNbFaceNodes;
	  m_gradB(YY, d) += nd*nodalState[ByID]*ovNbFaceNodes;
	  m_gradB(ZZ, d) += nd*nodalState[BzID]*ovNbFaceNodes;
	} 
      }
    }
    m_gradB /= volumes[elemID];
    
    // m_gradB(i,d) = d(B_i)/d(x_d)
    m_curlB[XX] = m_gradB(ZZ,YY) - m_gradB(YY,ZZ);
    m_curlB[YY] = m_gradB(XX,ZZ) - m_gradB(ZZ,XX);
    m_curlB[ZZ] = m_gradB(YY,XX) - m_gradB(XX,YY);
    
    const CFreal ovMuO = 1./PhysicalConsts::VacuumPermeability(); 
    m_J = m_curlB*ovMuO; 
    
    // here fill in Lorentz force
    MathFunctions::crossProd(m_J, m_B, m_LorentzForce);
    
    const CFuint estart = elemID*dim;
    for (CFuint d = 0; d < dim; ++d) {
      source[m_velIDs[d]] = m_LorentzForce[d];
      // store the Lorentz force in a data socket to be able to plot it
      LorentzForce[estart+d] = m_LorentzForce[d];
    }
    
    const CFuint TID = m_velIDs[dim-1]+1;
    CFreal pdim = eulerTerm->getPressureFromState(this->m_pdataArray[PTERM::P]);
    cf_assert(pdim > 0.01);
    CFreal Tdim = this->m_pdataArray[PTERM::T];
    cf_assert(Tdim > 0.01);
    CFreal* tVec = CFNULL; // AL: this will need to change for multi-temperature
    const CFreal sigma = m_library->sigma(Tdim, pdim, tVec);
    const CFreal J2 = MathFunctions::innerProd(m_J, m_J)/sigma; //Joule heating term
    source[TID] = MathFunctions::innerProd(m_u, m_LorentzForce) + J2;
    
    source *= volumes[elemID];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet
    
  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

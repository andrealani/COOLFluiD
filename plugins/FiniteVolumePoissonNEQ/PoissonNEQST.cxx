#include "FiniteVolumePoissonNEQ/PoissonNEQST.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalConsts.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/MultiScalarTerm.hh"
#include "FiniteVolumePoissonNEQ/FiniteVolumePoissonNEQ.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumePoissonNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<PoissonNEQST, 
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumePoissonNEQModule> 
poissonNEQSTFVMCCProvider("PoissonNEQST");

//////////////////////////////////////////////////////////////////////////////

PoissonNEQST::PoissonNEQST(const std::string& name) :
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

PoissonNEQST::~PoissonNEQST()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
PoissonNEQST::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result = 
    ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_LorentzForce);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void PoissonNEQST::setup()
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
      
void PoissonNEQST::computeSource
(Framework::GeometricEntity *const element, RealVector& source, RealMatrix& jacobian)
{  
  
  const EquationSubSysDescriptor& eqSSD = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint nbEqs = eqSSD.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // this is needed for coupling
  if (eqSSD.getEqSS() == 1 || nbEqs == totalNbEqs) {
    CFLogDebugMin( "PoissonNEQST::computeSource()" << "\n");
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint elemID = element->getID();
    DataHandle<CFreal> LorentzForce = socket_LorentzForce.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
    DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
    DataHandle<RealVector> nstates = this->_sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
    DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();

    SafePtr<MultiScalarTerm<EulerTerm > > eulerTerm = PhysicalModelStack::getActive()->getImplementor()->
      getConvectiveTerm().d_castTo<MultiScalarTerm<EulerTerm > >();
    
    State *const currState = element->getState(0);
    typedef EulerTerm PTERM;
    
    // compute physical data array, for which entries are defined in MultiScalarTerm<NavierStokes::EulerTerm>
    // EulerTerm:
    // {RHO=0, P=1, H=2, E=3, A=4, T=5, V=6, VX=7, VY=8, VZ=9, GAMMA=10, XP=11, YP=12, ZP=13}
    // MultiScalarTerm:
    // 0) mass fractions (ZP+1 -> ZP+nbSpecies)
    // 1) vibrational/electronic energies per unit volume (ZP+nbSpecies+1 ->...)
    this->getMethodData().getUpdateVar()->computePhysicalData(*currState, this->m_pdataArray);

    const CFuint nbSpecies = eulerTerm->getNbScalarVars(0);
    const CFuint nbTv = eulerTerm->getNbScalarVars(1); // number of Tv(s) and/or Te 
    
    // computation of cell-centered value of electrical conductivity
    CFreal pdim = this->m_pdataArray[PTERM::P];
    cf_assert(pdim > 0.01);
    CFreal Tdim = this->m_pdataArray[PTERM::T];
    cf_assert(Tdim > 0.01);
    CFreal* tVec = &(*currState)[nbSpecies+dim+1]; // array pointing to the Tv(s)/Te
    const CFreal sigma = m_library->sigma(Tdim, pdim, tVec);
    CFLog(DEBUG_MAX, "PoissonNEQST::computeSource() => sigma = " << sigma << "\n");
    
    // velocities (3 components needed for cross product, even in 2D)
    for (CFuint d = 0; d < DIM_3D; ++d) {
      m_u[d] = this->m_pdataArray[PTERM::VX+d];
    }
    
    // make sure that uz = 0 in 2D
    if (dim == DIM_2D) {m_u[ZZ] = 0.;}

    const CFuint startID = elemID*totalNbEqs;
    const CFuint phiID   = totalNbEqs-1;
    
    // gradients of phi computed from cell-centered gradients (ux, uy, uz) used for Least Square Reconstruction
    // uxi stores d/dx for each component of each state:
    // (drho0/dxi, drho1/dxi ... dTe/dxi)_S0, (drho0/dxi, drho1/dxi... dTe/dxi)_S1, ... (drho0/dxi, drho1/dxi... dTe/dxi)_SN-1    
    const CFreal dPhidX = m_ux[startID+phiID]; 
    const CFreal dPhidY = m_uy[startID+phiID];


    // Vatsalya

    
    
    /* The following is just some example for a different model for which Bx,By,Bz are state variables
       you might find something useful in it... or not
       
    // compute the gradients of B by applying Green Gauss in the cell volume
    m_gradB = 0.;

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
    const CFreal J2 = MathFunctions::innerProd(m_J, m_J)/sigma; //Joule heating term
    source[TID] = MathFunctions::innerProd(m_u, m_LorentzForce) + J2;
    */



    
    // only at the end we multiply for the cell volume
    source *= volumes[elemID];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumePoissonNEQ
    
  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

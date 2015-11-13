#include "ArcJetPhiST.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;
  
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ArcJetPhiST, 
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumeArcJetModule> 
arcJetPhiSTFVMCCProvider("ArcJetPhiST");
      
//////////////////////////////////////////////////////////////////////////////

void ArcJetPhiST::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFreal> >("Bfield", "Three components of magnetic field.");
  options.addConfigOption< CFreal >("ElectrodeX", "X component of the electrode position.");
  options.addConfigOption< CFreal >("ElectrodeRadius", "Radius of the electrode.");
  options.addConfigOption< CFreal >("ImposedCurrent", "Value for the imposed current");
  options.addConfigOption< bool , Config::DynamicOption<> >("DisableCurrent","Disables the current in the energy source");
}
      
//////////////////////////////////////////////////////////////////////////////

ArcJetPhiST::ArcJetPhiST(const std::string& name) :
  ComputeSourceTermFVMCC(name),
//   socket_CellID("CellID"), 
//   socket_kappa("kappa"),
//   socket_Ttable("Ttable"),
//   socket_Ptable("Ptable"),
  socket_volumes("volumes"),
  socket_normals("normals"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_isOutward("isOutward"),
  socket_LorentzForce("LorentzForce"),
  socket_Jx("Jx"),  
  socket_Jy("Jy"),
  socket_Jz("Jz"),
  m_library(CFNULL),
  m_normal(),
  m_curlB(),
  m_LorentzForce(),
  m_B(),
  m_u(),
  m_uf(),
  m_Bf(),
  m_xf(),
  m_UxB(),
  m_J(),
  m_gradPhi()
{
  addConfigOptionsTo(this);
  
  m_Bfield = vector<CFreal>();
  setParameter("Bfield",&m_Bfield);
  
  m_electrodeX = 0.;
  setParameter("ElectrodeX",&m_electrodeX); 
  
  m_electrodeRadius = 0.;
  setParameter("ElectrodeRadius",&m_electrodeRadius);

  m_imposedI = 0.;
  this->setParameter("ImposedCurrent", &m_imposedI);
  
  m_disableCurr = 0;
  this->setParameter("DisableCurrent", &m_disableCurr);
}

//////////////////////////////////////////////////////////////////////////////

ArcJetPhiST::~ArcJetPhiST()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ArcJetPhiST::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
//   result.push_back(&socket_CellID);
//   result.push_back(&socket_kappa);
//   result.push_back(&socket_Ttable);
//   result.push_back(&socket_Ptable);
  result.push_back(&socket_volumes);
  result.push_back(&socket_states);
  result.push_back(&socket_normals);
  result.push_back(&socket_nstates);
  result.push_back(&socket_isOutward);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
ArcJetPhiST::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result = 
    ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_LorentzForce);
  result.push_back(&socket_Jx);
  result.push_back(&socket_Jy);
  result.push_back(&socket_Jz);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ArcJetPhiST::setup()
{
  using namespace COOLFluiD::Framework;

  ComputeSourceTermFVMCC::setup();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert (m_library.isNotNull());
  
  const CFuint nbCells = socket_volumes.getDataHandle().size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  socket_LorentzForce.getDataHandle().resize(nbCells*dim);
  socket_Jx.getDataHandle().resize(nbCells);
  socket_Jy.getDataHandle().resize(nbCells);
  socket_Jz.getDataHandle().resize(nbCells);
  
  DataHandle<CFreal> Jx = socket_Jx.getDataHandle();
  DataHandle<CFreal> Jy = socket_Jy.getDataHandle();
  DataHandle<CFreal> Jz = socket_Jz.getDataHandle();
  Jx.resize(nbCells);
  Jx = 0;
  Jy.resize(nbCells);
  Jy = 0;
  Jz.resize(nbCells);
  Jz = 0;
  
  m_normal.resize(DIM_3D, 0.);
  m_curlB.resize(DIM_3D, 0.);
  m_LorentzForce.resize(DIM_3D, 0.);
  m_B.resize(DIM_3D, 0.);
  m_u.resize(DIM_3D, 0.);
  m_uf.resize(DIM_3D, 0.);
  m_Bf.resize(DIM_3D, 0.);
  m_xf.resize(DIM_3D, 0.);
  m_UxB.resize(DIM_3D, 0.);
  m_J.resize(DIM_3D, 0.);
  m_gradPhi.resize(DIM_3D, 0.);
  
  if (m_Bfield.size() < DIM_3D) {
    m_Bfield.resize(DIM_3D, 0.);
    cf_assert(m_Bfield.size() == DIM_3D);
  }
  
  // store the B field in a RealVector for later use
  for (CFuint d = 0; d < DIM_3D; ++d) {
    m_B[d] = m_Bfield[d];
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void ArcJetPhiST::computeSource
(Framework::GeometricEntity *const element, RealVector& source, RealMatrix& jacobian)
{  
  const EquationSubSysDescriptor& eqSSD = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint nbEqs = eqSSD.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(nbEqs <= totalNbEqs);
  cf_assert(nbEqs > 0);
  
//   DataHandle<CFreal> cellID = socket_CellID.getDataHandle();  
//   std::cout <<"ArcJetPhiST::computeSource => cellID" << cellID[0] <<"\n";
  
  // this is needed for coupling
  if (eqSSD.getEqSS() == 1 || nbEqs == totalNbEqs) {
    CFLog(DEBUG_MIN, "ArcJetPhiST::computeSource() START\n");
    const CFuint elemID = element->getID();
    DataHandle<CFreal> LorentzForce = socket_LorentzForce.getDataHandle();
    DataHandle<CFreal> Jx = socket_Jx.getDataHandle();
    DataHandle<CFreal> Jy = socket_Jy.getDataHandle();
    DataHandle<CFreal> Jz = socket_Jz.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
    DataHandle<CFreal> normals = socket_normals.getDataHandle();
    DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
    DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
    const CFreal ovVolume = 1./volumes[elemID];
    
    SafePtr<EulerTerm> eulerTerm = PhysicalModelStack::getActive()->getImplementor()->
      getConvectiveTerm().d_castTo<EulerTerm>();
    
    // please note the reference &
    const State& currState = *element->getState(0);
    this->getMethodData().getUpdateVar()->computePhysicalData(currState, this->m_pdataArray);
    
    // semplification for now, no support for multi-species or turbulence
    typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> PTERM;
    
    for (CFuint d = 0; d < DIM_3D; ++d) {
      m_u[d] = this->m_pdataArray[PTERM::VX+d];
    }
    
    // make sure that uz = 0 in 2D
    if (dim == DIM_2D) {m_u[ZZ] = 0.;}
    
    // AL: change this for NEQ
    const CFuint phiID = totalNbEqs-1;
    CFuint TID = phiID-1;
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    if ((totalNbEqs == 7 && dim == DIM_3D) || (totalNbEqs == 6 && dim == DIM_2D)) {
      TID = totalNbEqs - 3;
      cf_assert(TID == phiID-2);
    }
    
    const CFuint pID = 0;
    
    // compute the gradients of B by applying Green Gauss in the cell volume
    m_gradPhi = 0.;
    CFreal sigmaUxB = 0.;
    
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
      m_uf = 0.;
      m_xf = 0.;
      CFreal avTdim = 0.;
      CFreal avpdim = 0.;
      CFreal* tVec = CFNULL;
      
      for (CFuint n = 0; n < nbFaceNodes; ++n) {
	const CFuint nodeID = face->getNode(n)->getLocalID();
	const RealVector& nodalState = nstates[nodeID];
	// consider all 3D components here even in 2D
	for (CFuint d = 0; d < DIM_3D; ++d) {
	  const CFreal nd = m_normal[d];
	  m_gradPhi[d] += nd*nodalState[phiID]*ovNbFaceNodes;
	  m_uf[d]      += nodalState[1+d]; // change here for NEQ
	}
	
	m_xf   += *face->getNode(n);
	avTdim += nodalState[TID]; 
	avpdim += nodalState[pID];
      }
      
      m_uf   *= ovNbFaceNodes;
      m_xf   *= ovNbFaceNodes;
      avTdim *= ovNbFaceNodes; 
      avpdim *= ovNbFaceNodes;
      
      // AL: compute m_Bf in function of the midFace point
      //     Rough assumptions here ...
      computeB(m_xf[XX] - m_electrodeX, m_Bf);
      
      avpdim = eulerTerm->getPressureFromState(avpdim);
      cf_assert(avpdim > 0.01);
      const CFreal sigma = m_library->sigma(avTdim, avpdim, tVec);
      //      const CFreal sigma = 490; 	//WATCH OUT: only for debugging purposes
      //       cout <<"ArcJetPhiST::computeSource => sigma =" << sigma <<"\n";
      //       cout <<"ArcJetPhiST::computeSource => avTdim =" << avTdim <<"\n";
      //       cout <<"ArcJetPhiST::computeSource => avpdim =" << avpdim <<"\n";
      MathFunctions::crossProd(m_uf, m_Bf, m_UxB);
      sigmaUxB += sigma*MathFunctions::innerProd(m_UxB, m_normal);
    }
    
    m_gradPhi *= ovVolume;
    
    //cout << "m_gradPhi[OLD] = " << m_gradPhi << endl; 
    
    // AL: here we overwrite m_gradPhi with the value coming from the LS reconstruction
    const CFuint gradPhiID = elemID*totalNbEqs + phiID;
    m_gradPhi[XX] = this->m_ux[gradPhiID];
    m_gradPhi[YY] = this->m_uy[gradPhiID];
    m_gradPhi[ZZ] = this->m_uz[gradPhiID];
    
    //cout << "m_gradPhi[NEW] = " << m_gradPhi << endl << endl; 
    
    // source term contribution to 
    source[phiID] = -sigmaUxB*ovVolume; // the whole array will be multiplied by volume later
    
    //
    // source terms for the momentum and energy equations
    //
    
    CFreal pdim = eulerTerm->getPressureFromState(this->m_pdataArray[PTERM::P]);
    cf_assert(pdim > 0.01);
    CFreal Tdim = this->m_pdataArray[PTERM::T];
    cf_assert(Tdim > 0.01);
    CFreal* tVec = CFNULL; // AL: this will need to change for multi-temperature
    const CFreal sigma = m_library->sigma(Tdim, pdim, tVec);
    //    const CFreal sigma = 490; 	//WATCH OUT: only for debugging purposes
    
    
    
    computeB(element->getState(0)->getCoordinates()[XX] - m_electrodeX, m_B); 
    MathFunctions::crossProd(m_u, m_B, m_UxB);
    m_J = sigma*(-m_gradPhi + m_UxB); 
    
    // DEBUGGING: limiter for the current
    //CFreal Jmax = 1.15*10/(3.14159*0.04*0.04); //Max Value for 10A
    //CFreal Jmax   = 3000;  
    //CFreal Jmin   = -Jmax; //Min Value for 10A
    //CFreal xPos   = element->getState(0)->getCoordinates()[XX];
    //CFreal xLim   = 0.09;
    //if (xPos > xLim){
      //for (CFuint d = 0; d < DIM_3D; ++d) {
	//if(m_J[d] < Jmin) {m_J[d] = Jmin;}
	//if(m_J[d] > Jmax) {m_J[d] = Jmax;}
      //}
    //}
    // here fill in Lorentz force
    MathFunctions::crossProd(m_J, m_B, m_LorentzForce);
    
    const CFuint estart = elemID*dim;
    for (CFuint d = 0; d < dim; ++d) {
      source[m_velIDs[d]] = m_LorentzForce[d];
      
      // store the Lorentz force in a data socket to be able to plot it
      LorentzForce[estart+d] = m_LorentzForce[d];
    }
    
    //cout<<"ArcJetST before entering";
    
    const CFreal J2 = MathFunctions::innerProd(m_J, m_J)/sigma; //Joule heating term
    if(m_disableCurr == false){
      source[TID] = MathFunctions::innerProd(m_u, m_LorentzForce) + J2;
      //cout<<"Current enable";
    }
    else{
      source[TID] = 0;
      //cout<<"Current disabled";
    }
    
    source *= volumes[elemID];
    CFLog(DEBUG_MIN, "ArcJetPhiST::computeSource() END\n");
    
    // writing the data sockets to plot the electric current
    if (!getMethodData().isPerturb()) { // Condition to avoid writing when the source is perturbed
      if(currState.isParUpdatable()) { // Condition to write only once the partition cells
      //cout << "Storing the current" <<!getMethodData().isPerturb()<<"\n" ;
	Jx[elemID] = m_J[0]; //source[TID];/////This is a test
	Jy[elemID] = m_J[1];  //sigma; ///This is a test
	Jz[elemID] = m_J[2]; //m_J[2]; ///This is a test
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet
    
  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

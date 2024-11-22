#include "FiniteVolumePoissonNEQ/PoissonCNEQST.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalConsts.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/MultiScalarTerm.hh"
#include "FiniteVolumePoissonNEQ/FiniteVolumePoissonNEQ.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/BaseDataSocketSink.hh"
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

MethodStrategyProvider<PoissonCNEQST, 
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumePoissonNEQModule> 
PoissonCNEQSTFVMCCProvider("PoissonCNEQST");

//////////////////////////////////////////////////////////////////////////////

PoissonCNEQST::PoissonCNEQST(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_LorentzForce("LorentzForce"),
  socket_Bfield("Bfield"), // new
  socket_BfieldFaces("BfieldFaces"), //new
  socket_J("J"),  
  socket_sigma("sigma"),
  socket_E("E"), 
  m_library(CFNULL),
  m_normal(),
  m_curlB(),
  m_LorentzForce(),
  m_B(),
  m_u(),
  m_uf(),
  m_Bf(),
  m_J(),
  m_uB(),
  m_uB_f(),
  m_E(),
  m_gradB()
{
}

//////////////////////////////////////////////////////////////////////////////

PoissonCNEQST::~PoissonCNEQST()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
PoissonCNEQST::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result = 
    ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_LorentzForce);
   result.push_back(&socket_sigma);
  result.push_back(&socket_J);
  result.push_back(&socket_E);
  return result;
}
///* // check line 42
std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
PoissonCNEQST::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
    ComputeSourceTermFVMCC::needsSockets();
  result.push_back(&socket_Bfield);
  result.push_back(&socket_BfieldFaces);
  return result;
}
//*/
//////////////////////////////////////////////////////////////////////////////

void PoissonCNEQST::setup()
{
  using namespace COOLFluiD::Framework;

  ComputeSourceTermFVMCC::setup();
  // init a size of dim3*nbcells for source
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert (m_library.isNotNull());
  
  const CFuint nbCells = socket_volumes.getDataHandle().size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  socket_LorentzForce.getDataHandle().resize(nbCells*dim);
  
  socket_sigma.getDataHandle().resize(nbCells);
 
  socket_J.getDataHandle().resize(nbCells*3);
  
  socket_E.getDataHandle().resize(nbCells*3);

  m_curlB.resize(DIM_3D, 0.);
  m_LorentzForce.resize(DIM_3D, 0.);
  m_uB.resize(DIM_3D, 0.);
  m_B.resize(DIM_3D, 0.);
  m_u.resize(DIM_3D, 0.);
  m_J.resize(DIM_3D, 0.);
  m_E.resize(DIM_3D, 0.);
  m_normal.resize(DIM_3D, 0.);
  m_gradB.resize(DIM_3D, DIM_3D, 0.);
  m_uf.resize(DIM_3D, 0.);
  m_Bf.resize(DIM_3D, 0.);
  m_uB_f.resize(DIM_3D, 0.);
}

//////////////////////////////////////////////////////////////////////////////
      
void PoissonCNEQST::computeSource
(Framework::GeometricEntity *const element, RealVector& source, RealMatrix& jacobian)
{  
  const EquationSubSysDescriptor& eqSSD = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint nbEqs = eqSSD.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // this is needed for coupling
  if (eqSSD.getEqSS() == 1 || nbEqs == totalNbEqs) {
      CFLogDebugMin( "PoissonCNEQST::computeSource()" << "\n");
      const CFuint dim = PhysicalModelStack::getActive()->getDim();
      const CFuint elemID = element->getID();

      DataHandle<CFreal> E_ = socket_E.getDataHandle();
      DataHandle<CFreal> sigma_ = socket_sigma.getDataHandle();
      DataHandle<CFreal> J_ = socket_J.getDataHandle();
      DataHandle<CFreal> LorentzForce = socket_LorentzForce.getDataHandle();
      DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
      DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
      DataHandle<RealVector> nstates = this->_sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
      DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();
      DataHandle<CFreal> Bfield_cc = socket_Bfield.getDataHandle(); // cell center
      DataHandle<CFreal> Bfield_fc = socket_BfieldFaces.getDataHandle(); // face center
      
      /**************Cell center calculation *************************/
      const CFreal ovVolume = 1./volumes[elemID];
      SafePtr<MultiScalarTerm<EulerTerm > > eulerTerm = PhysicalModelStack::getActive()->getImplementor()->
      getConvectiveTerm().d_castTo<MultiScalarTerm<EulerTerm > >();
    
      State *const currState = element->getState(0);
      typedef EulerTerm PTERM;
      
      this->getMethodData().getUpdateVar()->computePhysicalData(*currState, this->m_pdataArray);

      const CFuint nbSpecies = eulerTerm->getNbScalarVars(0);
      const CFuint nbTv = eulerTerm->getNbScalarVars(1); // number of Tv(s) and/or Te 
      
      // computation of cell-centered value of electrical conductivity
      CFreal pdim = eulerTerm->getPressureFromState(this->m_pdataArray[PTERM::P]); //this->m_pdataArray[PTERM::P];

      cf_assert(pdim > 0.01);
      CFreal Tdim = this->m_pdataArray[PTERM::T];

      cf_assert(Tdim > 0.01);
      CFreal* tVec = CFNULL; //&(*currState)[nbSpecies+dim+1]; // array pointing to the Tv(s)/Te // commented to run CNEQ
      
      const CFuint startID = elemID*totalNbEqs;
      const CFuint idxB = elemID*DIM_3D;

      const CFuint uID = nbSpecies;     // ID of x-momentum equation
      const CFuint vID = nbSpecies+1;   // ID of y-momentum equation
      const CFuint eID = nbSpecies+dim; // ID of total energy equation

      const CFuint phiID   = totalNbEqs-1; // ID of phi, which is the last equation we solve

      // velocities and B (3 components needed for cross product, even in 2D)
      for (CFuint d = 0; d < DIM_3D; ++d) {
        m_u[d] = this->m_pdataArray[PTERM::VX+d];
        m_B[d] = Bfield_cc[idxB+d]; // this is reading through socket, is the id ok?
      }
    
      // make sure that uz = 0 in 2D
      if (dim == DIM_2D) {
        m_u[ZZ] = 0.0;
      }

      // Vatsalya
    
      const CFuint gradPhiID = startID + phiID;

      m_E[XX] = -this->m_ux[gradPhiID]; // electric field = grad phi
      m_E[YY] = -this->m_uy[gradPhiID]; // gradients of phi computed from cell-centered gradients (ux, uy, uz) used for Least Square Reconstruction
      m_E[ZZ] = -this->m_uz[gradPhiID];

      // Plato Devato 
      const CFreal sigma_vol = m_library->sigma(Tdim, pdim, tVec);
      
      CFLog(DEBUG_MAX, "PoissonCNEQST::computeSource() => sigma = " << sigma_vol << "\n");
      MathFunctions::crossProd(m_u, m_B, m_uB); // u x B vector of 3 size is done
      m_J = sigma_vol*( m_E + m_uB ); // current density using ohm's law

      // Electromagnetic Energy Deposition      
      const CFreal Ed = MathFunctions::innerProd(m_J, m_E);  
      const CFreal gamma_magnetic = 1.;
      const CFreal magn_coeff = gamma_magnetic/sigma_vol;

      // Joule Heting // not required in CNEQ
     // const CFreal JH = magn_coeff * MathFunctions::innerProd(m_J, m_J);

      // Lorentz Force Calculation
      MathFunctions::crossProd(m_J, m_B, m_LorentzForce);  

      // allocate right source term 
      const CFuint estart = elemID*dim;

      sigma_[elemID] = sigma_vol;

      for (CFuint d = 0; d < dim; ++d) {   

        source[uID+d] = m_LorentzForce[d]; //Lorentz force for momemntum equations
        LorentzForce[estart+d] = m_LorentzForce[d]; // store the Lorentz force in a data socket to be able to plot it. Can we store other source terms?? Will try once it works
        E_[estart+d] = m_E[d];
        J_[estart+d] = m_J[d];
       
      }
   
      source[eID]   = Ed ;//+ JH;     // source term for Total energy equation
       
  ////    source[eID+1] = JH;     // source term for electron energy equation, only for TCNEQ

      //---------------------- computation of source term for poisson equation at faces-------------------------//

      CFreal sigmaUxB = 0.;
 
      const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
      const CFuint nbFaces = faces.size();
      for (CFuint i = 0; i < nbFaces; ++i) {

        const GeometricEntity *const face = element->getNeighborGeo(i);
        const CFuint faceID = face->getID();
        const CFuint startID2 = faceID*PhysicalModelStack::getActive()->getDim();
        const CFuint nbFaceNodes = face->nbNodes();
        const CFreal ovNbFaceNodes = 1./(CFreal)nbFaceNodes;

        const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;

        //------------------ compute u_f and B_f at each face--------------------------//
        m_uf = 0. ;
        m_Bf = 0. ;

        for (CFuint d = 0; d < DIM_3D; ++d) {
          m_Bf[d] = Bfield_fc[startID2+d];
          m_normal[d] = normals[startID2+d]*factor;
        }

        CFreal Tdim2 = 0.;
        CFreal pdim2 =0. ;
        const CFreal rho_f = 0.0;
        CFreal* tVec3 = CFNULL; 

        for (CFuint n = 0; n < nbFaceNodes; ++n) {

          const CFuint nodeID = face->getNode(n)->getLocalID();
          const RealVector& nodalState = nstates[nodeID];
          // consider all 3D components here even in 2D
          
          for (CFuint d = 0; d < DIM_3D; ++d) {
            m_uf[d]     += nodalState[uID+d]; 
          }
           Tdim2 += nodalState[eID];    // total temperature at face
          CFreal rho_f = 0.; // can this be replaced by diff varset etc?? This can be calculated earlyon and stored, instead of calulcating again
          for (CFuint i = 0; i < nbSpecies; ++i) {
              rho_f += nodalState[i];
          }
           pdim2 += m_library->pressure(rho_f, Tdim2, tVec3); 
        }
        
        m_uf   *= ovNbFaceNodes;
        Tdim2  *= ovNbFaceNodes;
        pdim2  *= ovNbFaceNodes;
        const CFreal sigma_f = m_library->sigma(Tdim2, pdim2, tVec3); // sigma at face
        MathFunctions::crossProd(m_uf, m_Bf, m_uB_f);
        sigmaUxB += sigma_f*MathFunctions::innerProd(m_uB_f, m_normal);
      }
      
      source[phiID] = -sigmaUxB*ovVolume; 
      /*
      if(elemID ==7){
         cout<<"source[phiID] = "<<source[phiID]<< " volumes[elemID] = "<<volumes[elemID] << " pr = " <<source[phiID]*volumes[elemID] << "\n";
       }
      // only at the end we multiply for the cell volume
      */
      source *= volumes[elemID];//*0.000001;
      /*
      if(elemID ==7){
         cout<<"F_source[phiID] = "<< source[phiID]   << " source[T] = " << source[eID] << "\n";
         for(CFuint d = 0; d < dim; ++d){
            cout<< " d = " << d << " source[uID+d] = " << source[uID+d] << " m_LorentzForce[d] = "<<m_LorentzForce[d] <<"\n";
         }
       }
       */
      
     
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumePoissonNEQ
    
  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

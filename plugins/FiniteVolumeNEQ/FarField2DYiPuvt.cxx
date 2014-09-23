#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNEQ/FarField2DYiPuvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FarField2DYiPuvt, CellCenterFVMData, FiniteVolumeNEQModule>
farfieldEuler2DYiPuvtFVMCCProvider("FarFieldEuler2DYiPuvt");


//////////////////////////////////////////////////////////////////////////////

void FarField2DYiPuvt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Pout","pressure compressible"); 
  options.addConfigOption< std::vector<CFreal> > ("VTTv", "Inlet Velocity components, temperatures");
}

//////////////////////////////////////////////////////////////////////////////

FarField2DYiPuvt::FarField2DYiPuvt(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _library(CFNULL),
 _dataInnerState()

{
  addConfigOptionsTo(this);

  m_vTTv = std::vector<CFreal>();
  setParameter("VTTv",&m_vTTv);

  m_pressure = 0.0;
  setParameter("Pout",&m_pressure);
}

//////////////////////////////////////////////////////////////////////////////

FarField2DYiPuvt::~FarField2DYiPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void FarField2DYiPuvt::setGhostState(GeometricEntity *const face)
{
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(innerState, _dataInnerState);
  
  SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());
  
  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal v = _dataInnerState[EulerTerm::VY];
  const CFreal vn = u*nx + v*ny; 
  const CFreal aInnerState = _dataInnerState[EulerTerm::A];
  const CFreal machInner = vn / aInnerState;

if ((machInner < 1.0) && (machInner >= 0.0))  {
const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbSpecies = library->getNbSpecies();
    
  // extrapolate the partial pressure from inside
  const CFuint TID = nbSpecies + dim;
  
  State& ghostState = *face->getState(RIGHT);
  for (CFuint i = 0; i < m_vTTv.size(); ++i) {
    ghostState[nbSpecies+i] = 2.*m_vTTv[i] - innerState[nbSpecies+i];
  }
  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    ghostState[i] = innerState[i]*innerState[TID]/ghostState[TID];
  }
  
 }

if ((machInner > -1.0) && (machInner < 0.0)) {

  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
 const CFreal Rgas = library->getRgas();
  const CFuint nbSpecies = _varSet->getModel()->getNbScalarVars(0);
  const CFuint Tghost =  _dataInnerState[EulerTerm::T];
 

 library-> getMolarMasses(_mm);
  
      const CFuint firstSpecies = _varSet->getModel()->getFirstScalarVar(0);
      CFreal mass = 0;
      for (CFuint k = 0; k < nbSpecies; k++) {
        _Yi[k] = _dataInnerState[firstSpecies + k];
         mass += _Yi[k]/_mm[k];
        }
   
      //cout << "mass" << "" <<  mass << endl;
  
   const CFuint sum = mass*Rgas*Tghost;

  // first variable is the 
   const CFreal pGhost    = 2.*m_pressure  - pInnerState;
   const CFreal Rho_Ghost = pGhost /sum; 
               
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal yIGhost = _Yi[i];
    ghostState[i] = yIGhost*Rho_Ghost;
  }

  // following variables are exptrapolated
    ghostState[nbSpecies]    = _dataInnerState[EulerTerm::VX];
    ghostState[nbSpecies+1]  = _dataInnerState[EulerTerm::VY];
    ghostState[nbSpecies+2]  = _dataInnerState[EulerTerm::T];
    ghostState[nbSpecies+3]  = innerState[nbSpecies+3]; 
 }




}

//////////////////////////////////////////////////////////////////////////////

void FarField2DYiPuvt::setup()
{
   //cout << "mm" <<  endl;
   FVMCC_BC::setup();
  _varSet = getMethodData().getUpdateVar().d_castTo<MultiScalarVarSet<Euler2DVarSet> >();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _mm.resize(_varSet->getModel()->getNbScalarVars(0));
  _Yi.resize(_varSet->getModel()->getNbScalarVars(0));
   //m_pressure/=_varSet->getModel()->getPressRef();

 _library = Framework::PhysicalModelStack::getActive()->
    getImplementor()->getPhysicalPropertyLibrary
    <Framework::PhysicalChemicalLibrary>();
  
  m_ysIn.resize(_library->getNbSpecies());
  cf_assert(m_ysIn.size() > 0); 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

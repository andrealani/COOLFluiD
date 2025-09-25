#include "FluxReconstructionMHD/ManufacturedMHDSourceTerm.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "Framework/SubSystemStatus.hh"

#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathConsts.hh"

#include "MHD/MHDTerm.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

////////////////////////////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ManufacturedMHDSourceTerm, FluxReconstructionSolverData, FluxReconstructionMHDModule>
ManufacturedMHDSourceTermFRProvider("ManufacturedMHDSourceTerm");

///////////////////////////////////////////////////////////////////////////////////////////////////

void ManufacturedMHDSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("AddUpdateCoeff","Add the ST time step restriction.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void ManufacturedMHDSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

ManufacturedMHDSourceTerm::ManufacturedMHDSourceTerm(const std::string& name) :
  StdSourceTerm(name),
  socket_updateCoeff("updateCoeff"),
  m_order()
{ 
  addConfigOptionsTo(this);
  
  m_addUpdateCoeff = false;
  setParameter("AddUpdateCoeff",&m_addUpdateCoeff);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

ManufacturedMHDSourceTerm::~ManufacturedMHDSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void ManufacturedMHDSourceTerm::setup()
{
  StdSourceTerm::setup();
    
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();
  
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);
  

  const CFuint nbStates = updateCoeff.size();

}

//////////////////////////////////////////////////////////////////////////////////////////////////

void ManufacturedMHDSourceTerm::unsetup()
{
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  ManufacturedMHDSourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = StdSourceTerm::providesSockets();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
ManufacturedMHDSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedMHDSourceTerm::addSourceTerm(RealVector& resUpdates)
{       
  CFLog(VERBOSE, "ManufacturedMHDSourceTerm::addSourceTerm() => START\n");
  
  // set gradients
  const CFuint nbrStates = m_cellStates->size();
  const CFreal kappa = 0.017;
   
  
  for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
  {   
    const CFreal x = ((*m_cellStates)[iSol]->getCoordinates())[XX];
    const CFreal y = ((*m_cellStates)[iSol]->getCoordinates())[YY];
    const CFreal z = ((*m_cellStates)[iSol]->getCoordinates())[ZZ];   
    const CFreal r = sqrt(x*x+y*y+z*z);
    // density
    resUpdates[m_nbrEqs*iSol + 0] = 0.0;

    resUpdates[m_nbrEqs*iSol + 1] = 0.5*x*pow(r,-5./2.)*((1./r)-(5./(r*r))-kappa*z);

    resUpdates[m_nbrEqs*iSol + 2] = 0.5*y*pow(r,-5./2.)*((1./r)-(5./(r*r))-kappa*z);

    resUpdates[m_nbrEqs*iSol + 3] = 0.5*z*pow(r,-5./2.)*((1./r)-(5./(r*r))-kappa*z) + 5./2.*kappa*pow(r,-1./2.)*(1.+kappa*r*z) + kappa*pow(r,-1./2.);

    resUpdates[m_nbrEqs*iSol + 4] = 0.0;
    resUpdates[m_nbrEqs*iSol + 5] = 0.0;
    resUpdates[m_nbrEqs*iSol + 6] = 0.0;

    resUpdates[m_nbrEqs*iSol + 7] = 0.5/(r*r) + kappa*z*(3.5/r+2.*kappa*z) + 0.5*pow(kappa*r,2)*(7+5*kappa*r*z);


    resUpdates[m_nbrEqs*iSol + 8] = 0.0;

  }
  
  CFLog(VERBOSE, "ManufacturedMHDSourceTerm::addSourceTerm() => END\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////

/*void  ManufacturedMHDSourceTerm::getSToStateJacobian(const CFuint iState)
{
 
}*/

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

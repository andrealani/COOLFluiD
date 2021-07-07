#include "Common/PE.hh"
#include "Common/BadValueException.hh"
#include "Common/Stopwatch.hh"

#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PhysicalModel.hh"

#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolumeNEQ/ComputeFieldFromPotentialNEQ.hh"

#include <cmath>

/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeFieldFromPotentialNEQ,
		      CellCenterFVMData, FiniteVolumeNEQModule>
ComputeFieldFromPotentialNEQProvider("ComputeFieldFromPotentialNEQ");

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotentialNEQ::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >
    ("OtherNamespace", "Name of the other namespace (providing the potential).");
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeFieldFromPotentialNEQ::ComputeFieldFromPotentialNEQ(const std::string& name) :
  CellCenterFVMCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_pastStates("pastStates"),
  socket_nodes("nodes"),
  socket_otherUX("uX"),
  socket_otherUY("uY"),
  socket_otherUZ("uZ"),
  socket_Bfield("Bfield"),
  socket_otherStates("states"),
  m_applyProcessing(true)
{
  addConfigOptionsTo(this);
  
  m_otherNamespace = "";
  setParameter("OtherNamespace", &m_otherNamespace);
}

//////////////////////////////////////////////////////////////////////////////

ComputeFieldFromPotentialNEQ::~ComputeFieldFromPotentialNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeFieldFromPotentialNEQ::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_nodes);
  result.push_back(&socket_otherUX);
  result.push_back(&socket_otherUY);
  result.push_back(&socket_otherUZ);
  result.push_back(&socket_otherStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeFieldFromPotentialNEQ::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_Bfield);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotentialNEQ::setup()
{
  CFAUTOTRACE;

  CellCenterFVMCom::setup();
  
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();
  socket_Bfield.getDataHandle().resize(DIM_3D*nbStates);
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotentialNEQ::configure ( Config::ConfigArgs& args )
{
  CellCenterFVMCom::configure(args);

  cf_assert(m_otherNamespace != "");
  CFLog(VERBOSE, "ComputeFieldFromPotentialNEQ::configure() => m_otherNamespace = " <<
	m_otherNamespace << "\n");
  socket_otherUX.setDataSocketNamespace(m_otherNamespace);
  socket_otherUY.setDataSocketNamespace(m_otherNamespace);
  socket_otherUZ.setDataSocketNamespace(m_otherNamespace);
  socket_otherStates.setDataSocketNamespace(m_otherNamespace);
}
      
//////////////////////////////////////////////////////////////////////////////
 
void ComputeFieldFromPotentialNEQ::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "ComputeFieldFromPotentialNEQ::execute() => START\n");

  if (SubSystemStatusStack::getActive()->getNbIter() >= 1 && m_applyProcessing) {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
    Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getNamespace(m_otherNamespace);
    Common::SafePtr<SubSystemStatus> otherSubSystemStatus =
      SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

    DataHandle<CFreal> ux = socket_otherUX.getDataHandle();
    DataHandle<CFreal> uy = socket_otherUY.getDataHandle();
    DataHandle<CFreal> uz = socket_otherUZ.getDataHandle();
    DataHandle<CFreal> Bfield = socket_Bfield.getDataHandle();
    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
    
    const CFuint nbStates = states.size();
    CFLog(INFO, "ComputeFieldFromPotentialNEQ::execute() => transferring field\n");
    
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      const CFuint startIdx = iState*DIM_3D;
      Bfield[startIdx]   = ux[iState];
      Bfield[startIdx+1] = uy[iState];
      if (dim == DIM_3D) {
	Bfield[startIdx+2] = uz[iState];
      }
    }
    
    // AL: this has to be changed when doing unsteady cases!!!!!!
    // We need to implement a generic logic to apply the processing only when Poisson is solved
    m_applyProcessing = false;
    if (SubSystemStatusStack::getActive()->getNbIter()%this->getProcessRate()==0) {
      m_applyProcessing = true;
    }
    
    // update the pastStates
    /*DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();
    // Set Initial States to current states
    if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
      for(CFuint i = 0; i < states.size(); ++i) {
      *(pastStates[i]) = *(states[i]);
      }
      }*/
  } 
  
  CFLog(VERBOSE, "ComputeFieldFromPotentialNEQ::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////
 
void ComputeFieldFromPotentialNEQ::unsetup()
{
  CFAUTOTRACE;
  
  CellCenterFVMCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////


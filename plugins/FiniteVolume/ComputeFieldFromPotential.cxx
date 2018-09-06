#include "Common/PE.hh"
#include "Common/BadValueException.hh"

#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PhysicalModel.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/ComputeFieldFromPotential.hh"

/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeFieldFromPotential,
		      CellCenterFVMData, FiniteVolumeModule>
ComputeFieldFromPotentialProvider("ComputeFieldFromPotential");

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >
    ("VariableIDs", "IDs of the variable to be assigned to the newly computed field.");
  options.addConfigOption< string >
    ("OtherNamespace", "Name of the other namespace (providing the potential).");
  options.addConfigOption< CFreal >
    ("InterRadius",
     "Radius corresponding to the internal boundary between donor and current grids (<= 0 assumes one mesh).");
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeFieldFromPotential::ComputeFieldFromPotential(const std::string& name) :
  CellCenterFVMCom(name),
  socket_states("states"),
  socket_otherUX("uX"),
  socket_otherUY("uY"),
  socket_otherUZ("uZ"),
  socket_otherStates("states")
{
  addConfigOptionsTo(this);
  
  m_variableIDs = vector<CFuint>();
  setParameter("VariableIDs",&m_variableIDs);
  
  m_otherNamespace = "";
  setParameter("OtherNamespace", &m_otherNamespace);

  m_interRadius = -1.;
  setParameter("InterRadius", &m_interRadius);
}

//////////////////////////////////////////////////////////////////////////////

ComputeFieldFromPotential::~ComputeFieldFromPotential()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeFieldFromPotential::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_otherUX);
  result.push_back(&socket_otherUY);
  result.push_back(&socket_otherUZ);
  result.push_back(&socket_otherStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::setup()
{
  CFAUTOTRACE;

  CellCenterFVMCom::setup();

  if (m_variableIDs.size() == 0) {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    m_variableIDs.resize(dim);
    for (CFuint i = 0; i < dim; ++i) {
      m_variableIDs[i] = i;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::configure ( Config::ConfigArgs& args )
{
  CellCenterFVMCom::configure(args);

  cf_assert(m_otherNamespace != "");
  CFLog(VERBOSE, "ComputeFieldFromPotential::configure() => m_otherNamespace = " <<
	m_otherNamespace << "\n");
  socket_otherUX.setDataSocketNamespace(m_otherNamespace);
  socket_otherUY.setDataSocketNamespace(m_otherNamespace);
  socket_otherUZ.setDataSocketNamespace(m_otherNamespace);
  socket_otherStates.setDataSocketNamespace(m_otherNamespace);
}
      
//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::execute()
{
  CFAUTOTRACE;

  CFLog(INFO, "ComputeFieldFromPotential::execute() => START\n");
  
  if (SubSystemStatusStack::getActive()->getNbIter() >= 1) {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
    Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getNamespace(m_otherNamespace);
    Common::SafePtr<SubSystemStatus> otherSubSystemStatus =
      SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

    DataHandle<CFreal> ux = socket_otherUX.getDataHandle();
    DataHandle<CFreal> uy = socket_otherUY.getDataHandle();
    DataHandle<CFreal> uz = socket_otherUZ.getDataHandle();
    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
    cf_assert(ux.size() == states.size());
    cf_assert(uy.size() == states.size());
    if (dim == DIM_3D) {
      cf_assert(uz.size() == states.size());
    }
    
    const CFuint nbStates = states.size();
    cf_assert(dim >= DIM_2D);
    cf_assert(m_variableIDs.size() >= 2);
    const CFuint xVar = m_variableIDs[0];
    cf_assert(xVar < nbEqs);
    const CFuint yVar = m_variableIDs[1];
    cf_assert(yVar < nbEqs);
    const CFuint zVar = (dim == DIM_3D) ? m_variableIDs[2] : 0;
    cf_assert(zVar < nbEqs);

    if (m_interRadius <= 0.) {
      for (CFuint iState = 0; iState < nbStates; ++iState) {
	(*states[iState])[xVar] = ux[iState];
	(*states[iState])[yVar] = uy[iState];
	if (dim == DIM_3D) {
	  (*states[iState])[zVar] = uz[iState];
	}
      }
    }
    else {
      // those are the states correspondng to the smaller mesh
      DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle();
      
      // here follows algorithm to copy the B values for states
      // with R < m_interRadius and R > m_interRadius the algorithm will behave differently

      // State* s1 = states[cellID1];
      // Node& coord1 = s1->getCoordinates();
      // const CFreal radius1 = coord1.norm2();
      
      // State* s2 = otherStates[cellID2];
      // Node& coord2 = s2->getCoordinates();
      // const CFreal radius2 = coord2.norm2();
    }
  } 
  CFLog(INFO, "ComputeFieldFromPotential::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::unsetup()
{
  CFAUTOTRACE;
  
  CellCenterFVMCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////


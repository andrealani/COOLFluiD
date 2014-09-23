#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSVector.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/CopyFromTrsBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< CopyFromTrsBC,UFEMSolverData,UFEMPlugin > CopyFromTrsBCProvider("CopyFromTrsBC");

//////////////////////////////////////////////////////////////////////////////

void CopyFromTrsBC::defineConfigOptions(Config::OptionList& options)
{
   CFAUTOTRACE;

   options.addConfigOption< std::string >( "CopyFromTrs","From which trs to copy from (default is empty)." );
   options.addConfigOption< bool >( "CopyMatches","If number of states and states order matches (default is true)" );
}

//////////////////////////////////////////////////////////////////////////////

CopyFromTrsBC::CopyFromTrsBC(const std::string& name) :
  DirichletBC(name),
  m_copyFromTrs("")
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_copyMatches=true;
  setParameter( "CopyFromTrs",  &m_copyFromTrs);
  setParameter( "CopyMatches",  &m_copyMatches);

}

//////////////////////////////////////////////////////////////////////////////

void CopyFromTrsBC::setup()
{
  CFAUTOTRACE;
  DirichletBC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void CopyFromTrsBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  DirichletBC::configure(args);

  CFLog(INFO, getClassName() << ": CopyFromTrs: "      << m_copyFromTrs  << "\n");
  CFLog(INFO, getClassName() << ": CopyMatches: "      << m_copyMatches  << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void CopyFromTrsBC::executeOnTrs()
{
  CFAUTOTRACE;

  // when trs to trs
  if (m_copyFromTrs!=""){

    // states to modify
    DataHandle< State*, GLOBAL > states = socket_states.getDataHandle();

    // finding both trs
    vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
    SafePtr<TopologicalRegionSet> copyfromTrs=0;
    SafePtr<TopologicalRegionSet> currentTrs=getCurrentTRS();
    for (CFuint i = 0; i < trs.size(); ++i) 
      if (trs[i]->getName()==m_copyFromTrs)
        copyfromTrs=trs[i];
    if (copyfromTrs==0) Common::NoSuchValueException(FromHere(),"Trs not found with name: " + m_copyFromTrs);
    if (copyfromTrs->getStatesInTrs()!=currentTrs->getStatesInTrs()) Common::NoSuchValueException(FromHere(),"Current and copyfrom Trs differs in number of states");

    // number of equations and isupdated socket
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    DataHandle< bool > isUpdated = socket_isUpdated.getDataHandle();

    // loop on states and doing the copy
    Common::SafePtr< std::vector< CFuint > > currentStates = currentTrs->getStatesInTrs();
    Common::SafePtr< std::vector< CFuint > > copyfromStates = copyfromTrs->getStatesInTrs();
    std::vector< CFuint >::iterator icurrent;
    std::vector< CFuint >::iterator icopyfrom;
    CFuint istate=0;
    for (icurrent = currentStates->begin(),
         icopyfrom = copyfromStates->begin();
         icurrent != currentStates->end(),
         icopyfrom != copyfromStates->end();
         ++icurrent,
         ++icopyfrom,
         istate++){
      const CFuint nLocalID= *icurrent;
      State *currentState = states[*icurrent];
      const State *copyfromState = states[*icopyfrom];
      for (CFuint ieq=0; ieq<nbEqs; ieq++) if ((!isUpdated[nLocalID*nbEqs+ieq])&&(m_applyFlags[istate][ieq])) {
        (*currentState)[ieq]=(*copyfromState)[ieq];
      }
    }

  } else {

    Common::BadValueException(FromHere(),"Empty string in CopyFromTrs. Is it specified in the CFcase for this (" + getCurrentTRS()->getName() + ") boundary?");

  }

  DirichletBC::executeOnTrs();
}

//////////////////////////////////////////////////////////////////////////////

void CopyFromTrsBC::computeStateValuesCopyFromTrsBC(const Framework::State* currState)
{
  computeStateValuesBaseBC(currState);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD


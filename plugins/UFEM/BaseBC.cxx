#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSVector.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/BaseBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BaseBC,UFEMSolverData,UFEMPlugin > BaseBCProvider("BaseBC");

//////////////////////////////////////////////////////////////////////////////

void BaseBC::defineConfigOptions(Config::OptionList& options)
{
   CFAUTOTRACE;

   options.addConfigOption< std::vector< std::string > >( "Vars","Definition of the Variables (required)." );
   options.addConfigOption< std::vector< CFuint > >( "ApplyEqs","Apply the BC only to specified equations zero-based indexed (default all equations)." );
   options.addConfigOption< std::vector< std::string > >( "Def","Definition of the Functions (required)." );
}

//////////////////////////////////////////////////////////////////////////////

BaseBC::BaseBC(const std::string& name) :
  UFEMSolverCom(name),
  socket_interStates("interStates"),
  socket_states("states"),
  m_functions(0),
  m_vars(0),
  m_applyEqs(0)
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  setParameter( "Def",          &m_functions);
  setParameter( "Vars",         &m_vars);
  setParameter( "ApplyEqs",     &m_applyEqs);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  UFEMSolverCom::configure(args);

  // giving defaults for "Def" and "Vars"
  std::vector< std::string > varnames=getMethodData().getUpdateVar()->getVarNames();
  if (m_vars.size()==0) {
    Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getNamespace( getMethodData().getNamespace() );
    Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
    const CFuint nbdim = physModel->getDim();
    if ( nbdim > DIM_0D ) m_vars.push_back ( "x" );
    if ( nbdim > DIM_1D ) m_vars.push_back ( "y" );
    if ( nbdim > DIM_2D ) m_vars.push_back ( "z" );
    cf_assert ( nbdim <= DIM_3D );
    m_vars.push_back ( "t" );
    for (CFuint i=0; i<varnames.size(); i++) m_vars.push_back(varnames[i]);
  }
  if (m_functions.size()==0) {
    for (CFuint i=0; i<varnames.size(); i++) m_functions.push_back(varnames[i]);
  }

  // printing feedback
  std::ostringstream sstr;
  sstr.str(""); for (CFuint i=0; i<m_functions.size(); i++) sstr << m_functions[i] << " ";
  CFLog(INFO, getClassName() << ": Def: "              << sstr.str()     << "\n");
  sstr.str(""); for (CFuint i=0; i<m_vars.size(); i++)      sstr << m_vars[i]      << " ";
  CFLog(INFO, getClassName() << ": Vars: "             << sstr.str()     << "\n");
  sstr.str(""); for (CFuint i=0; i<m_applyEqs.size(); i++)  sstr << m_applyEqs[i]  << " ";
  CFLog(INFO, getClassName() << ": ApplyEqs: "         << sstr.str()     << "\n");

  // parser
  m_vFunction.setFunctions(m_functions);
  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

}

//////////////////////////////////////////////////////////////////////////////

void BaseBC::setup()
{
  CFAUTOTRACE;

  // first call parent method
  UFEMSolverCom::setup();

  // get nb of states, dim, eqs in trs and init m_applyFlags
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  std::string bndTrsName=getCurrentTRS()->getName();
  CFuint i; for (i = 0; i < trs.size(); ++i)
    if (trs[i]->getName() == bndTrsName)
      break;
  if (i==trs.size()) Common::NoSuchValueException(FromHere(),"Trs not found with name: "+bndTrsName);
  const CFuint nbTrsStates = trs[i]->getNbStatesInTrs();
  m_applyFlags.resize(nbTrsStates);
  for (CFuint i=0; i<nbTrsStates; i++) {
    m_applyFlags[i].resize(nbEqs);
    for (CFuint j=0; j<nbEqs; j++) m_applyFlags[i][j]=false;
  }

  // validate ApplyEqs config option and resize flags and vars
  m_parseVector.resize(nbDim+1+nbEqs);
  m_applyVars.resize(nbEqs);
  for (CFuint i=0; i<nbEqs; ++i){
    m_applyVars[i]=0.;
  }
  if (!m_applyEqs.size())
  {
    // apply to all PhysicalModel equations
    m_applyEqs.resize(nbEqs);
    for (CFuint i=0; i<nbEqs; ++i){
      m_applyEqs[i]=i;
    }
  }
  else
  {
    // validate the equations to apply the BC exist
    for (CFuint i=0; i<m_applyEqs.size(); ++i)
      if (m_applyEqs[i]>= nbEqs)
        throw BadValueException(FromHere(),
          "ApplyEqs refers to an equation that doesn't exist." );
  }
  for (CFuint i=0; i<nbTrsStates; i++) {
    for (CFuint j=0; j<m_applyEqs.size(); ++j){
    m_applyFlags[i][m_applyEqs[j]]=true;
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

void BaseBC::executeOnTrs()
{
}

//////////////////////////////////////////////////////////////////////////////

void BaseBC::computeStateValuesBaseBC(const State* currState)
{
  // if parse from CFcase
  const RealVector& temp = currState->getCoordinates();
  const CFuint nbDim=temp.size();
  const CFuint nbEqs=(*currState).size();
  for (CFuint i=0; i<nbDim; ++i) m_parseVector[i] = temp[i];
  m_parseVector[nbDim] = SubSystemStatusStack::getActive()->getCurrentTime();
  for (CFuint i=0; i<nbEqs; ++i) m_parseVector[nbDim+1+i] = (*currState)[i];
  m_vFunction.evaluate(m_parseVector,m_applyVars);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > BaseBC::needsSockets()
{
  CFAUTOTRACE;

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_interStates);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD


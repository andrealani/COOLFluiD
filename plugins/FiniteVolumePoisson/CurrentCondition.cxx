#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/PE.hh"
#include "Common/ParserException.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolumePoisson/CurrentCondition.hh"
#include "Poisson/PoissonConvTerm.hh"

#ifdef CF_HAVE_MPI
#include <mpi.h>
#endif


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CurrentCondition, CellCenterFVMData, FiniteVolumeModule> CurrentConditionFVMCCProvider("CurrentConditionFVMCC");

//////////////////////////////////////////////////////////////////////////////

void CurrentCondition::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >("ImposedCurrent", "Value for the imposed current");
  //options.addConfigOption< std::vector<std::string> >("Def","Definition of the Function.");
  //options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
}

//////////////////////////////////////////////////////////////////////////////

CurrentCondition::CurrentCondition(const std::string& name) :
  FVMCC_BC(name),
  socket_faceAreas("faceAreas"),
  m_library(CFNULL),
  m_sigmaIntegral(0.),
  _iter(),
  _resultI()
{
  addConfigOptionsTo(this);

  m_imposedI = 0.;
  setParameter("ImposedCurrent", &m_imposedI);

  //_function = std::vector<std::string>();
  //setParameter("Def",&_function);

  //_vars = std::vector<std::string>();
  //setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

CurrentCondition::~CurrentCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

void CurrentCondition::setup()
{

  FVMCC_BC::setup();

  //_vars.resize(PhysicalModelStack::getActive()->getDim() + 1);

  using namespace Framework;

  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  //cf_assert (m_library.isNotNull());

  _iter.resize(1);
  _resultI.resize(1);
}

//////////////////////////////////////////////////////////////////////////////

void CurrentCondition::unsetup()
{
  FVMCC_BC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CurrentCondition::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CurrentCondition::setGhostState(GeometricEntity *const face)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;

  //_vFunction.setFunctions(_function);
  //_vFunction.setVariables(_vars);
  //try {
    //_vFunction.parse();
  //}
  //catch (Common::ParserException& e) {
    //CFout << e.what() << "\n";
    //throw; // retrow the exception to signal the error to the user
  //}

  //_iter = SubSystemStatusStack::getActive()->getNbIter();

  //_vFunction.evaluate(_iter, _resultI);
  //m_imposedI = _resultI[0];

  // here it is assumed that the ArcJet induction equations are the last ones
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  //m_library->setComposition(Tdim, pdim, CFNULL);

  const CFreal sigma = 1.; //m_library->sigma(Tdim, pdim, tVec);
  cf_assert(sigma > 0.);
  cf_assert(m_sigmaIntegral > 0.);

  const CFreal dVz = -m_imposedI/m_sigmaIntegral;
  const CFreal dr = MathFunctions::getDistance(ghostState->getCoordinates(),
                                                   innerState->getCoordinates());
  (*ghostState)[0] = dVz*dr + (*innerState)[0]; // AAL: This line is the one imposing the current
}

//////////////////////////////////////////////////////////////////////////////

void CurrentCondition::preProcess()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;

  DataHandle<CFreal> faceAreas   = socket_faceAreas.getDataHandle();

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
    geoBuilder = this->getMethodData().getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(this->socket_states, this->socket_gstates, this->socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;
  SafePtr<TopologicalRegionSet> trs = this->getCurrentTRS();
  geoData.trs = trs;

  //const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // AL: this is not general, change for NEQ
  //const CFuint pID = 0;
  //const CFuint TID = nbEqs - 2;

  m_sigmaIntegral = 0.;
  CFreal localSigmaIntegral = 0.;
  const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
    CFLogDebugMed( "iFace = " << iFace << "\n");

    // build the GeometricEntity
    geoData.idx = iFace;

    GeometricEntity *const face = geoBuilder->buildGE();
//    CFreal pdim = eulerTerm->getPressureFromState((*face->getState(0))[pID]);
//    CFreal Tdim = (*face->getState(0))[TID];
//    CFreal* tVec = CFNULL;

//    m_library->setComposition(Tdim,pdim,CFNULL);

//    const CFreal sigma = m_library->sigma(Tdim, pdim, tVec);
    const CFreal sigma = 1; 				//WATCH OUT: only for debugging
    cf_assert(sigma > 0.);
    if (face->getState(0)->isParUpdatable()) {
      localSigmaIntegral += sigma*faceAreas[face->getID()];
    }

    // release the GeometricEntity
    geoBuilder->releaseGE();
  }

  m_sigmaIntegral = localSigmaIntegral;
#ifdef CF_HAVE_MPI
  const std::string nsp = this->getMethodData().getNamespace();
  MPI_Allreduce(&localSigmaIntegral, &m_sigmaIntegral, 1, MPI_DOUBLE, MPI_SUM,
        PE::GetPE().GetCommunicator(nsp));
#endif
  cf_assert(m_sigmaIntegral > 0.);

  CFLog(DEBUG_MIN, "CurrentCondition::sigmaIntegral = " << m_sigmaIntegral <<"\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > CurrentCondition::needsSockets()
{
  cout<<"CurrentCondition::needsSockets() \n";
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result  = FVMCC_BC::needsSockets();
  //std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_faceAreas);
  return result;
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

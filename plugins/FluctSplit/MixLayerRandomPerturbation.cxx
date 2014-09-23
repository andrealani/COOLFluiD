#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/SingleBehaviorFactory.hh"
//#include "Environment/DirPaths.hh"
//#include "Environment/FileHandlerInput.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/MixLayerRandomPerturbation.hh"
#include "NavierStokes/NavierStokes.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/SpaceMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Physics {

      namespace NavierStokes {

  //CFreal MixLayerRandomPerturbation::m_Velocity = 0.;
  //CFreal MixLayerRandomPerturbation::m_di = 0.;

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MixLayerRandomPerturbation,
                      DataProcessingData,
		      NavierStokesModule>
aMixLayerRandomPerturbationProvider("MixLayerRandomPerturbation");
//////////////////////////////////////////////////////////////////////////////

void MixLayerRandomPerturbation::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Velocity","Velocity.");
  options.addConfigOption< CFreal >("di","Initial vorticity thickness.");
  options.addConfigOption< CFreal >("order","Order relative to velocity.");
}

//////////////////////////////////////////////////////////////////////////////

MixLayerRandomPerturbation::MixLayerRandomPerturbation( const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  _solToUpdateVar(),
  _updateToSolVar()
{
  addConfigOptionsTo(this);

  //m_Velocity = 1.;
  setParameter("Velocity",&m_Velocity);
  //m_di = 1.;
  setParameter("di",&m_di);
  setParameter("order",&m_order);
}

//////////////////////////////////////////////////////////////////////////////

MixLayerRandomPerturbation::~MixLayerRandomPerturbation()
{
}

//////////////////////////////////////////////////////////////////////////////

void MixLayerRandomPerturbation::setup()
{
  DataProcessingCom::setup(); // first call setup of parent class
  
  SafePtr<SpaceMethodData> smdata = getMethodData().getCollaborator<SpaceMethod>()->
    getSpaceMethodData();
  cf_assert(smdata.isNotNull());
  
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getActive();
  const std::string provSolToUpdate = VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), smdata->getSolutionVarStr(), smdata->getUpdateVarStr());
  const std::string provUpdateToSol = VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), smdata->getUpdateVarStr(), smdata->getSolutionVarStr());
  
  _solToUpdateVar = Environment::Factory<VarSetTransformer>::getInstance().getProvider(provSolToUpdate)->
    create(physModel->getImplementor());
  
  _updateToSolVar = Environment::Factory<VarSetTransformer>::getInstance().getProvider(provUpdateToSol)->
    create(physModel->getImplementor());
  
}// setup()

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
MixLayerRandomPerturbation::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MixLayerRandomPerturbation::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  return result;
}


//////////////////////////////////////////////////////////////////////////////

void MixLayerRandomPerturbation::execute()
{
  //CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  CFreal nbstates = states.size();
  State consState; // state in conservative variables
  
  //CF_DEBUG_OBJ(nbstates);
  const CFuint k = SubSystemStatusStack::getActive()->getNbIter();
  //CF_DEBUG_OBJ(k);
  if (k == 0){
    CF_DEBUG_OBJ("Random perturbation introduced.");
    srand(time(NULL));
    for (CFuint i = 0; i < nbstates; ++i) {
      // transform from update to conservative variables 
      consState = *_updateToSolVar->transform(states[i]);
      Node& icoord = states[i]->getCoordinates();
      //CF_DEBUG_OBJ(m_di);
      //CF_DEBUG_OBJ(m_Velocity);
      if (icoord[0]!=0.0 && abs(icoord[0]-14.0*m_di)>1.0e-7){
	consState[1] = consState[1]+consState[0]*m_Velocity*m_order*exp(-pow(icoord[1]/m_di,2))*( ((rand()%1000)-500.0)/500.0);
	consState[2] =             consState[0]*m_Velocity*m_order*exp(-pow(icoord[1]/m_di,2))*( ((rand()%1000)-500.0)/500.0 + sin(0.8975/m_di*icoord[0]));
      }
      
      // transform from conservative to update variables 
      *states[i] = *_solToUpdateVar->transform(&consState);
    }
  }
}//execute

//////////////////////////////////////////////////////////////////////////////

void MixLayerRandomPerturbation::unsetup()
{
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void MixLayerRandomPerturbation::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




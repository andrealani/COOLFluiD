#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/SingleBehaviorFactory.hh"
//#include "Environment/DirPaths.hh"
//#include "Environment/FileHandlerInput.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/MixLayer3DRandomPerturbation.hh"
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

MethodCommandProvider<MixLayer3DRandomPerturbation,
                      DataProcessingData,
		      NavierStokesModule>
aMixLayer3dRandomPerturbationProvider("MixLayer3DRandomPerturbation");
//////////////////////////////////////////////////////////////////////////////

void MixLayer3DRandomPerturbation::defineConfigOptions(Config::OptionList& options)
{
  //  options.addConfigOption< CFreal >("Velocity","Velocity.");
  // options.addConfigOption< CFreal >("di","Initial vorticity thickness.");
  //options.addConfigOption< CFreal >("order","Order relative to velocity.");
}

//////////////////////////////////////////////////////////////////////////////

MixLayer3DRandomPerturbation::MixLayer3DRandomPerturbation( const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  _solToUpdateVar(),
  _updateToSolVar()
{
  addConfigOptionsTo(this);

  //m_Velocity = 1.;
  //  setParameter("Velocity",&m_Velocity);
  //m_di = 1.;
  // setParameter("di",&m_di);
  // setParameter("order",&m_order);
}

//////////////////////////////////////////////////////////////////////////////

MixLayer3DRandomPerturbation::~MixLayer3DRandomPerturbation()
{
}

//////////////////////////////////////////////////////////////////////////////

void MixLayer3DRandomPerturbation::setup()
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
MixLayer3DRandomPerturbation::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MixLayer3DRandomPerturbation::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  return result;
}


//////////////////////////////////////////////////////////////////////////////

void MixLayer3DRandomPerturbation::execute()
{
  //CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  CFreal nbstates = states.size();
  State consState; // state in conservative variables
  CFreal eps;
  CFreal s;
  //CF_DEBUG_OBJ(nbstates);
  const CFuint k = SubSystemStatusStack::getActive()->getNbIter();
  //CF_DEBUG_OBJ(k);
  if (k == 0){
    CF_DEBUG_OBJ("Random perturbation introduced.");
    srand(time(NULL));
    for (CFuint i = 0; i < nbstates; ++i) {
      //  CF_DEBUG_OBJ(i);
      // transform from update to conservative variables 
      eps = ((rand() % -1000 + 1000)+1)/1000.0;
      s = 1.0;
      consState = *_updateToSolVar->transform(states[i]);
      Node& icoord = states[i]->getCoordinates();
      consState[0] = 1.0*s*eps;
      consState[1] = 1.5*(1.0-(icoord[2]-1.0)*(icoord[2]-1.0))*s*eps*(1.0+s*eps);
      consState[2] = s*eps*s*eps;
      consState[3] = s*eps*s*eps;
      consState[4] = (287.046/0.4)*(1.0+(0.4/3.0)*0.72*0.25*2.25*(1.0-pow(icoord[2]-1.0,4)))+
				    0.5*1.5*1.5*(1.0-(icoord[2]-1.0)*(icoord[2]-1.0))*s*eps*(1.0+s*eps)*(1.0+s*eps);
      
     
      
      // transform from conservative to update variables 
      *states[i] = *_solToUpdateVar->transform(&consState);
    }
  }
}//execute

//////////////////////////////////////////////////////////////////////////////

void MixLayer3DRandomPerturbation::unsetup()
{
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void MixLayer3DRandomPerturbation::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




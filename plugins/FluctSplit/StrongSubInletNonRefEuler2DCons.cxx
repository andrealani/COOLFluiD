#include <numeric>
#include <algorithm>

#include "StrongSubInletNonRefEuler2DCons.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "MathTools/MatrixInverter.hh"
#include "Environment/ObjectProvider.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"
#include "Framework/CFL.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSubInletNonRefEuler2DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> strongSubInletNonRefEuler2DConsProvider("StrongSubInletNonRefEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletNonRefEuler2DCons::defineConfigOptions(Config::OptionList&
options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Variable names.");
  options.addConfigOption< std::vector<std::string> >("InFlow","Function defining the incoming flow.");
}

//////////////////////////////////////////////////////////////////////////////

StrongSubInletNonRefEuler2DCons::StrongSubInletNonRefEuler2DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_isBState("isBState"),
  _varSet(),
  m_nbEqs(),
  _r1(),
  _r2(),
  _r3()
{
  addConfigOptionsTo(this);

  m_vars_inflow.resize(0);

  m_vars_inflow = std::vector<std::string>();
  setParameter("Vars",&m_vars_inflow);

  m_function_inflow = vector<std::string>();
  setParameter("InFlow",&m_function_inflow);
}

//////////////////////////////////////////////////////////////////////////////

StrongSubInletNonRefEuler2DCons::~StrongSubInletNonRefEuler2DCons()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubInletNonRefEuler2DCons::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isBState);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletNonRefEuler2DCons::setup()
{
  FluctuationSplitCom::setup();

  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();

// create boundary nodal normals, which pointing outwards
  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  
  _r1.resize(PhysicalModelStack::getActive()->getNbEq());
  _r2.resize(PhysicalModelStack::getActive()->getNbEq());
  _r3.resize(PhysicalModelStack::getActive()->getNbEq());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
 obj.setDataSockets(socket_normals,socket_faceNeighCell);
 obj.create(getTrsList(), _bcNormals);
 
 // handling the inflow disturbances
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  inflow.resize(states.size());

  const CFuint nb_funs = m_function_inflow.size();

  for (CFuint i = 0; i < states.size(); ++i)
  {
    inflow[i].resize(m_nbEqs);
    for (CFuint j = 0; j < m_nbEqs; ++j)
      inflow[i][j] = 0.0;
  }
 
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletNonRefEuler2DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &(_bcNormals[getCurrentTrsID()]);

  DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<bool> isBState = socket_isBState.getDataHandle();

  Common::SafePtr< vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();

  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFuint TT   = dim;
  const CFreal time =
  SubSystemStatusStack::getActive()->getCurrentTimeDim();
  m_var_values[TT] = time;
  
  const CFreal gamma = _varSet->getModel()->getGamma();
  
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
   const CFuint stateID = (*statesIdx)[iState]; 
   Node& coord = states[stateID]->getCoordinates();
   for (CFuint iCoor = 0; iCoor < dim; ++iCoor) {
    m_var_values[iCoor] = coord[iCoor];
   }

   m_function_parser_inflow.evaluate(m_var_values, inflow[stateID]);
  }
  

// go through all the states involved in the boundary trs
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *state = states[stateID];

    
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    const CFreal nx = (*bcNormal)[0];
    const CFreal ny = (*bcNormal)[1];
    
    (*state)[0] = inflow[stateID][0];
    (*state)[1] = inflow[stateID][1];
    (*state)[2] = inflow[stateID][2];
    (*state)[3] = inflow[stateID][3];
    

    if (!isUpdated[stateID]) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
        _varSet->setEigenVect1(_r1, *state, *bcNormal);
        _varSet->setEigenVect2(_r2, *state, *bcNormal);
        _varSet->setEigenVect3(_r3, *state, *bcNormal);

        CFreal *const rhsStart = &rhs(stateID, 0, m_nbEqs);
	
	const RealVector& linearData = _varSet->getModel()->getPhysicalData();
        const CFreal c     = linearData[EulerTerm::A];

        CFreal drho = rhsStart[0];
        CFreal drho0u = rhsStart[1];
        CFreal drho0v = rhsStart[2];
        CFreal dp = (gamma-1.0)*(rhsStart[3] - 0.5* (rhsStart[1]*rhsStart[1]+rhsStart[2]*rhsStart[2])/rhsStart[0]);
	

        const CFreal beta1 = -(drho - dp/(c*c));
        const CFreal beta2 = -(ny*drho0u - nx*drho0v);
        const CFreal beta4 = -(-nx*drho0u - ny*drho0v + dp/c);

        _r1 *= beta1;
        _r2 *= beta2;
        _r3 *= beta4;
        _r3 += _r2 += _r1; // This is actually the whole boundary correction term

        for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
          rhs(stateID, iEq, m_nbEqs) += _r3[iEq];
        }
      }







    isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletNonRefEuler2DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  const std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
		 create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());
  
  cf_assert(_varSet.isNotNull());

    m_var_values.resize(m_vars_inflow.size());

  if(m_function_inflow.empty())
     throw BadValueException(FromHere(),"StrongSubInletNonRefEuler2DCons::setFuntion(): no incoming flow function provided.");

  // configure the expression for the mean flow
  m_function_parser_inflow.setFunctions(m_function_inflow);
  m_function_parser_inflow.setVariables(m_vars_inflow);

  try
  {
    m_function_parser_inflow.parse();
  }
  catch (Common::ParserException& e)
  {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletNonRefEuler2DCons::unsetup()
{
    inflow.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

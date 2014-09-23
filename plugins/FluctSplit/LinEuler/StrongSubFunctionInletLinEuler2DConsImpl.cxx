#include <numeric>
#include <algorithm>

#include "StrongSubFunctionInletLinEuler2DConsImpl.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitLinEuler.hh"
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
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSubFunctionInletLinEuler2DConsImpl, FluctuationSplitData, FluctSplitLinEulerModule> strongSubFunctionInletLinEuler2DConsImpProvider("StrongSubFunctionInletLinEuler2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

StrongSubFunctionInletLinEuler2DConsImpl::StrongSubFunctionInletLinEuler2DConsImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  _varSet(),
  m_nbEqs(),
  _jacobElem(4,4),
  _jacob(4,4),
  _jacobAll(),
  _in(4),
  _ira(4), 
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

void StrongSubFunctionInletLinEuler2DConsImpl::defineConfigOptions(Config::OptionList&
options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Variable names.");
  options.addConfigOption< std::vector<std::string> >("InFlow","Function defining the incoming flow.");
}

//////////////////////////////////////////////////////////////////////////////

StrongSubFunctionInletLinEuler2DConsImpl::~StrongSubFunctionInletLinEuler2DConsImpl()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubFunctionInletLinEuler2DConsImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubFunctionInletLinEuler2DConsImpl::setup()
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

//   executeOnTrs();
  
  
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubFunctionInletLinEuler2DConsImpl::executeOnTrs()
{
  
  Common::SafePtr<LinearSystemSolver> lss =
  getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  jacobMatrix->finalAssembly();  
  
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  
  _jacobAll.resize(statesIdx->size());
  for (CFuint i = 0; i < _jacobAll.size(); ++i) {
    _jacobAll[i].resize(m_nbEqs, m_nbEqs);
  }
  
  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();
  
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFuint TT   = dim;
  const CFreal time =
  SubSystemStatusStack::getActive()->getCurrentTimeDim();
  m_var_values[TT] = time;
  
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
   const CFuint stateID = (*statesIdx)[iState]; 
   Node& coord = states[stateID]->getCoordinates();
   for (CFuint iCoor = 0; iCoor < dim; ++iCoor) {
    m_var_values[iCoor] = coord[iCoor];
   }

   m_function_parser_inflow.evaluate(m_var_values, inflow[stateID]);
  }
      
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    const CFuint globalID = idxMapping.getColID(stateID)*m_nbEqs;
    
    State *state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    const CFreal nx = (*bcNormal)[0];
    const CFreal ny = (*bcNormal)[1];
    
    for(int i=0; i<m_nbEqs; i++)
      for(int j=0; j<m_nbEqs; j++)
	_jacob(i,j) = 0.0;
    
    (*state)[0] = inflow[stateID][0];
    (*state)[1] = inflow[stateID][1];
    (*state)[2] = inflow[stateID][2];
    (*state)[3] = inflow[stateID][3];
    

    if (!isUpdated[stateID] && state->isParUpdatable()) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
        _varSet->setEigenVect1(_r1, *state, *bcNormal);
        _varSet->setEigenVect2(_r2, *state, *bcNormal);
        _varSet->setEigenVect3(_r3, *state, *bcNormal);

        CFreal *const rhsStart = &rhs(stateID, 0, m_nbEqs);

        CFreal drho = rhsStart[0];
        CFreal drho0u = rhsStart[1];
        CFreal drho0v = rhsStart[2];
        CFreal dp = rhsStart[3];

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
        
      for (CFuint iVar = 0; iVar < m_nbEqs; ++iVar) {
      _ira[iVar] = globalID + iVar;
     }

	  // get the elements of the jacobian matrix
     jacobMatrix->getValues(m_nbEqs,
  			 &_ira[0],
  			 m_nbEqs,
  			 &_ira[0],
  			 &_jacobElem[0]);

    
     _jacob(0,3) = _jacobElem(0,3);
     _jacob(1,3) = _jacobElem(1,3);
     _jacob(2,3) = _jacobElem(2,3);
     _jacob(3,3) = _jacobElem(3,3);
  
     
     for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
       for (CFuint jEq = 0; jEq < m_nbEqs; ++jEq) {
       _jacobAll[iState](iEq, jEq) = _jacob(iEq, jEq);
       }
     }
     
      }
     }
    isUpdated[stateID] = true; // flagging is important!!!!!
    }
    
    
    jacobMatrix->finalAssembly();
    
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    const CFuint globalID = idxMapping.getColID(stateID)*m_nbEqs;
    State *const state = states[stateID];   
    
    if (state->isParUpdatable()) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
    
    for (CFuint iVar = 0; iVar < m_nbEqs; ++iVar) {
     _ira[iVar] = globalID + iVar;
    }
    
    
    for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < m_nbEqs; ++jEq) {
       _jacob(iEq, jEq) = _jacobAll[iState](iEq, jEq);
       }
     }
    
    jacobMatrix->setValues(m_nbEqs,
  			 &_ira[0],
  			 m_nbEqs,
  			 &_ira[0],
  			 &_jacob[0]);
        }
      }
    }
    
    jacobMatrix->finalAssembly();  
    
    
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubFunctionInletLinEuler2DConsImpl::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "LinEuler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::LinearizedEuler::LinEuler2DCons>());

  cf_assert(_varSet.isNotNull());
  
  
  m_var_values.resize(m_vars_inflow.size());

  if(m_function_inflow.empty())
     throw BadValueException(FromHere(),"StrongSubFunctionInletLinEuler2DConsImpl::setFuntion(): no incoming flow function provided.");

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
void StrongSubFunctionInletLinEuler2DConsImpl::unsetup()
{
  inflow.resize(0);
  _bcNormals.resize(0);
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

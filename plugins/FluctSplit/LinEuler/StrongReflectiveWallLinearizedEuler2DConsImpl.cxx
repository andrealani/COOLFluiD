#include "FluctSplitLinEuler.hh"
#include <numeric>
#include <algorithm>

#include "StrongReflectiveWallLinearizedEuler2DConsImpl.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
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

MethodCommandProvider<StrongReflectiveWallLinearizedEuler2DConsImpl, FluctuationSplitData, FluctSplitLinEulerModule> StrongReflectiveWallLinearizedEuler2DConsImplProvider("StrongReflectiveWallLinearizedEuler2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

StrongReflectiveWallLinearizedEuler2DConsImpl::StrongReflectiveWallLinearizedEuler2DConsImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  socket_bStatesNeighbors("bStatesNeighbors"),  
  _jacobElem(4,4),
  _jacob(4,4),
  _in(4),
  _ira(4),
  _varSet(),
  _bcNormals(),
  _r3(),
  _jacobAll()
{
}

//////////////////////////////////////////////////////////////////////////////

StrongReflectiveWallLinearizedEuler2DConsImpl::~StrongReflectiveWallLinearizedEuler2DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongReflectiveWallLinearizedEuler2DConsImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_bStatesNeighbors);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongReflectiveWallLinearizedEuler2DConsImpl::setup()
{

  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  CFuint nbE = PhysicalModelStack::getActive()->getNbEq();
  _r3.resize(nbE);

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);
  
}

//////////////////////////////////////////////////////////////////////////////

void StrongReflectiveWallLinearizedEuler2DConsImpl::executeOnTrs()
{
  Common::SafePtr<LinearSystemSolver> lss =
  getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  jacobMatrix->finalAssembly();

  CFuint m_nbEqs = PhysicalModelStack::getActive()->getNbEq();

  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< std::valarray<Framework::State*> > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle(); 
  
  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  
    _jacobAll.resize(statesIdx->size());
  for (CFuint i = 0; i < _jacobAll.size(); ++i) {
    _jacobAll[i].resize(m_nbEqs, m_nbEqs);
  }

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];
  const CFreal oneoverc = 1./c;

  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    const CFuint globalID = idxMapping.getColID(stateID)*m_nbEqs;
    
    State *const state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    const CFreal nx = (*bcNormal)[0];
    const CFreal ny = (*bcNormal)[1];
    
    for(int i=0; i<m_nbEqs; i++)
      for(int j=0; j<m_nbEqs; j
	++)
	_jacob(i,j) = 0.0;

    if (!isUpdated[stateID] && state->isParUpdatable()) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {

        _varSet->setEigenVect3(_r3, *state, *bcNormal);

        CFreal *const rhsStart = &rhs(stateID, 0, m_nbEqs);

        //CFreal drho = rhsStart[0]; not needed
        CFreal drho0u = rhsStart[1];
        CFreal drho0v = rhsStart[2];
        CFreal dp = rhsStart[3];

        // Calculate DW4 (corresponding to lambda_4 = Vn - c
        const CFreal DW4 = -nx*drho0u - ny*drho0v + dp/c;
        // Set DW3 = DW4 and solve for beta3
        const CFreal beta3 = DW4 - (nx*drho0u + ny*drho0v + dp/c);

        _r3 *= beta3; // this is actually beta3*r3 == boundary correction

        for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
          rhs(stateID, iEq, m_nbEqs) += _r3[iEq];
        }

      for (CFuint iVar = 0; iVar < m_nbEqs; ++iVar) {
      _ira[iVar] = globalID + iVar;
      _in[iVar] = globalID + iVar;
     }

	  // get the elements of the jacobian matrix
    jacobMatrix->getValues(m_nbEqs,
  			 &_ira[0],
  			 m_nbEqs,
  			 &_in[0],
  			 &_jacobElem[0]);

    
     CFreal A0 =  2*(nx*_jacobElem(0,1)+ny*_jacobElem(0,2)); //drho
     CFreal A1 =  2*(nx*_jacobElem(1,1)+ny*_jacobElem(1,2)); //drho0u
     CFreal A2 =  2*(nx*_jacobElem(2,1)+ny*_jacobElem(2,2)); //drho0v
     CFreal A3 =  2*(nx*_jacobElem(3,1)+ny*_jacobElem(3,2)); //dp
     

     _jacob(0,2) = 0.5*oneoverc*(A3+A2)+A0;
     _jacob(1,2) = -0.5*nx*(A3-A2)+ny*A1; 
     _jacob(2,2) = -0.5*ny*(A3-A2)-nx*A1; 
     _jacob(3,2) = 0.5*c*(A3+A2); 
    
     
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
      _in[iVar] = globalID + iVar;      
     }
    
    
    for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < m_nbEqs; ++jEq) {
       _jacob(iEq, jEq) = _jacobAll[iState](iEq, jEq);
       }
     }
    
    jacobMatrix->addValues(m_nbEqs,
  			 &_ira[0],
  			 m_nbEqs,
  			 &_in[0],
  			 &_jacob[0]);
    
        }
      }
    }
    
    jacobMatrix->finalAssembly();  
    
    
}

//////////////////////////////////////////////////////////////////////////////

void StrongReflectiveWallLinearizedEuler2DConsImpl::configure( Config::ConfigArgs& args)
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "LinEuler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::LinearizedEuler::LinEuler2DCons>());

  cf_assert(_varSet.isNotNull());

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

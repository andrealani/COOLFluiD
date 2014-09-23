#include "FluctSplitLinEuler.hh"
#include <numeric>

#include "StrongSubInletHedstrom2DConsImpl.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "MathTools/MathFunctions.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"


#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/LinEuler2DVarSet.hh"
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

MethodCommandProvider<StrongSubInletHedstrom2DConsImpl, FluctuationSplitData, FluctSplitLinEulerModule> StrongSubInletHedstrom2DConsImplProvider("StrongSubInletHedstrom2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

StrongSubInletHedstrom2DConsImpl::StrongSubInletHedstrom2DConsImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _bcNormals(),
  _jacobElem(4,4),
  _jacob(4,4),
  _jacobAll(),
  _in(4),
  _ira(4),
  _r1(),
  _r2(),
  _r3()
{
}

//////////////////////////////////////////////////////////////////////////////

StrongSubInletHedstrom2DConsImpl::~StrongSubInletHedstrom2DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubInletHedstrom2DConsImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletHedstrom2DConsImpl::setup()
{

  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  _r1.resize(PhysicalModelStack::getActive()->getNbEq());
  _r2.resize(PhysicalModelStack::getActive()->getNbEq());
  _r3.resize(PhysicalModelStack::getActive()->getNbEq());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletHedstrom2DConsImpl::executeOnTrs()
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

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  _jacobAll.resize(statesIdx->size());
  for (CFuint i = 0; i < _jacobAll.size(); ++i) {
    _jacobAll[i].resize(nbEqs, nbEqs);
  }
  
  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();  
    
  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];
  const CFreal oneoverc = 1./c;

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    const CFuint globalID = idxMapping.getColID(stateID)*nbEqs;
    State *const state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    const CFreal nx = (*bcNormal)[0];
    const CFreal ny = (*bcNormal)[1];
    
    for(int i=0; i<nbEqs; i++)
      for(int j=0; j<nbEqs; j++)
	_jacob(i,j) = 0.0;

    if (!isUpdated[stateID] && state->isParUpdatable()) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
        _varSet->setEigenVect1(_r1, *state, *bcNormal);
        _varSet->setEigenVect2(_r2, *state, *bcNormal);
        _varSet->setEigenVect3(_r3, *state, *bcNormal);

        CFreal *const rhsStart = &rhs(stateID, 0, nbEqs);

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

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) += _r3[iEq];
        }


      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      _ira[iVar] = globalID + iVar;
     }

	  // get the elements of the jacobian matrix
     jacobMatrix->getValues(nbEqs,
  			 &_ira[0],
  			 nbEqs,
  			 &_ira[0],
  			 &_jacobElem[0]);

    
     CFreal A0 =  (nx*_jacobElem(0,1)+ny*_jacobElem(0,2)+oneoverc*_jacobElem(0,3));
     CFreal A1 =  (nx*_jacobElem(1,1)+ny*_jacobElem(1,2)+oneoverc*_jacobElem(1,3));
     CFreal A2 =  (nx*_jacobElem(2,1)+ny*_jacobElem(2,2)+oneoverc*_jacobElem(2,3));
     CFreal A3 =  (nx*_jacobElem(3,1)+ny*_jacobElem(3,2)+oneoverc*_jacobElem(3,3));
     
     _jacob(0,0) = A0*oneoverc*0.5; 
     _jacob(0,1) = A1*oneoverc*0.5;  
     _jacob(0,2) = A2*oneoverc*0.5; 
     _jacob(0,3) = A3*oneoverc*0.5;
     
     _jacob(1,0) = A0*nx*0.5; 
     _jacob(1,1) = A1*nx*0.5;  
     _jacob(1,2) = A2*nx*0.5; 
     _jacob(1,3) = A3*nx*0.5;     
     
     _jacob(2,0) = A0*ny*0.5; 
     _jacob(2,1) = A1*ny*0.5;  
     _jacob(2,2) = A2*ny*0.5; 
     _jacob(2,3) = A3*ny*0.5;    

     _jacob(3,0) = A0*c*0.5; 
     _jacob(3,1) = A1*c*0.5;  
     _jacob(3,2) = A2*c*0.5; 
     _jacob(3,3) = A3*c*0.5;
     
  
     
     for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
       for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
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
    const CFuint globalID = idxMapping.getColID(stateID)*nbEqs;
    State *const state = states[stateID];   
    
    if (state->isParUpdatable()) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
     _ira[iVar] = globalID + iVar;
    }
    
    
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
       _jacob(iEq, jEq) = _jacobAll[iState](iEq, jEq);
       }
     }
    
    jacobMatrix->setValues(nbEqs,
  			 &_ira[0],
  			 nbEqs,
  			 &_ira[0],
  			 &_jacob[0]);
        }
      }
    }
    
    jacobMatrix->finalAssembly();  
    
    
}



//////////////////////////////////////////////////////////////////////////////

void StrongSubInletHedstrom2DConsImpl::configure( Config::ConfigArgs& args)
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

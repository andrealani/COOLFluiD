#include "UpwindBiasedMUSCLPolyRec.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/VarSetTransformer.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<UpwindBiasedMUSCLPolyRec, 
		       CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, 
		       FiniteVolumeModule> 
upwindBiasedMUSCLPolyRecProvider("UpwindBiasedMUSCL");

//////////////////////////////////////////////////////////////////////////////

UpwindBiasedMUSCLPolyRec::UpwindBiasedMUSCLPolyRec(const std::string& name) :
  FVMCC_PolyRec(name),
  _sixth(1./6.),
  _third(1./3.),
  _k1(1. + _sixth - _third),
  socket_stencil("stencil"),
  _upStates(4),
  _lStateBkp(),
  _rStateBkp()
{
}

//////////////////////////////////////////////////////////////////////////////

UpwindBiasedMUSCLPolyRec::~UpwindBiasedMUSCLPolyRec()
{
}

//////////////////////////////////////////////////////////////////////////////

void UpwindBiasedMUSCLPolyRec::configure ( Config::ConfigArgs& args )
{
  FVMCC_PolyRec::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpwindBiasedMUSCLPolyRec::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result =
    FVMCC_PolyRec::needsSockets();

  // Add the needed DataSocketSinks
  result.push_back(&socket_stencil);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UpwindBiasedMUSCLPolyRec::computeGradients()
{
  //  throw Common::NotImplementedException
  // (FromHere(), "UpwindBiasedMUSCLPolyRec::computeGradients()");
}

//////////////////////////////////////////////////////////////////////////////

void UpwindBiasedMUSCLPolyRec::extrapolateImpl(GeometricEntity* const face)
{
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  SafePtr<VarSetTransformer> updateToSolutionVS = 
    getMethodData().getUpdateToSolutionVecTrans();
  
  getMethodData().getUpdateVar()->setExtraData(true);
  
  cf_assert(_zeroGradient != CFNULL);
  
  // should we use conservative variables here???
  cf_assert(face->getID() < stencil.size());
  
  const vector<State*>& ss = stencil[face->getID()];  
  
  if (!isBoundaryFace()) {
    _upStates.resize(4);
    _upStates[0] = ss[0];
    _upStates[1] = ss[1];
    _upStates[2] = ss[2];
    _upStates[3] = ss[3];
    
    vector<State*>& cs = *updateToSolutionVS->transform(&_upStates);
    cf_assert(cs.size() == 4);
    getValues(LEFT)  = (*cs[0])*_k1 + (*cs[1])*_third - (*cs[2])*_sixth;
    getValues(RIGHT) = (*cs[0])*_third + (*cs[1])*_k1 - (*cs[3])*_sixth;
  }
  else { 
    _upStates.resize(3);
    _upStates[0] = ss[0];
    _upStates[1] = ss[1];
    _upStates[2] = ss[2];
    
    vector<State*>& cs = *updateToSolutionVS->transform(&_upStates);
    cf_assert(cs.size() == 4);
    
    getValues(LEFT)  = (*cs[0])*_k1 + (*cs[1])*_third - (*cs[2])*_sixth;
    getValues(RIGHT) =  getValues(LEFT);
    // (*cs[0])*_third + (*cs[1])*(2.*_third);
    //here use weighted reconstruction 
    //(*cs[0])*d1 + (*cs[1])*d2;
  }
  
  getBackupValues(LEFT)  = getValues(LEFT);
  getBackupValues(RIGHT) = getValues(RIGHT);
  
  // _lStateBkp = static_cast<RealVector>(*ss[0]);
  // _rStateBkp = static_cast<RealVector>(*ss[1]);
  
  getMethodData().getUpdateVar()->setExtraData(false);
}

//////////////////////////////////////////////////////////////////////////////

void UpwindBiasedMUSCLPolyRec::extrapolateImpl(GeometricEntity* const face,
					       CFuint iVar, CFuint leftOrRight)
{
  // the following is inefficient because of the variable transformation!!!
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  SafePtr<VarSetTransformer> updateToSolutionVS = 
    getMethodData().getUpdateToSolutionVecTrans();
  
  getMethodData().getUpdateVar()->setExtraData(true);
  
  // please note that you work with references !!!
  // don't forget the "&" !!!
  State& v  = getValues(leftOrRight);
  const vector<State*>& ss = stencil[face->getID()];  
  
  copyBackupValues();
  if (leftOrRight == LEFT) { 
    _upStates.resize(4);
    _upStates[0] = ss[0];
    _upStates[1] = ss[1];
    _upStates[2] = ss[2];
    _upStates[3] = ss[3];
    
    vector<State*>& cs = *updateToSolutionVS->transform(&_upStates);
    v[iVar] = (*cs[0])[iVar]*_k1 + (*cs[1])[iVar]*_third - (*cs[2])[iVar]*_sixth;
  }
  
  if (leftOrRight == RIGHT && !isBoundaryFace()) {
    _upStates.resize(4);
    _upStates[0] = ss[0];
    _upStates[1] = ss[1];
    _upStates[2] = ss[2];
    _upStates[3] = ss[3];
    
    vector<State*>& cs = *updateToSolutionVS->transform(&_upStates);
    v[iVar] = (*cs[0])[iVar]*_third + (*cs[1])[iVar]*_k1 - (*cs[3])[iVar]*_sixth;
  }
  
  if (leftOrRight == RIGHT && isBoundaryFace()) {
    _upStates.resize(3);
    _upStates[0] = ss[0];
    _upStates[1] = ss[1];
    _upStates[2] = ss[2];

    vector<State*>& cs = *updateToSolutionVS->transform(&_upStates);
    v[iVar] = (*cs[0])[iVar]*_k1 + (*cs[1])[iVar]*_third - (*cs[2])[iVar]*_sixth;
    //(*cs[0])[iVar]*_third + (*cs[1])[iVar]*(2.*_third);
  } 
  
  getMethodData().getUpdateVar()->setExtraData(false);
}
      
//////////////////////////////////////////////////////////////////////////////
      
void UpwindBiasedMUSCLPolyRec::updateWeights()
{
  //  FVMCC_PolyRec::updateWeights();
  cout << "UpwindBiasedMUSCLPolyRec::updateWeights()" << endl;
  throw Common::NotImplementedException (FromHere(),"UpwindBiasedMUSCLPolyRec::updateWeights()");
}

//////////////////////////////////////////////////////////////////////////////

void UpwindBiasedMUSCLPolyRec::setup()
{
  FVMCC_PolyRec::setup();
  
  // re-setup the variable transformer
  getMethodData().getUpdateToSolutionVecTrans()->setup(4);
  _lStateBkp.resize(PhysicalModelStack::getActive()->getNbEq());
  _rStateBkp.resize(PhysicalModelStack::getActive()->getNbEq());
}
      
//////////////////////////////////////////////////////////////////////////////

void UpwindBiasedMUSCLPolyRec::computeLimiters()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

#include "FiniteVolume/FiniteVolume.hh"
#include "SilverMullerCondition2DProjection.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SilverMullerCondition2DProjection, CellCenterFVMData, FiniteVolumeModule> silverMullerCondition2DProjectionFVMCCProvider("SilverMullerCondition2DProjectionFVMCC");

//////////////////////////////////////////////////////////////////////////////

SilverMullerCondition2DProjection::SilverMullerCondition2DProjection(const std::string& name) :
  SuperInlet(name),
  _variables()
{
}

//////////////////////////////////////////////////////////////////////////////

SilverMullerCondition2DProjection::~SilverMullerCondition2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

void SilverMullerCondition2DProjection::setup()
{

  SuperInlet::setup();

  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);

}

//////////////////////////////////////////////////////////////////////////////

void SilverMullerCondition2DProjection::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void SilverMullerCondition2DProjection::configure ( Config::ConfigArgs& args )
{
  SuperInlet::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void SilverMullerCondition2DProjection::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;
  
   const CFuint faceID = face->getID();
   const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
   DataHandle<CFreal> normals = socket_normals.getDataHandle();
   
   CFreal nx = normals[startID];
   CFreal ny = normals[startID + 1];
   CFreal nz = 0;  
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  for(CFuint iDim = 0; iDim < nbDim; iDim++)
  {
    _variables[iDim] = _bCoord[iDim];
  }
  _variables[PhysicalModelStack::getActive()->getDim()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  //Evaluate the function
  _vFunction.evaluate(_variables,*_input);

  //Set the state value
  if (_inputAdimensionalValues){
    *ghostState = *_inputToUpdateVar->transform(_input); 
    (*ghostState)[6] = -(*innerState)[6] + (((*ghostState)[0] - (*innerState)[0])*nx + ((*ghostState)[1] - (*innerState)[1])*ny + ((*ghostState)[2] - (*innerState)[2])*nz);    
    (*ghostState)[7] = -(*innerState)[7] + (((*ghostState)[3] - (*innerState)[3])*nx + ((*ghostState)[4] - (*innerState)[4])*ny + ((*ghostState)[5] - (*innerState)[5])*nz);     
  }
  else{
    *_dimState = *_inputToUpdateVar->transform(_input);
    _varSet->setAdimensionalValues(*_dimState, *ghostState);
    (*ghostState)[6] = -(*innerState)[6] + (((*ghostState)[0] - (*innerState)[0])*nx + ((*ghostState)[1] - (*innerState)[1])*ny + ((*ghostState)[2] - (*innerState)[2])*nz);    
    (*ghostState)[7] = -(*innerState)[7] + (((*ghostState)[3] - (*innerState)[3])*nx + ((*ghostState)[4] - (*innerState)[4])*ny + ((*ghostState)[5] - (*innerState)[5])*nz);     
   
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

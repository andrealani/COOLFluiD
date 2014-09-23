#include "FiniteVolume/FiniteVolume.hh"
#include "Irradiation2DProjectionDim.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Maxwell/Maxwell2DProjectionAdimVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Irradiation2DProjectionDim, CellCenterFVMData, FiniteVolumeModule> irradiation2DProjectionDimFVMCCProvider("Irradiation2DProjectionDimFVMCC");

//////////////////////////////////////////////////////////////////////////////

Irradiation2DProjectionDim::Irradiation2DProjectionDim(const std::string& name) :
  SuperInlet(name),
  _variables()
{
}

//////////////////////////////////////////////////////////////////////////////

Irradiation2DProjectionDim::~Irradiation2DProjectionDim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Irradiation2DProjectionDim::setup()
{

  SuperInlet::setup();

//   _variables.resize(PhysicalModelStack::getActive()->getDim() + 1); AAL: This is for unsteady solutions
  _variables.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void Irradiation2DProjectionDim::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void Irradiation2DProjectionDim::configure ( Config::ConfigArgs& args )
{
  SuperInlet::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void Irradiation2DProjectionDim::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
   
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = 0;  
  CFreal const c_e = 299792458;
  
  // coordinate of the boundary point
  
//   cout << "Irradiation2DProjectionDim::setGhostState 1 \n";
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  for(CFuint iDim = 0; iDim < nbDim; iDim++)
  {
    _variables[iDim] = _bCoord[iDim];
  }
  
  //Evaluate the function
  _vFunction.evaluate(_variables,*_input);
  
  *ghostState = *_inputToUpdateVar->transform(_input);

  /*if (_inputAdimensionalValues){
    cout << "Irradiation2DProjectionDim::setGhostState 2.a \n";        
    *ghostState = *_inputToUpdateVar->transform(_input);
    cout << "Irradiation2DProjectionDim::setGhostState 2.b \n";    
    (*ghostState)[7] = (*innerState)[7] + (1/c_e)*((*ghostState)[3] - (*innerState)[3])*nx + ((*ghostState)[4] - (*innerState)[4])*ny + ((*ghostState)[5] - (*innerState)[5])*nz;     
    (*ghostState)[2] = (*innerState)[2];
    
  }
  else{
    *ghostState = *_inputToUpdateVar->transform(_input);    
    //cout << "Irradiation2DProjectionDim::setGhostState 2.c \n";
    // *_dimState = *_inputToUpdateVar->transform(_input);
    _varSet->setAdimensionalValues(*_dimState, *ghostState);
    cout << "Irradiation2DProjectionDim::setGhostState 2.d \n";
    (*ghostState)[7] = (*innerState)[7] + (1/c_e)*((*ghostState)[3] - (*innerState)[3])*nx + ((*ghostState)[4] - (*innerState)[4])*ny + ((*ghostState)[5] - (*innerState)[5])*nz;        
    (*ghostState)[2] = (*innerState)[2];
  }*/
//   cout << "Irradiation2DProjectionDim::setGhostState 3 \n";
  
    *ghostState *= 2.;
    *ghostState -= *innerState;
    
    (*ghostState)[5] = (*innerState)[5]; //AAL: special boundary condition for the Ez. dEz/dn = 0

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

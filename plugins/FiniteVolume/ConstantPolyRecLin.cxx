#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/ConstantPolyRecLin.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ConstantPolyRecLin, CellCenterFVMData, PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule> constantPolyRecLinProvider("ConstantLin");

//////////////////////////////////////////////////////////////////////////////

ConstantPolyRecLin::ConstantPolyRecLin(const std::string& name) :
  PolyReconstructorLin(name),
  socket_linearizedStates("linearizedStates"),
  socket_linearizedGhostStates("linearizedGhostStates")
{
}

//////////////////////////////////////////////////////////////////////////////

ConstantPolyRecLin::~ConstantPolyRecLin()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConstantPolyRecLin::extrapolateImpl(GeometricEntity* const face)
{
 //  VolumeIntegratorImpl* const impl =
//     getMethodData().getVolumeIntegrator()->getSolutionIntegrator(face);
  
//   const CFuint nbPoints = 1; // enough
  
//   // for each quadrature point, set
//   // left state  = current state
//   // right state = neighbor state
//   vector<State*>& leftValues = getValues(LEFT);
//   vector<State*>& rightValues = getValues(RIGHT);
//   vector<RealVector>& backUpLeftValues = getBackupValues(LEFT);
//   vector<RealVector>& backUpRightValues = getBackupValues(RIGHT);
  
//   for (CFuint ip = 0; ip < nbPoints; ++ip) {
//     *leftValues[ip]   = *face->getState(LEFT);
//     *rightValues[ip]  = *face->getState(RIGHT);
    
//     backUpLeftValues[ip] = *leftValues[ip];
//     backUpRightValues[ip] = *rightValues[ip];
//   }
  
//   // for each quadrature point, set
//   // left state  = current state
//   // right state = neighbor state
//   vector<State*>& leftValues2 = getValues2(LEFT);
//   vector<State*>& rightValues2 = getValues2(RIGHT);
//   vector<RealVector>& backUpLeftValues2 = getBackupValues2(LEFT);
//   vector<RealVector>& backUpRightValues2 = getBackupValues2(RIGHT);

//   const bool isGhostLeft  = (*face->getState(LEFT)).isGhost();
//   const bool isGhostRight = (*face->getState(RIGHT)).isGhost();

//   State* leftState2(CFNULL);
//   State* rightState2(CFNULL);

//   if(isGhostLeft) leftState2 = (socket_linearizedGhostStates.getDataHandle())[(*face->getState(LEFT)).getLocalID()];
//   else leftState2 = (socket_linearizedStates.getDataHandle())[(*face->getState(LEFT)).getLocalID()];

//   if(isGhostRight) rightState2 = (socket_linearizedGhostStates.getDataHandle())[(*face->getState(RIGHT)).getLocalID()];
//   else rightState2 = (socket_linearizedStates.getDataHandle())[(*face->getState(RIGHT)).getLocalID()];

//   for (CFuint ip = 0; ip < nbPoints; ++ip) {
//     *leftValues2[ip]   = *leftState2;
//     *rightValues2[ip]  = *rightState2;

//     backUpLeftValues2[ip] = *leftValues2[ip];
//     backUpRightValues2[ip] = *rightValues2[ip];
//   }

}

//////////////////////////////////////////////////////////////////////////////

void ConstantPolyRecLin::extrapolateImpl(GeometricEntity* const face,
				      CFuint iVar, CFuint leftOrRight)
{
  // VolumeIntegratorImpl* const impl =
//     getMethodData().getVolumeIntegrator()->getSolutionIntegrator(face);

//   const std::valarray<CFreal>& coeff = impl->getCoeff();
//   const CFuint nbPoints = coeff.size();

//   // for each quadrature point, set
//   // reconstructed variable = same variable in the current state
//   vector<State*>& values = getValues(leftOrRight);
//   for (CFuint ip = 0; ip < nbPoints; ++ip) {
//     // set the new values
//     (*values[ip])[iVar]   = (*face->getState(leftOrRight))[iVar];
//   }

//   const bool isGhost = (*face->getState(leftOrRight)).isGhost();

//   State* state2(CFNULL);

//   if(isGhost) state2 = (socket_linearizedGhostStates.getDataHandle())[(*face->getState(leftOrRight)).getLocalID()];
//   else state2 = (socket_linearizedStates.getDataHandle())[(*face->getState(leftOrRight)).getLocalID()];

//   // for each quadrature point, set
//   // reconstructed variable = same variable in the current state
//   vector<State*>& values2 = getValues2(leftOrRight);
//   for (CFuint ip = 0; ip < nbPoints; ++ip) {
//     // set the new values
//     (*values2[ip])[iVar]   = (*state2)[iVar];
//   }


}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

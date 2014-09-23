// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshTools/ComputeShortestDistance.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

ComputeShortestDistance::ComputeShortestDistance() :
  _numericalJacob(new NumericalJacobian("NumericalJacobian")),
  _currentFace(CFNULL),
  _mappedCoord(),
  _coord(),
  _residual(),
  _otherResidual(),
  _jacobian(),
  _relaxation(1.),
  _tolerance(0.0000000001),
  _maxIter(100),
  _isSetup(false)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeShortestDistance::~ComputeShortestDistance()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShortestDistance::setDimension(const CFuint dim)
{
  CFAUTOTRACE;

  _coord.resize(PhysicalModelStack::getActive()->getDim());
  _mappedCoord.resize(dim);
  _residual.resize(dim);
  _otherResidual.resize(dim);
  _jacobian.resize(dim);

  /// Reference Values for the Numerical Jacobian
  RealVector refValues(1.,dim);
  _numericalJacob->setRefValues(refValues);

  _isSetup = true;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShortestDistance::compute(GeometricEntity* face, RealVector& coord, CFreal& minimumDistance)
{
  CFAUTOTRACE;

  cf_assert(_isSetup);

  _currentFace = face;
  _coord = coord;

  setupNewton();

  bool isAchieved = false;
  CFreal previousDistance = 0.;
  CFreal newDistance = 0.;

  for(CFuint k = 0; !isAchieved; ++k) {

     if(k>0) previousDistance = _residual[0];

     takeStep();

     updateSolution();

     newDistance = _residual[0];

     ///@todo find a different stop condition
     if((k>_maxIter) || ((std::abs(previousDistance - newDistance) < _tolerance) && (k>0)))
     {
       isAchieved = true;
       CFLogDebugMax("Finished after: " << k << " iterations \n");
       CFLogDebugMax("            dL: " << std::abs(previousDistance - newDistance) << "\n");
       CFLogDebugMax("      distance: " << newDistance << "\n");
     }
  }

  minimumDistance = computeResidual();
  CFLogDebugMax("===========>  minimumDistance: " << minimumDistance << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShortestDistance::setupNewton()
{
  /// Initialization
  _mappedCoord = 0.333;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ComputeShortestDistance::computeResidual()
{
  RealVector coord = _currentFace->computeCoordFromMappedCoord(_mappedCoord);

  CFreal d2 = 0;
  for(CFuint i=0;i<coord.size();i++)
  {
    d2 += (_coord[i]-coord[i])*(_coord[i]-coord[i]);
  }

  return sqrt(d2);
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShortestDistance::takeStep()
{

  for (CFuint iVar = 0; iVar < _mappedCoord.size(); ++iVar) {
    _residual[iVar] = computeResidual();
    _numericalJacob->perturb(iVar, _mappedCoord[iVar]);
    _otherResidual[iVar] = computeResidual();
    _numericalJacob->restore(_mappedCoord[iVar]);
  }

  // jacobian contribution (dR/dU)_k
  // compute (R[U + dU] - R[U])/eps
  _numericalJacob->computeDerivative(_residual,
                                     _otherResidual,
                                     _jacobian);
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShortestDistance::updateSolution()
{
  CFreal dU = 0.;
  CFreal maxi = 0.;
  CFuint nbNodes = _currentFace->getNbNodesSolutionShapeFunction();

  for(CFuint iVar=0;iVar<_mappedCoord.size();++iVar)
  {
    dU = -_residual[iVar] * _jacobian[iVar];
    _mappedCoord[iVar] += _relaxation*dU;

    ///Add constraints
    if(nbNodes==3) //triangle
      {
      cf_assert(_mappedCoord.size() == 2);
      if (iVar == 0) maxi = 1.-min(1.,max(0.,_mappedCoord[0]));
      else maxi = 1.-min(1.,max(0.,_mappedCoord[1]));

      _mappedCoord[iVar] = min(maxi,max(0.,_mappedCoord[iVar]));
      }
    else //line or quad
      _mappedCoord[iVar] = min(1.,max(-1.,_mappedCoord[iVar]));
  }

//  std::cout << _residual[0] << " " <<_jacobian[0] << " " << dU << std::endl;

//   for(CFuint iVar=0;iVar<_mappedCoord.size();++iVar)
//   {
//     _residual[iVar] = computeResidual();
//   }

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

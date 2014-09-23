// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre5LagrangeLine.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP3.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeLine<LagrangeShapeFunctionLineP0> >
volGaussLegendre5LagrangeLineP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeLine<LagrangeShapeFunctionLineP1> >
volGaussLegendre5LagrangeLineP1Provider;

// Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
//                        VolumeGaussLegendre5LagrangeLine<LagrangeShapeFunctionLineP2> >
// volGaussLegendre5LagrangeLineP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre5LagrangeLine<INTERPOLATOR>::setMappedCoordinates()
{
  this->_mappedCoord.resize(3);

  this->_mappedCoord[0].resize(1);
  this->_mappedCoord[1].resize(1);
  this->_mappedCoord[2].resize(1);
  
  this->_mappedCoord[0][0] = -0.77459667;
  this->_mappedCoord[1][0] = 0.0;
  this->_mappedCoord[2][0] = 0.77459667;
}
     
//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre5LagrangeLine<INTERPOLATOR>::setWeights()
{
  this->_coeff.resize(getNbQuadraturePoints());
  cf_assert(getNbQuadraturePoints() == 3);
  
  this->_coeff[0] = 0.555555555;
  this->_coeff[1] = 0.888888889;
  this->_coeff[2] = 0.555555555;  
}
     
//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

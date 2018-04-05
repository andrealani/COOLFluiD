// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NullIntegrableEntity.hh"
#include "Framework/ContourIntegrator.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullIntegrableEntity,
               IntegrableEntity,
               FrameworkLib>
nullIntegrableEntityProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullIntegrableEntity::NullIntegrableEntity() :
  IntegrableEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

NullIntegrableEntity::~NullIntegrableEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullIntegrableEntity::integrateContourInGeo(GeometricEntity* const cell, RealVector& result)
{
}

//////////////////////////////////////////////////////////////////////////////

void NullIntegrableEntity::integrateContourInGeo(GeometricEntity* const cell, RealMatrix& result)
{
}

//////////////////////////////////////////////////////////////////////////////

void NullIntegrableEntity::integrateVolumeInGeo(GeometricEntity* const cell,  RealVector& result)
{
}

//////////////////////////////////////////////////////////////////////////////

void NullIntegrableEntity::integrateVolumeInGeo(GeometricEntity* const cell,  RealMatrix& result)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "IntegrableEntity.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

IntegrableEntity::IntegrableEntity() :
  Common::OwnedObject()
{
}

//////////////////////////////////////////////////////////////////////////////

IntegrableEntity::~IntegrableEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void IntegrableEntity::integrateContourInGeo(GeometricEntity* const cell,
                                             RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"IntegrableEntity::integrateContourInGeo()");
}

//////////////////////////////////////////////////////////////////////////////

void IntegrableEntity::integrateVolumeInGeo(GeometricEntity* const cell,
                                            RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"IntegrableEntity::integrateVolumeInGeo()");
}

//////////////////////////////////////////////////////////////////////////////

void IntegrableEntity::integrateContourInGeo(GeometricEntity* const cell,
                                             RealMatrix& result)
{
  throw Common::NotImplementedException (FromHere(),"IntegrableEntity::integrateContourInGeo()");
}

//////////////////////////////////////////////////////////////////////////////

void IntegrableEntity::integrateVolumeInGeo(GeometricEntity* const cell,
                                            RealMatrix& result)
{
  throw Common::NotImplementedException (FromHere(),"IntegrableEntity::integrateVolumeInGeo()");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

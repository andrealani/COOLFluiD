// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/GeometricEntityRegister.hh"
#include "Framework/BaseGeometricEntityProvider.hh"
#include "Framework/GeometricEntityFactory.hh"
#include "Framework/InterpolatorRegister.hh"
#include "Config/DuplicateNameException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {


//////////////////////////////////////////////////////////////////////////////

GeometricEntityRegister& GeometricEntityRegister::getInstance()
{
  static GeometricEntityRegister instance;
  return instance;
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntityRegister::GeometricEntityRegister()
{
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntityRegister::~GeometricEntityRegister()
{
}

//////////////////////////////////////////////////////////////////////////////

CFuint GeometricEntityRegister::regist(const std::string& providerName)
{
  CFuint idx = 0;

  if(exists(providerName)) {
    // find which is the index and return it
    DatabaseType::iterator itr = _database.begin();
    for(;itr != _database.end();++itr) {
      if((*itr)->getName() == providerName) return idx;
      ++idx;
    }
  }
  else {
    // else add it to the database and return the index

    BaseGeometricEntityProvider* provider =
    GeometricEntityFactory::getProvider(providerName);

    cf_assert(!exists(provider));

    idx = _database.size();
    _database.push_back(provider);


    // register the solution interpolator
    InterpolatorID solID =
    InterpolatorRegister::getInstance().
      registInterpolator(provider->getSolutionShapeFunctionName(),
                        provider->getSolutionShapeFunctionType(),
                        provider->getSolutionShapeFunctionOrder(),
                        provider->getShape());

    // register the geometry interpolator
    InterpolatorID geoID =
    InterpolatorRegister::getInstance().
      registInterpolator(provider->getGeometryShapeFunctionName(),
                        provider->getGeometryShapeFunctionType(),
                        provider->getGeometryShapeFunctionOrder(),
                        provider->getShape());

    provider->setSolInterpolatorID(solID);
    provider->setGeomInterpolatorID(geoID);
   }

   return idx;
}

//////////////////////////////////////////////////////////////////////////////

BaseGeometricEntityProvider* GeometricEntityRegister::getProvider(const CFuint idx)
{
  cf_assert(idx < _database.size());
  return _database[idx];
}

//////////////////////////////////////////////////////////////////////////////

BaseGeometricEntityProvider* GeometricEntityRegister::getProvider(const std::string& providerName)
{
  DatabaseType::iterator itr = _database.begin();
  for(;itr != _database.end();++itr) {
    if((*itr)->getName() == providerName) return (*itr);
  }
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

bool GeometricEntityRegister::exists(const std::string& providerName)
{
  DatabaseType::iterator itr = _database.begin();
  for(;itr != _database.end();++itr) {
    if((*itr)->getName() == providerName) return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool GeometricEntityRegister::exists(BaseGeometricEntityProvider* provider)
{
  DatabaseType::iterator itr = _database.begin();
  for(;itr != _database.end();++itr) {
    if(*itr == provider) return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

void GeometricEntityRegister::clear()
{
  _database.clear();
}

//////////////////////////////////////////////////////////////////////////////


  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

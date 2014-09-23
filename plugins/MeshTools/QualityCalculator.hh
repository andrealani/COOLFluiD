// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_QualityCalculator_hh
#define COOLFluiD_MeshTools_QualityCalculator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ConcreteProvider.hh"
#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a base class
 * for a calculator of quality of elements
 *
 * @author Thomas Wuilbaut
 *
 */

class QualityCalculator : public Common::OwnedObject,
                          public Config::ConfigObject {
public:

  typedef Environment::ConcreteProvider<QualityCalculator,1> PROVIDER;
  typedef const std::string& ARG1;

  /**
   * Constructor
   */
  QualityCalculator(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~QualityCalculator();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ConfigObject::configure(args);

    // add here configuration, specific of this class
  }

  /**
   * Calculate the quality of an element
   * @param geoEnt the geometric entity
   */
  virtual CFreal computeQuality(Framework::GeometricEntity* geoEnt) = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "QualityCalculator";
  }

protected: // data

  // the quality of the element
  CFreal _quality;

  // the element
  Framework::GeometricEntity* _geoEnt;

}; // end of class QualityCalculator

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_QualityCalculator_hh

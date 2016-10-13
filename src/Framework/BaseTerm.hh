// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BaseTerm_hh
#define COOLFluiD_Framework_BaseTerm_hh




//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Common/NonCopyable.hh"
#include "Config/ConfigObject.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a basic physical term holding
/// some array data.
/// @author Andrea Lani
class Framework_API BaseTerm : public Common::NonCopyable<BaseTerm>,
			       public Config::ConfigObject {
  
 public:
  
  /// enumerator 
  enum {END=0};
  
  /// Constructor without arguments
  BaseTerm(const std::string& name) :
    ConfigObject(name),
    m_physicalData(),
    m_refPhysicalData(),
    m_startVar(0),
    m_currNbEqs(0)
  {
  }

  /// Default destructor
  virtual ~BaseTerm()
  {
  }

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    ConfigObject::configure(args);
  }

  /// Set physical data
  virtual void setupPhysicalData() = 0;

  /// Get the array based physical data
  RealVector& getPhysicalData()
  {
    return m_physicalData;
  }

  /// Get the reference array based physical data
  RealVector& getReferencePhysicalData()
  {
    return m_refPhysicalData;
  }

  /// Resize physical data
  void resizePhysicalData(RealVector& physicalData)
  {
    cf_assert(getDataSize() > 0);
    physicalData.resize(getDataSize());

    // by default the current number of equations is set to the total number
    m_currNbEqs = PhysicalModelStack::getActive()->getNbEq();
  }

protected:

  /// Physical data size
  virtual CFuint getDataSize() const {return 0;}
  
protected:

  /// array data
  RealVector m_physicalData;

  /// array reference data
  RealVector m_refPhysicalData;

  /// start variable ID
  CFuint m_startVar;

  /// current number of equations
  CFuint m_currNbEqs;

}; // end of class BaseTerm

//////////////////////////////////////////////////////////////////////////////

  } // namespace BaseTerm

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BaseTerm_hh

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_InterpolatorRegister_hh
#define COOLFluiD_Framework_InterpolatorRegister_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/InterpolatorProperties.hh"
#include "Common/NoSuchValueException.hh"
#include "Common/Trio.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class that has knowledge about
/// which interpolators are being used in the current SubSystem.
/// This differs from a Factory of providers as the interpolators
/// registered in this class are already in use as the ones on the Factory
/// are just to be potencially used.
/// This class is a Singleton.
/// @author Tiago Quintino
class Framework_API InterpolatorRegister {
private: // typedefs

  /// entry type in database
  typedef std::pair<std::string,
                    Common::Trio<CFPolyForm::Type,CFPolyOrder::Type,CFGeoShape::Type> > entry_type;

  /// database type
  typedef std::vector< entry_type > database_type;

private: // helper classes

  class Framework_API EqualName {
  public:

    EqualName(const std::string& name) : _name(name)
    {
    }

    /// Operator to check is name is equal
    bool operator()(const entry_type& elem)
    {
      return (elem.first == _name);
    }

    private: // data

    /// a refence to the name being compared
    const std::string& _name;

  }; // end class EqualName

  class Framework_API EqualProperties {
  public:

    EqualProperties(const CFPolyForm::Type&  interpolType,
                    const CFPolyOrder::Type& interpolOrder,
                    const CFGeoShape::Type&  shape)
                    : _prop(Common::make_Trio(interpolType,
                                             interpolOrder,
                                             shape))
    {
    }

    /// Operator to check if the properties of the
    /// interpolator are equal
    bool operator()(const entry_type& elem)
    {
      return (elem.second == _prop);
    }

    private: // data

    /// the trio of properties to be found
    Common::Trio<CFPolyForm::Type,CFPolyOrder::Type,CFGeoShape::Type> _prop;

  }; // end class EqualProperties

public: // methods

  /// Gets the instance of the Singleton
  static InterpolatorRegister& getInstance()
  {
    static InterpolatorRegister instance;
    return instance;
  }

  /// Cleans the database of registered interpolators
  void clearRegister()
  {
    _database.clear();
  }

  /// Gets the number of interpolators registered for this SubSystem
  CFuint getNbInterpolators() const
  {
    return _database.size();
  }

  /// Gets the InterpolatorProperties given a certain InterpolatorID present in the register.
  InterpolatorProperties
  getInterpolatorProperties(const InterpolatorID& id) const throw (Common::NoSuchValueException);

  /// Consults the register to get the InterpolatorID of a certain interpolator
  /// given its name.
  InterpolatorID getInterpolatorID(const std::string& name) const throw (Common::NoSuchValueException);

  /// Consults the register to get the InterpolatorID of a certain interpolator
  /// given its properties.
  InterpolatorID getInterpolatorID(
                           const CFPolyForm::Type&  interpolType,
                           const CFPolyOrder::Type& interpolOrder,
                           const CFGeoShape::Type&  shape) const throw (Common::NoSuchValueException);

  /// Regists an Interpolator to be used in the SubSystem.
  InterpolatorID registInterpolator(const std::string&    name,
                                     const CFPolyForm::Type&  interpolType,
                                     const CFPolyOrder::Type& interpolOrder,
                                     const CFGeoShape::Type&  shape);

private: // methods

  /// Constructor.
  /// Private to make it a Singleton.
  InterpolatorRegister()
  {
  }

  /// Destructor
  /// Private to make it a Singleton.
  ~InterpolatorRegister()
  {
  }

  /// Copy Constructor.
  /// Private to make it a Singleton.
  InterpolatorRegister(const InterpolatorRegister&);

  /// Copy Operator.
  /// Private to make it a Singleton.
  InterpolatorRegister& operator=(const InterpolatorRegister&);

private: // data

  // the register database
  database_type _database;

}; // end class InterpolatorRegister

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_InterpolatorRegister_hh

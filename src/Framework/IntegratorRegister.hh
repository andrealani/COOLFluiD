// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IntegratorRegister_hh
#define COOLFluiD_Framework_IntegratorRegister_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NoSuchValueException.hh"
#include "Common/Trio.hh"
#include "Common/CFMap3D.hh"
#include "Common/NonCopyable.hh"

#include "Common/CFLog.hh"

#include "Framework/IntegratorProperties.hh"
#include "Framework/VolumeIntegratorImpl.hh"
#include "Framework/ContourIntegratorImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides a class that has knowledge about
/// which Integrator's are available to be used in the current SubSystem.
/// This differs from a Factory of Provider's as the information
/// stored is just about the functionality, and not actually the objects.
/// This class is a Singleton.
/// @author Tiago Quintino
template <typename INTEGRATORIMPL>
class IntegratorRegister : public Common::NonCopyable<IntegratorRegister<INTEGRATORIMPL> > {
private: // typedefs

  /// entry type in database
  typedef Common::Trio<std::string,CFQuadrature::Type,CFPolyOrder::Type> EntryType;

  /// Database type
  typedef Common::CFMap3D<CFGeoShape::Type,CFPolyForm::Type,CFPolyOrder::Type,EntryType> DatabaseType;

public: // methods

  /// Gets the instance of the Singleton
  static IntegratorRegister& getInstance()
  {
    static IntegratorRegister instance;
    return instance;
  }

  /// Gets the IntegratorProperties given a certain IntegratorID present in the register.
  std::vector<std::string>
  getIntegratorsMatching(const CFGeoShape::Type& shape,
                         const CFPolyForm::Type& interpolType,
                         const CFPolyOrder::Type& interpolOrder);

  /// Regists an Integrator to be used in the SubSystem.
  void registIntegrator(const std::string&       name,
                        const CFIntegration::Type& integratType,
                        const CFQuadrature::Type&  quadratureType,
                        const CFPolyOrder::Type&       integratOrder,
                        const CFGeoShape::Type&        shape,
                        const CFPolyForm::Type&        interpolType,
                        const CFPolyOrder::Type&       interpolOrder);

  /// Access to the database for external quieries
  DatabaseType& getIntegratorDB () { return m_intdb; }

private: // methods

  /// Constructor.
  /// Private to make it a Singleton.
  IntegratorRegister() {}

  /// Destructor
  /// Private to make it a Singleton.
  ~IntegratorRegister() {}

private: // data

  // the register database
  DatabaseType m_intdb;

}; // end class IntegratorRegister

//////////////////////////////////////////////////////////////////////////////

template <typename INTEGRATORIMPL>
void IntegratorRegister<INTEGRATORIMPL>::registIntegrator(
                        const std::string&          name,
                        const CFIntegration::Type& integratType,
                        const CFQuadrature::Type&  quadratureType,
                        const CFPolyOrder::Type&       integratOrder,
                        const CFGeoShape::Type&        shape,
                        const CFPolyForm::Type&        interpolType,
                        const CFPolyOrder::Type&        interpolOrder)
{
  INTEGRATORIMPL::assertIntegratorType(integratType);
  m_intdb.insert(shape,interpolType,interpolOrder,Common::make_Trio(name,quadratureType,integratOrder));
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTEGRATORIMPL>
std::vector<std::string>
IntegratorRegister<INTEGRATORIMPL>::getIntegratorsMatching(
                       const CFGeoShape::Type& shape,
                       const CFPolyForm::Type& interpolType,
                       const CFPolyOrder::Type& interpolOrder)
{
  CFAUTOTRACE;

  typedef typename DatabaseType::MapIterator IteratorType;

  typedef typename std::pair<IteratorType,IteratorType> ResultType;

  ResultType searchResult = m_intdb.findBounds(shape,interpolType,interpolOrder);

  // list of names of IntegratorImpl's to return
  std::vector<std::string> result;

  IteratorType itr    = searchResult.first;
  IteratorType endItr = searchResult.second;
  for(; itr != endItr; ++itr) {
    result.push_back(itr->fourth.first);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IntegratorRegister_hh

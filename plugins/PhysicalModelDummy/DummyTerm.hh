// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_PhysicalModelDummy_DummyTerm_hh
#define COOLFluiD_PhysicalModelDummy_DummyTerm_hh

#include "Framework/BaseTerm.hh"

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a DummyTerm
class DummyTerm : public Framework::BaseTerm {

public:

  /// Constructor without arguments
  DummyTerm(const std::string& name);

  /// Default destructor
  virtual ~DummyTerm();

  /// Set physical data
  virtual void setupPhysicalData();

  /// Physical data size
  virtual CFuint getDataSize() const
  {
    /*
    return Framework::PhysicalModelStack::getActive()->getImplementor()->
      getNbEquations();
    */
    return m_var;
  }

  /// Set the variables names (number of equations per state)
  void setVarNames(std::vector< std::string >& varnames)
  {
    m_varnames = varnames;
    m_var      = varnames.size();
  }

  /// Get the variables names (size equal to number of equations per state)
  const std::vector< std::string >& getVarNames() const
  {
    return m_varnames;
  }

  /// Get the name
  static std::string getName()
  {
    return "DummyTerm";
  }


private:

  /// Number of equations per State
  CFuint m_var;

  /// Vector of equations' names
  std::vector< std::string > m_varnames;

}; // end of class DummyTerm

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

#endif // COOLFluiD_PhysicalModelDummy_DummyTerm_hh


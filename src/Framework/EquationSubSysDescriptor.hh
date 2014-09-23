// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_EquationSubSysDescriptor_hh
#define COOLFluiD_Framework_EquationSubSysDescriptor_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"

#include "MathTools/RealVector.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class holds some data describing the current equation subsystem
/// to be solved
/// @author Andrea Lani
class Framework_API EquationSubSysDescriptor : public Common::NonCopyable<EquationSubSysDescriptor> {
public:

  /// Constructor
  EquationSubSysDescriptor() :
    _iStartVar(0),
    _nbEqsSubSys(0),
    _iEqSubSys(0),
    _totalNbEqSS(1),
    _eqVarPatterns()
  {
  }

  /// Default destructor
  ~EquationSubSysDescriptor()
  {
  }

  /// Set the total number of equations subsystems
  void setTotalNbEqSS(CFuint totalNbEqSS)
  {
    _totalNbEqSS = totalNbEqSS;
  }

  /// Set some data describing the current equation subsystem to solve
  /// This is meant to be particularly useful if weak coupling is used
  /// @param iStartVar    first variable ID for the current equation subsystem
  /// @param nbEqsSubSys  number of equations of the current equation subsystem
  /// @param iEqSubSys    ID of the current equation subsystem
  void set(CFuint iStartVar, CFuint nbEqsSubSys, CFuint iEqSubSys)
  {
    _iStartVar   = iStartVar;
    _nbEqsSubSys = nbEqsSubSys;
    _iEqSubSys   = iEqSubSys;
  }

  /// Reset to default values some data describing the current equation
  /// subsystem to solve.
  /// This is meant to be particularly useful if weak coupling is used
  /// @param total number of equations
  void reset(CFuint totalNbEqs)
  {
    _iStartVar = _iEqSubSys = 0;
    _nbEqsSubSys = totalNbEqs;
  }

  /// Set the equation variable pattern for each equation subsystem
  void setEquationVarPatterns
  (const std::vector<std::vector<CFuint> >& eqVarPatterns)
  {
    _eqVarPatterns.resize(eqVarPatterns.size());
    for (CFuint i = 0; i < _eqVarPatterns.size(); ++i) {
      _eqVarPatterns[i].resize(eqVarPatterns[i].size());
      for (CFuint j = 0; j < _eqVarPatterns[i].size(); ++j) {
	_eqVarPatterns[i][j] = eqVarPatterns[i][j];
      }
    }
  }

  /// Get the ID of the first variable of the equation subsystem
  CFuint getStartVarSS() const {return _iStartVar;}

  /// Get the number of equations of the current equation subsystem
  CFuint getNbEqsSS() const {return _nbEqsSubSys;}

  /// Get the ID of the current equation subsystem
  CFuint getEqSS() const {return _iEqSubSys;}

  /// Get the total number of equation subsystems
  CFuint getTotalNbEqSS() const {return _totalNbEqSS;}

  /// Get the equation variable pattern for each equation subsystem
  const std::vector<std::vector<CFuint> >& getEqVarPatterns() const
  {
    return _eqVarPatterns;
  }

private:

  /// first variable ID for the current equation subsystem
  CFuint _iStartVar;

  // number of equations of the current equation subsystem
  CFuint _nbEqsSubSys;

  // ID of the current equation subsystem
  CFuint _iEqSubSys;

  // total number of equation subsystems
  CFuint _totalNbEqSS;

  // equation variable pattern for each equation subsystem
  std::vector<std::vector<CFuint> > _eqVarPatterns;

}; // end class EquationSubSysDescriptor

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_EquationSubSysDescriptor_hh

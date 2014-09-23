// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_EquationSetData_hh
#define COOLFluiD_Framework_EquationSetData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Framework.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class holds some data describing the current equation set
/// Each sub set of equations in a @see ConvectiveVarSet should ideally
/// be associated with an EquationSetData
/// @author Andrea Lani
class Framework_API EquationSetData {
public:

  /// Constructor
  EquationSetData() :
    _eqSetID(0),
    _startEqSetID(0),
    _eqSetVarIDs()
  {
  }

  /// Default destructor
  ~EquationSetData()
  {
  }

  /// Get the ID of the equation set
  CFuint getEqSetID() const {return _eqSetID;}

  /// Get the starting ID of the equation set
  CFuint getStartEqID() const {return _startEqSetID;}

  /// Get the array with the IDs of the variables belong to this
  /// equation set
  const std::vector<CFuint>&  getEqSetVarIDs() const
  {
    return _eqSetVarIDs;
  }

  /// Set up the equation set data
  void setup(CFuint eqSetID, CFuint startEqSetID, CFuint nbEqs)
  {
    CFLog(VERBOSE, "EquationSetData::setup() => (" << eqSetID << ", " << startEqSetID << ", " << nbEqs << "\n");
    
    _eqSetID = eqSetID;
    _startEqSetID = startEqSetID;
    if (nbEqs > 0) _eqSetVarIDs.resize(nbEqs);
    // set the IDs of the variables
    for (CFuint i = 0; i < nbEqs; ++i) {
      _eqSetVarIDs[i] = _startEqSetID + i;
    }
  }
  
  /// get the siae of this equation set data
  CFuint size() const {return _eqSetVarIDs.size();}
  
private:

  /// ID of the corresponding subset of equations
  CFuint _eqSetID;

  /// start ID of the corresponding subset of equations
  CFuint _startEqSetID;

  /// list of the IDs of the variables corresponding to the current
  /// subset of equations
  std::vector<CFuint> _eqSetVarIDs;

}; // end class EquationSetData

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_EquationSetData_hh

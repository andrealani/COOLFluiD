// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_BDF2Intermediate_hh
#define COOLFluiD_Numerics_NewtonMethod_BDF2Intermediate_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.

class BDF2Intermediate : public NewtonIteratorCom {
public:

  /// Constructor.
  explicit BDF2Intermediate(std::string name) : NewtonIteratorCom(name)
  {
  }

  /// Destructor.
  ~BDF2Intermediate()
  {
  }

  /// Execute Processing actions
  void execute();

}; // class BDF2Intermediate

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_BDF2Intermediate_hh


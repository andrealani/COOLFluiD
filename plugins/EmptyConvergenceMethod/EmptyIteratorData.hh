// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptyConvergenceMethod_EmptyIteratorData_hh
#define COOLFluiD_EmptyConvergenceMethod_EmptyIteratorData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/MethodCommand.hh"
#include "MathTools/RealVector.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethodData.hh"

#include "EmptyConvergenceMethod/EmptyConvergenceMethodAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace EmptyConvergenceMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a Data Object that is accessed by the different
  /// EmptyConvergenceMethodCom 's that compose the EmptyConvergenceMethod.
  /// @see EmptyConvergenceMethodCom
  /// @author Andrea Lani
class  EmptyConvergenceMethod_API EmptyIteratorData : public Framework::ConvergenceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  // static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  EmptyIteratorData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~EmptyIteratorData();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "EmptyIterator";
  }
  
}; // end of class EmptyIteratorData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for EmptyConvergenceMethod
typedef Framework::MethodCommand<EmptyIteratorData> EmptyIteratorCom;

/// Definition of a command provider for EmptyConvergenceMethod
typedef Framework::MethodCommand<EmptyIteratorData>::PROVIDER EmptyIteratorComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace EmptyConvergenceMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_EmptyConvergenceMethod_EmptyIteratorData_hh

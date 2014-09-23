// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOFluiD_Common_ProviderBase_hh
#define COOFluiD_Common_ProviderBase_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NonCopyable.hh"

#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// @brief Base class for provider types
/// @author Dries Kimpe
/// @author Tiago Quintino
class Common_API ProviderBase : public Common::NonCopyable<ProviderBase> {

public: // methods

  /// Constructor
  ProviderBase ();

  /// Virtual destructor
  virtual ~ProviderBase ();

  /// Free an instance created by this factory
  /// @param ptr pointer to be freed
  virtual void freeInstance ( void * ptr ) = 0;

  /// @return the name of this provider
  virtual std::string getProviderName () const = 0;

  /// @return the type of this provider
  virtual std::string getProviderType () const = 0;

}; // end ProviderBase

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Common_ProviderBase_hh

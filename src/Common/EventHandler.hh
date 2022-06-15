// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_EventHandler_hh
#define COOLFluiD_Common_EventHandler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Common/DynamicObject.hh"
#include "Common/OwnedObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Global Event Handler class
/// @author Tiago Quintino
class Common_API EventHandler :
    public Common::OwnedObject,
    public Common::DynamicObject,
    public Common::NonCopyable<EventHandler>
{
public: // methods

  /// Constructor private because is singleton
  EventHandler();

  /// Destructor private because is singleton
  ~EventHandler();

  /// Regists a signal on this EventHandler
  template < typename PTYPE, typename FTYPE >
  void addListener ( const std::string& sname, PTYPE* ptr, FTYPE pfunc, const std::string& desc = "" )
  {
#ifdef CF_HAVE_BOOST_1_76
    regist_signal ( sname , desc )->connect ( boost::bind ( pfunc, ptr, std::placeholders::_1 ) );
#else
    regist_signal ( sname , desc )->connect ( boost::bind ( pfunc, ptr, _1 ) ); 
#endif
  }

}; // class EventHandler

//////////////////////////////////////////////////////////////////////////////

} // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_EventHandler_hh

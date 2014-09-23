// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_EnumT_hh
#define COOLFluiD_Common_EnumT_hh

#include "Common/COOLFluiD.hh"
#include "Common/NonInstantiable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Common {

//////////////////////////////////////////////////////////////////////////////

template < typename EClass >
class EnumT : public Common::NonInstantiable< EnumT<EClass> >
{
  public: // typedefs

    typedef typename EClass::Type EnumType;

  protected:
    /// type for forward map
    typedef std::map < EnumType, std::string > FwdMap_t;
    /// type for reverse map
    typedef std::map < std::string, EnumType > BwdMap_t;

  public:

   /// conversion from enum to string
   static bool is_valid ( const EnumType& in )
   {
    FwdMap_t& all_fwd = EClass::Convert::all_fwd;
    return !( in == EClass::INVALID || all_fwd.find(in) == all_fwd.end() );
   }

   /// conversion from enum to string
  static std::string to_str ( const EnumType& in )
  {
    FwdMap_t& all_fwd = EClass::Convert::all_fwd;
    typename FwdMap_t::const_iterator fitr = all_fwd.find(in);
    return ( fitr == all_fwd.end() ) ? all_fwd[EClass::INVALID] : fitr->second;
  }

  /// conversion from int to enum then to string
  static std::string to_str ( const int& in )
  {
    FwdMap_t& all_fwd = EClass::Convert::all_fwd;
    EnumType ft = static_cast< EnumType >( in );
    typename FwdMap_t::const_iterator fitr = all_fwd.find(ft);
    return ( fitr == all_fwd.end() ) ? all_fwd[EClass::INVALID] : fitr->second;
  }

  /// conversion from int to enum (with validity check)
  static EnumType to_enum ( const int& in )
  {
    FwdMap_t& all_fwd = EClass::Convert::all_fwd;
    EnumType ft = static_cast< EnumType >( in );
    typename FwdMap_t::const_iterator fitr = all_fwd.find(ft);
    return ( fitr == all_fwd.end() ) ? EClass::INVALID : ft;
  }

  /// conversion from string to enum
  static EnumType to_enum ( const std::string& in )
  {
    BwdMap_t& all_rev = EClass::Convert::all_rev;
    typename BwdMap_t::const_iterator fitr = all_rev.find(in);
    return ( fitr == all_rev.end() ) ? EClass::INVALID : fitr->second;
  }

}; // class EnumT

//////////////////////////////////////////////////////////////////////////////

} // Common
} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_EnumT_hh

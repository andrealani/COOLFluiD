// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFGeoEnt_hh
#define COOLFluiD_CFGeoEnt_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/EnumT.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class Framework_API CFGeoEnt
{
  public:

  /// Enumeration of the types of Geometric Entities
  enum Type  { INVALID =  0,
               CELL    =  1,
               FACE    =  2,
               EDGE    =  3  };

  typedef Common::EnumT< CFGeoEnt > ConverterBase;

  struct Framework_API Convert : public ConverterBase
  {
    /// storage of the enum forward map
    static ConverterBase::FwdMap_t all_fwd;
    /// storage of the enum reverse map
    static ConverterBase::BwdMap_t all_rev;
  };

}; // class CFGeoEnt

///////////////////////////////////////////////////////////////////////////

Framework_API std::ostream& operator<< ( std::ostream& os, const CFGeoEnt::Type& in );
Framework_API std::istream& operator>> ( std::istream& is, CFGeoEnt::Type& in );

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFGeoEnt_hh

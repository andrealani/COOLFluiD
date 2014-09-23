// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFGeoShape_hh
#define COOLFluiD_CFGeoShape_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/EnumT.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class Framework_API CFGeoShape
{
  public:

  /// Enumeration of the Shapes recognized in COOLFluiD
  enum Type  { INVALID = -1,
               POINT   = 0,
               LINE    = 1,
               TRIAG   = 2,
               QUAD    = 3,
               TETRA   = 4,
               PYRAM   = 5,
               PRISM   = 6,
               HEXA    = 7  };

  typedef Common::EnumT< CFGeoShape > ConverterBase;

  struct Framework_API Convert : public ConverterBase
  {
    /// storage of the enum forward map
    static ConverterBase::FwdMap_t all_fwd;
    /// storage of the enum reverse map
    static ConverterBase::BwdMap_t all_rev;
  };

}; // class CFGeoShape

///////////////////////////////////////////////////////////////////////////

Framework_API std::ostream& operator<< ( std::ostream& os, const CFGeoShape::Type& in );
Framework_API std::istream& operator>> ( std::istream& is, CFGeoShape::Type& in );

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFGeoShape_hh

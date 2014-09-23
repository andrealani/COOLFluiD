// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFIntegration_hh
#define COOLFluiD_CFIntegration_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/EnumT.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class Framework_API CFIntegration
{
  public:

  /// Enumeration of the integration types recognized in COOLFluiD
  enum Type   { INVALID,
                CONTOUR,
                VOLUME    };

  typedef Common::EnumT< CFIntegration > ConverterBase;

  struct Framework_API Convert : public ConverterBase
  {
    /// storage of the enum forward map
    static ConverterBase::FwdMap_t all_fwd;
    /// storage of the enum reverse map
    static ConverterBase::BwdMap_t all_rev;
  };

}; // class CFIntegration

///////////////////////////////////////////////////////////////////////////

Framework_API std::ostream& operator<< ( std::ostream& os, const CFIntegration::Type& in );
Framework_API std::istream& operator>> ( std::istream& is, CFIntegration::Type& in );

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFIntegration_hh

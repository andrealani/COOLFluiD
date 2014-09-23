// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFQuadrature_hh
#define COOLFluiD_CFQuadrature_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/EnumT.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class Framework_API CFQuadrature
{
  public:

  /// Enumeration of the integration types recognized in COOLFluiD
  enum Type  { INVALID,
               GAUSSLEGENDRE,
               GAUSSLOBATO,
               NEWTONCOTES,
               SIMPSON,
               TRAPEZIUM,
               DGGAUSSLEGENDRE };

  typedef Common::EnumT< CFQuadrature > ConverterBase;

  struct Framework_API Convert : public ConverterBase
  {
    /// storage of the enum forward map
    static ConverterBase::FwdMap_t all_fwd;
    /// storage of the enum reverse map
    static ConverterBase::BwdMap_t all_rev;
  };

}; // class CFQuadrature

///////////////////////////////////////////////////////////////////////////

Framework_API std::ostream& operator<< ( std::ostream& os, const CFQuadrature::Type& in );
Framework_API std::istream& operator>> ( std::istream& is, CFQuadrature::Type& in );

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFQuadrature_hh

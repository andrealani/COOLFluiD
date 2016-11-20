// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFPolyForm_hh
#define COOLFluiD_CFPolyForm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/EnumT.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class Framework_API CFPolyForm
{
  public:

  /// Enumeration of the polynomial fornms supported
  enum Type  { INVALID            = 0,
               LAGRANGE           = 1,
               HERMITE            = 2,
               SERENDIPITY        = 3,
               SPECTRALFV         = 4,
               LEGENDRE           = 5,
               SPECTRALFD         = 6,
	       FLUXRECONSTRUCTION = 7
  };

  typedef Common::EnumT< CFPolyForm > ConverterBase;

  struct Framework_API Convert : public ConverterBase
  {
    /// storage of the enum forward map
    static ConverterBase::FwdMap_t all_fwd;
    /// storage of the enum reverse map
    static ConverterBase::BwdMap_t all_rev;
  };

}; // class CFPolyForm

///////////////////////////////////////////////////////////////////////////

Framework_API std::ostream& operator<< ( std::ostream& os, const CFPolyForm::Type& in );
Framework_API std::istream& operator>> ( std::istream& is, CFPolyForm::Type& in );

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFPolyForm_hh

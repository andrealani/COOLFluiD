// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFPolyOrder_hh
#define COOLFluiD_CFPolyOrder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/EnumT.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class Framework_API CFPolyOrder
{
  public:

  /// Enumeration of the polynomial orders recognized in COOLFluiD
  /// CFPolyOrder::MAXORDER stands for the last one, and should not be used except
  /// for checking boundness of the order.
  enum  Type { INVALID =-1,
               ORDER0  = 0,
               ORDER1  = 1,
               ORDER2  = 2,
               ORDER3  = 3,
               ORDER4  = 4,
               ORDER5  = 5,
               ORDER6  = 6,
               ORDER7  = 7,
               ORDER8  = 8,
               ORDER9  = 9,
               ORDER10 = 10,
               MAXORDER };

  typedef Common::EnumT< CFPolyOrder > ConverterBase;

  struct Framework_API Convert : public ConverterBase
  {
    /// storage of the enum forward map
    static ConverterBase::FwdMap_t all_fwd;
    /// storage of the enum reverse map
    static ConverterBase::BwdMap_t all_rev;
  };

}; // class CFPolyOrder

///////////////////////////////////////////////////////////////////////////

Framework_API std::ostream& operator<< ( std::ostream& os, const CFPolyOrder::Type& in );
Framework_API std::istream& operator>> ( std::istream& is, CFPolyOrder::Type& in );

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFPolyOrder_hh

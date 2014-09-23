// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_OptionMarkers_hh
#define COOLFluiD_Config_OptionMarkers_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/Config.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

struct Config_API EndMarkOption
{
  static void mark (Option * opt) {}
};

//////////////////////////////////////////////////////////////////////////////

template < typename ValidationType, typename NEXT = EndMarkOption >
struct ValidateOption
{
  static void mark (Option * opt)
  {
    opt->addValidation( new ValidationType(opt) );
    NEXT::mark(opt);
  }
};

//////////////////////////////////////////////////////////////////////////////

template < typename NEXT = EndMarkOption >
struct DynamicOption
{
  static void mark (Option * opt)
  {
    opt->makeOptionDynamic();
    NEXT::mark(opt);
  }
};

//////////////////////////////////////////////////////////////////////////////

template < typename NEXT = EndMarkOption >
struct BasicOption
{
  static void mark (Option * opt)
  {
    opt->makeOptionBasic();
    NEXT::mark(opt);
  }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_OptionMarkers_hh

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_PM_hh
#define COOLFluiD_Common_PM_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NonInstantiable.hh"

//////////////////////////////////////////////////////////////////////////////

/// This header defines and selects the correct parallelisation mode.
/// The PM_MPI / PM_SHM modes should never be used directly.
/// Instead, PM_CUR should always be used. This ensures that the
/// optimal parallel model will be used depending on the compilation mode.

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD  {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

class Common_API PM_MPI  : public Common::NonInstantiable<PM_MPI> {};

class Common_API PM_SHM  : public Common::NonInstantiable<PM_SHM> {};

class Common_API PM_SERIAL : public Common::NonInstantiable<PM_SERIAL> {};

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
  typedef PM_MPI    PM_CUR;
#else
  typedef PM_SERIAL   PM_CUR;
#endif // CF_HAVE_MPI

//////////////////////////////////////////////////////////////////////////////

    } // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_PM_hh

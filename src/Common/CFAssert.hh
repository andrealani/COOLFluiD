// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CFASSERT_HH
#define COOLFluiD_Common_CFASSERT_HH

/// @note This header should be included by including COOLFluiD.hh instead.
#include <cassert>

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

/// Manager of behavior of assertions
class Common_API AssertionManager : public Common::NonCopyable <AssertionManager> {
public:

  /// Constructor
  AssertionManager();

  /// Gets the instance of the manager
  static AssertionManager& getInstance ();

  /// If AssertionManager is not handling assertions and those are
  /// passed to the standard assert function
  /// Controlled by the build option CF_ENABLE_STDASSERT
  static bool notHandlingAssertions ()
  {
    #ifdef CF_ENABLE_STDASSERT
      return true;
    #else
      return false;
    #endif
  }

  /// Forward declaration of the function that inplements the always present assert
  static void do_assert (  bool condition,
                           const char * cond_str,
                           const char * file,
                           int line,
                           const char * func,
                           const char * desc = 0 );

  /// flag to  dynamically turn off assertions
  bool DoAssertions;
  /// assertions dump backtraces
  bool AssertionDumps;
  /// assertions throw exceptions
  bool AssertionThrows;

}; // class AssertionManager

//////////////////////////////////////////////////////////////////////////////

#ifndef CF_ENABLE_STDASSERT

/* using coolfluid assertion manager */

/// Assertion that always should be checked, even in optimized builds.
/// Use this for testing in non-speed-critical code.
#define cf_always_assert(a) \
     { if (!(a)) { ::COOLFluiD::AssertionManager::do_assert((a), #a, __FILE__, __LINE__, __FUNCTION__); } }

/// Assertion that always should be checked, even in optimized builds.
/// Use this for testing in non-speed-critical code.
/// Adds a description of the assertion
#define cf_always_assert_desc(msg,a) \
     { if (!(a)) { ::COOLFluiD::AssertionManager::do_assert((a), #a, __FILE__, __LINE__, __FUNCTION__, msg); } }

/// COOLFluiD assertion
/// Assertions are off if compiled with DNDEBUG
#ifndef NDEBUG
  #define cf_assert(a)        cf_always_assert((a))
  #define cf_assert_desc(m,a) cf_always_assert_desc(m,(a))
#else
  #define cf_assert(a)
  #define cf_assert_desc(m,a)
#endif

#else // CF_ENABLE_STDASSERT

/* using standard assertions */

/// COOLFluiD assertion
/// Assertions are off if compiled with DNDEBUG
#ifndef NDEBUG
#define   cf_assert(a)               assert(a)
#define cf_assert_desc(m,a)          assert(a)
#define cf_always_assert(a)          assert(a)
#define cf_always_assert_desc(msg,a) assert(a)
#else
  #define cf_assert(a)
  #define cf_assert_desc(m,a)
  #define cf_always_assert(a)
  #define cf_always_assert_desc(msg,a)
#endif

#endif // CF_ENABLE_STDASSERT

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CFASSERT_HH

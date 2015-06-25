// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GlobalReduceNone_hh
#define COOLFluiD_Framework_GlobalReduceNone_hh

#include "Common/PM.hh"
#include "Common/NonCopyable.hh"
#include "Framework/GlobalReduce.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Global reduce for PE SERIAL.
/// Resolves into returning GetLocalValue() for a call to GetGlobalValue()
template <typename PROVIDER>
class GlobalReduce<PROVIDER,Common::PM_SERIAL> : public Common::NonCopyable<GlobalReduce<PROVIDER,Common::PM_SERIAL> >
{
public:

  typedef typename PROVIDER::GR_RESULTTYPE RESULTTYPE;

  GlobalReduce (const PROVIDER & P);
  ~GlobalReduce ();
  RESULTTYPE GetGlobalValue () const;

private:
  const PROVIDER & _provider;
};

//////////////////////////////////////////////////////////////////////////////

template <typename PROVIDER>
inline GlobalReduce<PROVIDER,Common::PM_SERIAL>::GlobalReduce (const PROVIDER & S)
  : _provider(S)
{
}

template <typename PROVIDER>
inline GlobalReduce<PROVIDER,Common::PM_SERIAL>::~GlobalReduce ()
{
}

template <typename PROVIDER>
typename GlobalReduce<PROVIDER,Common::PM_SERIAL>::RESULTTYPE
inline GlobalReduce<PROVIDER,Common::PM_SERIAL>::GetGlobalValue () const
{
  return _provider.GR_GetLocalValue ();
}

//////////////////////////////////////////////////////////////////////////////

/// GlobalReduceOperation for PE PM_SERIAL:
///   for all operations, just copy Source to Dest
template <typename TAGCLASS, typename BASETYPE>
class GlobalReduceOperationHelper<TAGCLASS, BASETYPE, Common::PM_SERIAL>
{
public:
    static inline void DoReduce (BASETYPE * Source, BASETYPE * Dest,
      unsigned int Count)
    {
  for (unsigned int i=0; i<Count; ++i)
      *Dest++ = *Source++;
    }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

#endif // COOLFluiD_Framework_GlobalReduceNone_hh

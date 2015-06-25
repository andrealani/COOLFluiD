// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

// initialize the static data
PEInterface<>* PE::m_curr_PE = CFNULL;
bool           PE::m_is_init = false;
    
//////////////////////////////////////////////////////////////////////////////

bool PE::IsInitialised ()
{
  return m_is_init;
}

//////////////////////////////////////////////////////////////////////////////

void PE::InitPE (int* argc, char*** args)
{
  cf_assert (m_curr_PE == CFNULL);
  
  m_curr_PE = new PEInterface<> (argc, args);
  cf_assert (m_curr_PE != CFNULL);
  m_is_init = true;
}

//////////////////////////////////////////////////////////////////////////////

PEInterface<>& PE::GetPE ()
{
  cf_assert(m_is_init);
  cf_assert (m_curr_PE != CFNULL);
  return *m_curr_PE;
}

//////////////////////////////////////////////////////////////////////////////

void PE::DonePE ()
{
  cf_assert(m_curr_PE != CFNULL);

  // must be first to make sure all destructors dependent of PE acknowledge MPI is down
  m_is_init = false;
  deletePtr(m_curr_PE);
}
    
//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD


// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PEInterface.hh"

namespace COOLFluiD {
   namespace Common  {

//////////////////////////////////////////////////////////////////////////////

PEInterfaceBase::~PEInterfaceBase () {}

unsigned int
PEInterfaceBase::GetProcessorCount () const
{ throw Common::NotImplementedException (FromHere(),"GetProcessorCount()"); return 0; }

unsigned int
PEInterfaceBase::GetRank () const
{ throw Common::NotImplementedException (FromHere(),"PEInterface()::GetRank()"); return 0; }

void
PEInterfaceBase::setBarrier()
{ throw Common::NotImplementedException (FromHere(),"PEInterface()::setBarrier()"); }

bool
PEInterfaceBase::IsParallel () const
{ throw Common::NotImplementedException (FromHere(),"PEInterface()::IsParallel()"); return false; }

std::string
PEInterfaceBase::GetName () const
{ throw Common::NotImplementedException (FromHere(),"PEInterface::Getname");  return ""; }

bool
PEInterfaceBase::IsParallelCapable () const
{ throw Common::NotImplementedException (FromHere(),"IsParallelCapable"); return false;}

void
PEInterfaceBase::AdvanceCommunication ()
{ throw Common::NotImplementedException (FromHere(),"AdvanceCommunication"); }

//////////////////////////////////////////////////////////////////////////////

   }
}

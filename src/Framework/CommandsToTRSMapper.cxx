// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NullPointerException.hh"
#include "Common/MemFunArg.hh"
#include "Common/CFLog.hh"

#include "Framework/CommandsToTRSMapper.hh"
#include "Framework/NumericalCommand.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Framework/TrsNotFoundException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

CommandsToTRSMapper::CommandsToTRSMapper(
                      std::vector< Common::SafePtr<NumericalCommand> >& coms,
                      std::vector< Common::SafePtr<TopologicalRegionSet> >& trs) :
  m_coms(coms), m_trs(trs)
{
}

//////////////////////////////////////////////////////////////////////////////

CommandsToTRSMapper::~CommandsToTRSMapper()
{
}

//////////////////////////////////////////////////////////////////////////////

void CommandsToTRSMapper::mapComsToTrs() const
{
  // process all the commands by placing the TRS pointers in them
  std::for_each( m_coms.begin(),
                 m_coms.end(),
                 mem_fun_arg(*this,(&CommandsToTRSMapper::processCommand)) );
}

//////////////////////////////////////////////////////////////////////////////

void CommandsToTRSMapper::processCommand( Common::SafePtr<NumericalCommand> comPtr ) const
{
  vector< SafePtr<TopologicalRegionSet> > trsInCom;

  // safety check
  if ( comPtr.isNull() )
  {
    std::string msg;
    msg += "Pointer to a NumericalCommand\n";
    msg += "Did you forget to configure it !??";
    throw Common::NullPointerException (FromHere(), msg);
  }

  // get the names of the TRS that the command applies to
  const vector<std::string>& comTrsNames = comPtr->getTrsNames();

  trsInCom.reserve(comTrsNames.size());

  // loop over the names
  vector<std::string>::const_iterator trsNameItr;
  for (trsNameItr = comTrsNames.begin(); trsNameItr != comTrsNames.end(); ++trsNameItr) {

    // find the TRS with name
    bool nameFound = false;

    vector< SafePtr<TopologicalRegionSet> >::iterator trsItr;

    // loop over all TRS to check the name
    for (trsItr = m_trs.begin(); trsItr != m_trs.end(); ++trsItr)
    {
      std::string name = (*trsItr)->getName();

      CFLogDebugMin( "TRS name to find = " << name << "\n");

      if (name == *trsNameItr)
      {
        trsInCom.push_back(*trsItr);
        nameFound = true;
        break; // process the next name
      }
    }

    // if Trs name does not exist throw an exception
    if ( !nameFound ) {
      std::string msg;
      msg += "In command " + comPtr->getName() + " ";
      msg += "the TRS with name ";
      msg += *trsNameItr;
      msg += " was not found. Is the TRS";
      msg += " name in CFcase the same as in CFmesh or THOR SP file??\n";
      throw TrsNotFoundException (FromHere(),msg);
    }

  }
  comPtr->setTrsList(trsInCom);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD


// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_ParJFSolveSysGMRESR_hh
#define COOLFluiD_Numerics_Petsc_ParJFSolveSysGMRESR_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/StdParSolveSys.hh" // must come before any header

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a standard command to be executed during
  * solution of the linear system using Petsc library
  *
  * @author Jiri Simonek
  *
  */
class ParJFSolveSysGMRESR : public StdParSolveSys {
public:
	
  /**
   * Constructor.
   */
	explicit ParJFSolveSysGMRESR(const std::string& name);
	
  /**
   * Destructor.
   */
	~ParJFSolveSysGMRESR();
	
  /**
   * Execute Processing actions
   */
	void execute();
	
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
	void setup();
	
}; // class SolveSys

//////////////////////////////////////////////////////////////////////////////

} // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_ParJFSolveSysGMRESR_hh

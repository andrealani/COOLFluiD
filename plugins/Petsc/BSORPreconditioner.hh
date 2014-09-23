// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_BSORPreconditioner_hh
#define COOLFluiD_Numerics_Petsc_BSORPreconditioner_hh

//////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "Petsc/ShellPreconditioner.hh"
#include "Petsc/BSORPcContext.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a shell preconditioner object
 *
 * @author Jiri Simonek
 *
 */

class BSORPreconditioner : public ShellPreconditioner {
public:
 
  /**
   * Constructor
   */
  BSORPreconditioner(const std::string& name);

  /**
   * Default destructor
   */
  ~BSORPreconditioner();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set the preconditioner
   */
  virtual void setPreconditioner();
	
	virtual void computeBeforeSolving() {}
	virtual void computeAfterSolving() {}
	
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// BSOR context 
  BSORPcContext _pcc;

  /// Type of Gauss-Seidel relaxation
  //MatSORType _relaxType;

}; // end of class BSORPreconditioner
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_BSORPreconditioner_hh

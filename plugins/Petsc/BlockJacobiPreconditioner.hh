// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_BlockJacobiPreconditioner_hh
#define COOLFluiD_Numerics_Petsc_BlockJacobiPreconditioner_hh

//////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "Petsc/ShellPreconditioner.hh"
#include "Petsc/BlockJacobiPcJFContext.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace MathTools { 
    class MatrixInverter; 
  }
  
  namespace Petsc {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a shell preconditioner object
 *
 * @author Jiri Simonek
 *
 */

class BlockJacobiPreconditioner : public ShellPreconditioner {
public:

 /**
  * Constructor
  */
  BlockJacobiPreconditioner(const std::string& name);
  
  /**
   * Default destructor
   */
  ~BlockJacobiPreconditioner();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set the preconditioner
   */
  virtual void setPreconditioner();
  
  /**
   * Compute before solving the system
   */
  virtual void computeBeforeSolving();
  
  /**
   * Compute after solving the system
   */
  virtual void computeAfterSolving();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private:
  
  /// socket for the states
  Framework::DataSocketSink<CFreal> socket_diagMatrices;
  
  /// storage of the local updatable IDs or -1 (ghost) for all local states
  Framework::DataSocketSink<CFint> socket_upLocalIDsAll;
  
  /// BlockJacobi context 
  BlockJacobiPcJFContext _pcc;
  
  /// temporary data for holding the matrix inverter
  std::auto_ptr<MathTools::MatrixInverter> _inverter;
  
}; // end of class BlockJacobiPreconditioner
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_BlockJacobiPreconditioner_hh

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_TridiagPreconditioner_hh
#define COOLFluiD_Numerics_Petsc_TridiagPreconditioner_hh

//////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "Petsc/ShellPreconditioner.hh"
#include "Petsc/TridiagPcJFContext.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools { 
    class MatrixInverter; 
  }

  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a shell preconditioner object (Tridiagonal preconditioner)
 *
 * @author Jiri Simonek
 *
 */

class TridiagPreconditioner : public ShellPreconditioner {
public:

 /**
  * Constructor
  */
  TridiagPreconditioner(const std::string& name);

  /**
   * Default destructor
   */
  ~TridiagPreconditioner();

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

  /// socket for diagonal jacobians
  Framework::DataSocketSink<CFreal> socket_diagMatrices;

  /// socket for under-diagonal jacobians
  Framework::DataSocketSink<CFreal> socket_underDiagMatrices;

  /// socket for above-diagonal jacobians
  Framework::DataSocketSink<CFreal> socket_aboveDiagMatrices;

  /// storage of the local updatable IDs or -1 (ghost) for all local states
  Framework::DataSocketSink<CFint> socket_upLocalIDsAll;
	
	/// storage of booleans, which indicates if the cell is the first in a line (for DPLR line searching algorithm)
	Framework::DataSocketSink<bool> socket_dplrIsFirstInLine;

  /**
	 * DPLR IDs to local IDs converter
	 * storage of integers, which stores IDs of cells in the lines (connectivity).
	 * The begining of the next line is indicated by vector of booleans "socket_dplrIsFirstLine"
	 */
	Framework::DataSocketSink<CFint> socket_dplrToLocalIDs;

  /// Local IDs to DPLR IDs converter
	Framework::DataSocketSink<CFint> socket_localToDplrIDs;

  /// storage of integers, which stores the information to which line the cell belongs
	Framework::DataSocketSink<CFint> socket_dplrCellInLine;

  /// Tridiag context 
  TridiagPcJFContext _pcc;

  /// temporary data for holding the matrix inverter
  std::auto_ptr<MathTools::MatrixInverter> _inverter;

}; // end of class TridiagPreconditioner

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_TridiagPreconditioner_hh

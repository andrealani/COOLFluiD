#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtMisMPI_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtMisMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "mpi.h"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic outlet command
 *
 * @author Alessandro Sanna
 *
 *
 */
class SubOutletEuler2DPuvtMisMPI : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubOutletEuler2DPuvtMisMPI(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubOutletEuler2DPuvtMisMPI();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Set the preProcesses information of the whole outlet
   */
  virtual void preProcess();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

private: // data

  /// isoentropic Mach number
  CFreal M_is;
 
  /// total pressure at inlet
  CFreal Pt_inlet;

  ///my process
  CFuint _myP;

  ///number processes
  CFuint _nP;

  ///MPI_COMM_WORLD
  MPI_Comm _comm;

  /// sum of pressure values along outlet's faces belonging to the current trs
  double _P_sum; 

  /// sum of pressure values along outlet's faces of all trses
  double _Tot_P_sum;

  /// number of outlet faces belonging to the trs of the process
  unsigned int _n_adds; 

  /// total number of outlet faces 
  unsigned int _Tot_N_adds;

  /// average pressure along the whole outlet
  CFreal _P_average_inner;

  /// isentropic pressure along the outlet
  CFreal _P_is;

  /// damping factor
  CFreal _damping;

}; // end of class SubOutletEuler2DPuvtMisMPI

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtMisMPI_hh

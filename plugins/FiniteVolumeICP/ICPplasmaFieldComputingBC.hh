#ifndef COOLFluiD_Numerics_FiniteVolume_ICPplasmaFieldComputingBC_hh
#define COOLFluiD_Numerics_FiniteVolume_ICPplasmaFieldComputingBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "ICP/ICPReactionTerm.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {
    
    namespace FiniteVolumeICP {
          
//////////////////////////////////////////////////////////////////////////////
      
 /**
  * This class represents a command that applies the 
  * 2D EM field BCs ("integral boundary" approach) for the 
  * induction equation in ICP simulation.
  * 
  * @author Emanuele Sartori
  */
template <typename BASE, typename ST>
class ICPplasmaFieldComputingBC : public BASE {

//////////////////////////////////////////////////////////////////////////////
//                      Class Write Tecplot File                            //
//                                                                          //
//  This class implements induction equations BCs when using 2D EM field    //
//  model. Every ghost cell is coupled with each cell in the domain (or,    //
//  let say, to the current running into them).                             //
//                                                                          //
//  This computation is slow, so we try to reduce to the minimum the        //
//  computing effort: .preProcess computes J one time. At the end this      //
//  is done 4x4=16 times per iteration..... we should improve this.         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

public: 
  
  /**
   * Constructor
   */
  ICPplasmaFieldComputingBC(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ICPplasmaFieldComputingBC();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /**
   * compute el conductivity just one time per iteration per instance of this class.
   * It is possible to compute also J in cells; but Ep is changing during one single iteration.
   * TODO let's check how much (time) do we save using currentInCells_Re and currentInCells_Im
   */
  virtual void preProcess();

  /**
   * elliptic integrals (vector potential)
   */
  CFreal ellipticIntegralFirstKind(CFreal const& k);
  CFreal ellipticIntegralSecondKind(CFreal const& k);
  CFreal ellipticIntegralCombined(CFreal const& k);
  
private: //data

  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  /// socket for the vacuum electric field intensity (real and imaginary component)
  Framework::DataSocketSink<RealVector> socket_vacuumElFieldIntensity;
  /// socket for the elCondfield storage
  Framework::DataSocketSink<CFreal> socket_elCondField;
  /// socket for the current in cells (real and imaginary component)
  //  (used to write .plt file)
  Framework::DataSocketSink<RealVector> socket_currentInCells;

  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  /// source set
  Common::SafePtr<ST> m_srcTerm;
  
  /// mapping between ghost state and unique local ID
  Common::CFMap<Framework::State*, CFuint> m_mapGhostState2ID;
  
  /// real part of electric field array in ghost states
  RealVector m_EpR_inGhostCell_sum;
  
  /// imaginary part of electric field array in ghost states
  RealVector m_EpI_inGhostCell_sum;
  
  /// cell centers coordinates
  RealVector m_cellCentersCoord;

  /// current in cell (Re and Im parts)
  RealVector m_currentInCells;
  
  /// physical data array
  RealVector m_physicalData;
  
  /// number of states in processor
  std::vector<CFuint> m_nbStatesInProc;
  
}; // end of class ICPplasmaFieldComputingBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICPplasmaFieldComputingBC.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ICPInductionBC_hh

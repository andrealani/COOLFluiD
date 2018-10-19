#ifndef COOLFluiD_Numerics_FiniteVolumeICP_RMSJouleHeatSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolumeICP_RMSJouleHeatSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a time-averaged Joule heat source term
 *
 * @author Radek Honzatko
 *
 */
class RMSJouleHeatSourceTerm : public ComputeSourceTermFVMCC {

public:

   /**
    * Defines the Config Option's of this class
    * @param options a OptionList where to add the Option's
    **/
   static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  RMSJouleHeatSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RMSJouleHeatSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Set the variable set
   * @pre the input pointer is non const to allow dynamic_cast
   */
  virtual void setVarSet(Common::SafePtr<Framework::ConvectiveVarSet> varSet)
  {
  }
  
  /**
   * Compute the source term and jacobian
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
                             RealVector& source, RealMatrix& jacobian);
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:
  
   /// Compute the radiation loss term 
   virtual CFreal computeRadiationLoss(Framework::GeometricEntity *const element);

private: // data
  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_rmsJouleHeatSource;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  /// Euler physical data
  RealVector m_physicalData;

  /// flag to add joule heating to both energy equations (wrong)
  bool m_addToBothEnergyEquations;

  /// flag to calculate the Simple Radiation Heat Transfer Model
  bool m_computeArgonRadiativeHeatTransferLoss;

  CFreal m_ArRadC1;
  CFreal m_ArRadC2;
  CFreal m_ArRadC3;
  CFreal m_ArRadC4;
  CFreal m_ArRadC5;
  CFreal m_ArRadC6;

}; // end of class RMSJouleHeatSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RMSJouleHeatSourceTerm_hh

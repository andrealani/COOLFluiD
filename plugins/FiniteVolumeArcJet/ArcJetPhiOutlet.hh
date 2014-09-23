#ifndef COOLFluiD_Numerics_FiniteVolume_ArcJetPhiOutlet_hh
#define COOLFluiD_Numerics_FiniteVolume_ArcJetPhiOutlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/VectorialFunction.hh"
#include "Common/ParserException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {
    
    namespace FiniteVolumeArcJet {
          
//////////////////////////////////////////////////////////////////////////////
      
      /**
       * This class represents a command that applies the BC for the 
       * electric potential equation imposing the electric field
       * for a given current intensity
       * in ArcJet simulation
       * 
       * @author Andrea Lani
       * @author Alejandro Alvarez Laguna
       * @author Amrita Lonkar
       */
template <class BASE>
class ArcJetPhiOutlet : public BASE {

public: 
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  ArcJetPhiOutlet(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ArcJetPhiOutlet();
  
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
   * Set the preProcesses connectivity between faces belonging to different process
   */
  virtual void preProcess();
 
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private:
  
  /// storage of the face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
   
  /// integral of sigma*area over the exit face
  CFreal m_sigmaIntegral;
  
  /// imposed current
  CFreal m_imposedI;
  
  /// Limit in voltage given by the machine
  CFreal _machineLimit;
  
  /// CFL function of iteration
  Framework::VectorialFunction _vFunction;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// a string holding the function definition
  std::vector<std::string> _function;
  
  ///the iter number
  RealVector _iter;
  
  ///The value after the function
  RealVector _resultI;
  
  
}; // end of class ArcJetPhiOutlet

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetPhiOutlet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ArcJetPhiOutlet_hh

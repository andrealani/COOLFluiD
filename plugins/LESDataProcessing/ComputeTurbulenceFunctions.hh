#ifndef COOLFluiD_Numerics_LESDataProcessing_ComputeTurbulenceFunctions_hh
#define COOLFluiD_Numerics_LESDataProcessing_ComputeTurbulenceFunctions_hh

//////////////////////////////////////////////////////////////////////////////

#include "LESProcessingData.hh"
#include "Framework/SubSystemStatus.hh"
#include "LESProcessingCom.hh"

#include "LES/LESVarSet.hh"


//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {
  
  namespace Numerics {

    namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

/**
 * This command computes all the TurbulenceFunctions that are declared in
 * LESProcessingData. 
 *
 * They are calculated in a cell-centred way and can be outputted with the TecplotWriter module
 * Simulator.SubSystem.Tecplot.WriteSol = WriteSolutionBlockFV
 * Simulator.SubSystem.Tecplot.Data.DataHandleOutput.CCSocketNames = Qcriterion Vorticity SGSViscosity
 * Simulator.SubSystem.Tecplot.Data.DataHandleOutput.CCVariableNames = Q X-vorticity Y-vorticity Z-vorticity muT
 * 
 * @author Willem Deconinck
 *
 */
class ComputeTurbulenceFunctions : public LESProcessingCom {
public:
  
  typedef LES::LESVarSet LESVAR;
  
  /**
   * Constructor.
   */
  ComputeTurbulenceFunctions(const std::string& name) ;



  /**
   * Destructor.
   */
  virtual ~ComputeTurbulenceFunctions()
  {
  }

  static void defineConfigOptions(Config::OptionList& options);

  virtual void setup();
  
  virtual void unsetup();

  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  virtual void execute();

protected:
  
  /// arrray of gradients
  std::vector<RealVector> m_gradients;

private:

  CFuint m_processRate;
  
  std::vector<Common::SafePtr <TurbulenceFunction> > m_turbulenceFunctions;
    
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESDataProcessing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_ComputeTurbulenceFunctions_hh


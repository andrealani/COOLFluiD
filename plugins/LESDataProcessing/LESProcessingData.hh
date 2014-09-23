#ifndef COOLFluiD_Numerics_LESProcessing_LESProcessingData_hh
#define COOLFluiD_Numerics_LESProcessing_LESProcessingData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/MathConsts.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include <deque>
#include "Framework/VarSetTransformer.hh"
#include "LESProcessingCom.hh"
#include "GradientComputer.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/ConvectiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

  	namespace LESDataProcessing {
      
      // put strategy declarations here
      class GradientComputer;
      class TurbulenceFunction;

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * FilterCom 's that compose the LESProcessing.
   *
   * @see FilterCom
   *
   * @author Willem Deconinck
   */
class LESProcessingData : public Framework::DataProcessingData {

  
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  LESProcessingData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~LESProcessingData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up the member data
   */
   virtual void setup();
   
   /**
    * Unset up the member data
    */
    virtual void unsetup();
   

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "LESProcessingData";
  }

  virtual Common::SafePtr<Framework::VarSetTransformer> getUpdateToPrimTransformer()
  {
    return m_updateToPrimVar.getPtr();
  }
  
  /**
   * Get the gradient computer strategy
   * @return  reference to the gradient computer
   */
  Common::SafePtr<GradientComputer> getGradientComputer() const
  {
    cf_assert(m_gradientComputer.isNotNull());
    return m_gradientComputer.getPtr();
  }
  
  Common::SafePtr<Framework::DiffusiveVarSet> getLESVar();
  
  RealVector transformToPrim(RealVector& state);
  
  RealVector transformToPrimDim(const RealVector& state);
  
  CFreal getVolume(const CFuint& cellID);

  CFreal getVolumeAdim(const CFuint& cellID);
  
  std::vector<Common::SafePtr<TurbulenceFunction> > getTurbulenceFunctions()
  {
    std::vector<Common::SafePtr<TurbulenceFunction> > ptrs(m_turbulenceFunctions.size());
    for(CFuint i=0; i<m_turbulenceFunctions.size(); ++i) {
      ptrs[i] = m_turbulenceFunctions[i].getPtr();
    }
    return ptrs;
  }
  
  std::vector<Common::SafePtr<TurbulenceFunction> > getAverageTurbulenceFunctions()
  {
    std::vector<Common::SafePtr<TurbulenceFunction> > ptrs(m_averageTurbulenceFunctions.size());
    for(CFuint i=0; i<m_averageTurbulenceFunctions.size(); ++i) {
      ptrs[i] = m_averageTurbulenceFunctions[i].getPtr();
    }
    return ptrs;
  }

private:
  
  /// Euler Varset
  Common::SafePtr<Framework::ConvectiveVarSet> m_updateVar;
    
  /// State to dimensionalize updateVar
  RealVector m_dimState;
  
  /// Primitive state
  Framework::State* m_primState;
  
  /// Transformer from update to primitive Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_updateToPrimVar;
  
  /// Gradient Computer
  Common::SelfRegistPtr<GradientComputer> m_gradientComputer;
  std::string m_gradientComputerStr;
  
  /// Turbulence Functions
  std::vector<Common::SelfRegistPtr<TurbulenceFunction> > m_turbulenceFunctions;
  std::vector<std::string> m_turbulenceFunctionNames;
  
  /// Average Turbulence Functions
  std::vector<Common::SelfRegistPtr<TurbulenceFunction> > m_averageTurbulenceFunctions;
  std::vector<std::string> m_averageTurbulenceFunctionNames;
  
  
}; // end of class LESProcessingData

//////////////////////////////////////////////////////////////////////////////

		} // namespace LESProcessing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESProcessing_LESProcessingData_hh

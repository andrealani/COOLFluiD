#ifndef COOLFluiD_Numerics_ExplicitFilters_FilterData_hh
#define COOLFluiD_Numerics_ExplicitFilters_FilterData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/MathConsts.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include <deque>

#include "FilterStencil.hh"
#include "FilterWeight.hh"
#include "StencilComputer.hh"
#include "FilterStrategy.hh"
#include "CoordinateLinker.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  	namespace ExplicitFilters {
      
      class FilterStrategy;
      class StencilComputer;
      class CoordinateLinker;

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * FilterCom 's that compose the ExplicitFilters.
   *
   * @see FilterCom
   *
   * @author Willem Deconinck
   */
class FilterData : public Framework::DataProcessingData {

  
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  FilterData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~FilterData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up the member data
   */
   virtual void setup();

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "FilterData";
  }
  
public:
  /**
   * Get the stencil computer
   * @return  reference to the stencil computer
   */
  Common::SafePtr<StencilComputer> getStencilComputer() const
  {
    cf_assert(m_stencilComputer.isNotNull());
    return m_stencilComputer.getPtr();
  }
  
  /**
   * Get the filter strategy
   * @return  reference to the filter strategy
   */
  Common::SafePtr<FilterStrategy > getFilterStrategy() const
  {
    cf_assert(m_filterStrategy.isNotNull());
    return m_filterStrategy.getPtr();
  }
  
    /**
   * Get the filter strategy
   * @return  reference to the filter strategy
   */
  Common::SafePtr<CoordinateLinker > getCoordinateLinker() const
  {
    cf_assert(m_coordinateLinker.isNotNull());
    return m_coordinateLinker.getPtr();
  }
  
  /** 
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> >
    getGeoWithNodesBuilder() 
  {
   return &m_geoWithNodesBuilder;
  }

  /**
   * @return  The filter to grid ratio
   */
  CFreal getFilterGridRatio() const {
    return m_filterGridRatio;
  }
  
  /**
   * @return  The filter-stencil to grid ratio
   */
  CFreal getStencilGridRatio() const {
    return m_stencilGridRatio;
  }
  
  /**
   * @return A vector of cell-ID's to inspect closer
   */
  Common::SafePtr<std::vector<CFuint> > getInspectedCellIDs() {
    return &m_inspectedCellIDs;
  }

  /**
   * @return The order of the filter
   */
  CFuint getOrder() const {
    return m_order;
  }
  
  /**
   * @return Order of the target filter (analytical function)
   */
  CFuint getTargetFilterOrder() const {
    return m_targetFilterOrder;
  }
  
  /**
   * @return  The name of the prefix for transferFunction data files
   */
  std::string getTransferFunctionFileName() const {
    return m_transferFunctionFileName;
  }
  
  /**
   * @return A vector of ID's of cells that will be filtered (0=off, 1=on)
   */
  Common::SafePtr<std::deque<bool> > getFilterFlags() {
    return &m_filterFlag;
  }
   
  /**
   * @return A vector of ID's of cells that will be filtered (0=off, 1=on)
   */
  bool& getFilterFlag(const CFuint& idx) {
    return m_filterFlag[idx];
  }
  
  /**
   * @return A vector of stencils
   */
  Common::SafePtr<std::vector<FilterStencil> > getStencils() {
    return &m_stencil;
  }
   
  /**
   * @return A stencil object of the given index
   */
  Common::SafePtr<FilterStencil> getStencil(const CFuint& idx) {
    return (&m_stencil[idx]);
  }
   
   
  /**
   * @return A vector of weights
   */
  Common::SafePtr<std::vector<FilterWeight> > getWeights() {
    return &m_weight;
  }
  
  
  /**
   * @return A Weight object of the given index
   */
  Common::SafePtr<FilterWeight> getWeight(const CFuint& idx) {
    return (&m_weight[idx]);
  }
  
  
  bool outputDebug() {
    return m_outputDebug;
  }
  

private:
  
  /// Vector of filter stencils
  std::vector<FilterStencil> m_stencil;
  
  /// Vector of filter stencils
  std::vector<FilterWeight> m_weight;

  /// FilterStrategy
  Common::SelfRegistPtr<FilterStrategy>  m_filterStrategy;
  std::string m_filterTypeStr;
  
  /// StencilComputer strategy
  Common::SelfRegistPtr<StencilComputer> m_stencilComputer; 
  std::string m_stencilComputerStr;

  /// Coordinate linker strategy
  Common::SelfRegistPtr<CoordinateLinker>  m_coordinateLinker;
  std::string m_coordinateLinkerStr;

  /// builder of GeometricEntity's with Node's
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> m_geoWithNodesBuilder;

  /// Value of filter transfer function at filter cutoff
  CFreal m_Gcutoff;
  
  /// Filter to grid ratio
  CFreal m_filterGridRatio;
  
  /// Stencil size to grid size ratio
  CFreal m_stencilGridRatio;
  
  /// Order of the filter
  CFuint m_order;
  
  /// Order of the targetFilter for Least Squares optimization of some coefficients
  CFuint m_targetFilterOrder;
  
  /// Number of Rings used in the stencil of the filter
  CFuint m_nbRings;
  
  /// cells to inspect closer
  std::vector<CFuint> m_inspectedCellIDs;
  
  /// Base name of the output files for transfer functions
  std::string m_transferFunctionFileName;
  
  /// Flags if a cell should be filtered e.g. turn off (0) if bad transfer function
  std::deque<bool> m_filterFlag;
  
  /// Flag that tells if additional output has to be outputted in case of bad filters
  bool m_outputDebug;
  
}; // end of class ExplicitFiltersData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for Filter
typedef Framework::MethodCommand<FilterData> FilterCom;

/// Definition of a command provider for Filter
typedef Framework::MethodCommand<FilterData>::PROVIDER FilterComProvider;

//////////////////////////////////////////////////////////////////////////////

		} // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_FilterData_hh

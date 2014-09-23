#ifndef COOLFluiD_Numerics_ExplicitFilters_CoordinateLinker_hh
#define COOLFluiD_Numerics_ExplicitFilters_CoordinateLinker_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/MethodStrategy.hh"
#include "ExplicitFilters/FilterData.hh"
 
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {
		  
      class FilterData;
      class FilterStencil;

//////////////////////////////////////////////////////////////////////////////

/** 
 * This class offers a basic interface for coupling stencil element ID's with coordinates
 * 
 * @author Willem Deconinck
 */
class CoordinateLinker : public Framework::MethodStrategy<FilterData> {


public: // functions


  typedef Framework::BaseMethodStrategyProvider
  <FilterData,CoordinateLinker > PROVIDER;

  /** 
   * Constructor
   */
  CoordinateLinker(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~CoordinateLinker();

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "CoordinateLinker";
  }

  virtual Framework::Node& getCoordinates(const CFuint& idx, const bool isGhost = false) const = 0;
  virtual Framework::Node& getCoordinates(const Common::SafePtr<FilterStencil>& stencil, const CFuint& idx) const = 0;
  
  virtual Framework::State& getState(const CFuint& idx, const bool isGhost = false) const = 0;
  virtual Framework::State& getState(const Common::SafePtr<FilterStencil>& stencil, const CFuint& idx) const = 0;



}; // end of class CoordinateLinker


//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace ExplicitFilters

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_CoordinateLinker_hh

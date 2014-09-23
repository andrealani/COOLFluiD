#ifndef COOLFluiD_Numerics_ExplicitFilters_CoordinateLinkerFVM_hh
#define COOLFluiD_Numerics_ExplicitFilters_CoordinateLinkerFVM_hh

//////////////////////////////////////////////////////////////////////////////

#include "CoordinateLinker.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {
		  
//////////////////////////////////////////////////////////////////////////////

/** 
 * This class offers a basic interface for coupling stencil element ID's with coordinates
 * 
 * @author Willem Deconinck
 */
class CoordinateLinkerFVM : public CoordinateLinker {


public: // functions

  /** 
   * Constructor
   */
  CoordinateLinkerFVM(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~CoordinateLinkerFVM();

  /** 
   * Returns the DataSocket's that this numerical strategy needs 
   * as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
    needsSockets();
    
  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "CoordinateLinkerFVM";
  }


protected: // helper functions

  virtual Framework::Node& getCoordinates(const CFuint& idx, const bool isGhost) const;
  virtual Framework::Node& getCoordinates(const Common::SafePtr<FilterStencil>& stencil, const CFuint& idx) const;

  virtual Framework::State& getState(const CFuint& idx, const bool isGhost) const;
  virtual Framework::State& getState(const Common::SafePtr<FilterStencil>& stencil, const CFuint& idx)const;

  
private: // data

  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  Framework::DataSocketSink < Framework::State*> socket_gstates;


protected: // data



}; // end of class CoordinateLinkerFVM


//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace ExplicitFilters

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_CoordinateLinkerFVM_hh


#include "CoordinateLinkerFVM.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "ExplicitFilters.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<CoordinateLinkerFVM, 
                                  FilterData, 
                                  CoordinateLinker, 
                                  ExplicitFiltersModule> 
coordinateLinkerFVM("CoordinateLinkerFVM");
      

//////////////////////////////////////////////////////////////////////////////

CoordinateLinkerFVM::CoordinateLinkerFVM(const std::string& name) :
  CoordinateLinker(name),
  socket_states("states"),
  socket_gstates("gstates")
{
}

//////////////////////////////////////////////////////////////////////////////


CoordinateLinkerFVM::~CoordinateLinkerFVM()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > CoordinateLinkerFVM::needsSockets()
{
  // create socket sink for the stencil
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  result.push_back(&socket_states); 
  result.push_back(&socket_gstates);  
  return result;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
      
Framework::Node& CoordinateLinkerFVM::getCoordinates(const CFuint& idx, const bool isGhost) const
{
  if (isGhost) {
    Framework::DataHandle<Framework::State* > gstates = socket_gstates.getDataHandle();
    return gstates[idx]->getCoordinates();
  } else {
    Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
    return states[idx]->getCoordinates();  
  }
}
      
//////////////////////////////////////////////////////////////////////////////

Framework::Node& CoordinateLinkerFVM::getCoordinates(const Common::SafePtr<FilterStencil>& stencil, const CFuint& idx) const
{
  if (stencil->isGhost(idx)) {
    Framework::DataHandle<Framework::State* > gstates = socket_gstates.getDataHandle();
    return gstates[stencil->getElement(idx)]->getCoordinates();
  } else {
    Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
    return states[stencil->getElement(idx)]->getCoordinates();  
  }
}

      
//////////////////////////////////////////////////////////////////////////////

Framework::State& CoordinateLinkerFVM::getState(const CFuint& idx, const bool isGhost) const
{
  if (isGhost) {
    Framework::DataHandle<Framework::State* > gstates = socket_gstates.getDataHandle();
    return (*gstates[idx]);
  } else {
    Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
    return (*states[idx]);  
  }
}
      
//////////////////////////////////////////////////////////////////////////////

Framework::State& CoordinateLinkerFVM::getState(const Common::SafePtr<FilterStencil>& stencil, const CFuint& index) const
{
  CFuint idx = index;
  if (stencil->isGhost(idx)) {
    Framework::DataHandle<Framework::State* > gstates = socket_gstates.getDataHandle();
    return (*gstates[stencil->getElement(idx)]);
  } else {
    Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
    return (*states[stencil->getElement(idx)]);  
  }
}

//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace ExplicitFilters
		
  } // namespace Numerics

} // namespace COOLFluiD

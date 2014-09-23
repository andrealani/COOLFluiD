#include "ExplicitFilters/ExplicitFilters.hh"
#include "GetTransferFunction.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MathTools/SVDInverter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<GetTransferFunction, 
                      FilterData, 
                      ExplicitFiltersModule> 
GetTransferFunctionProvider("GetTransferFunction");

//////////////////////////////////////////////////////////////////////////////

void GetTransferFunction::execute()
{
  
  if((!(Framework::SubSystemStatusStack::getActive()->getNbIter() % m_processRate)) ||
    (Framework::SubSystemStatusStack::getActive()->getNbIter() == 0)) {
  
  
    // Compute and write the stencils for explicit filtering
    Common::SafePtr<StencilComputer> stencilComputer =
      getMethodData().getStencilComputer();
    stencilComputer->computeInspected();
    stencilComputer->outputStencil();
  
    // Write information of the filter
    Common::SafePtr<FilterStrategy> filter =
      getMethodData().getFilterStrategy();
    filter->getInfo();

    // Output the transferfunction of inspected cells
    filter->outputTransferFunction();
  
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > > GetTransferFunction::needsSockets()
{
  CFAUTOTRACE;
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

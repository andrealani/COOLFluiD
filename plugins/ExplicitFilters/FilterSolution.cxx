#include "ExplicitFilters/ExplicitFilters.hh"
#include "FilterSolution.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FilterSolution, 
                      FilterData, 
                      ExplicitFiltersModule> 
filterSolutionProvider("FilterSolution");

//////////////////////////////////////////////////////////////////////////////

void FilterSolution::execute()
{
  CFAUTOTRACE;
  
  if((!(Framework::SubSystemStatusStack::getActive()->getNbIter() % m_processRate)) ||
    (Framework::SubSystemStatusStack::getActive()->getNbIter() == 0)) {
    
    // Data handle of the solution states
    Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
    
    
    // Linker with the Space method
    Common::SafePtr<CoordinateLinker> coordinateLinker = getMethodData().getCoordinateLinker();

    // Other shortcut allocations and declarations
    CFuint nbCellsInStencil(0);
    const CFuint nbCells = getMethodData().getStencils()->size();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
    // Allocation for filtered states
    std::vector<Framework::State> filteredStates(nbCells);
    
    // Calculation of filtered states
    for(CFuint iStencil=0; iStencil<nbCells; ++iStencil) {
      if (getMethodData().getFilterFlag(iStencil)) {
        Common::SafePtr<FilterStencil> stencil = getMethodData().getStencil(iStencil);
        Common::SafePtr<FilterWeight > weights = getMethodData().getWeight(iStencil);
        if (weights->getNbElements() == 0)  CFLog(INFO, "no weights for cell " << iStencil << " \n");
        nbCellsInStencil = stencil->getNbElements();
        for(CFuint iCell=0; iCell<nbCellsInStencil; ++iCell) {
          for(CFuint iEq=0; iEq<nbEqs; ++iEq) {
            filteredStates[iStencil][iEq] += coordinateLinker->getState(stencil,iCell)[iEq] * weights->getWeight(iCell);
          }
        } 
      }
    }
  
    // Replace old states with filtered states
    for(CFuint iStencil=0; iStencil<nbCells; ++iStencil) {
      if (getMethodData().getFilterFlag(iStencil)) {
        (*states[iStencil]) = filteredStates[iStencil];
      }     
    }
  
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > > FilterSolution::needsSockets()
{
  CFAUTOTRACE;
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

#include "ExplicitFilters/ExplicitFilters.hh"
#include "FilterRHS.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FilterRhs, 
                      FilterData, 
                      ExplicitFiltersModule> 
FilterRhsProvider("FilterRhs");

//////////////////////////////////////////////////////////////////////////////

void FilterRhs::execute()
{
  CFAUTOTRACE;

  if((!(Framework::SubSystemStatusStack::getActive()->getNbIter() % m_processRate)) ||
    (Framework::SubSystemStatusStack::getActive()->getNbIter() == 0)) {

    // Data handle for the right hand side
    DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
    
    // Data handle of the solution states
    Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
        
    // Allocation of the filtered right hand side
    std::vector<CFreal> filteredRhs(rhs.size());
    
    // Other shortcut declarations and allocations
    const CFuint nbCells = getMethodData().getStencils()->size();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    CFuint nbCellsInStencil(0);
    CFuint stencilCellID(0);
    CFuint centreCellID(0);
    
    // Calculate the filtered right hand side
    for(CFuint iStencil=0; iStencil<nbCells; ++iStencil) {
      if (getMethodData().getFilterFlag(iStencil)) {
        Common::SafePtr<FilterStencil> stencil = getMethodData().getStencil(iStencil);
        Common::SafePtr<FilterWeight > weights = getMethodData().getWeight(iStencil);
        if (weights->getNbElements() == 0)  CFLog(INFO, "no weights for cell " << iStencil << " \n");
        nbCellsInStencil = stencil->getNbElements();
        centreCellID = stencil->getID();
        for(CFuint iCell=0; iCell<nbCellsInStencil; ++iCell) {
          stencilCellID = stencil->getElement(iCell);
          for(CFuint iEq=0; iEq<nbEqs; ++iEq) {
            filteredRhs[centreCellID+iEq*nbEqs] += rhs(stencilCellID,iEq,nbEqs) * weights->getWeight(iCell);
          }
        }
      }
    }
  
    // Correct the solution with the newly found filtered right hand side
    CFuint nbFilteredRHS = 0;
    for(CFuint iStencil=0; iStencil<nbCells; ++iStencil) {
      if (getMethodData().getFilterFlag(iStencil)) {
        bool filtered = false;
        for(CFuint iEq=0; iEq<nbEqs; ++iEq) {
          CFreal filteredRHS = filteredRhs[iStencil+iEq*nbEqs];
          CFreal RHS = rhs(iStencil,iEq,nbEqs);
          if (std::abs(RHS) > MathTools::MathConsts::CFrealEps()) {
            CFreal difference = 100*std::abs((RHS-filteredRHS)/(RHS));
            if (!MathTools::MathChecks::isNaN(filteredRHS) && difference < 50) {
              (*states[iStencil])[iEq] += filteredRHS - RHS;
              if (!filtered) filtered = true;
            }
          }
        }
        if (filtered) nbFilteredRHS++;
      }     
    }
    CFLog(INFO,"nbFilteredRHS = " << nbFilteredRHS << "/" << nbCells << "\n");
  
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > > FilterRhs::needsSockets()
{
  CFAUTOTRACE;
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

#include "FiniteVolume/FiniteVolume.hh"
#include "StdSetNodalStates.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetNodalStates,
                      CellCenterFVMData,
                      FiniteVolumeModule>
stdSetNodalStatesProvider("StdSetNodalStates");

//////////////////////////////////////////////////////////////////////////////

StdSetNodalStates::StdSetNodalStates(const std::string& name) :
  CellCenterFVMCom(name)
{
  addConfigOptionsTo(this);
  
  m_updateGradients = false;
  setParameter("updateGradients",&m_updateGradients);
}

   //////////////////////////////////////////////////////////////////////////////

StdSetNodalStates::~StdSetNodalStates()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSetNodalStates::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("updateGradients", "Flag telling whether to recompute the gradients."); 
}
      
//////////////////////////////////////////////////////////////////////////////

void StdSetNodalStates::execute()
{
   CFAUTOTRACE;

   // extrapolate the solution from cell centers to nodes
   getMethodData().getNodalStatesExtrapolator()->extrapolateInAllNodes();

   if (m_updateGradients) {
     // compute the gradients
     getMethodData().getPolyReconstructor()->computeGradients();
     
     CFLog(VERBOSE, "StdSetNodalStates::execute()\n");
   }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

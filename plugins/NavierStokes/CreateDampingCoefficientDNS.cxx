#include "Framework/MethodCommandProvider.hh"
//#include "LinEuler/LinearizedEuler.hh"
#include "CreateDampingCoefficientDNS.hh"
//#include "LinEuler/LinEulerTerm.hh"
//#include "Framework/SubSystemStatus.hh"
//#include "Common/BadValueException.hh"
#include "NavierStokes/NavierStokes.hh"
////////////////////////////////////////////////////////////////////////////
//

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

////////////////////////////////////////////////////////////////////////////
//

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

////////////////////////////////////////////////////////////////////////////
//

MethodCommandProvider<CreateDampingCoefficientDNS,
                      DataProcessingData,
                      NavierStokesModule>
aCreateDampingCoefficientdnsProvider("CreateDampingCoefficientDNS");

////////////////////////////////////////////////////////////////////////////
//

void CreateDampingCoefficientDNS::defineConfigOptions(Config::OptionList&
options)
{
  options.addConfigOption< CFreal >("nuMax","Maximum damping coefficient.");
  options.addConfigOption< std::vector<std::string> >("BoundaryTRS","Name of the outlet boundaries.");
  options.addConfigOption< CFreal >("r0","Radius of influence of the damping zone.");
  options.addConfigOption< CFreal >("beta","Exponent of the damping function.");

}

////////////////////////////////////////////////////////////////////////////
//

CreateDampingCoefficientDNS::CreateDampingCoefficientDNS(const std::string&
name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_dampingCoeff("dampingCoeff")
{
  addConfigOptionsTo(this);

  /// var names should be defined in CFcase.
  CF_DEBUG_POINT;
  /// Initialize the maximum value of the damping coefficient
  m_nuMax = 0.0;
  setParameter("nuMax",&m_nuMax);

  /// Radius of influence
  m_r0 = 1.0;
  setParameter("r0",&m_r0);

  /// Exponent of the damping function
  m_beta = 0.0;
  setParameter("beta",&m_beta);

  /// Name of the affected boundary
  m_BoundaryTRS = vector<std::string>();
  setParameter("BoundaryTRS",&m_BoundaryTRS);
}

////////////////////////////////////////////////////////////////////////////
//

CreateDampingCoefficientDNS::~CreateDampingCoefficientDNS()
{
}

////////////////////////////////////////////////////////////////////////////
//

std::vector<Common::SafePtr<BaseDataSocketSink> >
CreateDampingCoefficientDNS::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

////////////////////////////////////////////////////////////////////////////
//

std::vector<Common::SafePtr<BaseDataSocketSource> >
CreateDampingCoefficientDNS::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_dampingCoeff);
  return result;
}

////////////////////////////////////////////////////////////////////////////
//

void CreateDampingCoefficientDNS::setup()
{CF_DEBUG_POINT;
  DataHandle < Framework::State*, Framework::GLOBAL > states =
socket_states.getDataHandle();
  DataHandle<CFreal> dampingCoeff = socket_dampingCoeff.getDataHandle();
  dampingCoeff.resize(states.size());

  // Initialize the damping coefficient with zeros (This has to be doneexactly once, at the beginning during setup)
  for (CFuint iState = 0; iState < states.size(); ++iState)
  {
	dampingCoeff[iState] = 0.0;
  }
executeOnTrs();
}

////////////////////////////////////////////////////////////////////////////
//

void CreateDampingCoefficientDNS::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS(); //modify this later
  CFLogDebugMin( "CreateDampingCoefficient::executeOnTrs() called for TRS: "
<< trs->getName() << "\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states =
socket_states.getDataHandle();
  DataHandle<CFreal> dampingCoeff = socket_dampingCoeff.getDataHandle();

  // Get the size of the list containing all the TRSs
  vector<Common::SafePtr<TopologicalRegionSet> > trsList =
MeshDataStack::getActive()->getTrsList();
  const CFint nbTRSs = trsList.size();

  // Outer loop, cycle through all boundaries in m_BoundaryTRS
  for(CFuint iBoundaryTRS = 0; iBoundaryTRS < m_BoundaryTRS.size();
++iBoundaryTRS)
  {
    SafePtr<TopologicalRegionSet> boundaryTrs;
    for (CFint iTRS = 0; iTRS < nbTRSs; ++iTRS)
    {
      SafePtr<TopologicalRegionSet> currTrs = trsList[iTRS];
      if (currTrs->getName() == m_BoundaryTRS[iBoundaryTRS])
      {
        boundaryTrs = currTrs;
        break;
      }
    }

    // Now that I have the number of the current boundary, I can refer to it by calling "boundaryTrs"

    // This vector contains the global indices for every node(or state) of the TRS we are dealing with, which is "boundaryTrs"
    Common::SafePtr< vector<CFuint> > const statesIdxInTRS = boundaryTrs->getStatesInTrs();

    // Cycle through all the nodes of the given boundary, which happens to be "boundaryTrs"
    for (CFuint iBoundaryState = 0; iBoundaryState <boundaryTrs->getNbStatesInTrs(); ++iBoundaryState)
    {
      // Obtain coordinates for the given node
      Node& boundary_coord = states[(*statesIdxInTRS)[iBoundaryState]]->getCoordinates();
      CFreal x_b = boundary_coord[XX];
      CFreal y_b = boundary_coord[YY];
      //Cycle through all the nodes of the domain
      for (CFuint iState = 0; iState < states.size(); iState++) //is itrather ++iState?
      {
        // Obtain coordinates for the given node
        Node& coord = states[iState]->getCoordinates();
        const CFreal x = coord[XX];
        const CFreal y = coord[YY];
        const CFreal r = sqrt( (x-x_b)*(x-x_b)+(y - y_b)*(y - y_b));
        // If the node is inside the radius of influence...
        if (r < m_r0)
        {
          // Calculate (local) the damping coefficient
          const CFreal nu = m_nuMax * pow( ((m_r0 - r)/m_r0), m_beta);
          // Retrieve current value of the damping coefficient
          const CFreal nu_old = dampingCoeff[iState];
          // Update damping coefficient using the larger of the two values
          dampingCoeff[iState] = max(nu_old,nu);


        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////
//

void CreateDampingCoefficientDNS::unsetup()
{
  DataHandle<CFreal> dampingCoeff = socket_dampingCoeff.getDataHandle();

  dampingCoeff.resize(0);
}

////////////////////////////////////////////////////////////////////////////
//

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

////////////////////////////////////////////////////////////////////////////
//


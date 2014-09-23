#include "NavierStokesSATruncateK.hh"
#include "FiniteVolumeTurb/FiniteVolumeSA.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesSATruncateK,
                      DataProcessingData,
                      FiniteVolumeSAModule>
NavierStokesSATruncateKProvider("NavierStokesSATruncateK");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSATruncateK::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////
NavierStokesSATruncateK::NavierStokesSATruncateK(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states")
{
   addConfigOptionsTo(this);

}

//////////////////////////////////////////////////////////////////////////////

NavierStokesSATruncateK::~NavierStokesSATruncateK()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokesSATruncateK::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
void NavierStokesSATruncateK::execute()
{
  CFAUTOTRACE;

  const CFuint KIDX = 4;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  CFuint nbNegativesK = 0;
  CFuint nbNegativesOmega = 0;
  for (CFuint iState = 0; iState < states.size(); ++iState)
  {
	// truncate the K entry
	CFreal& K = (*states[iState])[KIDX];
        if(K < 0.) nbNegativesK++;
        K = std::max(K,0.000002);

        // truncate the K and Omega entries
        if(states.size() == 6){
          CFreal& Omega = (*states[iState])[KIDX+1];
          if(Omega < 0.) nbNegativesOmega++;
          K = std::max(K,0.000000166);
          Omega = std::max(Omega,121.);
        }
  }

  if(nbNegativesK > 0) cout << "**  TRUNCATED " << nbNegativesK << " states for K" << endl;
  if(nbNegativesOmega > 0) cout << "**  TRUNCATED " << nbNegativesOmega << " states for Omega" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSATruncateK::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

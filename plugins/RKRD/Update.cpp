#include "RKRD/RungeKuttaRD.hh"
#include "RKRD/Update.hpp"

#include "MathTools/RealVector.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////


using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Update, RKRDData, RKRDModule> UpdateProvider("Update");

//////////////////////////////////////////////////////////////////////////////

Update::Update(const std::string& name) : RKRDCom(name),
  socket_kstates("kstates"),
  socket_states("states")
{
}

Update::~Update()
{
}

void Update::execute()
{
  const CFuint nbeqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < RealMatrix > kstates = socket_kstates.getDataHandle();

  const CFuint nbstates = states.size();

  const CFuint om1 = getMethodData().getOrder() - 1 ;

  for ( CFuint s = 0; s < nbstates; ++s )
  {
    RealMatrix& kmat = kstates[s];
    for ( CFuint eq = 0; eq < nbeqs; ++eq )
      (*states[s])[eq] = kmat(eq,om1);
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > Update::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_kstates);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

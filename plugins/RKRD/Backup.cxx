#include "RKRD/RungeKuttaRD.hh"
#include "RKRD/Backup.hh"

#include "MathTools/RealVector.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Backup, RKRDData, RKRDModule> BackupProvider("Backup");

//////////////////////////////////////////////////////////////////////////////

Backup::Backup(const std::string& name) : RKRDCom(name),
  socket_kstates("kstates"),
  socket_states("states")
{
}

Backup::~Backup()
{
}

void Backup::execute()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  DataHandle<RealMatrix> kstates = socket_kstates.getDataHandle();

  CFuint order = getMethodData().getOrder();

  // put the states on all columns of the kstates U_k
  CFuint nbstates = states.size();
  for ( CFuint s = 0; s < nbstates ; ++s )
  {
    for ( CFuint k = 0; k < order ; ++k )
    {
      kstates[s].setColumn(*states[s], k);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > Backup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_kstates);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

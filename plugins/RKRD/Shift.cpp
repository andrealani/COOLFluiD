#include "RKRD/RungeKuttaRD.hh"
#include "RKRD/Shift.hpp"

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

MethodCommandProvider<Shift, RKRDData, RKRDModule> ShiftProvider("Shift");

//////////////////////////////////////////////////////////////////////////////

Shift::Shift(const std::string& name) : RKRDCom(name),
  socket_kstates("kstates"),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_median_areas("median_areas")
{
}

void Shift::execute()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealMatrix> dh_kstates = socket_kstates.getDataHandle();
  DataHandle<CFreal> rhs  = socket_rhs.getDataHandle();
  DataHandle<CFreal> median_areas = socket_median_areas.getDataHandle();

  const CFuint nbstates = states.size();
  const CFuint nbeqs = PhysicalModelStack::getActive()->getNbEq();

  const CFuint order = getMethodData().getOrder();
  const CFuint k     = getMethodData().K();

  const CFuint om1 = order-1;
  for ( CFuint s = 0; s < nbstates; ++s )
  {
    RealMatrix& kstates = dh_kstates[s];
    for ( CFuint eq = 0; eq < nbeqs; ++eq )
      kstates(eq,om1) = kstates(eq,0) +
                        rhs(s, eq, nbeqs) / median_areas[s];


    /// @todo need to be fixed
    // shifting goes here
    if ( k+2 < order )
    {
      CF_DEBUG_OBJ(k);

      // shift ( o-1-k) to (o-1-k-1) = (o-k-2)
      for ( CFuint sh = k; sh >= 0; --sh)
      {


        const CFuint shto = om1-k-sh-1;
        const CFuint shfr = om1-k-sh;

        CF_DEBUG_OBJ(shto);
        CF_DEBUG_OBJ(shfr);

        for ( CFuint eq = 0; eq < nbeqs; ++eq )
           kstates(eq,shto) = kstates(eq,shfr);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > Shift::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_kstates);
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_median_areas);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

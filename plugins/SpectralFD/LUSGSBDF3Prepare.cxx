#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/LUSGSBDF3Prepare.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LUSGSBDF3Prepare, SpectralFDMethodData, SpectralFDModule> lusgsbdf3PrepareProvider("LUSGSBDF3Prepare");

//////////////////////////////////////////////////////////////////////////////

LUSGSBDF3Prepare::LUSGSBDF3Prepare(const std::string& name) : LUSGSPrepare(name),
  m_3StepsTMSparams(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

LUSGSBDF3Prepare::~LUSGSBDF3Prepare()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSBDF3Prepare::execute()
{
  CFAUTOTRACE;

  // call the general LUSGSPrepare execute command 
  LUSGSPrepare::execute();

  // get subsystemstatus
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  // Derefence and resize the 3 steps time marching scheme parameters (BDF3 with variable time step)
  RealVector& params3StepsTMS = *m_3StepsTMSparams;
  const CFuint nbr3StepsTMSparam = 4;
  params3StepsTMS.resize(nbr3StepsTMSparam);

  // For the first and the second step, used the coefficients of the backward Euler scheme. The BDF3 is not self-starting
  if((subSysStatus->getNbIter() == 1) || (subSysStatus->getNbIter() == 2))
  {
    // coefficient aP1
    params3StepsTMS[0] = 1.0;

    // coefficient aP0PaP1
    params3StepsTMS[1] = 0.0;

    // coefficient aM1
    params3StepsTMS[2] = 0.0;

    // coefficient aM2
    params3StepsTMS[3] = 0.0;
  }
  else
  {
    // get time steps
    const CFreal newDT = subSysStatus->getDT();
    const CFreal oldDT = subSysStatus->getPreviousDT();
    const CFreal oldOldDT = subSysStatus->getPrevPrevDT();

    // The coefficients are multiplied by the newDT to take into account that
    // the diagonal values of the block jacobian are already divided by newDT
    // coefficient aP1
    params3StepsTMS[0] = newDT*(oldDT*oldOldDT + 2*newDT*oldOldDT + oldDT*oldDT + 4*newDT*oldDT + 3*newDT*newDT)/(newDT*(oldDT + newDT)*(oldOldDT + oldDT + newDT));

    // coefficient aP0PaP1
    params3StepsTMS[1] = params3StepsTMS[0] - newDT*(oldDT*oldOldDT + newDT*oldOldDT + oldDT*oldDT + 2*newDT*oldDT + newDT*newDT)/(newDT*oldDT*(oldOldDT + oldDT));

    // coefficient aM1
    params3StepsTMS[2] = newDT*(newDT*oldOldDT + newDT*oldDT + newDT*newDT)/(oldDT*(oldDT + newDT) * oldOldDT);

    // coefficient aM2
    params3StepsTMS[3] = -newDT*(newDT*oldDT + newDT*newDT)/(oldOldDT*(oldOldDT + oldDT)*(oldOldDT + oldDT + newDT));
   }
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSBDF3Prepare::setup()
{
  CFAUTOTRACE;

  // get the 3 steps time marching scheme parameters (BDF3 with variable time step)
  m_3StepsTMSparams = getMethodData().get3StepsTMSParams();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/BDF3TimeDiagBlockJacob.hh"
#include "SpectralFD/LUSGSBDF3Prepare.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BDF3TimeDiagBlockJacob,
                       SpectralFDMethodData,
                       SpectralFDModule>
BDF3_DiagBlockJacob("BDF3TimeDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

BDF3TimeDiagBlockJacob::BDF3TimeDiagBlockJacob(const std::string& name) :
    PseudoSteadyStdTimeDiagBlockJacob(name),
    m_3StepsTMSparams()
{
}

//////////////////////////////////////////////////////////////////////////////

BDF3TimeDiagBlockJacob::~BDF3TimeDiagBlockJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void BDF3TimeDiagBlockJacob::setup()
{
  PseudoSteadyStdTimeDiagBlockJacob::setup();

  // get the 3 steps time marching scheme parameters (BDF3 with variable time step)
  m_3StepsTMSparams = getMethodData().get3StepsTMSParams();
}

//////////////////////////////////////////////////////////////////////////////

void BDF3TimeDiagBlockJacob::addTimeResidual()
{
  // Derefence and resize the 3 steps time marching scheme parameters (BDF3 with variable time step)
  RealVector& params3StepsTMS = *m_3StepsTMSparams;

  // multiply diagonal values with the coefficient aP1 calculated in the LUSGSBDF3Prepare
  m_diagValues *= params3StepsTMS[0];

  // add time residual contribution BDF3
  const CFuint nbrSolPnts = m_cellStates->size();
  CFuint resIdx = 0;
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // add contribution to diagonal block jacobian
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resIdx)
    {
      (*m_currDiagMatrix)(resIdx,resIdx) += m_diagValues[iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
BDF3TimeDiagBlockJacob::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = PseudoSteadyStdTimeDiagBlockJacob::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

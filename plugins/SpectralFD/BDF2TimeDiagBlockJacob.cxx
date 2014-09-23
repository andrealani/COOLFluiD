#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/BDF2TimeDiagBlockJacob.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BDF2TimeDiagBlockJacob,
                       SpectralFDMethodData,
                       SpectralFDModule>
BDF2_DiagBlockJacob("BDF2TimeDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

BDF2TimeDiagBlockJacob::BDF2TimeDiagBlockJacob(const std::string& name) :
    PseudoSteadyStdTimeDiagBlockJacob(name)
{
}

//////////////////////////////////////////////////////////////////////////////

BDF2TimeDiagBlockJacob::~BDF2TimeDiagBlockJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void BDF2TimeDiagBlockJacob::setup()
{
  PseudoSteadyStdTimeDiagBlockJacob::setup();
}

//////////////////////////////////////////////////////////////////////////////

void BDF2TimeDiagBlockJacob::addTimeResidual()
{
  // get factor for the residual (in the Jacobian)
  const CFreal resFactor = getMethodData().getResFactor();
  const CFreal resFactorP1 = resFactor + 1.;

  // multiply diagonal values with the residual factor plus one
  m_diagValues *= resFactorP1;

  // add time residual contribution BDF2
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
BDF2TimeDiagBlockJacob::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = PseudoSteadyStdTimeDiagBlockJacob::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

#include "Environment/ObjectProvider.hh"
#include "MathTools/MathConsts.hh"
#include "Framework/MeshData.hh"
#include "LUSGSMethod/ComputeL2NormLUSGS.hh"
#include "LUSGSMethod/LUSGSMethod.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeL2NormLUSGS, ComputeNorm, LUSGSMethodModule, 1>
    computeL2NormLUSGSProvider("L2LUSGS");

//////////////////////////////////////////////////////////////////////////////

ComputeL2NormLUSGS::ComputeL2NormLUSGS(const std::string & name) :
  ComputeNormLUSGS(name),
  m_gr(*this),
  socket_rhsCurrStatesSet("rhsCurrStatesSet")
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeL2NormLUSGS::~ComputeL2NormLUSGS()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeL2NormLUSGS::GR_Combine (const CFreal & S1,
                                     const CFreal & S2,
                                     CFreal& D)
{
  D = S1 + S2;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ComputeL2NormLUSGS::GR_GetLocalValue () const
{
  return m_localResiduals[m_var_itr];
}

//////////////////////////////////////////////////////////////////////////////

void ComputeL2NormLUSGS::addStatesSetContribution()
{
  // get current states set index
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  const CFuint currStatesSetIdx = statesSetIdx[0];

  // get isStatesSetParUpdatable data handle
  DataHandle< bool > isStatesSetParUpdatable = socket_isStatesSetParUpdatable.getDataHandle();

  // if states set is parallel updatable, add contribution to residual norms
  if (isStatesSetParUpdatable[currStatesSetIdx])
  {
    // get state IDs in current states set
    DataHandle< vector< CFuint > > statesSetStateIDs = socket_statesSetStateIDs.getDataHandle();
    const vector< CFuint >& currStatesIDs = statesSetStateIDs[currStatesSetIdx];
    const CFuint nbStates = currStatesIDs.size();

    // get rhsCurrStatesSet data handle
    DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

    // loop over states in current states set
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {
      // loop over variables for which the residuals should be computed
      for (m_var_itr = 0; m_var_itr < m_residuals.size(); ++m_var_itr)
      {
        const CFreal tmp = rhsCurrStatesSet(iState,m_compute_var_id[m_var_itr],m_nbrEqs);
        m_localResiduals[m_var_itr] += tmp*tmp;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

RealVector ComputeL2NormLUSGS::compute ()
{
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  for(m_var_itr = 0; m_var_itr < m_residuals.size(); ++m_var_itr)
  {
    const CFreal globalValue = m_gr.GetGlobalValue (nsp);
    if(globalValue > 0.)
    {
      m_residuals[m_var_itr] = log10(sqrt(globalValue));
    }
    else
    {
      m_residuals[m_var_itr] = -MathTools::MathConsts::CFrealMax();
    }
  }

  return m_residuals;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeL2NormLUSGS::setup()
{
  CFAUTOTRACE;

  ComputeNormLUSGS::setup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
    ComputeL2NormLUSGS::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeNormLUSGS::needsSockets();

  result.push_back(&socket_rhsCurrStatesSet);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "LUSGSMethod/LUSGSIteratorData.hh"
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

MethodCommandProvider<NullMethodCommand<LUSGSIteratorData>,
                      LUSGSIteratorData, LUSGSMethodModule>
nullLUSGSIteratorComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void LUSGSIteratorData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Norm","L2 Norm of dU to reach to stop the (nonlinear) LU-SGS loop.");
   options.addConfigOption< bool >("PrintHistory","Print convergence history for each (nonlinear) LU-SGS Iterator step");
   options.addConfigOption< vector<CFuint> >("JacobFreezFreq","Number of time-steps to perform in the (nonlinear) LU-SGS iterator before to recompute the block Jacobian matrices.");
   options.addConfigOption< vector<CFuint> >("MaxSweepsPerStep","Maximum number of sweeps to perform in one LU-SGS step.");
}

//////////////////////////////////////////////////////////////////////////////

LUSGSIteratorData::LUSGSIteratorData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner),
    m_computeNormLUSGS(CFNULL),
    m_numJacob(CFNULL),
    m_achieved(false),
    m_forwardSweep(),
    m_stopSweep(),
    m_stopStatesLoop(),
    m_stopEqsLoop(),
    m_beforePertResComputation(),
    m_nbrStatesSets(),
    m_resAux(),
    m_withPivot()
{
  addConfigOptionsTo(this);

  m_jacobFreezFreq = vector<CFuint>();
  setParameter("JacobFreezFreq",&m_jacobFreezFreq);

  m_maxSweeps = vector<CFuint>();
  setParameter("MaxSweepsPerStep",&m_maxSweeps);

  m_maxNorm = -10.;
  setParameter("Norm",&m_maxNorm);

  m_printHistory = false;
  setParameter("PrintHistory",&m_printHistory);
}

//////////////////////////////////////////////////////////////////////////////

LUSGSIteratorData::~LUSGSIteratorData()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSIteratorData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);

  // if the frequency has not been specified, just resize
  // the corresponding vector and set it to 1
  if (m_jacobFreezFreq.size() == 0) {
    m_jacobFreezFreq.resize(1);
    m_jacobFreezFreq[0] = 1;
  }
  cf_assert(m_jacobFreezFreq.size() > 0);
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSIteratorData::setup()
{
  ConvergenceMethodData::setup();

  // dynamic cast the ComputeNorm object to a ComputeNormLUSGS class
  m_computeNormLUSGS = getNormComputer().d_castTo< ComputeNormLUSGS >();

  // create numerical Jacobian computer
  m_numJacob.reset(new NumericalJacobian("NumericalJacobian"));

  // set reference values in numerical Jacobian computer
  RealVector refValues = PhysicalModelStack::getActive()->getImplementor()->getRefStateValues();
  m_numJacob->setRefValues(refValues);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RusanovSchemeCSys.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RusanovSchemeCSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
rusanovSchemeCSysProvider("SysRusanovC");

//////////////////////////////////////////////////////////////////////////////

RusanovSchemeCSys::RusanovSchemeCSys(const std::string& name) :
  RDS_SplitterSys(name),
  m_alpha(),
  m_maxEV(),
  m_minEV(),
  m_sumUmin()
{
}

//////////////////////////////////////////////////////////////////////////////

RusanovSchemeCSys::~RusanovSchemeCSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void RusanovSchemeCSys::setup()
{
  RDS_SplitterSys::setup();

  m_sumUmin.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void RusanovSchemeCSys::distribute(vector<RealVector>& residual)
{
   cf_assert(_nbStatesInCell != 0);

   const RealVector& phiT = getMethodData().getDistributionData().phi;
   const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
   const CFreal invnstates = 1.0 / _nbStatesInCell;

   // calculate the alpha artificial diffusion coefficient
   // this should be the
   m_maxEV = m_eValues[0]->max();
   m_minEV = m_eValues[0]->min();
   for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
     m_maxEV = std::max(m_maxEV,m_eValues[iState]->max());
     m_minEV = std::min(m_minEV,m_eValues[iState]->min());
   }

   // maximum spectral radius of all matrices Kj for
   m_alpha = std::abs(m_maxEV - m_minEV);

   // calculate the residual distribution
   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

      // calculate the diffusion value
      m_sumUmin = (*tStates[iState] - *tStates[0]);
      for (CFuint jState = 1; jState < _nbStatesInCell; ++jState) {
        m_sumUmin += (*tStates[iState] - *tStates[jState]);
      }

      residual[iState] = m_alpha * m_sumUmin;
      residual[iState] += phiT;
      residual[iState] *= invnstates;
   }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

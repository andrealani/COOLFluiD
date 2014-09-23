#ifndef COOLFluiD_Numerics_ExplicitFilters_FilterWeight_hh
#define COOLFluiD_Numerics_ExplicitFilters_FilterWeight_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MathConsts.hh"
#include <deque>


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  	namespace ExplicitFilters {
      
//////////////////////////////////////////////////////////////////////////////

class FilterWeight
{
  public:
    FilterWeight() : m_mustCompute(true) {}
    bool                           mustCompute()const { return m_mustCompute; }
    void                           clear() { m_weights.resize(0); }
    void                           setCompute(const bool mustCompute) { m_mustCompute = mustCompute;}
    Common::SafePtr<RealVector >   getWeights() { return &m_weights; }
    void                           setWeights(const RealVector& weights) { m_weights.resize(weights.size()); m_weights = weights; } 
    CFreal                         getWeight(const CFuint& idx) const { return m_weights[idx];}
    CFuint                         getNbElements() const { return m_weights.size(); }
  private:
    bool         m_mustCompute;
    RealVector   m_weights;
};

//////////////////////////////////////////////////////////////////////////////

		} // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_FilterWeight_hh

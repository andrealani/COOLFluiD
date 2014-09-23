#ifndef COOLFluiD_FluctSplit_NullIntegral_hh
#define COOLFluiD_FluctSplit_NullIntegral_hh

#include "Common/StringOps.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

/// Trait class to does nothing when called to integrate the residual
/// @author Tiago Quintino
class NullResidual
{
  public:

  template < typename ELEMGEO, typename PHYSICS >
  struct Integrator
  {
    enum { NBQDPT = 0 };

    static std::string getClassName () { return "NullResidual"; };

    template < unsigned int SELEM, typename FSDATA >
    static inline void exec ( FSDATA& data )
    {
//       CFout << "NullResidual::Integrator::compute() does nothing\n" << CFendl;
    }
  };
};

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_NullIntegral_hh


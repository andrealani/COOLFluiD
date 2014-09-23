#ifndef COOLFluiD_FluctSplit_LRD_hh
#define COOLFluiD_FluctSplit_LRD_hh

#include "FluctSplit/MetaSchemes/NullIntegral.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

/// Trait class that describes the residual computation and distribution
/// based on a linearized residual computed from the k upwind parameters.
/// This algorithms does not use any integration.
/// @author Tiago Quintino
template < typename SPLITTER, typename INTEGRATION = NullResidual >
class LRD
{
  public: // typedefs
    typedef SPLITTER    Splitter_type;
    typedef INTEGRATION Integration_type;

  public: // functions

    static std::string getClassName () { return "LRD"; };

    /// Distributes the subelem residual using the splitter
    template < unsigned int SELEM, typename FSDATA >
    static inline void distribute_residual ( FSDATA& data )
    {
      data.splitter.template computeK <SELEM, FSDATA> (data);
      data.splitter.template distribute_linres <SELEM, FSDATA> (data);
    }

    /// Computes the subelem residual an stores it in the residual variable of SubElem.
    /// This implementation is empty because residual is computed from k upwind parameters
    /// during the distribute function of the scheme
    /// Distributes the in the subelem residual using the splitter
    template < typename FSDATA >
    static inline void compute_residual ( FSDATA& data )
    {
//       CFout << "compute_residual() does nothing in LRD\n" << CFendl;
    }
};

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_LRD_hh


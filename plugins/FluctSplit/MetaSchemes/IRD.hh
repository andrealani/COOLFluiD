#ifndef COOLFluiD_FluctSplit_IRD_hh
#define COOLFluiD_FluctSplit_IRD_hh

#include "FluctSplit/MetaSchemes/VolumeIntegral.hh"
#include "FluctSplit/MetaSchemes/ContourIntegral.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

/// Trait class that describes the residual computation and distribution
/// based on a integrated subelement residual computed at quadrature points.
/// This algorithm may use VolumeResidual integration  or ContourResidual integration.
/// @author Tiago Quintino
template < typename SPLITTER, typename INTEGRATION >
class IRD
{
  public:
    typedef SPLITTER    Splitter_type;
    typedef INTEGRATION Integration_type;

  public: // functions

    static std::string getClassName () { return "IRD"; };

    /// Distributes the subelem residual using the splitter
    template < unsigned int SELEM, typename FSDATA >
    static inline void distribute_residual ( FSDATA& data )
    {
      data.splitter.template computeK <SELEM, FSDATA> (data);
      data.splitter.template distribute_intres <SELEM, FSDATA> (data);
    }

    /// Computes the subelem residual an stores it in the residual variable of SubElem.
    /// This implementation is empty because residual is computed from k upwind parameters
    /// during the distribute function of the scheme
    template < typename FSDATA >
    static inline void compute_residual ( FSDATA& data )
    {
//       CFout << "compute_residual() integrates in IRD\n" << CFendl;
      FSDATA::Integrator_type::compute_all_quad_pts ( data );
      Common::Loop1<typename FSDATA::Integrator_type, FSDATA::NBSUBELEM>::run( data );
    }
};

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_IRD_hh


#ifndef COOLFluiD_FluctSplit_VolumeIntegral_hh
#define COOLFluiD_FluctSplit_VolumeIntegral_hh

#include "Common/StringOps.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

/// Trait class to compute the residual via volume integration
/// @author Tiago Quintino
template < unsigned int TORDER >
class VolumeResidual
{
  public:

  enum { ORDER = TORDER };

  template < typename ELEMGEO, typename PHYSICS >
  struct Integrator
  {
    /// @todo to implement the volume integration trait class which defines number of quadrature points
    enum { NBQDPT = 0 };

    static std::string getClassName () { return "VolumeResidual" + Common::StringOps::to_str((int) ORDER); };

    /// Computes integration data on all quadrature points
    template < typename FSDATA >
    static void compute_all_quad_pts ( FSDATA& data )
    {
//       CFout << "VolumeResidual::Integrator::compute_all_quad_pts()\n" << CFendl;
      // loop on all quadrature points
      // compute residual on each
    }

    /// Assemble the quadrature points contribution to
    /// the subelem residual
    template < unsigned int SELEM, typename FSDATA >
    static inline void exec ( FSDATA& data )
    {
//       CFout << "VolumeResidual::Integrator::exec() subelem [" << SELEM << "]\n" << CFendl;
      // on this subelem
      // loop on quad points that make subelem
      // sum contribution from all quads to subelem res
    }
  };
};

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_VolumeIntegral_hh


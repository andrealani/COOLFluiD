#ifndef COOLFluiD_FluctSplit_ContourIntegral_hh
#define COOLFluiD_FluctSplit_ContourIntegral_hh

#include "Common/StringOps.hh"
#include "Common/Meta/Loop.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

/// @todo Support for non isoparametric elements
template < typename SOLSF, unsigned int ORDER >
struct CIntegral
{
};

/// Trait class to compute the residual via contour integration
/// @author Tiago Quintino
template < unsigned int TORDER >
class ContourResidual
{
  public:

  enum { ORDER = TORDER };

  template < typename ELEMGEO, typename PHYSICS >
  struct Integrator
  {
    typedef CIntegral< typename ELEMGEO::SolSF_type, ORDER> Integral_type;

    enum { NBQDPT    = Integral_type::NBQDPT    };
    enum { NBSUBQDPT = Integral_type::NBSUBQDPT };

    static std::string getClassName () { return "ContourResidual" + Common::StringOps::to_str((int) ORDER); };

    template < typename INTEGRAL, unsigned int QDPT >
    struct PrintSF
    {
      template < unsigned int NODE >
      static void exec ()
      {
        typedef typename INTEGRAL::template QdPt<QDPT> Coord_t;
        typedef typename INTEGRAL::Elem::template ShapeF<NODE> ShapeF_t;
        PRINT_DEBUG ( CFout << "SQdPt SF " << NODE << " [" << ShapeF_t::template value< Coord_t >() << "]\n"; );
      }
    };

    template < typename INTEGRAL >
    struct PrintQD
    {
      template < unsigned int QDPT >
      static void exec ()
      {
        // print the weight
        PRINT_DEBUG ( CFout << "QdPt [" << QDPT << "] weight [" << INTEGRAL::template QdPt<QDPT>::weight() << "] xy [" << INTEGRAL::template QdPt<QDPT>::xi() << "," << INTEGRAL::template QdPt<QDPT>::eta() << "]\n"; );
        // print the values at each quadratrue poitn of each shape functions
        Common::Loop< PrintSF< INTEGRAL, QDPT >, INTEGRAL::NBSF >::run();
        PRINT_DEBUG ( CFout << "  -------  \n"; );
      }
    };

    /// Class to compute the physical model related data at each quadrature point
    /// @author Tiago Quintino
    template < typename INTEGRAL >
    struct ComputePhysicalData
    {
        /// Class to sum up the interpolation of the states.
        /// Specialization to close the loop.
        template < unsigned int QDPT >
        struct InterpolateStates
        {
          template < unsigned int NBSF, typename FSDATA >
          static void exec ( FSDATA& d )
          {
            typedef typename INTEGRAL::template QdPt<QDPT> QdPtCoord;
            typedef typename INTEGRAL::Elem::template ShapeF<NBSF> ShapeF;

            d.qdpts[QDPT].state += ShapeF::template value<QdPtCoord>() * *(d.states[NBSF]);
          }
        };

        /// Class to sum up the interpolation of the states.
        /// Specialization to close the loop.
        template < unsigned int QDPT >
        struct InterpolateNodes
        {
          template < unsigned int NBSF, typename FSDATA >
          static void exec ( FSDATA& d )
          {
            typedef typename INTEGRAL::template QdPt<QDPT> QdPtCoord;
            typedef typename INTEGRAL::Elem::template ShapeF<NBSF> ShapeF;
            d.qdpts[QDPT].node += ShapeF::template value<QdPtCoord>() * *(d.nodes[NBSF]);
          }
        };

      /// Function that gets called for each quadrature point
      template < unsigned int QDPT, typename FSDATA >
      static void exec ( FSDATA& data )
      {
//         CFout << "\n";
        // interpolate states @ quadrature points
        data.qdpts[QDPT].state = 0.0;
        Common::Loop1<InterpolateStates<QDPT>, INTEGRAL::NBSF>::run(data);
//         CFout << "\n";
        // interpolate nodes @ quadrature points
        data.qdpts[QDPT].node = 0.0;
        Common::Loop1<InterpolateNodes<QDPT>, INTEGRAL::NBSF>::run(data);
//         CFout << "\n";

        // compute the data at each quadrature point
        data.linear_varset->setExtraPhysicalVars(&data.qdpts[QDPT].extra_vars);
        data.linear_varset->computePhysicalData(data.qdpts[QDPT].state, data.qdpts[QDPT].physdata);

        // compute the flux at each quadrature point
        PHYSICS::PhysicalModel_type::computeFlux(data.qdpts[QDPT].state, data.qdpts[QDPT].physdata, data.qdpts[QDPT].flux);
      }
    };

    /// Class to sum each quadrature point flux contribution to each subelement
    /// @author Tiago Quintino
    template < typename INTEGRAL, unsigned int SELEM >
    struct SumQuadPointFlux
    {
      /// This function is called for each quadrature point
      template < unsigned int QDPT, typename FSDATA >
      static void exec ( FSDATA& data )
      {
        typedef typename INTEGRAL::template SubElem<SELEM> SE_t;
        typedef typename SE_t::template QdPt<QDPT> SEQD_t;
        typedef typename INTEGRAL::template QdPt<SEQD_t::ID> QdPt_t;

       // weight * face_jacobian * flux * normal
        data.sub_elem[SELEM].residual -=
          QdPt_t::weight() *
          SEQD_t::face_jacob(data) * ( data.qdpts[SEQD_t::ID].flux * SEQD_t::normal(data) );

        PRINT_DEBUG ( CFout << " *** normal elem [" << SELEM << "] " << SEQD_t::normal(data) << "\n"; );
        PRINT_DEBUG ( CFout << " *** weight      [" << SELEM << "] " << QdPt_t::weight()  << "\n"; );
        PRINT_DEBUG ( CFout << " *** face_jacob  [" << SELEM << "] " << SEQD_t::face_jacob(data)  << "\n"; );
        PRINT_DEBUG ( CFout << " *** normal area [" << SELEM << "] " << QdPt_t::weight() * SEQD_t::face_jacob(data) << "\n"; )
        PRINT_DEBUG ( RealVector fflux ( data.sub_elem[SELEM].residual );  fflux = QdPt_t::weight() * SEQD_t::face_jacob(data) * ( data.qdpts[SEQD_t::ID].flux * SEQD_t::normal(data) );  CFout << " *** FLUX [" << QDPT << "] " << fflux << "\n"; );
      }
    };

    /// Computes integration data on all quadrature points
    template < typename FSDATA >
    static void compute_all_quad_pts ( FSDATA& data )
    {
      PRINT_DEBUG ( CFout << "ContourResidual::Integrator::compute_all_quad_pts()\n" << CFendl; )
      // Common::Loop< PrintQD< Integral_type >, NBQDPT >::run();
      // loop on all quadrature points and compute flux on them
      Common::Loop1< ComputePhysicalData< Integral_type >, NBQDPT >::run(data);
    }

    /// Assemble the quadrature points contribution to
    /// the subelem residual
    template < unsigned int SELEM, typename FSDATA >
    static inline void exec ( FSDATA& data )
    {
      PRINT_DEBUG ( CFout << "ContourResidual::Integrator::exec() subelem [" << SELEM << "]\n" << CFendl; )
      // on this subelem
      // loop on quad points that make subelem
      // project them on normal and sum the flux
      data.sub_elem[SELEM].residual = 0.0;
      Common::Loop1< SumQuadPointFlux< Integral_type, SELEM >, NBSUBQDPT >::run(data);
      PRINT_DEBUG ( CFout << "+++++ phi [" << data.sub_elem[SELEM].residual << "]\n"; )
    }

  };
};

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_ContourIntegral_hh


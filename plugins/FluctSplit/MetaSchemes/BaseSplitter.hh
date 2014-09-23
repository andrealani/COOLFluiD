#ifndef COOLFluiD_FluctSplit_BaseSplitter_hh
#define COOLFluiD_FluctSplit_BaseSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// Base class for numerical schemes
/// @author Tiago Quintino
template < typename ELEMGEO , typename PHYSICS >
class BaseSplitter
{
    public: // typedefs

    enum { NBSTATES = ELEMGEO::NBSTATES };
    enum { NEQS     = PHYSICS::NEQS };

    public: // functions

    /// Computes the K upwind coefficients
    template < unsigned int SELEM, typename FSDATA > inline void computeK ( FSDATA& data );

    protected: // helper classes

      template < unsigned int SELEM, typename FSDATA >
      struct ComputeStateK
      {
        template < unsigned int STATE >
        static inline void exec ( FSDATA& data )
        {
          PRINT_DEBUG ( CFout << "ComputeStateK selem [" << SELEM << "] state ["<< STATE << "]\n" << CFendl; )

          std::vector<CFreal>& nodal_areas = data.sub_elem[SELEM].nodal_areas;
          std::vector<RealMatrix>& kplus   = data.kplus;
          std::vector<RealMatrix>& kmin    = data.kmin;
          std::vector<RealVector>& evalues = data.evalues;

          data.distribute_varset->splitJacobian( kplus[STATE], kmin[STATE], evalues[STATE], data.sub_elem[SELEM].nodal_normals[STATE] );

          const CFreal kcoeff = ( 1. / ELEMGEO::DIM ) * nodal_areas[STATE];

          kplus[STATE] *= kcoeff;
          kmin [STATE] *= kcoeff;
        }
      };

      template < unsigned int SELEM, typename FSDATA >
      struct ComputeUpdateCoeff
      {
        template < unsigned int STATE >
        static inline void exec ( FSDATA& data )
        {
//           CFout << "ComputeUpdateCoeff selem [" << SELEM << "] state ["<< STATE << "]\n" << CFendl;
          const CFreal max_evalue = std::max(0.0, data.evalues[STATE].max());
          const CFreal kcoeff = ( 1. / ELEMGEO::DIM ) * data.sub_elem[SELEM].nodal_areas[STATE];

          const CFuint local_id = data.sub_elem[SELEM].states[STATE]->getLocalID();
          data.updateCoeff[local_id] += kcoeff * max_evalue;
        }
      };
};

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO , typename PHYSICS >
template < unsigned int SELEM, typename FSDATA >
inline void BaseSplitter<ELEMGEO,PHYSICS>::computeK ( FSDATA& data )
{
  using namespace COOLFluiD::Framework;
  // loop on all the states and compute each state k parameters
  Common::Loop1<ComputeStateK<SELEM,FSDATA>, ELEMGEO::NBSUBSTATES>::run(data);


// skip update coeff in implicit numerical perturbations
/// @todo FSMHO: make it run with implicit
// if (!getMethodData().getDistributionData().isPerturb)
// {
   Common::Loop1<ComputeUpdateCoeff<SELEM,FSDATA>, ELEMGEO::NBSUBSTATES>::run(data);
// }

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_BaseSplitter_hh


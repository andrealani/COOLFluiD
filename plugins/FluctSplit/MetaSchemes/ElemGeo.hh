#ifndef COOLFluiD_FluctSplit_ElemGeo_hh
#define COOLFluiD_FluctSplit_ElemGeo_hh

#include <boost/static_assert.hpp>

#include "FluctSplit/MetaSchemes/ShapeFunctions.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

template <typename GEOSF , typename SOLSF>
struct ElemGeo
{
  public: // typedefs

    typedef GEOSF GeoSF_type;
    typedef SOLSF SolSF_type;

    // information about element
    enum { DIM      = GEOSF::DIM     };
    enum { SHAPE    = GEOSF::SHAPE   };
    enum { NBNODES  = GEOSF::NBNODES };
    enum { NBSTATES = SOLSF::NBNODES };
    enum { GEOORDER = GEOSF::ORDER   };
    enum { SOLORDER = SOLSF::ORDER   };

    // information about subelement
    enum { NBSUBELEM   = SOLSF::NBSUBELEM };
    enum { NBSUBSTATES = SOLSF::NBSUBNODES };

  /// Constructor
  ElemGeo ()
  {
    // some sanity checks
    BOOST_STATIC_ASSERT (GEOSF::DIM   == SOLSF::DIM);
    BOOST_STATIC_ASSERT (GEOSF::SHAPE == SOLSF::SHAPE);
  }

};

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_ElemGeo_hh

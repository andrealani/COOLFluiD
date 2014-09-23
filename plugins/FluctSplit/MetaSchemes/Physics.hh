#ifndef COOLFluiD_FluctSplit_Physics_hh
#define COOLFluiD_FluctSplit_Physics_hh

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

template < typename MODEL >
class Physics
{
  public: // typedfs

    enum { DIM      = MODEL::DIM      }; // dimensionality of the physics
    enum { NEQS     = MODEL::NEQS     }; // number of equations
    enum { PHYSDATA = MODEL::PHYSDATA }; // size of the physical data array
    enum { XVARS    = MODEL::XVARS    }; // size of the extra vars array

    typedef MODEL PhysicalModel_type;
};

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_Physics_hh

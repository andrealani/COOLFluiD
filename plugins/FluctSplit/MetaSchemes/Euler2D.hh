#ifndef COOLFluiD_FluctSplit_Euler2D_hh
#define COOLFluiD_FluctSplit_Euler2D_hh

#include "Common/StringOps.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

struct Euler2D
{
  enum {RHO=0, P=1, H=2, E=3, A=4, T=5, V=6, VX=7, VY=8, VZ=9, GAMMA=10};

  enum { DIM      = DIM_2D }; // dimensionality of the physics
  enum { NEQS     = 4      }; // dim momentum + continuity and energy
  enum { PHYSDATA = 11     }; // size of the physical data array
  enum { XVARS    = 0      }; // size of the extra vars array

  static std::string getClassName () { return "Euler2D"; };

  static void computeFlux ( const Framework::State& vars, const RealVector& physdata, RealMatrix& flux)
  {
    const CFreal u     = physdata[VX];
    const CFreal v     = physdata[VY];
    const CFreal rhoU  = physdata[RHO]*u;
    const CFreal rhoV  = physdata[RHO]*v;
    const CFreal rhoUU = rhoU*u;
    const CFreal rhoUV = rhoU*v;
    const CFreal rhoVV = rhoV*v;

    flux(0,XX) = rhoU;
    flux(0,YY) = rhoV;

    flux(1,XX) = physdata[P] + rhoUU;
    flux(1,YY) = rhoUV;

    flux(2,XX) = rhoUV;
    flux(2,YY) = physdata[P] + rhoVV;

    flux(3,XX) = rhoU*physdata[H];
    flux(3,YY) = rhoV*physdata[H];
  }

}; // end class Euler2D

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_Euler2D_hh

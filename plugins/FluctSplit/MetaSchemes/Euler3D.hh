#ifndef COOLFluiD_FluctSplit_Euler3D_hh
#define COOLFluiD_FluctSplit_Euler3D_hh

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

struct Euler3D
{
  enum {RHO=0, P=1, H=2, E=3, A=4, T=5, V=6, VX=7, VY=8, VZ=9, GAMMA=10};

  enum { DIM      = DIM_3D }; // dimensionality of the physics
  enum { NEQS     = 5      }; // dim momentum + continuity and energy
  enum { PHYSDATA = 11     }; // size of the physical data array
  enum { XVARS    = 0      }; // size of the extra vars array

  static std::string getClassName () { return "Euler3D"; };

  static void computeFlux ( const Framework::State& vars, const RealVector& physdata, RealMatrix& flux)
  {
    const CFreal u     = physdata[VX];
    const CFreal v     = physdata[VY];
    const CFreal w     = physdata[VZ];
    const CFreal rhoU  = physdata[RHO]*u;
    const CFreal rhoV  = physdata[RHO]*v;
    const CFreal rhoW  = physdata[RHO]*w;
    const CFreal rhoUU = rhoU*u;
    const CFreal rhoUV = rhoU*v;
    const CFreal rhoVV = rhoV*v;
    const CFreal rhoUW = rhoU*w;
    const CFreal rhoVW = rhoV*w;
    const CFreal rhoWW = rhoW*w;
    const CFreal press = physdata[P];
    const CFreal ht    = physdata[H];

    flux(0,XX) = rhoU;
    flux(0,YY) = rhoV;
    flux(0,ZZ) = rhoW;

    flux(1,XX) = rhoUU + press;
    flux(1,YY) = rhoUV;
    flux(1,ZZ) = rhoUW;

    flux(2,XX) = rhoUV;
    flux(2,YY) = rhoVV + press;
    flux(2,ZZ) = rhoVW;

    flux(3,XX) = rhoUW;
    flux(3,YY) = rhoVW;
    flux(3,ZZ) = rhoWW + press;

    flux(4,XX) = rhoU*ht;
    flux(4,YY) = rhoV*ht;
    flux(4,ZZ) = rhoW*ht;
  }

}; // end class Euler3D

///////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit
} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_Euler3D_hh

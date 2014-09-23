#ifndef COOLFluiD_UFEM_AssemblyData_hh
#define COOLFluiD_UFEM_AssemblyData_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { namespace Framework { class GeometricEntity; } }

namespace COOLFluiD {
    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

struct AssemblyData
{

  /// constructor
  AssemblyData ( const CFuint& blocksize ) :
  T   (blocksize, blocksize),
  A   (blocksize, blocksize),
  dTu (blocksize, blocksize),
  dAu (blocksize, blocksize),
  db  (blocksize, blocksize),
  b   (blocksize)
  {}

  /// zero the data
  void clean ()
  {
    CFAUTOTRACE;
    T   = 0.;
    A   = 0.;
    dTu = 0.;
    dAu = 0.;
    db  = 0.;
    b   = 0.;
  }

  RealMatrix T;
  RealMatrix A;
  RealMatrix dTu;
  RealMatrix dAu;
  RealMatrix db;
  RealVector b;

};

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_AssemblyData_hh


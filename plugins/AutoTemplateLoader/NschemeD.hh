#ifndef COOLFluiD_Benchmark_NschemeD_hh
#define COOLFluiD_Benchmark_NschemeD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Benchmark {

//////////////////////////////////////////////////////////////////////////////

/// This class benchmarks the NschemeD
/// @author Tiago Quintino
template < unsigned int NEQS, unsigned int NSTATES >
class NschemeD
{
  public:

  typedef RealMatrix Matrix_t;
  typedef RealVector Vector_t;

  enum { E_NEQS = NEQS };
  enum { E_NSTATES = NSTATES };

  static void resize ( Matrix_t& mat , CFuint m, CFuint n) { mat.resize(m,n); for( CFuint e = 0; e < NEQS; ++e) mat(e,e) = 1.0; }
  static void resize ( Vector_t& vec , CFuint m) { vec.resize(m); }
  static void init ( Matrix_t& mat ) {
    for( CFuint i = 0; i < NEQS; ++i)
      for( CFuint j = 0; j < NEQS; ++j)
        mat(i,j) = 1.0+i*j;
  }
  static void init ( Vector_t& vec ) {
    for( CFuint i = 0; i < NEQS; ++i)
      vec[i] = i;
  }

  NschemeD();
  
  ~NschemeD();
  
  void distribute ( const std::vector<Matrix_t*>& kplus,
                    const std::vector<Vector_t*>& states,
                    const RealVector&  res,
                          std::vector<Vector_t*>& st_res );


private:
    
    Vector_t skplus_u;
    Matrix_t skplus;
    Matrix_t inv_k;
    Vector_t u_in;

    CFuint nbstates;
    CFuint nbeqs;

    MathTools::MatrixInverterT<NEQS> inverter;

}; // end of class NschemeD

//////////////////////////////////////////////////////////////////////////////

  } // namespace Benchmark

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

# include "NschemeD.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Benchmark_NschemeD_hh

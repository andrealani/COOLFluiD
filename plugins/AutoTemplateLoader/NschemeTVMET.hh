#ifndef COOLFluiD_Benchmark_NschemeTVMET_hh
#define COOLFluiD_Benchmark_NschemeTVMET_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"

#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>

#include "MatrixInverterTVMET.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Benchmark {

//////////////////////////////////////////////////////////////////////////////

/// This class benchmarks the NschemeTVMET
/// @author Tiago Quintino
template < unsigned int NEQS, unsigned int NSTATES >
class NschemeTVMET
{
  public:

  typedef tvmet::Matrix<CFreal,NEQS,NEQS> Matrix_t;
  typedef tvmet::Vector<CFreal,NEQS>      Vector_t;

  enum { E_NEQS = NEQS };
  enum { E_NSTATES = NSTATES };

  static void resize ( Matrix_t& mat , CFuint m, CFuint n) { }
  static void resize ( Vector_t& vec , CFuint m) { }
  static void init ( Matrix_t& mat ) {
    for( CFuint i = 0; i < NEQS; ++i)
      for( CFuint j = 0; j < NEQS; ++j)
        mat(i,j) = 1.0+i*j;
  }
  static void init ( Vector_t& vec ) {
    for( CFuint i = 0; i < NEQS; ++i)
      vec(i) = i;
  }

  NschemeTVMET();
  
  ~NschemeTVMET();
  
  void distribute ( const std::vector<Matrix_t*>& kplus,
                    const std::vector<Vector_t*>& states,
                    const Vector_t&  res,
                          std::vector<Vector_t*>& st_res );


private:
    
    Vector_t skplus_u;
    Matrix_t skplus;
    Matrix_t inv_k;
    Vector_t u_in;

    MathTools::MatrixInverterTVMET<NEQS> inverter;

}; // end of class NschemeTVMET

//////////////////////////////////////////////////////////////////////////////

  } // namespace Benchmark

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

# include "NschemeTVMET.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Benchmark_NschemeTVMET_hh

#ifndef COOLFluiD_MathTools_LeastSquaresSolver_hh
#define COOLFluiD_MathTools_LeastSquaresSolver_hh

#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/LUInverter.hh"

namespace COOLFluiD{
  
  namespace MathTools{

class LeastSquaresSolver{
public:

  //default constructor with number of paramteres to estimate
  LeastSquaresSolver(CFuint nbParameters);

  //Solve the system
 void solve(RealMatrix &coefficients, RealVector &weights, 
           RealVector  &dependentVar, RealVector &solution);

private:
  
 CFuint m_nbParameters;
 RealMatrix m_A;
 RealMatrix m_invA;
 
 RealVector m_b;

 LUInverter m_inverter;

};


LeastSquaresSolver::LeastSquaresSolver(CFuint nbParameters) : 
  m_nbParameters(nbParameters), 
  m_A(nbParameters,nbParameters), 
  m_invA(nbParameters,nbParameters),
  m_b(nbParameters),
  m_inverter(nbParameters)
{

}

void LeastSquaresSolver::solve(RealMatrix &coefficients, RealVector &weights, 
                               RealVector &dependentVar, RealVector &solution)
{
  CFuint nbPoints = coefficients.nbRows();
  cf_assert(coefficients.nbCols() == m_nbParameters);


  cf_assert(           nbPoints >= m_nbParameters );
  cf_assert(    solution.size() == m_nbParameters ); 
  cf_assert(     weights.size() == nbPoints       );
  cf_assert(dependentVar.size() == nbPoints       );
 
  // Assemble matrix A
  // A_jk = sum_i..nbPoints coefficients_ij  weights_ii  coefficients_ik 
  for( CFuint j = 0; j<m_nbParameters; ++j  )
  {
    for( CFuint k = 0; k<m_nbParameters; ++k )
    {
      m_A(j,k) = 0.;
      for(CFuint i = 0; i <nbPoints; ++i)
      {
        m_A(j,k) += coefficients(i,j) * weights[i] * coefficients(i,k);
      }
    }
  }

  // Assemble vector b
  // b_j = sum_i...nbPoints coefficients_ij weights_ii dependentVar_i
  for(CFuint j = 0; j < m_nbParameters; ++j)
  {
    m_b[j] = 0.;
    for(CFuint i = 0; i < nbPoints; ++i)
    {
      m_b[j] += coefficients(i,j) * weights[i] * dependentVar[i];
    }
  }
 
  //Solve Ax = b
 
  //try 1: use LU decomposition and perform x = inverse(A)*b 
  m_inverter.invert(m_A, m_invA);
  
  for(CFuint i = 0; i< m_nbParameters; ++i)
  {
    solution[i] = 0;
    for(CFuint j = 0; j<m_nbParameters; ++j)
    {
      solution[i] += m_invA(i,j) * m_b[j];
    }
  }

}


//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_LeastSquaresSolver


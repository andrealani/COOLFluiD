#include <iomanip>
#include "ReconstructionFilter2D.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "ExplicitFilters.hh"

//#include "Environment/Factory.hh"

//#include "Framework/VolumeIntegrator.hh"
//#include "Framework/MeshData.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<ReconstructionFilter2D, 
                                  FilterData, 
                                  FilterStrategy, 
                                  ExplicitFiltersModule> 
leastSquaresReconstructionFilter("ReconstructionFilter2D");

// //////////////////////////////////////////////////////////////////////////////
// 
// void ReconstructionFilter2D::defineConfigOptions(Config::OptionList& options)
// {
//  options.addConfigOption< CFreal >("RelaxationFactor","Factor to add weight to the centre cell");
//  options.addConfigOption< CFreal >("WeightingFactor","Weighting factor to influence filter width");
// }
		
//////////////////////////////////////////////////////////////////////////////
		
ReconstructionFilter2D::ReconstructionFilter2D(const std::string& name) :
  ReconstructionFilter(name)
{
}

//////////////////////////////////////////////////////////////////////////////


ReconstructionFilter2D::~ReconstructionFilter2D()
{
}

//////////////////////////////////////////////////////////////////////////////


void ReconstructionFilter2D::setup()
{
  CFAUTOTRACE;
  ReconstructionFilter::setup();
}

//////////////////////////////////////////////////////////////////////////////

std::string ReconstructionFilter2D::reconstructionEquation(const CFuint order) const
{    
  std::string equation = "\\phi_i = \\overline{\\phi}_0";
  for (CFuint k=1; k<=order; k++) {
    // for all terms of degree k
      
     for (CFuint n2=0; n2<=k; n2++) {
      for (CFuint n1=0; n1<=k; n1++) {
        if (n1+n2 == k) {
          std::stringstream term, numerator, denominator;
          if (k==1) {
            term << " + " ;
            numerator << "\\partial \\phi_0";
          }
          else {
            term << " + \\frac{" << binomialCoefficient(n1,n2) << "}{" << k << "!} ";
            numerator << "\\partial^" << k << "\\phi_0";
          }
          if (n1>0) {
            if (n1==1) {
              term << " \\Delta x_{0i}";                              
              denominator << " \\partial x";
            }    
            else {
              term << " \\Delta x_{0i}^" << n1;
              denominator << " \\partial x^" << n1;
            }
          }
          if (n2>0) {
            if (n2==1) {
              term << " \\Delta y_{0i}";
              denominator << " \\partial y";                          
            }                      
            else {
              term << " \\Delta y_{0i}^" << n2;
              denominator << " \\partial y^" << n2;                              
            }
          }
          
          // assemble the term                        
          term << "\\frac{" << numerator.str() << "}{" << denominator.str() << "}";
          
          // add term to the equation
          equation += term.str();
          
          // clear the stringstreams
          term.str(std::string());
          numerator.str(std::string());
          denominator.str(std::string());
        
        }
      }
    } 
  }
  return equation;
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter2D::calculateSystemMatrix(RealMatrix& A)
{
  CFuint nbCells = m_stencil->getNbElements();
  Framework::Node& centreCell = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,0);
	CFuint j=0;
  CFreal dx, dy, coefficient;
  CFuint order = getMethodData().getOrder();
	
  for (CFuint k=0; k<=order; k++) {
  	// for all terms of degree k
    for (CFuint n2=0; n2<=k; n2++) {
      for (CFuint n1=0; n1<=k; n1++) {
        if (n1+n2 == k) {
           
          // fill this column
          coefficient = binomialCoefficient(n1,n2)/CFreal(MathTools::MathFunctions::faculty(k));
          for (CFuint i=0; i<nbCells; i++) {
            Framework::Node& iCell = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,i);
            dx = iCell[XX] - centreCell[XX];
            dy = iCell[YY] - centreCell[YY];
            A(i,j) = coefficient * pow(dx,n1)*pow(dy,n2);               
          }
          j++; // go to next column
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal ReconstructionFilter2D::binomialCoefficient(const CFuint n1, const CFuint n2) const
{
  //  (n1 + n2)!                       (n1 + n2 + n3 + ... + nk)!   
  //  ----------         multinomial:  ---------------------------  
  //   n1!  n2!                          n1!  n2!  n3!  ...  nk!    
  
  return ( CFreal(MathFunctions::faculty(n1+n2))/CFreal( MathFunctions::faculty(n1) * MathFunctions::faculty(n2)) );
}

//////////////////////////////////////////////////////////////////////////////

	  } // namespace ExplicitFilters
		
  } // namespace Numerics

} // namespace COOLFluiD

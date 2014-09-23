
#include "ReconstructionFilter3D.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "ExplicitFilters.hh"


//#include "Environment/Factory.hh"

//#include "Framework/VolumeIntegrator.hh"
//#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<ReconstructionFilter3D, 
                                  FilterData, 
                                  FilterStrategy, 
                                  ExplicitFiltersModule> 
leastSquaresReconstructionFilter3D("ReconstructionFilter3D");

//////////////////////////////////////////////////////////////////////////////

// void ReconstructionFilter3D::defineConfigOptions(Config::OptionList& options)
// {
//  options.addConfigOption< CFreal >("RelaxationFactor","Factor to add weight to the centre cell");
//  options.addConfigOption< CFreal >("WeightingFactor","Weighting factor to influence filter width");
// }
		
//////////////////////////////////////////////////////////////////////////////
		
ReconstructionFilter3D::ReconstructionFilter3D(const std::string& name) :
  ReconstructionFilter(name)
{
  //   addConfigOptionsTo(this);
  // _relaxationFactor = 0.;
  // setParameter("RelaxationFactor",&_relaxationFactor);
  // _weightingFactor = 1.;
  // setParameter("WeightingFactor",&_weightingFactor);
}

//////////////////////////////////////////////////////////////////////////////


ReconstructionFilter3D::~ReconstructionFilter3D()
{
	;
}

//////////////////////////////////////////////////////////////////////////////


void ReconstructionFilter3D::setup()
{
  CFAUTOTRACE;
  ReconstructionFilter::setup();
}

//////////////////////////////////////////////////////////////////////////////

std::string ReconstructionFilter3D::reconstructionEquation(const CFuint order) const
{    
  std::string equation = "\\phi_i = \\overline{\\phi}_0";
  for (CFuint k=1; k<=order; k++) {
    // for all terms of degree k
      
    for (CFuint n3=0; n3<=k; n3++) {
      for (CFuint n2=0; n2<=k; n2++) {
        for (CFuint n1=0; n1<=k; n1++) {
          if (n1+n2+n3 == k) {
            std::stringstream term, numerator, denominator;
            if (k==1) {
              term << " + " ;
              numerator << "\\partial \\phi_0";
            }
            else {
              term << " + \\frac{" << trinomialCoefficient(n1,n2,n3) << "}{" << k << "!} ";
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
            if (n3>0) {
              if (n3==1) {
                term << " \\Delta z_{0i}";
                denominator << " \\partial z";
              }
              else {
                term << " \\Delta z_{0i}^" << n3;
                denominator << " \\partial z^" << n3;
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
  }
  return equation;
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter3D::calculateSystemMatrix(RealMatrix& A)
{
  CFuint nbCells = m_stencil->getNbElements();
  Framework::Node& centreCell = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,0);
  CFuint j=0;
  CFreal dx, dy,dz, coefficient;
  CFuint order = getMethodData().getOrder();
	
  for (CFuint k=0; k<=order; k++) {
  	// for all terms of degree k
    for(CFuint n3=0; n3<=k; n3++) {
      for (CFuint n2=0; n2<=k; n2++) {
        for (CFuint n1=0; n1<=k; n1++) {
          if (n1+n2+n3 == k) {
           
            // fill this column
            coefficient = trinomialCoefficient(n1,n2,n3)/CFreal(MathTools::MathFunctions::faculty(k));
            for (CFuint i=1; i<nbCells; i++) { // don't include first cell (i=0)
              Framework::Node& iCell = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,i);
              dx = iCell[XX] - centreCell[XX];
              dy = iCell[YY] - centreCell[YY];
              dz = iCell[ZZ] - centreCell[ZZ];
              A(i-1,j) = coefficient * pow(dx,n1)*pow(dy,n2)*pow(dz,n3);               
            }
            j++; // go to next column
          }
        }
      }
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

CFreal ReconstructionFilter3D::trinomialCoefficient(const CFuint n1, const CFuint n2, const CFuint n3) const
{
  //  (n1 + n2 + n3)!                       (n1 + n2 + n3 + ... + nk)!   
  //  ---------------         multinomial:  ---------------------------  
  //   n1!  n2!  n3!                          n1!  n2!  n3!  ...  nk!    
  // return ( CFreal(ReconstructionFilter::faculty(n1+n2+n3))/CFreal( ReconstructionFilter::faculty(n1) * ReconstructionFilter::faculty(n2) * ReconstructionFilter::faculty(n3) ) );
    return ( CFreal(MathTools::MathFunctions::faculty(n1+n2+n3))/CFreal( MathTools::MathFunctions::faculty(n1) * MathTools::MathFunctions::faculty(n2) * MathTools::MathFunctions::faculty(n3) ) );
}

//////////////////////////////////////////////////////////////////////////////

	  } // namespace ExplicitFilters
		
  } // namespace Numerics

} // namespace COOLFluiD

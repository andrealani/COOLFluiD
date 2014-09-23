
#include "InterpolationFilter.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "ExplicitFilters.hh"
#include "Framework/PhysicalModel.hh"
#include "MathTools/SVDInverter.hh"

#include "nnls.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

		namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<InterpolationFilter, 
                                  FilterData, 
                                  FilterStrategy, 
                                  ExplicitFiltersModule> 
interpolationFilter("InterpolationFilter");

//////////////////////////////////////////////////////////////////////////////

void InterpolationFilter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("NonNegativeLeastSquares","Flag for non negative least squares (default=true)");
	options.addConfigOption< CFreal >("CentreWeight","Factor to add weight to the centre cell (default=0.5)");
	options.addConfigOption< CFreal >("RelaxationFactor","Factor to add weight to the centre cell (default=0.0)");
}

//////////////////////////////////////////////////////////////////////////////

InterpolationFilter::InterpolationFilter(const std::string& name) :
  FilterStrategy(name),
  socket_triangles("filterTriangles")
{
  addConfigOptionsTo(this);
  m_nnls = true;
  setParameter("NonNegativeLeastSquares",&m_nnls);
	m_centreWeight = 0.5;
	setParameter("CentreWeight",&m_centreWeight);
	m_relaxationFactor = -1.;
	setParameter("RelaxationFactor",&m_relaxationFactor);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > InterpolationFilter::needsSockets()
{
  // create socket sink for the stencil
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FilterStrategy::needsSockets();
  result.push_back(&socket_triangles);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

InterpolationFilter::~InterpolationFilter()
{

}

//////////////////////////////////////////////////////////////////////////////


void InterpolationFilter::setup()
{
  CFAUTOTRACE;
  FilterStrategy::setup();
}

//////////////////////////////////////////////////////////////////////////////

void InterpolationFilter::configure(Config::ConfigArgs& args)
{
  FilterStrategy::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void InterpolationFilter::getInfo()
{
  CFLog(INFO, "Explicit Filter Information: \n" <<
              "---------------------------- \n" );
  CFLog(INFO, "Type: " << getClassName() << " \n");
}

//////////////////////////////////////////////////////////////////////////////

void InterpolationFilter::addCentreWeight(const CFuint& centreCellID)
{
  if (m_relaxationFactor < 0) 
    return;
  
  
  CFreal w0_orig = m_weights[0];
  m_weights /= (1.-w0_orig);
  m_weights[0] = 0;
    
    
    
  // The maximum wave number in the stencil
  // RealVector kmax = calculateMaximalWaveNumber(centreCellID);

  CFreal kvalue = calculateMinimumTransferFunction();
  RealVector kmin(kvalue,2);
  kmin *= sqrt(2.)/2.;
  // //kmin *= 0.8;
  // CFLog(INFO, "kmin = " << kmin << " \t");
  // CFLog(INFO, "kmax = " << kmax << " \n");
  //
  // CFreal G_max = integrateTransferFunction(centreCellID,kmin,kmax);

  // The current real value of the transfer function at k=kmax
  // CFreal G_max = real(transferFunction(centreCellID,kmax));
  CFreal G_max = real(transferFunction(kmin)) ;

  CFreal a,b,m,s,fa,fb,fm,f,w0;
  CFuint counter;
  a=-100., b=100.;
  fa = (a+(1.-a)*G_max) - m_relaxationFactor;
  fb = (b+(1.-b)*G_max) - m_relaxationFactor;
  counter = 0;
  do {
    m = a + (b-a)/2.;
    fm = (m+(1.-m)*G_max) - m_relaxationFactor;

    if (fa<0.)    s = -1.;
    else          s = 1.;
    w0 = m + (m-a)*s*fm/sqrt(fm*fm - fa*fb);

    f = (w0+(1.-w0)*G_max) - m_relaxationFactor;

    if (fa*fm < 0.) { b = m ; fb = fm; } else { a = m ; fa = fm; };
    if (fa*f  < 0.) { b = w0; fb = f ; } else { a = w0; fa = f ; };

    counter++;
    if(counter >= 100) {
      w0 = 0.;
      break;
    }
  } while ( std::abs(f) >= 0.0001 );
  m_weights *= (1.-w0);
  m_weights[0] = w0;

}

//////////////////////////////////////////////////////////////////////////////

void InterpolationFilter::calculateWeights()
{
  CFuint centreCellID = getStencilID();
  CFLog(INFO, "\nCalculating weights for cell " << centreCellID << " \n");

  Framework::DataHandle<std::vector<Triangle> > triangles =
    socket_triangles.getDataHandle();    
    
  CFuint nbBasicFilters = 0;//triangles[centreCellID].size();
  
  // Calculate nodal weights for each basic filter
  std::vector<RealVector> nodeWeights(nbBasicFilters);
  MathTools::MatrixInverter* inverter = MathTools::MatrixInverter::create(3,false);
  RealVector X0 = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,0);
  for(CFuint i=0; i<nbBasicFilters; ++i) {
    // calculate matrix A
    RealMatrix A(3,3), invA(3,3);
    A(0,0) = 1.;
    A(1,0) = 1.;
    A(2,0) = 1.;
    //A(0,1) = triangles[centreCellID][i][0]->getCoordinates()[XX] - X0[XX];
    //A(1,1) = triangles[centreCellID][i][1]->getCoordinates()[XX] - X0[XX];
    //A(2,1) = triangles[centreCellID][i][2]->getCoordinates()[XX] - X0[XX];
    //A(0,2) = triangles[centreCellID][i][0]->getCoordinates()[YY] - X0[YY];
    //A(1,2) = triangles[centreCellID][i][1]->getCoordinates()[YY] - X0[YY];
    //A(2,2) = triangles[centreCellID][i][2]->getCoordinates()[YY] - X0[YY];
    
    // calculate b
    RealVector b(3);
    b[0] = 1.;
    b[1] = 0.;
    b[2] = 0.;
    
    // inverse A
    inverter->invert(A,invA);
    
    // node weights
    nodeWeights[i].resize(3);
    
    // nodeWeights[i] = b*invA;    // this doesn't work!!! --> do it by hand
    for(CFuint m=0; m<3; ++m) { //columns
      nodeWeights[i][m] = 0.;
      for(CFuint n=0; n<3; ++n) { //rows
        nodeWeights[i][m] += b[n] * invA(n,m);
      }
    }
  }
  
  // Calculate basic filter weights
  RealMatrix B(4,nbBasicFilters);
  RealVector d(4);
  CFreal equalityConstraintWeight = 100;
  for(CFuint i=0; i<nbBasicFilters; ++i) {
    B(0,i) = calculateTriangleFilterMoment(2, 0, triangles[centreCellID][i], nodeWeights[i]);
    B(1,i) = calculateTriangleFilterMoment(0, 2, triangles[centreCellID][i], nodeWeights[i]);
    B(2,i) = calculateTriangleFilterMoment(1, 1, triangles[centreCellID][i], nodeWeights[i]);  
    
    // Equality constraints LHS
    B(3,i) = (1.) * equalityConstraintWeight;
  }

  d[0] = pow(getFilterWidth(),2)/12.;
  d[1] = pow(getFilterWidth(),2)/12.;
  d[2] = 0.0; //pow(getFilterWidth(),2)/4.;  // or zero ??
  
  // Equality constraints RHS
  d[3] = (1.-m_centreWeight) * equalityConstraintWeight;
  
  RealVector beta(nbBasicFilters);
  if (m_nnls) {
    MathTools::NonNegativeLeastSquares leastSquaresSolver;
    leastSquaresSolver.solve(B,d,beta);    
  }
  else {    
    MathTools::SVDInverter leastSquaresSolver(B);
    try {      
      leastSquaresSolver.solve(d,beta);  
    }
    catch(MathTools::SVDException& e) {
      CFLog(VERBOSE, e.what() << "\n");
      throw;
    }
  }


  CFLog(INFO, "beta = " << beta << " \n");
  CFLog(INFO, "Satisfaction: \n");
  for(CFuint i=0; i<nbBasicFilters; ++i) {
    CFLog(INFO, "eq"<<i<<":\t");
    CFLog(INFO, MathTools::MathFunctions::innerProd(B.getRow<RealVector>(i),beta) << " = "<< d[i]  <<" \n");
  }
  
  CFuint nbCells = m_stencil->getNbElements();
  m_weights.resize(nbCells);
  m_weights[0] = m_centreWeight;
  CFuint count=1;
  for(CFuint i=0; i<nbBasicFilters; ++i) {
    m_weights[count] = beta[i]*nodeWeights[i][0]; count++;
    m_weights[count] = beta[i]*nodeWeights[i][1]; count++;
    m_weights[count] = beta[i]*nodeWeights[i][2]; count++;
  }
  
  CFLog(INFO, "before liftoff \n");
  CFLog(INFO, "M(2,0) = " << calculateFilterMoment(2,0) << " \n");
  CFLog(INFO, "M(0,2) = " << calculateFilterMoment(0,2) << " \n");
  CFLog(INFO, "M(1,1) = " << calculateFilterMoment(1,1) << " \n");
  
  
  addCentreWeight(centreCellID);
  CFLog(INFO, "after liftoff \n");
  CFLog(INFO, "M(2,0) = " << calculateFilterMoment(2,0) << " \n");
  CFLog(INFO, "M(0,2) = " << calculateFilterMoment(0,2) << " \n");
  CFLog(INFO, "M(1,1) = " << calculateFilterMoment(1,1) << " \n");
  

  CFLog(INFO, "FGR_reference = " << calculateMaximalWaveNumber()/calculateFilterWaveNumberRootFinding(XYdiagonal) << " \n");
  CFLog(INFO, "calculateFilterGridRatioTransferFunctionMoment20() = "   
  << calculateFilterGridRatioTransferFunctionMoment20() << " \n");
  CFLog(INFO, "calculateFilterGridRatioTransferFunctionMoment11() = " 
  << calculateFilterGridRatioTransferFunctionMoment11() << " \n");
  CFLog(INFO, "calculateFilterGridRatioM2() = " << calculateFilterGridRatioFilterMoment2() << " \n");
  CFLog(INFO, "calculateFilterWaveNumberRootFinding(XYdiagonal) = " << calculateFilterWaveNumberRootFinding(XYdiagonal) << " \n");
  CFLog(INFO, "calculateMinimumTransferFunction() = " << calculateMinimumTransferFunction() << " \n");
}

//////////////////////////////////////////////////////////////////////////////

CFreal InterpolationFilter::calculateTriangleFilterMoment(const CFuint& q, const CFuint& r, Triangle& triangle, const RealVector& nodeWeights)
{
  CFreal M(0.);
//  const RealVector X0 = triangle.getCentroid()->getCoordinates();
//  RealVector dX(X0.size());
//  for(CFuint i=0; i<3; ++i) {
//    dX = triangle[i]->getCoordinates() - X0;
//    M += nodeWeights[i] * pow(dX[XX],q) * pow(dX[YY],r) ;
//  }
// 
  return M;
}

//////////////////////////////////////////////////////////////////////////////

	  } // end of namespace ExplicitFilters

  } // end of namespace Numerics

} // end of namespace COOLFluiD

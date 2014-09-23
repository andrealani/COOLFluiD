
#include "ReconstructionFilter.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "ExplicitFilters.hh"
#include "Framework/PhysicalModel.hh"
#include "MathTools/SVDInverter.hh"
#include "MathTools/FindMinimum.hh"
#include <boost/math/tools/minima.hpp>
#include "nnls.hh"
#include "Quad2D.hh"

// includes for file manipulations
#include <iomanip>
#include "Environment/DirPaths.hh"
#include "Framework/PathAppender.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Common/StringOps.hh"

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

void ReconstructionFilter::defineConfigOptions(Config::OptionList& options)
{
	options.addConfigOption< CFreal >("WeightingFactor","Weighting factor to influence filter width");
  options.addConfigOption< bool >("Weighting","Defines if weighting is used or not in the reconstruction");
  options.addConfigOption< bool >("NonNegativeLeastSquares","Flag for non negative least squares (default=true)");
	options.addConfigOption< CFreal >("CentreWeight","Factor to add weight to the centre cell (default=0.5)");
}

//////////////////////////////////////////////////////////////////////////////

ReconstructionFilter::ReconstructionFilter(const std::string& name) :
  FilterStrategy(name)
{
  addConfigOptionsTo(this);
	m_relaxationFactor = 0.;
	setParameter("RelaxationFactor",&m_relaxationFactor);
	m_weightingFactor = 1.;
	setParameter("WeightingFactor",&m_weightingFactor);
  m_weighting = true;
  setParameter("Weighting",&m_weighting);
  m_nnls = true;
  setParameter("NonNegativeLeastSquares",&m_nnls);
	m_centreWeight = 0.5;
	setParameter("CentreWeight",&m_centreWeight);
}

//////////////////////////////////////////////////////////////////////////////


ReconstructionFilter::~ReconstructionFilter()
{
	;
}

//////////////////////////////////////////////////////////////////////////////


void ReconstructionFilter::setup()
{
  CFAUTOTRACE;
  FilterStrategy::setup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > ReconstructionFilter::needsSockets()
{
  // create socket sink for the stencil
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FilterStrategy::needsSockets();
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter::configure ( Config::ConfigArgs& args)
{
  FilterStrategy::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter::getInfo()
{
  CFLog(INFO, "Explicit Filter Information: \n" <<
              "---------------------------- \n" );
  CFLog(INFO, "Type: " << getClassName() << " \n");
  CFLog(INFO, "Order: " << getMethodData().getOrder() << " \n");
  CFLog(INFO, "Number of unknowns: " << calculateNbUnknowns(getMethodData().getOrder()) << " \n");
  CFLog(INFO, "Reconstruction Equation (LaTeX source): \n" << reconstructionEquation(getMethodData().getOrder()) << "\n");


}

//////////////////////////////////////////////////////////////////////////////

CFuint ReconstructionFilter::calculateNbUnknowns(const CFuint& order)
{
	CFuint dim(Framework::PhysicalModelStack::getActive()->getDim());
	CFuint prod(1);
	for (CFuint i=1; i<=dim; ++i) {
		prod *= order+i;
	}
	prod /= MathTools::MathFunctions::faculty(dim);

	return prod;
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter::calculateWeightsMatrix(RealMatrix& W)
{

	CFuint nbCells = m_stencil->getNbElements();

	/// @todo make the diagonal matrix into a vector
  cf_assert(W.nbRows() == nbCells);
  for (CFuint i=0; i<W.nbRows(); ++i)
  {
    W(i,i) = 1.;
  }

  if (m_weighting == true) {
    CFreal Delta = getFilterWidth();
    Delta *= m_weightingFactor;
    for (CFuint i=0; i<W.nbRows(); ++i) {
      RealVector DR = RealVector(getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,i) - getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,0));
      W(i,i) = sqrt(6./(MathTools::MathConsts::CFrealPi()*Delta))*exp(- ( MathTools::MathFunctions::innerProd(DR,DR)/(Delta*Delta)) ) ;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter::addRelaxationFactor()
{
  if (m_relaxationFactor < 0) 
    return;
    
  CFreal w0_orig = m_weights[0];
  
  m_weights /= (1.-w0_orig);
  m_weights[0] = 0.;
  
    
    
  // The maximum wave number in the stencil
  // RealVector kmax = calculateMaximalWaveNumberVector();
	CFuint dim(Framework::PhysicalModelStack::getActive()->getDim());
  CFreal kvalue = calculateMinimumTransferFunction();
  
  RealVector kmin(kvalue,dim);
  kmin *= sqrt(static_cast<CFreal>(dim))/static_cast<CFreal>(dim);
  // //kmin *= 0.8;
  // CFLog(INFO, "kmin = " << kmin << " \t");
  // CFLog(INFO, "kmax = " << kmax << " \n");
  //
  // CFreal G_max = integrateTransferFunction(centreCellID,kmin,kmax);

  // The current real value of the transfer function at k=kmax
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
    w0 = m + (m-a)*s*fm/sqrt(std::abs(fm*fm - fa*fb));
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

void ReconstructionFilter::calculateWeights()
{
   if (m_weighting == true){
    // Least Squares Reconstruction where the weighting factor is being optimized
    // to obtain the correct filter grid ratio (several least squares evaluations)
    try {
      optimizedLeastSquaresReconstruction();
      getMethodData().getWeight(getStencilID())->setWeights(m_weights);
      getMethodData().getWeight(getStencilID())->setCompute(false);
      if (!isGoodTransferfunction()) {
        throw ("Bad transferfunction for cell "+Common::StringOps::to_str(getStencilID()));
      }
    }
    catch(std::string str) {
      CFLog(DEBUG_MED, "Exception raised: " << str << " \n");
      m_errorMessagesDuringWeightCalculation += str + " \n";
      getMethodData().getStencil(getStencilID())->setCompute(true);
      getMethodData().getWeight(getStencilID())->setCompute(true);
    }
  } else {
    // Use user-defined weight and see which filter grid ratio is obtained after
    m_weighting = true;
    leastSquaresReconstruction();
    getMethodData().getWeight(getStencilID())->setWeights(m_weights);
    CFLog(DEBUG_MED, "calculateFilterGridRatio(XYdiagonal) = " << calculateFilterGridRatio(Diagonal) << " \n");
    // getFilterGridRatioVSweightingFactor();
    m_weighting = false;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter::leastSquaresReconstruction()
{
  CFuint nbCells = getStencil()->getNbElements();
  m_weights.resize(nbCells);
  CFuint order = getMethodData().getOrder();
	CFuint nbUnknowns = calculateNbUnknowns(order);
	
  
    // calculate matrix A
    RealMatrix A(nbCells,nbUnknowns);
    calculateSystemMatrix(A);

    // calculate matrix W
    RealMatrix W(nbCells,nbCells);
    calculateWeightsMatrix(W);

    // calculate W*A
    RealMatrix WA(nbCells,nbUnknowns);
    WA = W*A;
  
    // calculate pseudo inverse of (W*A)
    MathTools::SVDInverter pseudoInverter(WA);
    pseudoInverter.invert(A); // put inverse in A
  
    // arrive to the weighted pseudo inverse of A
    A = A*W;

    // node weights are first row of this pseudo inverse
    m_weights = A.getRow<RealVector>(0);
    
    // add weight to centre cell to improve filter transfer function
    addRelaxationFactor();
  
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter::getFilterGridRatioVSweightingFactor() 
{
  boost::filesystem::path file;
  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle;
  std::string base = boost::filesystem::basename("FGRvsWeight");
  std::string centreCellIDstr = Common::StringOps::to_str(getStencilID());
  
  CFLog(DEBUG_MED, "Writing file FGRvsWeight-" << getStencilID() << ".dat \n");
  m_errorMessagesDuringWeightCalculation += "Wrote file FGRvsWeight-" + centreCellIDstr + ".dat \n";
  
  /* ---------------------- data file ------------------ */
  std::string data_file = base + ".dat";
  file = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(data_file);
  file = Framework::PathAppender::getInstance().appendParallel( file );
  file = Framework::PathAppender::getInstance().appendCustom( file, centreCellIDstr  );


  fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  std::ofstream& fout = fhandle->open(file);

  CFreal originalWeightingFactor = m_weightingFactor;
  bool originalWeighting = m_weighting;
  RealVector originalWeights(m_weights.size());
  originalWeights = m_weights;
  m_weighting = true;
  
  RealVector kmax = calculateMaximalWaveNumberVector();
  RealVector kmax_x(kmax.size());
  kmax_x[XX] = kmax[XX];
  for (CFreal D=0.5; D<=5; D+=0.1) {
    m_weightingFactor = D;
    leastSquaresReconstruction();
    CFreal transferFunction_kmax = real(transferFunction(kmax));
    CFreal transferFunction_kmax_x = real(transferFunction(kmax_x));
    
    fout << D << "\t";
    // if ((std::abs(transferFunction_kmax_x)>1.0) || (std::abs(transferFunction_kmax)>1.0) || MathChecks::isNaN(transferFunction_kmax)){
    if (!isGoodTransferfunction()){
      fout << 1 << " \n";
    } 
    else {
      fout << calculateFilterGridRatio(Diagonal) << " \n";
    }
  }
  fhandle->close();
  m_weightingFactor = originalWeightingFactor;
  m_weighting = originalWeighting;
  m_weights = originalWeights;
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter::outputAdditional()
{
  getFilterGridRatioVSweightingFactor();
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionFilter::optimizedLeastSquaresReconstruction() 
{  
  CFreal FGR = getMethodData().getFilterGridRatio();
  CFreal a,b,m,s,fa,fb,fm,fp,p,u,l;
  CFuint counter;
  RealVector kmax = calculateMaximalWaveNumberVector();

  CFreal origRelaxationFactor = m_relaxationFactor;
  if (m_relaxationFactor > 0.0)
   m_relaxationFactor = 0.0;
  a=0.5, b=3.0, l=0.5, u=5.0;
  
  m_weightingFactor = a;
  leastSquaresReconstruction();
  fa = calculateFilterGridRatio(Diagonal) - FGR;  
  while(fa==1-FGR || !isGoodTransferfunction()) {
    a+=0.1;
    b+=0.1;
    l=a;
    u+=0.1;
    m_weightingFactor = a;
    leastSquaresReconstruction();
    fa = calculateFilterGridRatio(Diagonal) - FGR;
    if (a<3.0){
      CFLog(DEBUG_MED, "a = " << a << " \n");
      break;
    }
  }
  
  m_weightingFactor = b;
  leastSquaresReconstruction();
  fb = calculateFilterGridRatio(Diagonal) - FGR;
  while(fb==1-FGR) {
    b-=0.1;
    u=b;
    m_weightingFactor = b;
    leastSquaresReconstruction();
    fb = calculateFilterGridRatio(Diagonal) - FGR;
    if (b<a){
      CFLog(DEBUG_MED, "b = " << b << " \n");
      break;
    }
  }
  
  CFreal astart = a, bstart=b;
  counter = 0;
  do {
    m = a + (b-a)/2.;
    m_weightingFactor = m;
    leastSquaresReconstruction();
    fm = calculateFilterGridRatio(Diagonal) - FGR;
    
    if (fa<0.)    s = -1.;
    else          s = 1.;
    p = m + (m-a)*s*fm/sqrt(std::abs(fm*fm - fa*fb));
    p = p<u?p:u;
    p = p>l?p:l;
    m_weightingFactor = p;
    leastSquaresReconstruction();
    fp = calculateFilterGridRatio(Diagonal) - FGR;
    
    if (fa*fm < 0.) { b = m; fb = fm; } else { a = m; fa = fm; };
    if (fa*fp < 0.) { b = p; fb = fp; } else { a = p; fa = fp; };

    counter++;
    if(counter >= 5) {
      m_relaxationFactor = origRelaxationFactor;
      throw ("Really bad transferfunction or could not find optimized Filter Grid Ratio for cell "+Common::StringOps::to_str(getStencilID())+" within domain ["+Common::StringOps::to_str<CFreal>(astart)+","+Common::StringOps::to_str<CFreal>(bstart)+"]");
      break;
    }

  } while ( std::abs(fp/FGR) >= 0.01 );
  m_weightingFactor = p;
  m_relaxationFactor = origRelaxationFactor;
  addRelaxationFactor();
}

//////////////////////////////////////////////////////////////////////////////

	  } // end of namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

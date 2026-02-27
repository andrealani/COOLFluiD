
#include "FilterStrategy.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

// includes for file manipulations
#include <iomanip>
#include "Environment/DirPaths.hh"
#include "Framework/PathAppender.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Common/StringOps.hh"


#include "MathTools/FindMinimum.hh"
#include "MathTools/MathChecks.hh"
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#ifdef CF_HAVE_BOOST_1_85
#define BOOST_TIMER_ENABLE_DEPRECATED
#endif
#include <boost/progress.hpp>

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

void FilterStrategy::defineConfigOptions(Config::OptionList& options)
{
  // Configuration option definitions here	
  options.addConfigOption< CFreal >("TransferFunctionPower","Power of the output of the transfer function");
  options.addConfigOption< CFuint >("MaxNbStencilRecomputations","Maximum number of stencil recomputations");
  options.addConfigOption< CFreal >("AllowedDeviation","The allowed deviation from 0 of the transferfunction at the maximal wavenumber");
  options.addConfigOption< CFreal >("AllowedOvershoot","The allowed overshoot of the transferfunction at half the filter wavenumber");
  options.addConfigOption< CFreal >("RelaxationFactor","Factor to add weight to the centre cell");

}

//////////////////////////////////////////////////////////////////////////////

FilterStrategy::FilterStrategy(const std::string& name) :
  Framework::MethodStrategy<FilterData>(name),
  m_allWeightsCalculated(false)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.
  m_transferFunctionPower = 1.;
  setParameter("TransferFunctionPower",&m_transferFunctionPower);
  m_maxNbStencilRecomputations = 5;
  setParameter("MaxNbStencilRecomputations",&m_maxNbStencilRecomputations);
  m_allowedDeviation = 0.65;
  setParameter("AllowedDeviation",&m_allowedDeviation);
  m_allowedOvershoot = 1.05;
  setParameter("AllowedOvershoot",&m_allowedOvershoot);
  m_relaxationFactor = 0.;
	setParameter("RelaxationFactor",&m_relaxationFactor);
  
}

//////////////////////////////////////////////////////////////////////////////


FilterStrategy::~FilterStrategy()
{
	;
}

//////////////////////////////////////////////////////////////////////////////

void FilterStrategy::setup()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "setting up FilterStrategy \n");
 

}

//////////////////////////////////////////////////////////////////////////////

void FilterStrategy::configure ( Config::ConfigArgs& args )
{
  Framework::MethodStrategy<FilterData>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > FilterStrategy::needsSockets()
{
  // create socket sink for the stencil
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > FilterStrategy::providesSockets()
{
  // create socket sink for the filter weights
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::getFilterWidth()
{
  return getStencil()->getCellWidth()*getMethodData().getFilterGridRatio();
}

//////////////////////////////////////////////////////////////////////////////

CFuint FilterStrategy::getNbWeightsToCompute()
{
  {
    CFuint counter = 0;
    CFuint nbStencils = getMethodData().getWeights()->size();
    for (CFuint iStencil=0; iStencil<nbStencils; iStencil++) {
      if (getMethodData().getWeight(iStencil)->mustCompute()) {
        counter++;
      }
    }
    return counter;
  }
}

//////////////////////////////////////////////////////////////////////////////

      
void FilterStrategy::calculateAllWeights()
{
  // only calculate in case they are not already calculated
  if (!m_allWeightsCalculated) {
        
    
    // First pass to calculate all weights    
    const CFuint nbStencils = getMethodData().getStencils()->size();
    boost::progress_display progress(nbStencils);
    for(CFuint iStencil=0; iStencil<nbStencils; iStencil++, ++progress) {
      setStencil(getMethodData().getStencil(iStencil));
      calculateWeights();
    }  
    
    CFLog(INFO, "\n");
    
    
    
    
    //////////////////////////
    
    
    // Datahandles and declarations
    CFuint iRecompute = 0;
    CFreal enlargementFactor = 1.15;
    CFuint nbFiltersToRecompute = getNbWeightsToCompute();
    while (nbFiltersToRecompute && iRecompute < m_maxNbStencilRecomputations) {
      CFLog(INFO, "Recomputing " << nbFiltersToRecompute << " stencils and weights...");
      
      
      // Enlarge radius for all stencils
      for(CFuint iStencil=0; iStencil<nbStencils; iStencil++){
        if (getMethodData().getStencil(iStencil)->mustCompute()) {
          getMethodData().getStencil(iStencil)->enlargeRadiusWithFactor(enlargementFactor);
        }
      }
      
      // Compute stencils with enlarged radius
      getMethodData().getStencilComputer()->compute();
      
      
      
      CFLog(INFO, "Recomputing " << nbFiltersToRecompute << " weights...");
      // Compute weights with new stencils
      m_errorMessagesDuringWeightCalculation = "";
      boost::progress_display* progressRecompute(CFNULL);
      if (nbFiltersToRecompute > 200) {
        progressRecompute = new boost::progress_display(nbFiltersToRecompute);
      }
      else {
        CFLog(INFO, " \n");
      }
      // Enlarge radius for all stencils
      for(CFuint iStencil=0; iStencil<nbStencils; iStencil++){
        if (getMethodData().getWeight(iStencil)->mustCompute()) {
          if (progressRecompute != CFNULL) ++(*progressRecompute);
          setStencil(getMethodData().getStencil(iStencil));
          calculateWeights();              
        }
      }
      delete progressRecompute;
      
      nbFiltersToRecompute = getNbWeightsToCompute();
      iRecompute++;
    }
    
    
    // Output additional things for bad filters
    if (iRecompute == m_maxNbStencilRecomputations && nbFiltersToRecompute) {
      CFLog(INFO, nbFiltersToRecompute << " stencils are not good after " << m_maxNbStencilRecomputations << " recomputations. \n");

      for(CFuint iStencil=0; iStencil<nbStencils; iStencil++){
        if (getMethodData().getStencil(iStencil)->mustCompute() || getMethodData().getWeight(iStencil)->mustCompute()) {
          if (getMethodData().outputDebug()) {
            // Output some extra data and plots for this filter
            setStencil(getMethodData().getStencil(iStencil));
            outputAdditional();
          }
          // getMethodData().getStencil(iStencil)->clear();
          // getMethodData().getWeight(iStencil)->clear();
          getMethodData().getFilterFlag(iStencil)=false;
        }
      }
    }
    
    
    if (m_errorMessagesDuringWeightCalculation.size() != 0) {
      
      CFLog(INFO, "Errors during the weight calculation: \n" 
                  << m_errorMessagesDuringWeightCalculation);            
    }
    
    m_allWeightsCalculated = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void FilterStrategy::calculateCentreCellWeight()
{
  // only calculate in case they are not already calculated
  if (!m_allWeightsCalculated) {
        
    
    // First pass to calculate weight
    m_errorMessagesDuringWeightCalculation.clear();
    calculateWeights();
    
    CFuint iRecompute = 0;
    CFreal enlargementFactor = 1.15;
    while (getMethodData().getWeight(getStencilID())->mustCompute() && iRecompute < m_maxNbStencilRecomputations) {
      CFLog(INFO, "Recomputing stencil and weight...");
      
      // Enlarge radius for stencil
      if (getMethodData().getStencil(getStencilID())->mustCompute()) {
        getMethodData().getStencil(getStencilID())->enlargeRadiusWithFactor(enlargementFactor);
      }
      
      // Compute stencil with enlarged radius
      getMethodData().getStencilComputer()->computeStencil(getStencilID());
      
      // Compute weights with new stencil
      m_errorMessagesDuringWeightCalculation.clear();
      calculateWeights();
      iRecompute++;
    }
    
    
    // Output additional things for bad filters
    if (iRecompute == m_maxNbStencilRecomputations && getMethodData().getWeight(getStencilID())->mustCompute()) {
      CFLog(INFO, "transfer function is not good after " << m_maxNbStencilRecomputations << " recomputations. \n");
      if (getMethodData().outputDebug()) {
        // Output some extra data and plots for this filter
        outputAdditional();
      }
      // getMethodData().getStencil(getStencilID())->clear();
      // getMethodData().getWeight(getStencilID())->clear();
      getMethodData().getFilterFlag(getStencilID())=false;
    }
    
    
    if (m_errorMessagesDuringWeightCalculation.size() != 0) {
      
      CFLog(INFO, "Errors during the weight calculation: \n" 
                  << m_errorMessagesDuringWeightCalculation);            
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

RealVector FilterStrategy::getRealVectorInDirection(direction dir)
{
  CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
  RealVector v(dim);
  
  switch (dim) {
    case DIM_2D:
      switch (dir) {
        case Xdirection :
          v[XX] = 1;
          v[YY] = 0;
          return v;
        case Ydirection :
          v[XX] = 0;
          v[YY] = 1;
          return v;
        case Diagonal :
        case XYdiagonal :
          v[XX] = 1./sqrt(2.);
          v[YY] = 1./sqrt(2.);
          return v;
        default:
          v[XX] = 1;
          v[YY] = 0;
          return v;
      }
      break;
      
    case DIM_3D:
      switch (dir) {
        case Xdirection :
          v[XX] = 1;
          v[YY] = 0;
          v[ZZ] = 0;
          return v;
        case Ydirection :
          v[XX] = 0;
          v[YY] = 1;
          v[ZZ] = 0;
          return v;
        case Zdirection :
          v[XX] = 0;
          v[YY] = 0;
          v[ZZ] = 1;
          return v;
        case XYdiagonal :
          v[XX] = 1./sqrt(2.);
          v[YY] = 1./sqrt(2.);
          v[ZZ] = 0;
          return v;
        case YZdiagonal :
          v[XX] = 0;
          v[YY] = 1./sqrt(2.);
          v[ZZ] = 1./sqrt(2.);
          return v;
        case XZdiagonal :
          v[XX] = 1./sqrt(2.);
          v[YY] = 0;
          v[ZZ] = 1./sqrt(2.);
          return v;
        case Diagonal :
        case XYZdiagonal :
          v[XX] = 1./sqrt(3.);
          v[YY] = 1./sqrt(3.);
          v[ZZ] = 1./sqrt(3.);
          return v;
      }
      break;
  }
  
  // Should not get here
  return v;
  
}

//////////////////////////////////////////////////////////////////////////////

RealVector FilterStrategy::calculateMaximalWaveNumberVector()
{
  return RealVector(calculateMaximalWaveNumber(),Framework::PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateMaximalWaveNumber()
{
  // kmax = (pi / mesh_size)
  return MathTools::MathConsts::CFrealPi()/m_stencil->getCellWidth();
}

//////////////////////////////////////////////////////////////////////////////

CFcomplex FilterStrategy::transferFunction(RealVector& k) 
{

  CFcomplex G(0,0);
  Framework::Node& X0 = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,0);

  CFuint nbCells = m_stencil->getNbElements();
  for(CFuint iCell=0; iCell<nbCells; ++iCell) {
      Framework::Node& X = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,iCell);;
      G += m_weights[iCell] * exp(-MathTools::MathConsts::CFcomplexI()*MathTools::MathFunctions::innerProd(k,RealVector(X-X0)));
  }
  return G;
}

//////////////////////////////////////////////////////////////////////////////

CFcomplex FilterStrategy::transferFunction(const RealVector& k, const std::vector<Framework::State*>& stencil, const RealVector& weights)
{
  Framework::Node& X0 = m_centreCell->getCoordinates();
  RealVector dX(X0.size());
  CFcomplex G(0,0);
  CFuint nbCells = stencil.size();
  for(CFuint iCell=0; iCell<nbCells; ++iCell) {
    dX = stencil[iCell]->getCoordinates() - X0;
    G += weights[iCell] * exp(-MathTools::MathConsts::CFcomplexI()*MathTools::MathFunctions::innerProd(k,dX));        
  }
  return G;
}

//////////////////////////////////////////////////////////////////////////////

bool FilterStrategy::isGoodTransferfunction()
{
  
  // Initializations
  bool isGood = true;
  std::vector<RealVector> vec;
  std::vector<CFreal> val;
  RealVector kmax = calculateMaximalWaveNumberVector();
  RealVector overshootLocation = kmax/(2.0*getMethodData().getFilterGridRatio());
  CFuint dim=Framework::PhysicalModelStack::getActive()->getDim();
  CFreal tf;
  RealVector tmp;
  
  vec.reserve(10);
  val.reserve(10);
  
  // Definition of vectors and values under which the absolute
  // value of the transferfunction must lie.
  
  // XY(Z)-diagonal end
  tmp.resize(dim);  tmp = kmax;
  vec.push_back(tmp);
  val.push_back(m_allowedDeviation);
    
  // X-direction end
  tmp.resize(dim,0);  tmp[XX] = kmax[XX];
  vec.push_back(tmp);
  val.push_back(m_allowedDeviation);
  
  // Y-direction end
  tmp.resize(dim,0);  tmp[YY] = kmax[YY];
  vec.push_back(tmp);
  val.push_back(m_allowedDeviation);
  
  // X(-Y)-diagonal end
  tmp = kmax;  tmp[YY] = -tmp[YY];
  vec.push_back(tmp);
  val.push_back(m_allowedDeviation);
  
  // XY(Z)-diagonal overshoot
  tmp = overshootLocation;
  vec.push_back(tmp);
  val.push_back(m_allowedOvershoot);
  
  // X-direction overshoot
  tmp.resize(dim,0);  tmp[XX] = overshootLocation[XX];
  vec.push_back(tmp);
  val.push_back(m_allowedOvershoot);
  
  // Y-direction overshoot
  tmp.resize(dim,0);  tmp[YY] = overshootLocation[YY];
  vec.push_back(tmp);
  val.push_back(m_allowedOvershoot);
  
    
  if (dim > 2) {

    // XY-diagonal end
    tmp.resize(dim,0);  tmp[XX] = kmax[XX];   tmp[YY] = kmax[YY];
    vec.push_back(tmp);
    val.push_back(m_allowedDeviation);
    
    // XZ-diagonal end
    tmp.resize(dim,0);  tmp[XX] = kmax[XX];   tmp[ZZ] = kmax[ZZ];
    vec.push_back(tmp);
    val.push_back(m_allowedDeviation);
    
    // YZ-diagonal end
    tmp.resize(dim,0);  tmp[YY] = kmax[YY];   tmp[ZZ] = kmax[ZZ];
    vec.push_back(tmp);
    val.push_back(m_allowedDeviation);
    
  }
    
  // Check the transferfunctions
  for(CFuint i=0; i<val.size(); ++i) {
    tf = real(transferFunction(vec[i]));
    isGood = isGood && !MathChecks::isNaN(tf) && (std::abs(tf-m_relaxationFactor)-val[i]*(1.-m_relaxationFactor)<0.0);
    if (!isGood) {
      return false;
    } 
  }
  return isGood;
}

//////////////////////////////////////////////////////////////////////////////

void FilterStrategy::outputTransferFunction()
{

  CFuint N=50;
    
    
  CFuint dim(Framework::PhysicalModelStack::getActive()->getDim());
	
  if (dim==DIM_2D) {
    RealVector** K = new RealVector* [N];
    CFcomplex** G = new CFcomplex* [N];
    for(CFuint i=0; i<N; ++i) {
      K[i] = new RealVector [N];
      G[i] = new CFcomplex [N];
      for(CFuint j=0; j<N; ++j) {
        K[i][j].resize(2);
      }
    }

  	// Loop over inspectedCellIDs
  	Common::SafePtr<std::vector<CFuint> > inspectedCellIDs = getMethodData().getInspectedCellIDs();
    CFuint nbInspectedCellIDs = inspectedCellIDs->size();
    for (CFuint iInspected=0; iInspected<nbInspectedCellIDs; ++iInspected) {    
      CFuint centreCellID=(*inspectedCellIDs)[iInspected];
      CFLog(INFO, "\n" << "inspected ID = " << centreCellID << "\n");
      setStencil(getMethodData().getStencil(centreCellID));
      // calculate these filter weights again, just to be sure.
      calculateCentreCellWeight();

      RealVector kmax = calculateMaximalWaveNumberVector();

      /* -------------- Output gnuplot ----------------- */
      RealVector k_11(N);
      RealVector k_10(N);
      RealVector G_11(N);
      RealVector G_10(N);
      RealVector G_im_11(N);
      RealVector G_im_10(N);

      RealVector K_11(2);
      RealVector K_10(2);

      for (CFuint i=0; i<N; i++) {
          K_11 = CFreal(i)*kmax/(N-1.0);
          k_11[i] = K_11.norm2();
          K_10[XX] = i*kmax[XX]/(N-1.0);     K_10[YY] = 0.0;
          k_10[i] = K_10.norm2();

          G_11[i]=real(transferFunction(K_11));
          G_10[i]=real(transferFunction(K_10));
          G_im_11[i]=real(transferFunction(K_11));
          G_im_10[i]=real(transferFunction(K_10));

      }

      outputTransferFunctionGnuplot(k_11,G_11,k_10,G_10,centreCellID);

      /* ------------------- Output Tecplot ------------------- */

      for(CFuint i=0; i<N; ++i) {
        for(CFuint j=0; j<N; ++j) {
          K[i][j][XX] = i*kmax[XX]/(N-1.);
          K[i][j][YY] = j*kmax[YY]/(N-1.);
          G[i][j] = transferFunction(K[i][j]);
        }
      }
      outputTransferFunctionTecplot(K,G,N,centreCellID);


      std::string filename = getMethodData().getTransferFunctionFileName();
      boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename);
      std::string centreCellIDstr = Common::StringOps::to_str(centreCellID);
      file = Framework::PathAppender::getInstance().appendCustom( file , centreCellIDstr);
      
      /* ------------------ Output additional ------------------ */
      outputAdditional();
    }

    for(CFuint i=0; i<N; ++i) {
      delete[] K[i];
      delete[] G[i];
    }
  }
  else if (dim==DIM_3D) {
    CFLog(INFO, "Output for 3D not yet implemented \n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void FilterStrategy::outputTransferFunctionGnuplot(RealVector& k_11, RealVector& G_11, RealVector& k_10, RealVector& G_10, CFuint& centreCellID)
{ 

  // preparation of the output
  
  boost::filesystem::path file;
  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle;
#ifdef CF_HAVE_BOOST_1_85
  std::string base = boost::filesystem::path(getMethodData().getTransferFunctionFileName()).stem().string();
#else
  std::string base = boost::filesystem::basename(getMethodData().getTransferFunctionFileName());
#endif
  std::string centreCellIDstr = Common::StringOps::to_str(centreCellID);
  
  
  /* ---------------------- gnuplot data file ------------------ */
  std::string data_file = base + ".dat";
  file = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(data_file);
  file = Framework::PathAppender::getInstance().appendParallel( file );
  file = Framework::PathAppender::getInstance().appendCustom( file, centreCellIDstr  );


  fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  std::ofstream& fout = fhandle->open(file);


  // writing some comments 
  fout 
  << "#\n"
  << "# Transferfunction for cell " << centreCellID << "\n"
  << "#\n"
  << "# FGR_target = " << getMethodData().getFilterGridRatio() << "\n"
  << "# FGR_xy = " << calculateFilterGridRatio(XYdiagonal) << "\n"
  << "# FGR_x = " << calculateFilterGridRatio(Xdirection) << "\n"
  << "#\n"
  << "#   k_xy                   G_target             G_xy               k_x                   G_x\n"
  << "#   -----------------------------------------------------------------------------------------------------\n";
  
  // writing data  
  CFreal kmax = calculateMaximalWaveNumber();

  for (CFuint i=0; i<k_11.size(); i++) {
      fout
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << k_11[i]
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << pow(targetTransferFunction(k_11[i]*MathTools::MathConsts::CFrealPi()/kmax),m_transferFunctionPower) 
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << pow(G_11[i],m_transferFunctionPower) 
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << k_10[i]  
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << pow(G_10[i],m_transferFunctionPower) 
      << "\n";
  } 

  //closing the file
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void FilterStrategy::outputTransferFunctionTecplot(RealVector** K, CFcomplex** G, const CFuint& N, const CFuint& centreCellID) 
{
  // preparation of the output
    
  boost::filesystem::path file;
  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle;
#ifdef CF_HAVE_BOOST_1_85
  std::string base = boost::filesystem::path(getMethodData().getTransferFunctionFileName()).stem().string();
#else
  std::string base = boost::filesystem::basename(getMethodData().getTransferFunctionFileName());
#endif 
  std::string centreCellIDstr = Common::StringOps::to_str(centreCellID);
  
  std::string data_file = base + ".plt";
  file = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(data_file);
  file = Framework::PathAppender::getInstance().appendParallel( file );
  file = Framework::PathAppender::getInstance().appendCustom( file, centreCellIDstr  );

  fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  std::ofstream& fout = fhandle->open(file);


  // writing data
        
  fout << "TITLE = \"" << "Transfer function, "
  << "\"" << "\n"
  << "VARIABLES = "
  << "\"kx\" \\ \n"
  << "\"ky\" \\ \n"
  << "\"G real\" \\ \n"
  << "\"G imag\" \\ \n"
  << "\"G abs\" \\ \n";
  
  fout << "ZONE T =  \"Block Number = 0" 
  << "\" \\ \n"
  << "I = " << N << " \\ \n"
  << "J = " << N << " \\ \n"
  << "DATAPACKING = POINT \n";
  
  for (CFuint j=0; j<N; j++) {
    for (CFuint i=0; i<N ; i++) {
      CFcomplex Gpow = pow(G[i][j],m_transferFunctionPower);
      fout
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << K[i][j]
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << real(Gpow) 
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << imag(Gpow)
      << " " << std::setw(20) << std::fixed << std::setprecision(12) << abs(Gpow)
      << "\n";
    } 
  } 

  //closing the file
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::targetTransferFunction(const RealVector& k)
{
  return targetTransferFunction(k.norm2());
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::targetTransferFunction(const CFreal& k)
{
  CFreal G_cutoff = 0.5;
  CFreal m = getMethodData().getTargetFilterOrder();
  CFreal a = -(2.*m)*log(G_cutoff)*pow(MathTools::MathConsts::CFrealPi()/getMethodData().getFilterGridRatio(),-2.*m);
  return exp(-a/(2.*m)*pow(k,2.*m));
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::targetTransferFunction(const CFreal& kx, const CFreal& ky)
{
  CFreal k = sqrt(kx*kx + ky*ky);
  return targetTransferFunction(k);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterGridRatioTransferFunctionMoment2()
{
  CFreal DeltaX = m_stencil->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;
  return pow( 3.*calculateTransferFunctionMoment2() * pow(D[XX]*D[YY],3) 
          / (pow(MathConsts::CFrealPi(),4) * (pow(D[XX],2) + pow(D[YY],2))) ,-1./4.);
          
  // If Moment is calculated with polar integration
  // return pow( 2.*calculateTransferFunctionMoment2(centreCellID) * pow(D[XX],4) 
  //         / pow(MathConsts::CFrealPi(),5) ,-1./4.);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterGridRatioTransferFunctionMoment20()
{
  CFreal DeltaX = m_stencil->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;
  return pow( 3.*calculateTransferFunctionMoment20() * pow(D[XX],3)*D[YY] / pow(MathConsts::CFrealPi(),4)  ,-1./4.);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterGridRatioTransferFunctionMoment11()
{
  CFreal DeltaX = m_stencil->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;
  return pow( 4.*calculateTransferFunctionMoment11() * pow(D[XX],2)*pow(D[YY],2) / pow(MathConsts::CFrealPi(),4)  ,-1./4.);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterGridRatioFilterMoment2()
{
  CFreal Delta = m_stencil->getCellWidth();
  CFreal M2 = calculateFilterMoment(2,0);
  return sqrt(12.*M2)/Delta;
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterMoment(const CFuint& q, const CFuint& r)
{
  Framework::Node& X0 = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,getStencilID());
  RealVector dX(X0.size());
  CFreal M(0.);
  CFuint nbCells = m_stencil->getNbElements();
  for(CFuint iCell=0; iCell<nbCells; ++iCell) {
      dX = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,iCell) - X0;
      M += pow(dX[XX],q) * pow(dX[YY],r) * m_weights[iCell];
  }
  return M;
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterMoment(const CFuint& q, const CFuint& r, std::vector<Framework::State*>& stencil, RealVector& weights)
{
  Framework::Node& X0 = m_centreCell->getCoordinates();
  RealVector dX(X0.size());
  CFreal M(0.);
  CFuint nbCells = stencil.size();
  for(CFuint iCell=0; iCell<nbCells; ++iCell) {
      dX = stencil[iCell]->getCoordinates() - X0;
      M += pow(dX[XX],q) * pow(dX[YY],r) * weights[iCell];
  }
  return M;
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateTransferFunctionMoment20()
{
  Framework::Node& X0 = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,getStencilID());
  CFreal DeltaX = m_stencil->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;

  CFreal M(0.);
  CFuint nbCells = m_stencil->getNbElements();
  for(CFuint iCell=1; iCell<nbCells; ++iCell) {
      Framework::Node& X = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,iCell);
      RealVector dX = X-X0;

      M += m_weights[iCell] * 2.*sin(MathConsts::CFrealPi()*dX[YY]/(2.*D[YY])) / (pow(D[XX],2)*pow(dX[XX],3)*dX[YY]) *
        (
          -4.*cos(MathConsts::CFrealPi()/2.*(dX[XX]/D[XX]+dX[YY]/D[YY]))*sin(MathConsts::CFrealPi()/2.*dX[XX]/D[XX])*pow(D[XX],2.)
          + 2.*MathConsts::CFrealPi()*cos(MathConsts::CFrealPi()*(dX[XX]/D[XX]+dX[YY]/(2.*D[YY])))*D[XX]*dX[XX]
          + pow(MathConsts::CFrealPi(),2)*sin(MathConsts::CFrealPi()*(dX[XX]/D[XX]+dX[YY]/(2.*D[YY])))*pow(dX[XX],2)
        );
      // CFLog(INFO, "Mloop = " << M << " \t");
      // CFLog(INFO, "dX = " << dX << " \n");
  }

  M += m_weights[0]*pow(MathConsts::CFrealPi(),4)/(3.*pow(D[XX],3)*D[YY]);
  // CFLog(INFO, "Mfinal = " << M << " \n");
  return (M);
}

//////////////////////////////////////////////////////////////////////////////

/**
 * For stencils without centre cell!!!!!!
 */
CFreal FilterStrategy::calculateTransferFunctionMoment20(const std::vector<Framework::State*>& stencil, const RealVector& weights)
{
  Framework::Node& X0 = m_centreCell->getCoordinates();
  CFreal DeltaX = getMethodData().getStencil(m_centreCell->getLocalID())->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;

  CFreal M(0.);
  CFuint nbCells = stencil.size();
  for(CFuint iCell=0; iCell<nbCells; ++iCell) {
      Framework::Node& X = stencil[iCell]->getCoordinates();
      RealVector dX = X-X0;

      M += weights[iCell] * 2.*sin(MathConsts::CFrealPi()*dX[YY]/(2.*D[YY])) / (pow(D[XX],2)*pow(dX[XX],3)*dX[YY]) *
        (
          -4.*cos(MathConsts::CFrealPi()/2.*(dX[XX]/D[XX]+dX[YY]/D[YY]))*sin(MathConsts::CFrealPi()/2.*dX[XX]/D[XX])*pow(D[XX],2.)
          + 2.*MathConsts::CFrealPi()*cos(MathConsts::CFrealPi()*(dX[XX]/D[XX]+dX[YY]/(2.*D[YY])))*D[XX]*dX[XX]
          + pow(MathConsts::CFrealPi(),2)*sin(MathConsts::CFrealPi()*(dX[XX]/D[XX]+dX[YY]/(2.*D[YY])))*pow(dX[XX],2)
        );
      // CFLog(INFO, "Mloop = " << M << " \t");
      // CFLog(INFO, "dX = " << dX << " \n");
  }

  // CFLog(INFO, "Mfinal = " << M << " \n");
  return (M);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateTransferFunctionMoment11()
{
  Framework::Node& X0 = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,getStencilID());
  CFreal DeltaX = m_stencil->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;

  // RealVector K(2);
  // K[XX] = MathConsts::CFrealPi()/D[XX];
  // K[YY] = MathConsts::CFrealPi()/D[YY];

  CFreal kvalue = calculateMinimumTransferFunction();
  RealVector K(kvalue,2);
  // K *= sqrt(2.)/2.;

  CFcomplex M(0.,0.);
  CFuint nbCells = m_stencil->getNbElements();
  for(CFuint iCell=1; iCell<nbCells; ++iCell) {
      Framework::Node& X = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,iCell);
      RealVector dX = X-X0;

      M += m_weights[iCell] / (dX[XX]*dX[XX] * dX[YY]*dX[YY])
            * exp(-MathConsts::CFcomplexI()*MathFunctions::innerProd(K,dX))
            * (-1. + exp(MathConsts::CFcomplexI()*K[XX]*dX[XX]) - MathConsts::CFcomplexI()*K[XX]*dX[XX])
            * (-1. + exp(MathConsts::CFcomplexI()*K[YY]*dX[YY]) - MathConsts::CFcomplexI()*K[YY]*dX[YY]);
  }

  M += m_weights[0]/4. * (K[XX]*K[XX] * K[YY]*K[YY]);
  // CFLog(INFO, "Mfinal = " << M << " \n");
  return real(M);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateTransferFunctionMoment11(const std::vector<Framework::State*>& stencil, const RealVector& weights)
{
  Framework::Node& X0 = m_centreCell->getCoordinates();
  CFreal DeltaX = getMethodData().getStencil(m_centreCell->getLocalID())->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;

  // RealVector K(2);
  // K[XX] = MathConsts::CFrealPi()/D[XX];
  // K[YY] = MathConsts::CFrealPi()/D[YY];

  CFreal kvalue = calculateMinimumTransferFunction();
  RealVector K(kvalue,2);
  // K *= sqrt(2.)/2.;

  CFcomplex M(0.,0.);
  CFuint nbCells = stencil.size();
  for(CFuint iCell=1; iCell<nbCells; ++iCell) {
      Framework::Node& X = stencil[iCell]->getCoordinates();
      RealVector dX = X-X0;

      M += weights[iCell] / (dX[XX]*dX[XX] * dX[YY]*dX[YY])
            * exp(-MathConsts::CFcomplexI()*MathFunctions::innerProd(K,dX))
            * (-1. + exp(MathConsts::CFcomplexI()*K[XX]*dX[XX]) - MathConsts::CFcomplexI()*K[XX]*dX[XX])
            * (-1. + exp(MathConsts::CFcomplexI()*K[YY]*dX[YY]) - MathConsts::CFcomplexI()*K[YY]*dX[YY]);
  }

  M += weights[0]/4. * (K[XX]*K[XX] * K[YY]*K[YY]);
  // CFLog(INFO, "Mfinal = " << M << " \n");
  return real(M);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateTransferFunctionMoment2()
{
  Framework::Node& X0 = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,getStencilID());
  CFreal DeltaX = m_stencil->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;

  // RealVector K(2);
  // K[XX] = MathConsts::CFrealPi()/D[XX];
  // K[YY] = MathConsts::CFrealPi()/D[YY];

  CFreal kvalue = calculateMinimumTransferFunction();
  RealVector K(kvalue,2);
  // K *= sqrt(2.)/2.;

  CFcomplex M(0.,0.);
  CFcomplex I = MathConsts::CFcomplexI();
  CFuint nbCells = m_stencil->getNbElements();
  for(CFuint iCell=1; iCell<nbCells; ++iCell) {
      Framework::Node& X = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,iCell);
      RealVector dX = X-X0;

      M += (1./pow(dX[XX]*dX[YY],3)) *
        (
      	  (
      		  m_weights[iCell]*
      		  (
      		    2.*(-1. + exp(I*K[XX]*dX[XX]))*(-1. + exp(I*K[YY]*dX[YY]))*pow(dX[YY],2) 
      		    - 2.*I*(-1. + exp(I*K[YY]*dX[YY]))*K[XX]*dX[XX]*pow(dX[YY],2) 
      		    + pow(dX[XX],2)*
      		    (
      		      - 2.*I*(-1. + exp(I*K[XX]*dX[XX]))*K[YY]*dX[YY] 
      		      + (-1. + exp(I*K[XX]*dX[XX]))*pow(K[YY]*dX[YY],2) 
      		      + (-1. + exp(I*K[YY]*dX[YY]))*(-2. + 2.*exp(I*K[XX]*dX[XX]) + pow(K[XX]*dX[YY],2))
      		    )
      		  )
      		)
      		/exp(I*(K[XX]*dX[XX] + K[YY]*dX[YY]))
          );
  }

  M += (1./3.)*m_weights[0]*K[XX]*K[YY]*(pow(K[XX],2) + pow(K[YY],2));
  
  // CFLog(INFO, "Mfinal = " << M << " \n");
  return real(M);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterGridRatio(const direction dir) 
{
  return calculateMaximalWaveNumber() / calculateFilterWaveNumberRootFinding(dir);
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterWaveNumberRootFinding(const direction dir)
{
  using namespace boost::math::policies;
  
  
  typedef policy<
     domain_error<throw_on_error>,
     pole_error<errno_on_error>,
     overflow_error<errno_on_error>,
     evaluation_error<errno_on_error> 
  > myPolicy; 
  myPolicy pol;
  Gfunc_cutoff G(this,dir);
  Tol tol;
  CFreal location, lower, upper, kmax;
  boost::uintmax_t max_iter=50;
  kmax = calculateMaximalWaveNumber();
  lower = kmax/getMethodData().getFilterGridRatio()/3.;
  // upper = sqrt(Framework::PhysicalModelStack::getActive()->getDim())*kmax;
  upper = kmax/1.3;
  // CFLog(INFO, "lower = " << lower << " \n");
  // CFLog(INFO, "upper = " << upper << " \n");
  try {
    location = boost::math::tools::toms748_solve(G, lower, upper, tol, max_iter, pol).first;
  }
  catch(std::domain_error& err) {
    // CFLog(INFO, "Exception caught in function ExplicitFilters::FilterStrategy::calculateFilterWaveNumberRootFinding: \n");
    // CFLog(INFO, err.what() << "\n");
    location = kmax;
  }
  
  return location;
}



  // CFuint dim(Framework::PhysicalModelStack::getActive()->getDim());
  // 
  // CFreal G_value = 0.5;
  // 
  // RealVector k_dir = getRealVectorInDirection(dir);
  // 
  // CFreal a,b,m,s,fa,fb,fm,fp,p;
  // RealVector ka(dim), kb(dim), km(dim), kp(dim);
  // CFuint counter;
  // a=0., b=0.5*calculateMaximalWaveNumber();
  // ka = a*k_dir;
  // kb = b*k_dir;
  // fa = real(transferFunction(ka)) - G_value ;
  // fb = real(transferFunction(kb)) - G_value ;
  // counter = 0;
  // do {
  //   m = a + (b-a)/2.;
  //   km = m*k_dir;
  //   fm = real(transferFunction(km)) - G_value ;
  // 
  //   if (fa<0.)    s = -1.;
  //   else          s = 1.;
  //   p = m + (m-a)*s*fm/(sqrt(fabs(fm*fm - fa*fb))+MathTools::MathConsts::CFrealEps());
  //   kp = p*k_dir;
  //   fp = real(transferFunction(kp)) - G_value ;
  // 
  //   if (fa*fm < 0.) { b = m; fb = fm; } else { a = m; fa = fm; };
  //   if (fa*fp < 0.) { b = p; fb = fp; } else { a = p; fa = fp; };
  // 
  //   counter++;
  //   if(counter >= 100) {
  //     CFLog(INFO, "maximum iteration reached to find filter cutoff wave number \n");
  //     break;
  //   }
  // } while( std::abs(fp) >= 0.001 );
  // 
  // CFreal k_value = p;
  // return k_value;



//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateFilterWaveNumberRootFinding()
{

  CFuint dim(Framework::PhysicalModelStack::getActive()->getDim());

  CFreal G_value = 0.5;

  RealVector k_dir = getRealVectorInDirection(Diagonal);

  CFreal a,b,m,s,fa,fb,fm,fp,p;
  RealVector ka(dim), kb(dim), km(dim), kp(dim);
  CFuint counter;
  a=0., b=0.5*calculateMaximalWaveNumber();
  ka = a*k_dir;
  kb = b*k_dir;
  fa = real(transferFunction(ka)) - G_value ;
  fb = real(transferFunction(kb)) - G_value ;
  counter = 0;
  do {
    m = a + (b-a)/2.;
    km = m*k_dir;
    fm = real(transferFunction(km)) - G_value ;

    if (fa<0.)    s = -1.;
    else          s = 1.;
    p = m + (m-a)*s*fm/(sqrt(fabs(fm*fm - fa*fb))+MathTools::MathConsts::CFrealEps());
    kp = p*k_dir;
    fp = real(transferFunction(kp)) - G_value ;

    if (fa*fm < 0.) { b = m; fb = fm; } else { a = m; fa = fm; };
    if (fa*fp < 0.) { b = p; fb = fp; } else { a = p; fa = fp; };

    counter++;
    if(counter >= 100) {
      CFLog(INFO, "maximum iteration reached to find filter cutoff wave number \n");
      break;
    }
  } while( std::abs(fp) >= 0.001 );

  CFreal k_value = p;
  // return k_value*k_dir;
  return k_value;

}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::calculateMinimumTransferFunction()
{
  // CFreal maxStep = calculateMaximalWaveNumber(centreCellID).norm2()/8;
  // CFreal k1(0.0);
  // CFreal k2(maxStep/10.);
  // Gfunc G(centreCellID,this);
  // 
  // MathTools::Brent findMinimum;
  // findMinimum.setMaxStep(maxStep);
  // findMinimum.bracket(k1,k2,G);
  // return findMinimum.minimize(G);
  
  CFreal kmax = calculateMaximalWaveNumber();
  Gfunc G(this);
  
  // CFreal kmin = kmax/getMethodData().getFilterGridRatio()/3;
  CFreal kmin = 0.0;
  CFuint zones = 4;
  CFreal step = (kmax-kmin)/CFreal(zones);
  CFreal location, lower, upper;
  
  for(CFuint i=0; i<zones; ++i) {
    lower = kmin + i*step;
    upper = kmin + (i+1)*step;
    
    try {
      location = boost::math::tools::brent_find_minima(G, lower, upper, 32).first;
    }
    catch(std::domain_error& err) {
       CFLog(INFO, "Exception caught in function ExplicitFilters::FilterStrategy::calculateMinimumTransferFunction: \n");
       CFLog(INFO, err.what() << "\n");
      location = kmax;
    }
    if (location < upper && location > kmin + 1e-5) {
      // minimum is found within these bounds
      // CFLog(INFO, "location minimum = " << location << " \n");
      break;
    }
  }
  return location;
}

//////////////////////////////////////////////////////////////////////////////

CFreal FilterStrategy::integrateTransferFunction(const RealVector &kmin, const RealVector& kmax)
{
  Framework::Node& X0 = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,getStencilID());
  CFreal DeltaX = m_stencil->getCellWidth();
  CFreal DeltaY = DeltaX;
  RealVector D(2);
  D[XX] = DeltaX;     D[YY] = DeltaY;

  CFreal I = m_weights[0]*(kmax[XX]-kmin[XX])*(kmax[YY]-kmin[YY]);

  CFuint nbCells = m_stencil->getNbElements();
  for(CFuint iCell=1; iCell<nbCells; ++iCell) {
    Framework::Node& X = getMethodData().getCoordinateLinker()->getCoordinates(m_stencil,iCell);
    RealVector dX = X-X0;

    I += m_weights[iCell] / (dX[XX]*dX[YY]) *
      (
        - cos(MathFunctions::innerProd(kmax,dX)) - cos(MathFunctions::innerProd(kmin,dX))
        + cos(kmin[XX]*dX[XX] + kmax[YY]*dX[YY]) + cos(kmax[XX]*dX[XX] + kmin[YY]*dX[YY])
      );
  }
  I /= (kmax[XX]-kmin[XX])*(kmax[YY]-kmin[YY]);

  return I;
}

//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace ExplicitFilters
		
  } // namespace Numerics

} // namespace COOLFluiD

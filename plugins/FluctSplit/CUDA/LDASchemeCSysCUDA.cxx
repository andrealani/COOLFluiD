#include "FluctSplit/CUDA/LDASchemeCSysCUDA.hh"
#include "FluctSplit/CUDA/LDAC_CUDA.hh"
#include "FluctSplit/CUDA/FluctSplitCUDA.hh"

#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MeshData.hh"
#include "Common/CUDA/CudaEnv.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LDASchemeCSysCUDA,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitCUDAModule>
ldacSchemeSysCUDAProvider("SysLDACCUDA");
      
//////////////////////////////////////////////////////////////////////////////

LDASchemeCSysCUDA::LDASchemeCSysCUDA(const std::string& name) :
  RDS_SplitterSys(name),
  _sumKplus(),
  _invK(),
  _uTemp(),
  dev_a(NULL),
  dev_b(NULL),
  dev_c(NULL),
  dev_d(NULL)
{
}

//////////////////////////////////////////////////////////////////////////////

LDASchemeCSysCUDA::~LDASchemeCSysCUDA()
{
}


//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSysCUDA::setup()
{
  RDS_SplitterSys::setup();

  _uTemp.resize(_nbEquations);
  _sumKplus.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _beta.resize(_nbEquations, _nbEquations);
  _k.resize(_nbEquations, _nbEquations);

  const CFuint matSize = _sumKplus.size(); 
  CudaEnv::allocDev(dev_a, matSize);
  CudaEnv::allocDev(dev_b, matSize);
  CudaEnv::allocDev(dev_c, matSize);
  CudaEnv::allocDev(dev_d, matSize);
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSysCUDA::unsetup()
{
  CudaEnv::free(dev_a);
  CudaEnv::free(dev_b);
  CudaEnv::free(dev_c);
  CudaEnv::free(dev_d);
  
  RDS_SplitterSys::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSysCUDA::distribute(vector<RealVector>& residual)
{
  DistributionData& ddata = getMethodData().getDistributionData();

  const RealVector& phiT = ddata.phi;
  
  addToKplusCUDA(_sumKplus.size(), dev_a, dev_b, dev_c, dev_d, 
		 &(*_kPlus[0])[0], &(*_kPlus[1])[0], &(*_kPlus[2])[0], &_sumKplus[0]);

/*  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState)
  {
    _sumKplus  += *_kPlus[iState];
  }*/
 
   
  _inverter->invert(_sumKplus, _invK);
  
  _uTemp = _invK * phiT;
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
  {
    residual[iState] = (*_kPlus[iState]) * _uTemp;

    if (ddata.computeBetas)
    {
      (*ddata.currBetaMat)[iState] = (*_kPlus[iState]) * _invK;
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSysCUDA::distributePart(vector<RealVector>& residual)
{
/*  RealSliceVector::setSize(_nbEquations);

  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _sumKplus  += *_kPlus[iState];
  }
  
  _inverter->invert(_sumKplus, _invK);
  
  RealVector& phi = getMethodData().getDistributionData().phi;
  _uTemp.slice(0) = _invK * phi.slice(_firstVarID);
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState].slice(_firstVarID) =
      (*_kPlus[iState]) * _uTemp.slice(0);
  }*/
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSysCUDA::computePicardJacob(vector<RealMatrix*>& jacob)
{
  // carefull with the signs !!!
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    // beta coefficient
    _beta = (*_kPlus[iState])*_invK;
    const CFuint nStart = iState*_nbStatesInCell;

    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      RealMatrix *const block = jacob[nStart + jState];

      _k =  (*_kPlus[jState]) + (*_kMin[jState]);
      (*block) = _beta*_k;
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#ifndef COOLFluiD_Numerics_FluctSplit_BSchemeBase_hh
#define COOLFluiD_Numerics_FluctSplit_BSchemeBase_hh

//////////////////////////////////////////////////////////////////////////////

#include "Splitter.hh"
#include "DistributionData.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the common data and functions for all B schemes
/// @author Tiago Quintino
template < typename BASE_SPLITTER >
class BSchemeBase : public BASE_SPLITTER {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit BSchemeBase(const std::string& name);

  /// Destructor
  virtual ~BSchemeBase();

  /// Configure this object with user defined parameters
  /// @param args arguments for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup this object with data depending on the mesh
  virtual void setup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
  providesSockets();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

protected: // functions

  /// Setup this object with data depending on the mesh
  virtual void storeThetas();

protected: // data

  /// socket with theta clending coefficient
  Framework::DataSocketSink<CFreal> socket_thetas;

  /// remembers what was the last cell computed
  CFuint m_cellid;

  /// remembers what was the last cell computed
  CFuint m_subcellid;

  /// max number of subcells
  CFuint m_max_nbsubcells;

  /// store the thetas
  bool m_store_thetas;

  RealVector m_uTemp;
  
  RealVector m_theta;

  std::vector<RealVector> m_phiN;

  RealVector m_sumPhiN;

  /// min theta
  CFreal m_min_theta;
  
  /// freeze thetas
  CFuint m_freezeTheta;

  /// 
  bool m_addExtraDiss;
  
  /// flag telling if to compute a first order jacobian
  bool m_firstOrderJacob;
  
}; // end of class BSchemeBase

//////////////////////////////////////////////////////////////////////////////

template < typename BASE_SPLITTER >
void BSchemeBase<BASE_SPLITTER>::defineConfigOptions(Config::OptionList& options)
{
   options.template addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
   options.template addConfigOption< bool >  ("StoreThetas","Store the thetas for visualization");
   options.template addConfigOption< CFreal, Config::DynamicOption<> >("MinTheta","Minimum theta, used to keep a minimum diffusion");
   options.template addConfigOption< CFuint, Config::DynamicOption<> >("FreezeTheta","If =1, freeze the theta coefficients to allow convergence");
   options.template addConfigOption< bool >("FirstOrderJacob","Compute a first order jacobian.");

   options.template addConfigOption< bool >("ExtraDiss","Include and isotropic dissipative term.");

}

//////////////////////////////////////////////////////////////////////////////

template < typename BASE_SPLITTER >
BSchemeBase<BASE_SPLITTER>::BSchemeBase(const std::string& name) :
  BASE_SPLITTER(name),
  socket_thetas("thetas",false),
  m_cellid(0),
  m_subcellid(0),
  m_uTemp(),
  m_theta(),
  m_phiN(),
  m_sumPhiN()
{
  this->addConfigOptionsTo(this);
  
  m_max_nbsubcells = 1;
  this->setParameter("MaxNbSubElems",&m_max_nbsubcells);
  
  m_store_thetas = false;
  this->setParameter("StoreThetas",&m_store_thetas);

  m_min_theta = 0.;
  this->setParameter("MinTheta",&m_min_theta);
  
  m_freezeTheta = 0;
  this->setParameter("FreezeTheta",&m_freezeTheta);
  
  m_firstOrderJacob = false;
  this->setParameter("FirstOrderJacob",&m_firstOrderJacob);

  m_addExtraDiss = true;
  this->setParameter("ExtraDiss",&m_addExtraDiss);

}

//////////////////////////////////////////////////////////////////////////////

template < typename BASE_SPLITTER >
BSchemeBase<BASE_SPLITTER>::~BSchemeBase()
{
}

//////////////////////////////////////////////////////////////////////////////

template < typename BASE_SPLITTER >
void BSchemeBase<BASE_SPLITTER>::configure ( Config::ConfigArgs& args )
{
  BASE_SPLITTER::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

template < typename BASE_SPLITTER >
void BSchemeBase<BASE_SPLITTER>::setup()
{
  BASE_SPLITTER::setup();
  
  const CFuint nbEqs = BASE_SPLITTER::_nbEquations; // AL: this is needed for sys+scalar schemes
  // Framework::PhysicalModelStack::getActive()->getNbEq();
  
  m_theta.resize(nbEqs);
  m_uTemp.resize(nbEqs);
  m_sumPhiN.resize(nbEqs);
  
  const CFuint maxNbStatesInCell = BASE_SPLITTER::_kPlus.size();

//   cf_assert ( _kPlus.size() == MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell() );
  
  m_phiN.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; iState++)
  {
    m_phiN[iState].resize(nbEqs);
  }

  if ( m_store_thetas && !socket_thetas.isConnected())
    throw Common::BadValueException (FromHere(),"User required storing of thetas but socket is not connected");
}

//////////////////////////////////////////////////////////////////////////////

template < typename BASE_SPLITTER >
inline void BSchemeBase<BASE_SPLITTER>::storeThetas()
{
  cf_assert(m_store_thetas);

  DistributionData& distdata = BASE_SPLITTER::getMethodData().getDistributionData();

  if (distdata.isPerturb) return; // skip if is being perturbed
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  
  // const CFuint nbEqs = BASE_SPLITTER::_nbEquations;
  // find the sub element index
  if (m_cellid != distdata.cellID)
  {
    m_cellid = distdata.cellID;
    m_subcellid = 0;
  }
  else
  {
    ++m_subcellid;
  }
  
  if (m_freezeTheta == 0) {    
    Framework::DataHandle< CFreal > thetas = socket_thetas.getDataHandle();
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      thetas( m_cellid*m_max_nbsubcells + m_subcellid, iEq, nbEqs) = m_theta[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

template < typename BASE_SPLITTER >
std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
BSchemeBase<BASE_SPLITTER>::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result = BASE_SPLITTER::providesSockets();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

template < typename BASE_SPLITTER >
std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
BSchemeBase<BASE_SPLITTER>::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = BASE_SPLITTER::needsSockets();
  result.push_back(&socket_thetas);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BSchemeBase_hh

#include "BasicCatalycityModel.hh"
#include "Catalycity.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Catalycity {
    
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<BasicCatalycityModel, 
			    CatalycityModel, 
			    CatalycityModule, 1>
basicCatalycityModelProv("Basic");
    
//////////////////////////////////////////////////////////////////////////////

void BasicCatalycityModel::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("ModelID","ID identifying the basic model to use."); 
  options.addConfigOption< vector<CFreal> >("CatalycityCoeff","Catalycity coefficients ranging from 0. to 1.0 one per each species");
}
    
//////////////////////////////////////////////////////////////////////////////

BasicCatalycityModel::BasicCatalycityModel(const std::string& name) : 
  CatalycityModel(name),
  m_mmasses(),
  m_mconcentrations()
{
  addConfigOptionsTo(this);
  
  m_modelID = 0;
  setParameter("ModelID",&m_modelID);
  
  m_catalycity = vector<CFreal>();
  setParameter("CatalycityCoeff",&m_catalycity);
}
    
//////////////////////////////////////////////////////////////////////////////
    
BasicCatalycityModel::~BasicCatalycityModel()
{
}
    
//////////////////////////////////////////////////////////////////////////////
 
void BasicCatalycityModel::setup()
{
  CatalycityModel::setup();
  
  const CFuint nbSpecies = m_library->getNbSpecies();
  m_mmasses.resize(nbSpecies);
  m_library->getMolarMasses(m_mmasses);
  
  m_mconcentrations.resize(nbSpecies);
  
  if (m_catalycity.size() != nbSpecies) {
    throw BadValueException(FromHere(), "CatalycityCoeff size != number of species");
  }
}

//////////////////////////////////////////////////////////////////////////////
    
void BasicCatalycityModel::computeMassProduction(const CFreal temp, 
						 const CFreal rho, 
						 const RealVector& ys, 
						 RealVector& wallMassProduction)
{  
  //	initialization
  CFreal reactantconsumption = 0.0;
  wallMassProduction = 0.;
  
  // [x_i] = rho_i/M_i = rho*y_i/M_i = rho*x_i/M
  m_mconcentrations = rho*(ys/m_mmasses);
  
  switch (m_modelID) {
  case 0:
    //air 5
    // Oxygen recombination only, mutation order
    catalycityMoleProduction(m_mconcentrations[0],m_mmasses[0],temp,m_catalycity[0], reactantconsumption);
    wallMassProduction[2] = -0.5*reactantconsumption;
    wallMassProduction[0] = reactantconsumption;
    break;
  case 1:
    //air 5
    // Nitrogen recombination only, mutation order
    catalycityMoleProduction(m_mconcentrations[1],m_mmasses[1],temp,m_catalycity[0], reactantconsumption);
    wallMassProduction[4] = -0.5*reactantconsumption;
    wallMassProduction[1] = reactantconsumption;
    break;
  case 2:
    //air 5
    // Oxygen and Nitrogen recombination only using equal catalycity, mutation order
    catalycityMoleProduction(m_mconcentrations[0],m_mmasses[0],temp, m_catalycity[0], reactantconsumption);
    wallMassProduction[2] = -0.5*reactantconsumption;
    wallMassProduction[0] = reactantconsumption;
    
    catalycityMoleProduction(m_mconcentrations[1],m_mmasses[1],temp, m_catalycity[0], reactantconsumption);
    wallMassProduction[4] = -0.5*reactantconsumption;
    wallMassProduction[1] = reactantconsumption;
    break;
    
  default:
    throw BadValueException(FromHere(), "ModelID does NOT exist");
    break;
  }
  
  m_library->setSpeciesFractions(ys); // this should not be necessary, MUST BE CHECKED
  wallMassProduction *= m_library->getMMass();
}

//////////////////////////////////////////////////////////////////////////////
  
} // namespace Catalycity
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

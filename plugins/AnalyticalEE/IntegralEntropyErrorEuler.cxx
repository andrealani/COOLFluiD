#include "Environment/DirPaths.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"

#include "AnalyticalEE/AnalyticalEE.hh"
#include "AnalyticalEE/IntegralEntropyErrorEuler.hh"

#include "Framework/MeshData.hh"
#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
// using namespace COOLFluiD::Physics::NavierStokes;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<IntegralEntropyErrorEuler, AnalyticEEData, AnalyticalEEModule> computeIntegralEntropyErrorErrorProvider("IntegralEntropyErrorEuler");

//////////////////////////////////////////////////////////////////////////////

void IntegralEntropyErrorEuler::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Functions","Definition of the Functions.");
  options.addConfigOption< std::string >("OutputFile","Name of Output File to write the Error Norm.");
  options.addConfigOption< CFreal >("ReferenceEntropy", "Reference value from which the entropy deviation is computed.");
}

//////////////////////////////////////////////////////////////////////////////

IntegralEntropyErrorEuler::IntegralEntropyErrorEuler(std::string name) : ComputeIntegralError(name), m_physicalData()
// ,
//   socket_nodes("nodes"),
//   socket_states("states"),
//   m_L2(0),
//   m_result(0)
{
  m_file = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);

  m_vars = std::vector<std::string>();
  setParameter("Vars",&m_vars);

  m_functions = std::vector<std::string>();
  setParameter("Functions",&m_functions);

  m_nameOutputFile = "IntegralErrorNorm.plt";
  setParameter("OutputFile",&m_nameOutputFile);

  m_refEntropy = 1.0;
  setParameter("ReferenceEntropy",&m_refEntropy);

}

//////////////////////////////////////////////////////////////////////////////

IntegralEntropyErrorEuler::~IntegralEntropyErrorEuler()
{
}

//////////////////////////////////////////////////////////////////////////////

void IntegralEntropyErrorEuler::setup()
{
  CFAUTOTRACE;
  // first call parent method (setup of abstract base class)
  ComputeIntegralError::setup();

  //we won't need the following two variables:
  m_exact.resize(0); 
  m_result.resize(0);

  m_L2.resize(1);
  prepareOutputFile();

  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  updateVarSet->setup();

  m_qdState = new Framework::State();
}

//////////////////////////////////////////////////////////////////////////////

void IntegralEntropyErrorEuler::unsetup()
{
  CFAUTOTRACE;

  delete m_qdState;

  // last call parent method
  AnalyticEECom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void IntegralEntropyErrorEuler::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ComputeIntegralError::configure(args);

  m_vFunction.setFunctions(m_functions);


  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
IntegralEntropyErrorEuler::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = ComputeIntegralError::providesSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
IntegralEntropyErrorEuler::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = ComputeIntegralError::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void IntegralEntropyErrorEuler::execute()
{
  CFAUTOTRACE;

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");

  //prepares to loop over cells by getting the GeometricEntityPool
  StdTrsGeoBuilder::GeoData& geoData = m_stdTrsGeoBuilder.getDataGE();
  geoData.trs = trs;

  m_L2[0] = 0.;

  /// loop over element types
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    /// loop over elements
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      // build the GeometricEntity
      geoData.idx = cellIdx;

      GeometricEntity *const cell = m_stdTrsGeoBuilder.buildGE();

      integralErrorInCell(cell,m_error);

      m_L2[0] += m_error[0];


      m_stdTrsGeoBuilder.releaseGE();

    } //Loop over elements

  }//Loop over element types


    m_L2[0] = std::sqrt(m_L2[0]);
    //m_result[0] = std::log(m_L2[0]);

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  openFile(true);

  cf_assert(m_file->isopen());
  ofstream& out = m_file->get();  

  out << iter << " " ;
    out << std::log(m_L2[0]) << " "<< m_L2[0] << " ";

  out <<"\n\n";

  m_file->close();
}

//////////////////////////////////////////////////////////////////////////////

void IntegralEntropyErrorEuler::openFile(bool append)
{
  CFAUTOTRACE;

  cf_assert(!m_file->isopen());

  boost::filesystem::path fpath = Environment::DirPaths::getInstance().getResultsDir() / 
  boost::filesystem::path(m_nameOutputFile);
#ifdef CF_HAVE_BOOST_1_85
  fpath.replace_extension(".plt");
#else
  boost::filesystem::change_extension(fpath,".plt");
#endif
 
  fpath = Framework::PathAppender::getInstance().appendParallel( fpath );

  if (append)
  {
    m_file->open(fpath,ios::app);
  }
  else
  {
    m_file->open(fpath);
  }

  m_file->get().precision(16);
}

//////////////////////////////////////////////////////////////////////////////

void IntegralEntropyErrorEuler::prepareOutputFile()
{
  CFAUTOTRACE;

  openFile(false);

  cf_assert(m_file->isopen());
  ofstream& out = m_file->get();

  out << "TITLE  =  Integrated Analytical Error Norm of Entropy Error" << "\n";
  out << "VARIABLES = Iter "; 
    out << "LogL2_Entropy L2_Entropy ";
  out <<"\n";

  m_file->close();
}

//////////////////////////////////////////////////////////////////////////////

CFreal IntegralEntropyErrorEuler::JacobianInTriagP2(const std::vector<Framework::Node*>* triagnodes, const RealVector& mapped_coord)
{
  const CFreal x0 = (*(*triagnodes)[0])[XX];
  const CFreal x1 = (*(*triagnodes)[1])[XX];
  const CFreal x2 = (*(*triagnodes)[2])[XX];
  const CFreal x3 = (*(*triagnodes)[3])[XX];
  const CFreal x4 = (*(*triagnodes)[4])[XX];
  const CFreal x5 = (*(*triagnodes)[5])[XX];

  const CFreal y0 = (*(*triagnodes)[0])[YY];
  const CFreal y1 = (*(*triagnodes)[1])[YY];
  const CFreal y2 = (*(*triagnodes)[2])[YY];
  const CFreal y3 = (*(*triagnodes)[3])[YY];
  const CFreal y4 = (*(*triagnodes)[4])[YY];
  const CFreal y5 = (*(*triagnodes)[5])[YY];

  const CFreal dN0dxi = -3 + 4*mapped_coord[YY]+4*mapped_coord[XX];
  const CFreal dN1dxi = 4*mapped_coord[XX]-1.0;
  const CFreal dN2dxi = 0.0;
  const CFreal dN3dxi = 4 - 8*mapped_coord[XX]-4*mapped_coord[YY];
  const CFreal dN4dxi = 4*mapped_coord[YY];
  const CFreal dN5dxi = - 4 * mapped_coord[YY];

  const CFreal dN0deta = -3 + 4*mapped_coord[YY]+4*mapped_coord[XX];
  const CFreal dN1deta = 0.0;
  const CFreal dN2deta = 4*mapped_coord[YY]-1.0;
  const CFreal dN3deta = - 4 * mapped_coord[XX];
  const CFreal dN4deta = 4*mapped_coord[XX];
  const CFreal dN5deta = 4 - 8*mapped_coord[YY]-4*mapped_coord[XX];

  const CFreal dxdxi = dN0dxi*x0 + dN1dxi*x1 + dN2dxi*x2 + dN3dxi*x3 + dN4dxi*x4 + dN5dxi*x5;
  const CFreal dydxi = dN0dxi*y0 + dN1dxi*y1 + dN2dxi*y2 + dN3dxi*y3 + dN4dxi*y4 + dN5dxi*y5;

  const CFreal dxdeta = dN0deta*x0 + dN1deta*x1 + dN2deta*x2 + dN3deta*x3 + dN4deta*x4 + dN5deta*x5;
  const CFreal dydeta = dN0deta*y0 + dN1deta*y1 + dN2deta*y2 + dN3deta*y3 + dN4deta*y4 + dN5deta*y5;

  return std::abs(dxdxi*dydeta - dxdeta*dydxi);

}


//////////////////////////////////////////////////////////////////////////////


CFreal IntegralEntropyErrorEuler::JacobianInTriagP3(const std::vector<Framework::Node*>* triagnodes, const RealVector& mapped_coord)
{

  const CFreal x0 = (*(*triagnodes)[0])[XX];
  const CFreal x1 = (*(*triagnodes)[1])[XX];
  const CFreal x2 = (*(*triagnodes)[2])[XX];
  const CFreal x3 = (*(*triagnodes)[3])[XX];
  const CFreal x4 = (*(*triagnodes)[4])[XX];
  const CFreal x5 = (*(*triagnodes)[5])[XX];
  const CFreal x6 = (*(*triagnodes)[6])[XX];
  const CFreal x7 = (*(*triagnodes)[7])[XX];
  const CFreal x8 = (*(*triagnodes)[8])[XX];
  const CFreal x9 = (*(*triagnodes)[9])[XX];

  const CFreal y0 = (*(*triagnodes)[0])[YY];
  const CFreal y1 = (*(*triagnodes)[1])[YY];
  const CFreal y2 = (*(*triagnodes)[2])[YY];
  const CFreal y3 = (*(*triagnodes)[3])[YY];
  const CFreal y4 = (*(*triagnodes)[4])[YY];
  const CFreal y5 = (*(*triagnodes)[5])[YY];
  const CFreal y6 = (*(*triagnodes)[6])[YY];
  const CFreal y7 = (*(*triagnodes)[7])[YY];
  const CFreal y8 = (*(*triagnodes)[8])[YY];
  const CFreal y9 = (*(*triagnodes)[9])[YY];

  const CFreal L0 = 1.0-mapped_coord[XX]-mapped_coord[YY];
  const CFreal L1 = mapped_coord[XX];
  const CFreal L2 = mapped_coord[YY];

  const CFreal dN0dxi = -0.5*(9.0*L0*(2.0*L0-1.0) + (3.0*L0-1.0)*(3.0*L0-2.0));
  const CFreal dN0deta = -0.5*(9.0*L0*(2.0*L0-1.0) + (3.0*L0-1.0)*(3.0*L0-2.0));

  const CFreal dN1dxi = 0.5*(9.0*L1*(2.0*L1-1.0) + (3.0*L1-1.0)*(3.0*L1-2.0));
  const CFreal dN1deta = 0.0;

  const CFreal dN2dxi = 0.0;
  const CFreal dN2deta = 0.5*(9.0*L2*(2.0*L2-1.0) + (3.0*L2-1.0)*(3.0*L2-2.0));

  const CFreal dN3dxi = 4.5*(3.0*L0*L0 - 6.0*L0*L1 + L1 - L0);
  const CFreal dN3deta = 4.5*(L1 - 6.0*L0*L1);

  const CFreal dN4dxi = 4.5*(-3.0*L1*L1 + 6.0*L0*L1 + L1 - L0);
  const CFreal dN4deta = 4.5*(L1 - 3.0*L1*L1);

  const CFreal dN5dxi = 4.5*(6.0*L1*L2 - L2);
  const CFreal dN5deta = 4.5*(3.0*L1*L1 - L1);

  const CFreal dN6dxi = 4.5*(3.0*L2*L2 - L2);
  const CFreal dN6deta = 4.5*(6.0*L1*L2 - L1);

  const CFreal dN7dxi = 4.5*(L2 - 3.0*L2*L2);
  const CFreal dN7deta = 4.5*(-3.0*L2*L2 + 6.0*L0*L2 - L0 + L2);

  const CFreal dN8dxi = 4.5*(L2 - 6.0*L0*L2);
  const CFreal dN8deta = 4.5*(3.0*L0*L0 - 6.0*L0*L2 + L2 - L0);

  const CFreal dN9dxi = 27.0*L2*(L0-L1);
  const CFreal dN9deta = 27.0*L1*(L0-L2);


  const CFreal dxdxi = dN0dxi*x0 + dN1dxi*x1 + dN2dxi*x2 + dN3dxi*x3 + dN4dxi*x4 + dN5dxi*x5 +\
		       dN6dxi*x6 + dN7dxi*x7 + dN8dxi*x8 + dN9dxi*x9;
  const CFreal dydxi = dN0dxi*y0 + dN1dxi*y1 + dN2dxi*y2 + dN3dxi*y3 + dN4dxi*y4 + dN5dxi*y5 +\
		       dN6dxi*y6 + dN7dxi*y7 + dN8dxi*y8 + dN9dxi*y9;

  const CFreal dxdeta = dN0deta*x0 + dN1deta*x1 + dN2deta*x2 + dN3deta*x3 + dN4deta*x4 + dN5deta*x5 +\
		      	dN6deta*x6 + dN7deta*x7 + dN8deta*x8 + dN9deta*x9;
  const CFreal dydeta = dN0deta*y0 + dN1deta*y1 + dN2deta*y2 + dN3deta*y3 + dN4deta*y4 + dN5deta*y5 +\
			dN6deta*y6 + dN7deta*y7 + dN8deta*y8 + dN9deta*y9;

  return std::abs(dxdxi*dydeta - dxdeta*dydxi);

}







//////////////////////////////////////////////////////////////////////////////

void IntegralEntropyErrorEuler::integralErrorInCell(Framework::GeometricEntity *const cell, RealVector& error)
{
  vector<Node*>*  cellNodes  = cell->getNodes();
  vector<State*>* cellStates = cell->getStates();


  const CFuint nbNodes = cellNodes->size();


  const CFuint nbStates = cellStates->size();
  RealVector shapeFunctionValues(nbStates);

//   m_exact = 0.0;
  error[0] = 0.0;

    SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

  ///Numerical quadrature: hardcoded for triangles
  for(CFuint iQdPt = 0; iQdPt < nbQdPts; ++iQdPt) {

  m_mappedCoord[XX] = xQdPtsTriag[iQdPt];
  m_mappedCoord[YY] = yQdPtsTriag[iQdPt];
  m_numerical = 0.0;

  CFreal J = 1.0;

    ///In case the of P2P2 triangle: nonlinear jacobian
    if(nbNodes == 6)
    {
	J = JacobianInTriagP2(cellNodes, m_mappedCoord);
    } 
    else if(nbNodes == 10)
    {
	J = JacobianInTriagP3(cellNodes, m_mappedCoord);
    } //nonlinear jacobian

  shapeFunctionValues = cell->computeShapeFunctionAtMappedCoord(m_mappedCoord);

  for(CFuint istate=0;istate < nbStates; ++istate)
  {
    m_numerical = m_numerical + shapeFunctionValues[istate] * ( *((*cellStates)[istate]->getData()) );
  }

    for(CFuint i=0;i<m_nbEqs; ++i) (*m_qdState)[i] = m_numerical[i];

    updateVarSet->setDimensionalValues(*m_qdState,m_numerical);

    const CFreal gamma = 1.4;
//     CFreal gamma = getModel()->getGamma();
    const CFreal rho = m_numerical[0];
    const CFreal p = ( gamma - 1.0 ) * (m_numerical[3] - 0.5 * (m_numerical[1]*m_numerical[1] + m_numerical[2]*m_numerical[2])/rho);

    ///Use the physical definition of entropy:
    //const CFreal s = std::log(p) - gamma * std::log(rho);  
    //error[0] += qdWeightsTriag[iQdPt] * (s - m_refEntropy) * (s - m_refEntropy) * J;

    ///Compute the entropy error the same way as Krivodonova:
    const CFreal s = p / std::pow(rho,gamma);  
    error[0] += qdWeightsTriag[iQdPt] * (s / m_refEntropy - 1) * (s / m_refEntropy - 1) * J;

//     error[0] += qdWeightsTriag[iQdPt] * J; //This would just compute the surface of the element

  } //Loop over quadrature points

  if(nbNodes == 3) {
  const CFreal volume = cell->computeVolume();
  error[0] = std::abs(volume) * error[0]; }

  else 
  error[0] = 0.5 * error[0]; //0.5 is the area of P2P2 triangle in reference space


}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

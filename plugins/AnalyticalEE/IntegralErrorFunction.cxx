#include "Environment/DirPaths.hh"

#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"

#include "AnalyticalEE/AnalyticalEE.hh"
#include "AnalyticalEE/IntegralErrorFunction.hh"

#include "Framework/MeshData.hh"
#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<IntegralErrorFunction, AnalyticEEData, AnalyticalEEModule> computeIntegralErrorFunctionProvider("IntegralErrorFunction");

//////////////////////////////////////////////////////////////////////////////

void IntegralErrorFunction::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Functions","Definition of the Functions.");
  options.addConfigOption< std::string >("OutputFile","Name of Output File to write the Error Norm.");
}

//////////////////////////////////////////////////////////////////////////////

IntegralErrorFunction::IntegralErrorFunction(std::string name) : ComputeIntegralError(name)
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

//   m_vars[0] = "x";
//   m_vars[1] = "y";

  m_functions = std::vector<std::string>();
  setParameter("Functions",&m_functions);

  m_nameOutputFile = "IntegralErrorNorm.plt";
  setParameter("OutputFile",&m_nameOutputFile);

}

//////////////////////////////////////////////////////////////////////////////

IntegralErrorFunction::~IntegralErrorFunction()
{
}

//////////////////////////////////////////////////////////////////////////////

void IntegralErrorFunction::setup()
{
  CFAUTOTRACE;

  // first call parent method (setup of abstract base class)
  ComputeIntegralError::setup();

//   const CFuint nbeqs = PhysicalModelStack::getActive()->getNbEq();
//   m_exact.resize(nbeqs);
//   m_L2.resize(nbeqs);
//   m_result.resize(nbeqs);
//   prepareOutputFile();

}

//////////////////////////////////////////////////////////////////////////////

void IntegralErrorFunction::unsetup()
{
  CFAUTOTRACE;

//   m_exact.resize(0);
//   m_L2.resize(0);
//   m_result.resize(0);
//   m_file->close();

  // last call parent method
  AnalyticEECom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void IntegralErrorFunction::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ComputeIntegralError::configure(args);

  m_vFunction.setFunctions(m_functions);


//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
//   const CFuint dim = 2;
// 
//   m_vars = std::vector<std::string>();
//   m_vars.resize(dim);
//   m_vars[0] = "x";
//   m_vars[1] = "y";
// 
//   if(dim == DIM_3D) m_vars[2] = "z";
// 
//   CFout << "These are the variables:\n";
//   CFout << m_vars[0] << "\n";
//   CFout << m_vars[1] << "\n";
//   CFout << m_vars[2] << "\n";


  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }


//   CFout << "Compute discrete error:\n";
//   CFout << "Functions: " << m_functions[0] << "\n";
//   CFout << "Variables:\n";
//   CFout << "\t\t" << m_vars[0] << "\n";
//   CFout << "\t\t" << m_vars[1] << "\n";
//   CFout << "Size of variables' field: " << m_vars.size() << "\n";
//   CFout << "Output file: " << m_nameOutputFile << "\n";
//   CF_DEBUG_EXIT;

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
IntegralErrorFunction::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
IntegralErrorFunction::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void IntegralErrorFunction::execute()
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

  m_L2 = 0.;

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

      for(CFuint iEq =0; iEq < m_nbEqs; iEq ++)
        m_L2[iEq] += m_error[iEq];


      m_stdTrsGeoBuilder.releaseGE();

    } //Loop over elements

  }//Loop over element types


//   CF_DEBUG_EXIT;


  for(CFuint iEq =0; iEq < m_L2.size(); iEq ++){
    m_L2[iEq] = std::sqrt(m_L2[iEq]);
    m_result[iEq] = std::log(m_L2[iEq]);
  }

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  openFile(true);

  cf_assert(m_file->isopen());
  ofstream& out = m_file->get();  

  out << iter << " " ;
  for(CFuint iEq =0; iEq < m_L2.size(); iEq ++){
    out << m_result[iEq] << " "<< m_L2[iEq] << " ";
  }
  out <<"\n";
/*

for(CFuint iEq =0; iEq < m_L2.size(); iEq ++){
     cout <<  m_result[iEq] << "\t" << m_L2[iEq] << "\t";
} */
  out <<"\n";

  m_file->close();
}

//////////////////////////////////////////////////////////////////////////////

void IntegralErrorFunction::openFile(bool append)
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

void IntegralErrorFunction::prepareOutputFile()
{
  CFAUTOTRACE;

  openFile(false);

  cf_assert(m_file->isopen());
  ofstream& out = m_file->get();

  out << "TITLE  =  Integrated Analytical Error Norm" << "\n";
  out << "VARIABLES = Iter "; 
  for (CFuint iEq = 0; iEq < m_L2.size(); ++iEq)
    out << "LogL2_"<<iEq<< " L2_" <<iEq<< " ";

  out <<"\n";

  m_file->close();
}

//////////////////////////////////////////////////////////////////////////////


void IntegralErrorFunction::integralErrorInCell(Framework::GeometricEntity *const cell, RealVector& error)
{
//   vector<Node*>* cellNodes = cell->getNodes();
  vector<State*>* cellStates = cell->getStates();

  const CFuint nbStates = cellStates->size();
  RealVector shapeFunctionValues(nbStates);

  m_exact = 0.0;
  error = 0.0;

  ///Numerical quadrature: hardcoded for triangles
  for(CFuint iQdPt = 0; iQdPt < nbQdPts; ++iQdPt) {

  m_mappedCoord[XX] = xQdPtsTriag[iQdPt];
  m_mappedCoord[YY] = yQdPtsTriag[iQdPt];
  m_numerical = 0.0;

  shapeFunctionValues = cell->computeShapeFunctionAtMappedCoord(m_mappedCoord);

  for(CFuint istate=0;istate < nbStates; ++istate)
  {
    m_numerical = m_numerical + shapeFunctionValues[istate] * ( *((*cellStates)[istate]->getData()) );
  }

  ///To evaluate the exact solution, we have to find coordinates in physical space!
  m_physicalCoord = cell->computeCoordFromMappedCoord(m_mappedCoord);
  m_vFunction.evaluate(m_physicalCoord,m_exact);

  for(CFuint ieq = 0; ieq < m_nbEqs; ++ieq)
  {
    error[ieq] += qdWeightsTriag[iQdPt] * (m_numerical[ieq] - m_exact[ieq]) * (m_numerical[ieq] - m_exact[ieq]);
  }

  } //Loop over quadrature points

  const CFreal volume = cell->computeVolume();
  error = volume * error;


}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"

#include "FluxReconstructionMethod/InitStateAddVar.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitStateAddVar, FluxReconstructionSolverData, FluxReconstructionModule>
InitStateAddVarProvider("InitStateAddVar");

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("InitVars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("InitDef","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

InitStateAddVar::InitStateAddVar(const std::string& name) :
  StdInitState(name),
  m_tmpFun(),
  m_tmpVars()
{
  addConfigOptionsTo(this);
  m_initFunctions = std::vector<std::string>();
  setParameter("InitDef",&m_initFunctions);

  m_initVars = std::vector<std::string>();
  setParameter("InitVars",&m_initVars);
}

//////////////////////////////////////////////////////////////////////////////

InitStateAddVar::~InitStateAddVar()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::setup()
{

  CFAUTOTRACE;
  StdInitState::setup();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_tmpFun.resize(m_initFunctions.size());
  cf_assert(m_tmpFun.size() > 0);

  m_tmpVars.resize(dim + m_tmpFun.size());
  cf_assert(m_tmpVars.size() > 0); 
}

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::unsetup()
{
  CFAUTOTRACE;

  StdInitState::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdInitState::configure(args);

  m_vInitFunction.setFunctions(m_initFunctions);
  m_vInitFunction.setVariables(m_initVars);
  try {
    m_vInitFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

}

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::executeOnTrs()
{
  CFAUTOTRACE;

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  CFLogDebugMin("InitStateAddVar::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  //geoBuilder -> setup();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  const CFuint nbAddVars = m_initFunctions.size();
  
  cf_assert(m_tmpVars.size() == dim + nbAddVars);

  // loop over elements
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    // local coordinates of initialization points and matrix used for initialization
    SafePtr< vector< RealVector > > locInitPntCoords = frLocalData[iElemType]->getSolPntsLocalCoords();
    SafePtr< RealMatrix >           initTransfMatrix = frLocalData[iElemType]->getInitTransfMatrix();

    // get number of states
    const CFuint nbrStates = frLocalData[iElemType]->getNbrOfSolPnts();

    // loop over elements
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      //c
      // build the GeometricEntity
      geoData.idx = cellIdx;
      GeometricEntity *const cell = geoBuilder->buildGE();

      // get the states
      vector<State*>* solPntStates = cell->getStates();
      cf_assert(solPntStates->size() == nbrStates);

      // loop over initialization points
      for (CFuint iPnt = 0; iPnt < nbrStates; ++iPnt)
      {
        // compute initialization node global coordinate
        m_initPntCoords = cell->computeCoordFromMappedCoord((*locInitPntCoords)[iPnt]);
        
        // first evaluate expressions for additional variables f1, f2, ...
        // as functions of (x,y,z)
        m_vInitFunction.evaluate(m_initPntCoords, m_tmpFun);

        // set the first dim components of the input vars for the vFunction
        m_tmpVars.slice(0,dim) = m_initPntCoords.slice(0,dim);

        // set the following additional variables with the result of the
        // already evaluated expressions
        m_tmpVars.slice(dim, nbAddVars) = m_tmpFun.slice(0,nbAddVars);

        // evaluate the state variables as functions of (x,y,z, f1, f2, ...)
        m_vFunction.evaluate(m_tmpVars, *m_inputState);
	
        // evaluate the function at the state coordinate
        //m_vFunction.evaluate(m_initPntCoords,*m_inputState);
        
        *m_initPntsStates[iPnt] = *m_inputToUpdateVar->transform(m_inputState);
	m_varSet->setAdimensionalValues(*m_initPntsStates[iPnt], *(*solPntStates)[iPnt]);
      }

//       // transform initialization points solutions to solution point solutions
//       for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
//       {
//         // variable for a dimensional state
//         State dimState(RealVector(0.0,m_nbrEqs));
// 
//         for (CFuint iPnt = 0; iPnt < nbrStates; ++iPnt)
//         {
//           dimState += (*initTransfMatrix)(iSol,iPnt)*(*m_initPntsStates[iPnt]);
//         }
// 
//         // adimensionalize the value if needed and store
//         m_varSet->setAdimensionalValues(dimState, *(*solPntStates)[iSol]);
//       }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

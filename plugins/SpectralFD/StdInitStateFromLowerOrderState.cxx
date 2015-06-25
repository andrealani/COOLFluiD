#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"

#include "Common/BadValueException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/FileFormatException.hh"

#include "SpectralFD/HexaSpectralFDElementData.hh"
#include "SpectralFD/QuadSpectralFDElementData.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/StdInitStateFromLowerOrderState.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdInitStateFromLowerOrderState, SpectralFDMethodData, SpectralFDModule>
StdInitStateFromLowerOrderProvider("StdInitStateFromLowerOrderState");

//////////////////////////////////////////////////////////////////////////////

void StdInitStateFromLowerOrderState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SourceSolutionFileName","Name of the file containing the lower order solution.");
}

//////////////////////////////////////////////////////////////////////////////

StdInitStateFromLowerOrderState::StdInitStateFromLowerOrderState(const std::string& name) :
  SpectralFDMethodCom(name),
  m_varSet(CFNULL),
  m_initPntsStates(),
  m_inputState(),
  m_nbrEqs(),
  m_dim(),
  m_initPntCoords(),
  m_srcSolCFmeshFileName(),
  m_srcSolPolyOrder(),
  m_srcSolSDData(CFNULL),
  m_solTransMatr(),
  m_cellSrcStates()
{
  addConfigOptionsTo(this);

  m_srcSolCFmeshFileName = "";
  setParameter("SourceSolutionFileName",&m_srcSolCFmeshFileName);
}

//////////////////////////////////////////////////////////////////////////////

StdInitStateFromLowerOrderState::~StdInitStateFromLowerOrderState()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    StdInitStateFromLowerOrderState::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdInitStateFromLowerOrderState::setup()
{

  CFAUTOTRACE;
  SpectralFDMethodCom::setup();

  // get number of physical variables
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

  // resize variables
  m_initPntCoords.resize(m_dim);

  // set maxNbStatesInCell
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_initPntsStates.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState)
  {
    m_initPntsStates[iState] = new State();
  }

  // create state for input variables
  m_inputState = new State();

  // get updateVarSet
  m_varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void StdInitStateFromLowerOrderState::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iState = 0; iState < m_initPntsStates.size(); ++iState)
  {
    deletePtr(m_initPntsStates[iState]);
  }
  m_initPntsStates.resize(0);

  deletePtr(m_inputState);

  SpectralFDMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void StdInitStateFromLowerOrderState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  SpectralFDMethodCom::configure(args);

  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
}

//////////////////////////////////////////////////////////////////////////////

void StdInitStateFromLowerOrderState::executeOnTrs()
{
  CFAUTOTRACE;

  // Path where to read files from
  const std::string append = Environment::DirPaths::getInstance().getWorkingDir().string() + "/";
  const std::string filename = append + m_srcSolCFmeshFileName;

  // open the CFmesh file containing the solution to initialize from
  ifstream srcSolCFmeshFile(filename.c_str());
  if (!srcSolCFmeshFile.is_open())
  {
    throw FileFormatException (FromHere(),"Cannot open " + m_srcSolCFmeshFileName + "...");
  }

  // read the header and perform some file checks
  readCFmeshFileHeaderData(srcSolCFmeshFile);

  // create solution transformation matrix
  computeSolutionTransformationMatrix();

  // skip ahead to states information
  std::string line;
  vector<std::string> words;
  do
  {
    getline(srcSolCFmeshFile,line);
  } while (line != "!LIST_STATE 1");

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  CFLogDebugMin("StdInitStateFromLowerOrderState::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  // loop over elements
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    // get number of states
    const CFuint nbrStatesCur = sdLocalData[iElemType]->getNbrOfSolPnts();
    cf_assert(nbrStatesCur == m_solTransMatr.nbRows());
    const CFuint nbrStatesSrc = m_solTransMatr.nbCols();

    // loop over elements
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      // build the GeometricEntity
      geoData.idx = cellIdx;
      GeometricEntity *const cell = geoBuilder->buildGE();

      // get the states
      vector<State*>* solPntStates = cell->getStates();
      cf_assert(solPntStates->size() == nbrStatesCur);

      // read source solution from file
      for (CFuint iSrc = 0; iSrc < nbrStatesSrc; ++iSrc)
      {
        getline(srcSolCFmeshFile,line);
        words = Common::StringOps::getWords(line);
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_cellSrcStates[iSrc][iEq] = Common::StringOps::from_str<CFreal>(words[iEq]);
        }
      }

      // transform to current solution
      for (CFuint iCur = 0; iCur < nbrStatesCur; ++iCur)
      {
        // variable for a dimensional state
        State dimState(RealVector(0.0,m_nbrEqs));

        for (CFuint iSrc = 0; iSrc < nbrStatesSrc; ++iSrc)
        {
          dimState += m_solTransMatr(iCur,iSrc)*m_cellSrcStates[iSrc];
        }

        // adimensionalize the value if needed and store
        m_varSet->setAdimensionalValues(dimState, *(*solPntStates)[iCur]);
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
  getline(srcSolCFmeshFile,line);

  // close the CFmesh file containing the solution to initialize from
  srcSolCFmeshFile.close();
}

//////////////////////////////////////////////////////////////////////////////

void StdInitStateFromLowerOrderState::readCFmeshFileHeaderData(ifstream& inputFile)
{
  // strings to store lines and words in
  std::string line;
  vector<std::string> words;

  // read all the lines in the file until !END tag, and take action depending on tag
  do
  {
    getline(inputFile,line);

    if (line[0] == '!')
    {
      words = Common::StringOps::getWords(line);
      cf_assert(words.size() > 0);
      if (words[0] == "!CFMESH_FORMAT_VERSION")
      {
        if (words.size() != 2)
        {
          throw FileFormatException (FromHere(),"Cannot read CFmesh file format version...");
        }
        try
        {
          if (words[1] != "1.3")
          {
            throw FileFormatException
                (FromHere(), "StdInitStateFromLowerOrderState is written to be used with CFmesh file format 1.3, format of given CFmesh file is " + words[1]);
          }
        }
        catch (FileFormatException e)
        {
        }
      }
      else if (words[0] == "!NB_DIM")
      {
        if (words.size() != 2)
        {
          throw FileFormatException (FromHere(),"Cannot read dimensionality...");
        }
        if (Common::StringOps::from_str<CFuint>(words[1]) != m_dim)
        {
          throw BadValueException (FromHere(),"Wrong dimensionality in given CFmesh file: " + words[1]);
        }
      }
      else if (words[0] == "!NB_EQ")
      {
        if (words.size() != 2)
        {
          throw FileFormatException (FromHere(),"Cannot read number of equations...");
        }
        if (Common::StringOps::from_str<CFuint>(words[1]) != m_nbrEqs)
        {
          throw BadValueException (FromHere(),"Wrong number of equations in given CFmesh file: " + words[1]);
        }
      }
      else if (words[0] == "!NB_ELEM")
      {
        if (words.size() != 2)
        {
          throw FileFormatException (FromHere(),"Cannot read number of elements...");
        }
        if (Common::StringOps::from_str<CFuint>(words[1]) != (*MeshDataStack::getActive()->getElementTypeData())[0].getNbElems())
        {
          throw BadValueException (FromHere(),"Wrong number of elements in given CFmesh file: " + words[1]);
        }
      }
      else if (words[0] == "!NB_ELEM_TYPES")
      {
        if (words.size() != 2)
        {
          throw FileFormatException (FromHere(),"Cannot read number of element types...");
        }
        if (Common::StringOps::from_str<CFuint>(words[1]) != 1)
        {
          throw BadValueException (FromHere(),"Wrong number of element types in given CFmesh file: " + words[1]);
        }
      }
      else if (words[0] == "!SOL_POLYORDER")
      {
        if (words.size() != 2)
        {
          throw FileFormatException (FromHere(),"Cannot read solution polynomial order...");
        }
        m_srcSolPolyOrder = Common::StringOps::from_str<CFuint>(words[1]);
      }
      else if (words[0] == "!ELEM_TYPES")
      {
        if (words.size() != 2)
        {
          throw FileFormatException (FromHere(),"Cannot read element types...");
        }
        if (!(words[1] == "Quad" && m_dim == 2) && !(words[1] == "Hexa" && m_dim == 3))
        {
          throw BadValueException (FromHere(),"Unsupported element type in given CFmesh file: " + words[1]);
        }
      }
      else if (words[0] == "!LIST_STATE")
      {
        if (words.size() != 2)
        {
          throw FileFormatException (FromHere(),"Cannot read states...");
        }
        if (words[1] != "1")
        {
          throw BadValueException (FromHere(),"States data is not present in given CFmesh file");
        }
      }
    }
  } while (line != "!END");

  // reset get pointer to beginning of file
  inputFile.seekg(0, ios::beg);
}

//////////////////////////////////////////////////////////////////////////////

void StdInitStateFromLowerOrderState::computeSolutionTransformationMatrix()
{
  // create SpectralFDElementData object corresponding to solution polynomial order in CFmesh file
  if (m_dim == 2)
  {
    m_srcSolSDData = new QuadSpectralFDElementData(static_cast<CFPolyOrder::Type>(m_srcSolPolyOrder));
  }
  if (m_dim == 3)
  {
    m_srcSolSDData = new HexaSpectralFDElementData(static_cast<CFPolyOrder::Type>(m_srcSolPolyOrder));
  }

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // number of states in input and current cells
  const CFuint nbrCellStatesSrc = m_srcSolSDData->getNbrOfSolPnts();
  const CFuint nbrCellStatesCur = sdLocalData[0]  ->getNbrOfSolPnts();

  // compute solution transformation matrix
  m_solTransMatr.resize(nbrCellStatesCur,nbrCellStatesSrc);
  if (nbrCellStatesSrc <= nbrCellStatesCur) // input is lower order
  {
    // get high-order solution point mapped coordinates
    SafePtr< vector< RealVector > > highOrderSolPntsCoords = sdLocalData[0]->getSolPntsLocalCoords();

    // evaluate low-order solution polynomials at high-order solution points
    vector< vector< CFreal > > solPolyVals = m_srcSolSDData->getSolPolyValsAtNode(*highOrderSolPntsCoords);

    // copy to matrix
    for (CFuint iSrc = 0; iSrc < nbrCellStatesSrc; ++iSrc)
    {
      for (CFuint iCur = 0; iCur < nbrCellStatesCur; ++iCur)
      {
        m_solTransMatr(iCur,iSrc) = solPolyVals[iCur][iSrc];
      }
    }
  }
  else
  {
    throw Common::NotImplementedException (FromHere(),"Initialization of low order solution from higher order solution is not yet implemented (should compute transformation matrix using a least squares approach)");
  }
  
  // resize m_cellSrcStates 
  // AL: here there was a potential bug
  m_cellSrcStates.resize(nbrCellStatesSrc,RealVector(0., m_nbrEqs)); 
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

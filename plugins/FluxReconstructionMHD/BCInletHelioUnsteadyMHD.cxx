#include "Framework/MethodStrategyProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Common/OldLookupTable.hh"
#include "Common/NotImplementedException.hh"


#include "MHD/MHD3DProjectionVarSet.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/BCInletHelioUnsteadyMHD.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include <fstream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <regex>
#include <algorithm>

// Filesystem includes will be handled via Environment/DirPaths which already includes boost filesystem

#include "Framework/MeshData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/TrsNotFoundException.hh"
#include "MathTools/MathChecks.hh"
#include "Common/SwapEmpty.hh"
#include "Common/FilesystemException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCInletHelioUnsteadyMHD,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMHDModule >
  BCInletHelioMHDProvider("SuperInletHelioMHD");

Framework::MethodStrategyProvider<
    BCInletHelioUnsteadyMHD,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMHDModule >
  BCInletHelioUnsteadyMHDProvider("SuperInletHelioUnsteadyMHD");

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<std::string> >
   ("TRSName","Name of the TRSs on which values must be prescribed");

  options.addConfigOption<CFreal>
    ("RotationAngle", "Rotation angle (in degrees) to apply to the wall distribution coordinates");

  options.addConfigOption<std::vector<CFuint> >
    ("RotationCoordIDs", "IDs (must be 2) of the coordinates lying in the rotation plane");
  
  options.addConfigOption<CFuint>
    ("NbClosestPoints", "Number of closest points for surface interpolation"); 
  
  options.addConfigOption<vector<CFint> >
    ("ExtractCoordXYIDs", "IDs corresponding to the x,y coordinate (z=0) for which plane is extracted");
  
  options.addConfigOption< std::vector<CFuint> >
    ("InitialSolutionIDs", "IDs of initial solution components that will be used as BC value.");

  options.addConfigOption< vector<std::string> >
    ("FileNameTw", "Names of the files with the given boundary distribution (multiple files for time interpolation)");

  options.addConfigOption< vector<CFreal> >
    ("FileNameTime", "Time corresponding to each file with the given boundary distribution");
    
  // Automatic file generation options
  options.addConfigOption<std::string>
    ("FilePattern", "Pattern for automatic file generation (e.g., 'corona_%04d_%02dh%02dm%02ds_surface.dat'). Use with AutoGenerate=true.");
    
  options.addConfigOption<bool>
    ("AutoGenerate", "Enable automatic file generation based on pattern and time range");
    
  options.addConfigOption<std::string>
    ("StartDateTime", "Starting date-time in format 'YYYY-MM-DD_HH:MM:SS' (e.g., '2019-07-04_02:32:54')");
    
  options.addConfigOption<std::string>
    ("EndDateTime", "Ending date-time in format 'YYYY-MM-DD_HH:MM:SS' (e.g., '2019-07-04_08:32:54')");
    
  options.addConfigOption<CFreal>
    ("DataTimeStep", "Time step between files in non-dimensional units (user converts from physical time)");
    
  options.addConfigOption<CFreal>
    ("StartTime", "Starting simulation time (non-dimensional) corresponding to StartDateTime");
    
  options.addConfigOption<CFreal>
    ("EndTime", "Ending simulation time (non-dimensional) corresponding to EndDateTime");  
}

//////////////////////////////////////////////////////////////////////////////

BCInletHelioUnsteadyMHD::BCInletHelioUnsteadyMHD(const std::string& name) :
  BCStateComputer(name),
  m_extractCoordZID(-1),
  m_flxPntsLocalCoords(),
  m_thisTRS(),
  m_flxLocalCoords(), 
  m_flxPntCoords(),
  m_nbrFaceFlxPnts(),
  m_nbrFaceFlxPntsMax(),
  m_dim(),      
  m_orient(),         
  _faceBuilder(),          
  m_globalToLocalTRSFaceID(),
  m_flxPntRhos(),
  m_flxPntUs(),
  m_flxPntVs(),
  m_flxPntWs(),
  m_flxPntBxs(),
  m_flxPntBys(),
  m_flxPntBzs(),
  m_flxPntPs(),
  m_nbGeoEnts(),
  m_nbEqs(),
  m_initialSolutionMap(),
  m_intCell(CFNULL),
  m_cellStates(CFNULL),
  m_currFace(CFNULL),
  m_cellStatesFlxPnt(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_faceFlxPntConn(CFNULL),
  m_fileNameTw(),
  m_fileNameTime(),
  m_allSurfaces(),
  m_surfaceAtTime(),
  m_useTimeInterpolation(false),
  m_filePattern(""),
  m_autoGenerate(false),
  m_startDateTime(""),
  m_endDateTime(""),
  m_dataTimeStep(1.0),
  m_startTime(0.0),
  m_endTime(100.0),
  m_lastInterpolationTime(-1.0)
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_trsName = vector<std::string>();
  this->setParameter("TRSName",&m_trsName);

  m_angle = 0.0;
  this->setParameter("RotationAngle",&m_angle);

  m_xvec = std::vector<CFuint>();
  this->setParameter("RotationCoordIDs",&m_xvec);
  
  m_nbClosestPoints = 0;
  this->setParameter("NbClosestPoints",&m_nbClosestPoints);
  
  m_extractCoordXYID = vector<CFint>();
  this->setParameter("ExtractCoordXYIDs",&m_extractCoordXYID);

  m_initialSolutionIDs = std::vector<CFuint>();
  setParameter("InitialSolutionIDs",&m_initialSolutionIDs);

  m_fileNameTw = vector<std::string>();
  this->setParameter("FileNameTw",&m_fileNameTw);

  m_fileNameTime = vector<CFreal>();
  this->setParameter("FileNameTime",&m_fileNameTime);

  // Automatic file generation parameters
  m_filePattern = "";
  this->setParameter("FilePattern",&m_filePattern);
  
  m_autoGenerate = false;
  this->setParameter("AutoGenerate",&m_autoGenerate);
  
  m_startDateTime = "";
  this->setParameter("StartDateTime",&m_startDateTime);
  
  m_endDateTime = "";
  this->setParameter("EndDateTime",&m_endDateTime);
  
  m_dataTimeStep = 1.0; // Default 1 non-dimensional time unit
  this->setParameter("DataTimeStep",&m_dataTimeStep);
  
  m_startTime = 0.0;
  this->setParameter("StartTime",&m_startTime);
  
  m_endTime = 100.0;
  this->setParameter("EndTime",&m_endTime);

  m_useTimeInterpolation = false;
}

//////////////////////////////////////////////////////////////////////////////

BCInletHelioUnsteadyMHD::~BCInletHelioUnsteadyMHD()
{
  CFAUTOTRACE;
}


//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::preProcess()
{  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  if (m_useTimeInterpolation) {
    // Update boundary data through time interpolation
    interpolateSurfaceDataInTime();
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  // Dynamic time interpolation for time-evolving boundary conditions
  if (m_useTimeInterpolation) {
    const CFreal currTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
    // Only update if time has changed (for efficiency)
    if (std::abs(currTime - m_lastInterpolationTime) > 1e-12) {
      interpolateSurfaceDataAtTime(currTime);
      m_lastInterpolationTime = currTime;
    }
  }
  
  // number of states
  //const CFuint nbrStates = ghostStates.size();
  //cf_assert(nbrStates == intStates.size());
  //cf_assert(nbrStates == normals.size());
  
  CFuint nbrStates = coords.size();//m_nbrFaceFlxPnts;  
  
  // Get the localFaceID from the map, knowing the faceGlobalID
  const CFuint faceLocalID = m_globalToLocalTRSFaceID.find(m_face->getID());
  
  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {              
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);


    // We compute the local index of the flux point
    const CFuint localFlxPntID = faceLocalID*m_nbrFaceFlxPntsMax + iState; 

    const CFreal xI_dimless = coords[iState][XX];
    const CFreal yI_dimless = coords[iState][YY];
    const CFreal zI_dimless = coords[iState][ZZ];
    const CFreal rI_dimless = (coords[iState]).norm2();
    const CFreal rhoI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless);
    const CFreal thetaBoundary_dimless = std::atan2(std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless),zI_dimless);

    CFreal rhoBoundary = m_flxPntRhos[localFlxPntID];
    CFreal uBoundary = m_flxPntUs[localFlxPntID];
    CFreal vBoundary = m_flxPntVs[localFlxPntID];
    CFreal wBoundary = m_flxPntWs[localFlxPntID];
    CFreal BxBoundary = m_flxPntBxs[localFlxPntID];
    CFreal ByBoundary = m_flxPntBys[localFlxPntID];
    CFreal BzBoundary = m_flxPntBzs[localFlxPntID];
    CFreal pBoundary = m_flxPntPs[localFlxPntID];


    const CFreal BxI_dimless = intState[4];
    const CFreal ByI_dimless = intState[5];
    const CFreal BzI_dimless = intState[6];


    ghostState[4] = 2.0*BxBoundary - BxI_dimless;
    ghostState[5] = 2.0*ByBoundary - ByI_dimless;
    ghostState[6] = 2.0*BzBoundary - BzI_dimless;

    ghostState[1] = 2.0*uBoundary - intState[1];
    ghostState[2] = 2.0*vBoundary - intState[2];
    ghostState[3] = 2.0*wBoundary - intState[3];

    //rho
    ghostState[0] = 2.0*rhoBoundary - intState[0];
    //p
    ghostState[7] = 2.0*pBoundary - intState[7];
    //phi
    ghostState[8] = intState[8];
    
       
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  //const CFuint nbrStateGrads = intGrads.size();
  //cf_assert(nbrStateGrads == ghostGrads.size());
  //cf_assert(nbrStateGrads == normals.size()); 
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::generateFileListFromPattern()
{
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Framework;
  
  CFLog(INFO, "BCInletHelioUnsteadyMHD::generateFileListFromPattern() => START\n");
  
  if (!m_autoGenerate || m_filePattern.empty()) {
    return;
  }
  
  // Extract directory from pattern (everything before the last '/')
  std::string directory = ".";  // Default to current directory
  std::string filenamePattern = m_filePattern;
  
  size_t lastSlash = m_filePattern.find_last_of("/\\");
  if (lastSlash != std::string::npos) {
    directory = m_filePattern.substr(0, lastSlash);
    filenamePattern = m_filePattern.substr(lastSlash + 1);
  }
  
  CFLog(INFO, "Scanning directory: " << directory << "\n");
  CFLog(INFO, "Looking for pattern: " << filenamePattern << "\n");
  
  // Helper function to parse date-time from filename
  std::vector<int> (*extractDateTimeFromFilename)(const std::string&) = 
    [](const std::string& filename) -> std::vector<int> {
    // Try pattern 1: YYYY-MM-DD_HHhMMmSSs (full format)
    std::regex dateTimeRegex1(R"((\d{4})-(\d{2})-(\d{2})_(\d{2})h(\d{2})m(\d{2})s)");
    std::smatch matches;
    
    if (std::regex_search(filename, matches, dateTimeRegex1)) {
      return {
        std::stoi(matches[1].str()), // year
        std::stoi(matches[2].str()), // month  
        std::stoi(matches[3].str()), // day
        std::stoi(matches[4].str()), // hour
        std::stoi(matches[5].str()), // minute
        std::stoi(matches[6].str())  // second
      };
    }
    
    // Try pattern 2: YYYY-MM-DD_HHh (hour only format)
    std::regex dateTimeRegex2(R"((\d{4})-(\d{2})-(\d{2})_(\d{2})h)");
    if (std::regex_search(filename, matches, dateTimeRegex2)) {
      return {
        std::stoi(matches[1].str()), // year
        std::stoi(matches[2].str()), // month  
        std::stoi(matches[3].str()), // day
        std::stoi(matches[4].str()), // hour
        0,                           // minute (default to 0)
        0                            // second (default to 0)
      };
    }
    
    return {}; // Empty if no match
  };
  
  // Helper function to convert datetime to seconds since epoch (simplified)
  long long (*dateTimeToSeconds)(int, int, int, int, int, int) = 
    [](int year, int month, int day, int hour, int min, int sec) -> long long {
    // Simple conversion - assumes each month has 30 days (good enough for time differences)
    return sec + min*60 + hour*3600 + day*86400LL + month*2592000LL + year*31536000LL;
  };
  
  // Helper function to parse date-time string
  std::vector<int> (*parseDateTime)(const std::string&) = 
    [](const std::string& dateTimeStr) -> std::vector<int> {
    // Parse "YYYY-MM-DD_HH:MM:SS"
    std::stringstream ss(dateTimeStr);
    std::string datePart, timePart;
    
    std::getline(ss, datePart, '_');
    std::getline(ss, timePart);
    
    // Parse date part "YYYY-MM-DD"
    std::stringstream dateStream(datePart);
    std::string yearStr, monthStr, dayStr;
    std::getline(dateStream, yearStr, '-');
    std::getline(dateStream, monthStr, '-');
    std::getline(dateStream, dayStr);
    
    // Parse time part "HH:MM:SS"
    std::stringstream timeStream(timePart);
    std::string hourStr, minStr, secStr;
    std::getline(timeStream, hourStr, ':');
    std::getline(timeStream, minStr, ':');
    std::getline(timeStream, secStr);
    
    return {std::stoi(yearStr), std::stoi(monthStr), std::stoi(dayStr),
            std::stoi(hourStr), std::stoi(minStr), std::stoi(secStr)};
  };
  
  // Get start and end time bounds
  std::vector<int> startParts = parseDateTime(m_startDateTime);
  std::vector<int> endParts = parseDateTime(m_endDateTime);
  
  long long startSeconds = dateTimeToSeconds(startParts[0], startParts[1], startParts[2], 
                                           startParts[3], startParts[4], startParts[5]);
  long long endSeconds = dateTimeToSeconds(endParts[0], endParts[1], endParts[2],
                                         endParts[3], endParts[4], endParts[5]);
  
  CFLog(INFO, "Time range: " << m_startDateTime << " to " << m_endDateTime << "\n");
  
  // Use Environment::DirPaths to work with directories (already has boost filesystem)
  Environment::DirPaths& dirPaths = Environment::DirPaths::getInstance();
  boost::filesystem::path dirPath;
  
  if (directory == ".") {
    dirPath = dirPaths.getWorkingDir();
  } else {
    dirPath = boost::filesystem::path(directory);
    if (!dirPath.is_absolute()) {
      dirPath = dirPaths.getWorkingDir() / dirPath;
    }
  }
  
  // Convert pattern to regex (replace *** with datetime pattern)  
  std::string regexPattern = filenamePattern;
  size_t pos = regexPattern.find("***");
  if (pos != std::string::npos) {
    // Support both formats: YYYY-MM-DD_HHhMMmSSs OR YYYY-MM-DD_HHh
    regexPattern.replace(pos, 3, R"(\d{4}-\d{2}-\d{2}_\d{2}h(?:\d{2}m\d{2}s)?)");
  }
  
  // Scan directory for matching files
  std::vector<std::pair<std::string, long long>> validFiles; // filename, seconds
  
  try {
    if (!boost::filesystem::exists(dirPath) || !boost::filesystem::is_directory(dirPath)) {
      throw Common::FilesystemException(FromHere(), "Directory not found: " + dirPath.string());
    }
    
    std::regex fileRegex(regexPattern);
    boost::filesystem::directory_iterator endIter;
    
    for (boost::filesystem::directory_iterator iter(dirPath); iter != endIter; ++iter) {
      if (boost::filesystem::is_regular_file(iter->status())) {
        std::string filename = iter->path().filename().string();
        
        // Check if filename matches our pattern
        if (std::regex_match(filename, fileRegex)) {
          // Extract datetime from filename
          std::vector<int> fileParts = extractDateTimeFromFilename(filename);
          if (!fileParts.empty()) {
            long long fileSeconds = dateTimeToSeconds(fileParts[0], fileParts[1], fileParts[2],
                                                    fileParts[3], fileParts[4], fileParts[5]);              // Check if file is within our time range
            if (fileSeconds >= startSeconds && fileSeconds <= endSeconds) {
              // Create path relative to base directory for compatibility with readSurfaceData
              boost::filesystem::path basePath = Environment::DirPaths::getInstance().getBaseDir();
              boost::filesystem::path fullPath = iter->path();
              std::string relativePath;
              
              // Convert to relative path for readSurfaceData compatibility
              std::string baseStr = basePath.string();
              std::string fullStr = fullPath.string();
              
              if (fullStr.find(baseStr) == 0) {
                // File is under base directory, make it relative
                relativePath = fullStr.substr(baseStr.length());
                // Remove leading slash if present
                if (!relativePath.empty() && (relativePath[0] == '/' || relativePath[0] == '\\')) {
                  relativePath = relativePath.substr(1);
                }
              } else {
                // File is outside base directory, use full absolute path for readSurfaceData
                relativePath = fullStr;
              }
              
              validFiles.push_back({relativePath, fileSeconds});
              CFLog(VERBOSE, "Found matching file: " << filename << " at time " << fileSeconds << "\n");
            }
          }
        }
      }
    }
    
  } catch (const boost::filesystem::filesystem_error& ex) {
    throw Common::FilesystemException(FromHere(), 
      "Error scanning directory " + dirPath.string() + ": " + ex.what());
  }
  
  // Sort files by time using explicit comparison function
  struct CompareFilesByTime {
    bool operator()(const std::pair<std::string, long long>& a, 
                   const std::pair<std::string, long long>& b) const {
      return a.second < b.second;
    }
  };
  
  std::sort(validFiles.begin(), validFiles.end(), CompareFilesByTime());
  
  // Generate file list and time mapping
  m_fileNameTw.clear();
  m_fileNameTime.clear();
  
  if (validFiles.empty()) {
    throw Common::BadValueException(FromHere(), 
      "No files found matching pattern " + filenamePattern + " in directory " + dirPath.string());
  }
  
  // Create time mapping based on DataTimeStep
  for (size_t i = 0; i < validFiles.size(); ++i) {
    m_fileNameTw.push_back(validFiles[i].first);
    
    // Calculate non-dimensional time
    CFreal simTime = m_startTime + i * m_dataTimeStep;
    m_fileNameTime.push_back(simTime);
    
    CFLog(VERBOSE, "File " << i << ": " << validFiles[i].first 
                  << " -> sim time " << simTime << "\n");
  }
  
  CFLog(INFO, "Found and sorted " << validFiles.size() << " files\n");
  CFLog(INFO, "Time range: " << m_startTime << " to " << (m_startTime + (validFiles.size()-1) * m_dataTimeStep) << " (non-dimensional)\n");
  CFLog(INFO, "BCInletHelioUnsteadyMHD::generateFileListFromPattern() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();
  m_needsSpatCoord = true;
  
  // Generate file list automatically if enabled
  generateFileListFromPattern();
  
  m_dim = PhysicalModelStack::getActive()->getDim();
  if (m_nbClosestPoints == 0) {
    m_nbClosestPoints = m_dim + (m_dim - 2); // 2 in 2D, 4 in 3D
  }
  
  if (m_extractCoordXYID.size() > 0) 
  {
    if ((m_extractCoordXYID[0] == 0 && m_extractCoordXYID[1] == 1) ||
	  (m_extractCoordXYID[0] == 1 && m_extractCoordXYID[1] == 0)) 
    {
      m_extractCoordZID = 2;
    }
      
    if ((m_extractCoordXYID[0] == 0 && m_extractCoordXYID[1] == 2) ||
	  (m_extractCoordXYID[0] == 2 && m_extractCoordXYID[1] == 0)) 
    {
      m_extractCoordZID = 1;
    } 
      
    if ((m_extractCoordXYID[0] == 1 && m_extractCoordXYID[1] == 2) ||
	  (m_extractCoordXYID[0] == 2 && m_extractCoordXYID[1] == 1)) 
    {
      m_extractCoordZID = 0;
    }
  }
    
  if (std::abs(m_angle) > 0.0) 
  {
    // conversion to radiants
    m_angle *= MathConsts::CFrealPi()/180.;
    cf_always_assert(m_xvec.size() == 2);
  }
  
  vector<SurfaceData*> surfaces;
  
  // Check if we have multiple files for time interpolation
  if (m_fileNameTw.size() > 1) {
    m_useTimeInterpolation = true;
    
    // Set default times if not provided
    if (m_fileNameTime.size() == 0) {
      m_fileNameTime.resize(m_fileNameTw.size());
      for (CFuint i = 0; i < m_fileNameTw.size(); ++i) {
        m_fileNameTime[i] = static_cast<CFreal>(i);
      }
    }
    
    cf_assert(m_fileNameTw.size() == m_fileNameTime.size());
    
    // Read all surface files for time interpolation
    m_allSurfaces.resize(m_fileNameTw.size());
    for (CFuint i = 0; i < m_fileNameTw.size(); ++i) {
      CFLog(INFO, "BCInletHelioUnsteadyMHD: setup: readSurfaceData file " << i << " => START\n");
      readSurfaceData(m_allSurfaces[i], m_fileNameTw[i]);
      CFLog(INFO, "BCInletHelioUnsteadyMHD: setup: readSurfaceData file " << i << " => END\n");
    }
    
    // Initialize current surface data structure
    m_surfaceAtTime.resize(m_allSurfaces[0].size());
    for (CFuint is = 0; is < m_allSurfaces[0].size(); ++is) {
      m_surfaceAtTime[is] = new SurfaceData();
      // Copy structure from first surface
      const SurfaceData& refSurf = *m_allSurfaces[0][is];
      m_surfaceAtTime[is]->xyz = refSurf.xyz;
      m_surfaceAtTime[is]->rho.resize(refSurf.rho.size());
      m_surfaceAtTime[is]->u.resize(refSurf.u.size());
      m_surfaceAtTime[is]->v.resize(refSurf.v.size());
      m_surfaceAtTime[is]->w.resize(refSurf.w.size());
      m_surfaceAtTime[is]->Bx.resize(refSurf.Bx.size());
      m_surfaceAtTime[is]->By.resize(refSurf.By.size());
      m_surfaceAtTime[is]->Bz.resize(refSurf.Bz.size());
      m_surfaceAtTime[is]->p.resize(refSurf.p.size());
    }
    
    // Copy first file data to current surface
    surfaces = m_surfaceAtTime;
  } else if (m_fileNameTw.size() == 1) {
    // Single file case
    m_useTimeInterpolation = false;
    CFLog(INFO, "BCInletHelioUnsteadyMHD: setup: readSurfaceData => START\n");
    readSurfaceData(surfaces, m_fileNameTw[0]);
    CFLog(INFO, "BCInletHelioUnsteadyMHD: setup: readSurfaceData => END\n");
  } else {
    throw Common::BadValueException(FromHere(), "No files specified in FileNameTw");
  }
  
  const CFuint nbSurf = surfaces.size();
  cf_assert(nbSurf >= 1);
  
  RealVector tmpNode(m_dim);
  
  // Get face builders from method data
  _faceBuilder = getMethodData().getSecondFaceBuilder();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);

  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  SafePtr< std::vector< std::vector< RealVector > > >              faceFlxPntsLocalCoordsPerType = frLocalData[0]->getFaceFlxPntsLocalCoordsPerType();
  m_nbrFaceFlxPnts = m_flxLocalCoords->size();
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();
  
  // Set max number of face flux points
  if (elemShape == CFGeoShape::PRISM)  // (Max number of face flx pnts)
  {
    m_nbrFaceFlxPntsMax = (order+1)*(order+1);
  }
  else
  {
    m_nbrFaceFlxPntsMax = m_nbrFaceFlxPnts;
  }
  
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();

  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsMax; ++iFlx)
  {
    m_cellStatesFlxPnt.push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsMax; ++iFlx)
  {
    m_cellStatesFlxPnt[iFlx]->setLocalID(iFlx);
  }

  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
  }
  
  m_flxPntsLocalCoords.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
    m_flxPntsLocalCoords[iFlx].resize(m_dim);
  }
  
  // get TRS list
  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();
  
  const CFuint nbTRSs = trsList.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    CFLog(INFO, "iTRS "<< iTRS << " \n");
    if (m_trsName[0]==trsList[iTRS]->getName() ){
       m_thisTRS = trsList[iTRS];
       CFLog(INFO, "Matching BC "<<m_trsName[0]<<" with "<<m_thisTRS->getName() << "\n");
    }
  }
  
  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = m_thisTRS;
  faceData.isBoundary = true;

  // We get the number of faces in this TRS
  m_nbGeoEnts = m_thisTRS->getLocalNbGeoEnts();  //getNbTRs();

  // We get the number of equations
  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  std::vector<FlxPntStruct> bndFlxPnts;
  bndFlxPnts.reserve(m_nbGeoEnts*m_nbrFaceFlxPnts); // oversized to avoid frequent memory reallocation

  FlxPntStruct thisFlxPnt;
  
  createFaceOrientationStartIndexes();
  
  map< std::string , vector< vector< CFuint > > >&
      bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
    vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[m_thisTRS->getName()];
    
    cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1; 

  const CFuint nbTRs = m_thisTRS->getNbTRs();

  // loop over TRs
  CFuint localFaceID = 0;
  bool isFirstOrientation = true;
  
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    for (m_orient = 0; m_orient < nbOrients; ++m_orient){
       // start and stop index of the faces with this orientation
       const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
       const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

       // loop over faces with this orientation
       for (CFuint iFace = startFaceIdx; iFace < stopFaceIdx; ++iFace){

          /* Face setup */
          faceData.idx = iFace;
          GeometricEntity *const face = _faceBuilder->buildGE();
          const CFuint faceGlobalID = face->getID();
          m_globalToLocalTRSFaceID.insert(faceGlobalID,localFaceID);
          
          if (elemShape == CFGeoShape::PRISM)
          {
            const CFGeoShape::Type geo = face->getShape();

            if (geo == CFGeoShape::TRIAG) // triag face
            {
              m_nbrFaceFlxPnts=((*faceFlxPntsLocalCoordsPerType)[0]).size();
              (*m_flxLocalCoords)=((*faceFlxPntsLocalCoordsPerType)[0]);
            }
            else  // quad face
            {
              m_nbrFaceFlxPnts=((*faceFlxPntsLocalCoordsPerType)[1]).size();
              (*m_flxLocalCoords)=((*faceFlxPntsLocalCoordsPerType)[1]);
            }
          }

          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
            // Compute coordinates
             m_flxPntCoords[iFlx] = face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
          }

          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
             /* Load this flux point in a structure */
             thisFlxPnt.setCentreCoordinates(m_flxPntCoords[iFlx]);
             thisFlxPnt.setGlobalFaceID(faceGlobalID);
             thisFlxPnt.setLocalFaceID(localFaceID);
             thisFlxPnt.setLocalFluxID(iFlx);
             thisFlxPnt.setOrient(m_orient);
             
             bndFlxPnts.push_back(thisFlxPnt);
             
             /* release geometry */
             _faceBuilder->releaseGE();
          } // End of the flux points loop
          localFaceID++;
          
       } // End of the faces loop
       
       if (startFaceIdx!=stopFaceIdx){
          isFirstOrientation=false;
       }
    } // End of the orientation loop
  } // End of the TR loop
  
  CFuint nbBndFlxPnts = bndFlxPnts.size();
  CFuint nbBndFlxPntsMax = m_nbrFaceFlxPntsMax*localFaceID;

  ClosestPointData closestPoint;
  closestPoint.surfaceIDs.resize(m_nbClosestPoints); 
  closestPoint.pointsIDs.resize(m_nbClosestPoints); 
  closestPoint.r.resize(m_nbClosestPoints); 
  
  m_flxPntRhos.resize(nbBndFlxPntsMax);
  m_flxPntUs.resize(nbBndFlxPntsMax);
  m_flxPntVs.resize(nbBndFlxPntsMax);
  m_flxPntWs.resize(nbBndFlxPntsMax);
  m_flxPntBxs.resize(nbBndFlxPntsMax);
  m_flxPntBys.resize(nbBndFlxPntsMax);
  m_flxPntBzs.resize(nbBndFlxPntsMax);
  m_flxPntPs.resize(nbBndFlxPntsMax);
  
  // Perform initial spatial extrapolation
  if (m_useTimeInterpolation) {
    // For time interpolation, copy first file to current surface and extrapolate
    for (CFuint is = 0; is < m_allSurfaces[0].size(); ++is) {
      const SurfaceData& refSurf = *m_allSurfaces[0][is];
      SurfaceData& currentSf = *m_surfaceAtTime[is];
      
      for (CFuint ip = 0; ip < refSurf.rho.size(); ++ip) {
        currentSf.rho[ip] = refSurf.rho[ip];
        currentSf.u[ip] = refSurf.u[ip];
        currentSf.v[ip] = refSurf.v[ip];
        currentSf.w[ip] = refSurf.w[ip];
        currentSf.Bx[ip] = refSurf.Bx[ip];
        currentSf.By[ip] = refSurf.By[ip];
        currentSf.Bz[ip] = refSurf.Bz[ip];
        currentSf.p[ip] = refSurf.p[ip];
      }
    }
    extrapolateSpatialData();
  } else {
    // Single file case - perform traditional spatial interpolation
    performInitialSpatialInterpolation(surfaces, bndFlxPnts);
  }
  
  // cleanup the memory for single file case
  if (!m_useTimeInterpolation) {
    for (CFuint is = 0; is < nbSurf; ++is) {
      deletePtr(surfaces[is]);
    }
  }
  
  
  m_tempStates.resize(m_nbrFaceFlxPntsMax);
  
  // number of equations
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  //for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsMax; ++iFlx)
  {
    m_tempStates[iFlx].resize(nbEqs); 
  }
}

//////////////////////////////////////////////////////////////////////////////


void BCInletHelioUnsteadyMHD::createFaceOrientationStartIndexes()
{
  CFAUTOTRACE;
   //MPI_Barrier(_comm);//the processing of a an individual processor is paused temporarily till all the processors compute till some certain values

   //MPI_Allgather(&m_nbGeoEnts, 1, MPIStructDef::getMPIType(&m_nbGeoEnts), 
	//&_nbFacesPerProcess[0], 1, MPIStructDef::getMPIType(&m_nbGeoEnts), _comm);

  CFLog(INFO, " BCInletHelioUnsteadyMHD::createFaceOrientationStartIndexes()\n");

  // START INDEXES FOR INNER CFGeoEnt::CELLS
  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get connectivity where the face start indexes are stored
  SafePtr< ConnectivityTable< CFuint > > innerFacesStartIdxsConnTable
    = MeshDataStack::getActive()->getConnectivity("innerFacesStartIdxs");

  // get number of orientations
  const CFuint nbrFaceOrientsPlus1 = innerFacesStartIdxsConnTable->nbRows();

  // resize innerFacesStartIdxs
  innerFacesStartIdxs.resize(nbrFaceOrientsPlus1);

  // put indexes in innerFacesStartIdxs
  for (CFuint iIdx = 0; iIdx < nbrFaceOrientsPlus1; ++iIdx)
  {
    innerFacesStartIdxs[iIdx] = (*innerFacesStartIdxsConnTable)(iIdx,0);
  }

  // START INDEXES FOR BOUNDARY CFGeoEnt::CELLS
  // get bndFacesStartIdxs from FluxReconstructionSolverData
  map< std::string , vector< vector< CFuint > > >& bndFacesStartIdxs = getMethodData().getBndFacesStartIdxs();

  // get TRS list
  vector<Common::SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  // number of TRSs
  const CFuint nbTRSs = trsList.size();

  // loop over boundary TRSs
  CFuint iBndTRS = 0;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    const std::string trsName = trsList[iTRS]->getName();

    if (trsList[iTRS]->hasTag("boundary"))
    {

      // get connectivity where the face start indexes are stored
      SafePtr< ConnectivityTable< CFuint > > boundaryFacesStartIdxsConnTable
        = MeshDataStack::getActive()->getConnectivity(trsName+"boundaryFacesStartIdxs");

      // get number of TRs in this TRS
      const CFuint nbTRs = trsList[iTRS]->getNbTRs();
      // number of boundary face orientations + 1
      const CFuint nbBndFaceOrientP1 = boundaryFacesStartIdxsConnTable->nbRows()/nbTRs;

      // array over TRs and startIdxs
      vector< vector< CFuint > > trStartIdxs(nbTRs);

      // loop over TRs
      CFuint iStartIdx = 0;
      for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
      {
        // resize trStartIdxs[iTR]
        trStartIdxs[iTR].resize(nbBndFaceOrientP1);

        // loop over start indexes
        for (CFuint iIdx = 0; iIdx < nbBndFaceOrientP1; ++iIdx, ++iStartIdx)
        {
          trStartIdxs[iTR][iIdx] = (*boundaryFacesStartIdxsConnTable)(iStartIdx,0);
        }
      }
      cf_assert(iStartIdx == boundaryFacesStartIdxsConnTable->nbRows());

      // store trStartIdxs in boundary TRS map
      bndFacesStartIdxs[trsName] = trStartIdxs;

      // increment boundary TRS counter
      ++iBndTRS;
    } 
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::interpolateSurfaceDataInTime()
{
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  if (!m_useTimeInterpolation) return;
  
  const CFreal currTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  interpolateSurfaceDataAtTime(currTime);
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::interpolateSurfaceDataAtTime(CFreal targetTime)
{
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  if (!m_useTimeInterpolation) return;
  
  // Add safety check for minimum number of files
  if (m_fileNameTime.size() == 0) {
    CFLog(ERROR, "No files available for time interpolation\n");
    return;
  }
  
  if (m_allSurfaces.size() == 0) {
    CFLog(ERROR, "No surface data loaded for time interpolation\n");
    return;
  }
  
  const CFreal maxTimeDim = SubSystemStatusStack::getActive()->getMaxTimeDim();
  const CFreal maxDT = SubSystemStatusStack::getActive()->getDT();
  
  CFLog(VERBOSE, "BCInletHelioUnsteadyMHD::interpolateSurfaceDataAtTime() => targetTime=" << targetTime << ", maxTime=" << maxTimeDim << "\n");
  CFLog(VERBOSE, "Available files: " << m_fileNameTime.size() << ", surfaces loaded: " << m_allSurfaces.size() << "\n");
  
  // time corresponding to the first input file used for the interpolation
  const CFreal startTime = m_fileNameTime[0];
  // time corresponding to the last input file used for the interpolation
  const CFreal endTime   = m_fileNameTime.back();
  
  CFLog(VERBOSE, "Time range: [" << startTime << ", " << endTime << "]\n");
  
  // Find the correct time interval for cubic Hermite interpolation
  CFuint idx_i = 0, idx_ip1 = 0;
  CFreal leftTime = m_fileNameTime[0];
  CFreal rightTime = m_fileNameTime.back();
  
  // Find the interval [i, i+1] containing targetTime
  bool intervalFound = false;
  for (CFuint i = 0; i < m_fileNameTime.size() - 1; ++i) {
    if ((targetTime - m_fileNameTime[i]) * (m_fileNameTime[i+1] - targetTime) >= 0.0) {
      idx_i = i;
      idx_ip1 = i + 1;
      leftTime = m_fileNameTime[i];
      rightTime = m_fileNameTime[i+1];
      intervalFound = true;
      break;
    }
  }
  
  // Handle boundary cases
  if (!intervalFound) {
    if (targetTime <= m_fileNameTime[0]) {
      idx_i = idx_ip1 = 0;
      leftTime = rightTime = m_fileNameTime[0];
      CFLog(VERBOSE, "Using first file (time <= startTime)\n");
    } else if (targetTime >= m_fileNameTime.back()) {
      idx_i = idx_ip1 = m_fileNameTime.size() - 1;
      leftTime = rightTime = m_fileNameTime.back();
      CFLog(VERBOSE, "Using last file (time >= endTime)\n");
    } else {
      CFLog(ERROR, "Failed to find time interval for targetTime=" << targetTime << "\n");
      return;
    }
  }
  
  CFLog(VERBOSE, "Using files idx_i=" << idx_i << " idx_ip1=" << idx_ip1 << " for targetTime=" << targetTime << "\n");
  
  // Additional safety check: ensure indices are valid
  if (idx_i >= m_allSurfaces.size() || idx_ip1 >= m_allSurfaces.size()) {
    CFLog(ERROR, "Invalid surface indices: idx_i=" << idx_i << " idx_ip1=" << idx_ip1 
                 << " size=" << m_allSurfaces.size() << "\n");
    return;
  }
  
  // Cubic Hermite interpolation coefficients
  CFreal Temp = rightTime - leftTime;
  CFreal Coefi1, Coefi2, time_norm;
  
  // Handle division by zero when times are equal
  if (std::abs(Temp) < 1e-12) {
    Coefi1 = 1.0;
    Coefi2 = 0.0;
    time_norm = 0.0;
  } else {
    Coefi1 = (rightTime - targetTime) / Temp;  // Linear coefficient for left point
    Coefi2 = (targetTime - leftTime) / Temp;   // Linear coefficient for right point
    time_norm = Coefi2;  // Normalized time parameter
  }
  
  // Hermite basis functions
  CFreal h00 = 2.0 * pow(time_norm, 3.0) - 3.0 * pow(time_norm, 2.0) + 1.0;
  CFreal h10 = pow(time_norm, 3.0) - 2.0 * pow(time_norm, 2.0) + time_norm;
  CFreal h01 = -2.0 * pow(time_norm, 3.0) + 3.0 * pow(time_norm, 2.0);
  CFreal h11 = pow(time_norm, 3.0) - pow(time_norm, 2.0);
  
  CFLog(VERBOSE, "Cubic Hermite interpolating between files " << idx_i << " (t=" << leftTime << ") and " 
                << idx_ip1 << " (t=" << rightTime << ") with time_norm=" << time_norm << "\n");
  CFLog(VERBOSE, "Hermite basis: h00=" << h00 << ", h10=" << h10 << ", h01=" << h01 << ", h11=" << h11 << "\n");
  
  // Get surface datasets for interpolation (need i-1, i, i+1, i+2 for cubic Hermite)
  const std::vector<SurfaceData*>& surfaces_i = m_allSurfaces[idx_i];
  const std::vector<SurfaceData*>& surfaces_ip1 = m_allSurfaces[idx_ip1];
  
  // Add safety checks for basic consistency
  if (surfaces_i.size() != surfaces_ip1.size()) {
    CFLog(ERROR, "Surface size mismatch: " << surfaces_i.size() << " vs " << surfaces_ip1.size() << "\n");
    return;
  }
  
  if (surfaces_i.size() != m_surfaceAtTime.size()) {
    CFLog(ERROR, "Surface size mismatch: " << surfaces_i.size() << " vs " << m_surfaceAtTime.size() << "\n");
    return;
  }
  
  cf_assert(surfaces_i.size() == surfaces_ip1.size());
  cf_assert(surfaces_i.size() == m_surfaceAtTime.size());
  
  // Handle case where we can do cubic Hermite with derivative estimation
  bool useCubicHermite = (idx_i > 0 && idx_ip1 < m_fileNameTime.size() - 1 && intervalFound);
  
  // Add safety checks for array bounds
  if (useCubicHermite) {
    // Verify array bounds before accessing
    if (idx_i - 1 >= m_allSurfaces.size() || idx_ip1 + 1 >= m_allSurfaces.size()) {
      CFLog(WARN, "Array bounds issue in cubic Hermite: idx_i-1=" << (idx_i-1) 
                  << ", idx_ip1+1=" << (idx_ip1+1) << ", size=" << m_allSurfaces.size() << "\n");
      useCubicHermite = false;
    }
  }
  
  if (useCubicHermite) {
    // Get additional datasets for derivative estimation (i-1 and i+2)
    const std::vector<SurfaceData*>& surfaces_im1 = m_allSurfaces[idx_i - 1];
    const std::vector<SurfaceData*>& surfaces_ip2 = m_allSurfaces[idx_ip1 + 1];
    
    CFreal Temp1 = m_fileNameTime[idx_i] - m_fileNameTime[idx_i - 1];
    CFreal Temp2 = m_fileNameTime[idx_ip1 + 1] - m_fileNameTime[idx_ip1];
    
    // Add safety checks for division by zero in derivative estimation
    if (std::abs(Temp1) < 1e-12 || std::abs(Temp2) < 1e-12) {
      CFLog(WARN, "Division by zero in derivative estimation: Temp1=" << Temp1 << ", Temp2=" << Temp2 << "\n");
      useCubicHermite = false;
    }
    
    CFLog(VERBOSE, "Using cubic Hermite with derivative estimation\n");
    
    // Update m_surfaceAtTime with cubic Hermite interpolation
    for (CFuint is = 0; is < surfaces_i.size(); ++is) {
      const SurfaceData& sf_im1 = *surfaces_im1[is];  // i-1
      const SurfaceData& sf_i = *surfaces_i[is];      // i
      const SurfaceData& sf_ip1 = *surfaces_ip1[is];  // i+1
      const SurfaceData& sf_ip2 = *surfaces_ip2[is];  // i+2
      SurfaceData& currentSf = *m_surfaceAtTime[is];
      
      cf_assert(sf_i.rho.size() == sf_ip1.rho.size());
      cf_assert(sf_i.rho.size() == currentSf.rho.size());
      
      // Cubic Hermite interpolation with derivative estimation for each variable at each point
      for (CFuint ip = 0; ip < sf_i.rho.size(); ++ip) {
        // Derivative estimation at point i (central difference where possible)
        CFreal deri1_rho = 0.5 * ((sf_i.rho[ip] - sf_im1.rho[ip]) / Temp1 + (sf_ip1.rho[ip] - sf_i.rho[ip]) / Temp);
        CFreal deri1_u = 0.5 * ((sf_i.u[ip] - sf_im1.u[ip]) / Temp1 + (sf_ip1.u[ip] - sf_i.u[ip]) / Temp);
        CFreal deri1_v = 0.5 * ((sf_i.v[ip] - sf_im1.v[ip]) / Temp1 + (sf_ip1.v[ip] - sf_i.v[ip]) / Temp);
        CFreal deri1_w = 0.5 * ((sf_i.w[ip] - sf_im1.w[ip]) / Temp1 + (sf_ip1.w[ip] - sf_i.w[ip]) / Temp);
        CFreal deri1_Bx = 0.5 * ((sf_i.Bx[ip] - sf_im1.Bx[ip]) / Temp1 + (sf_ip1.Bx[ip] - sf_i.Bx[ip]) / Temp);
        CFreal deri1_By = 0.5 * ((sf_i.By[ip] - sf_im1.By[ip]) / Temp1 + (sf_ip1.By[ip] - sf_i.By[ip]) / Temp);
        CFreal deri1_Bz = 0.5 * ((sf_i.Bz[ip] - sf_im1.Bz[ip]) / Temp1 + (sf_ip1.Bz[ip] - sf_i.Bz[ip]) / Temp);
        CFreal deri1_p = 0.5 * ((sf_i.p[ip] - sf_im1.p[ip]) / Temp1 + (sf_ip1.p[ip] - sf_i.p[ip]) / Temp);
        
        // Derivative estimation at point i+1
        CFreal deri2_rho = 0.5 * ((sf_ip2.rho[ip] - sf_ip1.rho[ip]) / Temp2 + (sf_ip1.rho[ip] - sf_i.rho[ip]) / Temp);
        CFreal deri2_u = 0.5 * ((sf_ip2.u[ip] - sf_ip1.u[ip]) / Temp2 + (sf_ip1.u[ip] - sf_i.u[ip]) / Temp);
        CFreal deri2_v = 0.5 * ((sf_ip2.v[ip] - sf_ip1.v[ip]) / Temp2 + (sf_ip1.v[ip] - sf_i.v[ip]) / Temp);
        CFreal deri2_w = 0.5 * ((sf_ip2.w[ip] - sf_ip1.w[ip]) / Temp2 + (sf_ip1.w[ip] - sf_i.w[ip]) / Temp);
        CFreal deri2_Bx = 0.5 * ((sf_ip2.Bx[ip] - sf_ip1.Bx[ip]) / Temp2 + (sf_ip1.Bx[ip] - sf_i.Bx[ip]) / Temp);
        CFreal deri2_By = 0.5 * ((sf_ip2.By[ip] - sf_ip1.By[ip]) / Temp2 + (sf_ip1.By[ip] - sf_i.By[ip]) / Temp);
        CFreal deri2_Bz = 0.5 * ((sf_ip2.Bz[ip] - sf_ip1.Bz[ip]) / Temp2 + (sf_ip1.Bz[ip] - sf_i.Bz[ip]) / Temp);
        CFreal deri2_p = 0.5 * ((sf_ip2.p[ip] - sf_ip1.p[ip]) / Temp2 + (sf_ip1.p[ip] - sf_i.p[ip]) / Temp);
        
        // Cubic Hermite interpolation formula: f(t) = f0*h00 + f0'*dt*h10 + f1*h01 + f1'*dt*h11
        currentSf.rho[ip] = sf_i.rho[ip] * h00 + deri1_rho * Temp * h10 + sf_ip1.rho[ip] * h01 + deri2_rho * Temp * h11;
        currentSf.u[ip] = sf_i.u[ip] * h00 + deri1_u * Temp * h10 + sf_ip1.u[ip] * h01 + deri2_u * Temp * h11;
        currentSf.v[ip] = sf_i.v[ip] * h00 + deri1_v * Temp * h10 + sf_ip1.v[ip] * h01 + deri2_v * Temp * h11;
        currentSf.w[ip] = sf_i.w[ip] * h00 + deri1_w * Temp * h10 + sf_ip1.w[ip] * h01 + deri2_w * Temp * h11;
        currentSf.Bx[ip] = sf_i.Bx[ip] * h00 + deri1_Bx * Temp * h10 + sf_ip1.Bx[ip] * h01 + deri2_Bx * Temp * h11;
        currentSf.By[ip] = sf_i.By[ip] * h00 + deri1_By * Temp * h10 + sf_ip1.By[ip] * h01 + deri2_By * Temp * h11;
        currentSf.Bz[ip] = sf_i.Bz[ip] * h00 + deri1_Bz * Temp * h10 + sf_ip1.Bz[ip] * h01 + deri2_Bz * Temp * h11;
        currentSf.p[ip] = sf_i.p[ip] * h00 + deri1_p * Temp * h10 + sf_ip1.p[ip] * h01 + deri2_p * Temp * h11;
      }
    }
  } else {
    // Fall back to linear interpolation when cubic Hermite is not possible
    CFLog(VERBOSE, "Using linear interpolation (cubic Hermite not available)\n");
    
    // Update m_surfaceAtTime with linear interpolation
    for (CFuint is = 0; is < surfaces_i.size(); ++is) {
      const SurfaceData& sf_i = *surfaces_i[is];
      const SurfaceData& sf_ip1 = *surfaces_ip1[is];
      SurfaceData& currentSf = *m_surfaceAtTime[is];
      
      cf_assert(sf_i.rho.size() == sf_ip1.rho.size());
      cf_assert(sf_i.rho.size() == currentSf.rho.size());
      
      // Linear interpolation for each variable at each point
      for (CFuint ip = 0; ip < sf_i.rho.size(); ++ip) {
        currentSf.rho[ip] = Coefi1 * sf_i.rho[ip] + Coefi2 * sf_ip1.rho[ip];
        currentSf.u[ip] = Coefi1 * sf_i.u[ip] + Coefi2 * sf_ip1.u[ip];
        currentSf.v[ip] = Coefi1 * sf_i.v[ip] + Coefi2 * sf_ip1.v[ip];
        currentSf.w[ip] = Coefi1 * sf_i.w[ip] + Coefi2 * sf_ip1.w[ip];
        currentSf.Bx[ip] = Coefi1 * sf_i.Bx[ip] + Coefi2 * sf_ip1.Bx[ip];
        currentSf.By[ip] = Coefi1 * sf_i.By[ip] + Coefi2 * sf_ip1.By[ip];
        currentSf.Bz[ip] = Coefi1 * sf_i.Bz[ip] + Coefi2 * sf_ip1.Bz[ip];
        currentSf.p[ip] = Coefi1 * sf_i.p[ip] + Coefi2 * sf_ip1.p[ip];
      }
    }
  }
  
  CFLog(VERBOSE, "Time interpolation completed (" << (useCubicHermite ? "cubic Hermite" : "linear") 
                << "), updating boundary flux point values\n");
  
  // After time interpolation, we need to update the boundary flux point arrays
  // that are used in computeGhostStates. The spatial mapping remains the same,
  // but we need to re-extract the interpolated values from m_surfaceAtTime.
  updateBoundaryFluxPointValues();
  
  CFLog(VERBOSE, "BCInletHelioUnsteadyMHD::interpolateSurfaceDataAtTime() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::extrapolateSpatialData()
{
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  const CFuint nbSurf = m_surfaceAtTime.size();
  cf_assert(nbSurf >= 1);
  
  RealVector tmpNode(m_dim);
  
  ClosestPointData closestPoint;
  closestPoint.surfaceIDs.resize(m_nbClosestPoints); 
  closestPoint.pointsIDs.resize(m_nbClosestPoints); 
  closestPoint.r.resize(m_nbClosestPoints);
  
  // Get TRS list to rebuild flux points
  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trsList.size();
  
  // Find our TRS
  SafePtr<TopologicalRegionSet> thisTRS;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    if (m_trsName[0] == trsList[iTRS]->getName()) {
      thisTRS = trsList[iTRS];
      break;
    }
  }
  
  // Rebuild boundary flux points structure
  std::vector<FlxPntStruct> bndFlxPnts;
  
  // Get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");
  
  // Get face builder data
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = thisTRS;
  faceData.isBoundary = true;
  
  map< std::string , vector< vector< CFuint > > >&
      bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[thisTRS->getName()];
  
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1; 
  const CFuint nbTRs = thisTRS->getNbTRs();
  
  CFuint localFaceID = 0;
  
  // Rebuild flux points for current surface data
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
    for (CFuint orient = 0; orient < nbOrients; ++orient) {
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][orient+1];
      
      for (CFuint iFace = startFaceIdx; iFace < stopFaceIdx; ++iFace) {
        faceData.idx = iFace;
        GeometricEntity *const face = _faceBuilder->buildGE();
        
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx) {
          m_flxPntCoords[iFlx] = face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
          
          FlxPntStruct thisFlxPnt;
          thisFlxPnt.setCentreCoordinates(m_flxPntCoords[iFlx]);
          thisFlxPnt.setLocalFaceID(localFaceID);
          thisFlxPnt.setLocalFluxID(iFlx);
          bndFlxPnts.push_back(thisFlxPnt);
        }
        
        _faceBuilder->releaseGE();
        localFaceID++;
      }
    }
  }
  
  // Now do the spatial interpolation with current surface data
  const CFuint nbBnd = bndFlxPnts.size();
  for(CFuint iFlx = 0; iFlx < nbBnd; iFlx++) {
    FlxPntStruct& currFlxPnt = bndFlxPnts[iFlx];
    closestPoint.reset();
    
    bool flagOut = false;
    for (CFuint is = 0; is < nbSurf && (!flagOut); ++is) {
      const SurfaceData& sf = *m_surfaceAtTime[is];
      const CFuint nbPoints = sf.rho.size();
      for (CFuint ip = 0; ip < nbPoints; ++ip) {
        sf.xyz.putRow(ip, tmpNode);
        const CFreal distance = MathFunctions::getDistance(currFlxPnt.getCentreCoordinates(), tmpNode);
        
        if (distance < 1.0e-8) {
          // Exact match - use this point directly
          m_flxPntRhos[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.rho[ip];
          m_flxPntUs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.u[ip];
          m_flxPntVs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.v[ip];
          m_flxPntWs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.w[ip];
          m_flxPntBxs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.Bx[ip];
          m_flxPntBys[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.By[ip];
          m_flxPntBzs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.Bz[ip];
          m_flxPntPs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.p[ip];
          flagOut = true;
          break;
        } else {
          // Build closest points stencil
          CFint counter = -1;
          for (CFuint n = 0; n < m_nbClosestPoints; ++n) {
            if (distance < closestPoint.r[n]) {
              counter++;
            }
          }
          
          if (counter >= 0) {
            for (CFuint i = 0; i < static_cast<CFuint>(counter); ++i) {
              closestPoint.regressionFromTo(i+1, i);
            }
            closestPoint.surfaceIDs[counter] = is;
            closestPoint.pointsIDs[counter] = ip;
            closestPoint.r[counter] = distance;
          }
        }
      }
    }
    
    if (!flagOut) {
      // Weighted interpolation from closest points
      CFreal matchingRho = 0.0, matchingU = 0.0, matchingV = 0.0, matchingW = 0.0;
      CFreal matchingBx = 0.0, matchingBy = 0.0, matchingBz = 0.0, matchingP = 0.0;
      CFreal sumWeights = 0.0;
      
      for (CFuint n = 0; n < m_nbClosestPoints; ++n) {
        const CFuint idxs = closestPoint.surfaceIDs[n];
        const SurfaceData& sf = *m_surfaceAtTime[idxs];
        const CFreal weight = 1./closestPoint.r[n];
        sumWeights += weight;
        const CFuint idxp = closestPoint.pointsIDs[n];
        
        matchingRho += weight * sf.rho[idxp];
        matchingU += weight * sf.u[idxp];
        matchingV += weight * sf.v[idxp];
        matchingW += weight * sf.w[idxp];
        matchingBx += weight * sf.Bx[idxp];
        matchingBy += weight * sf.By[idxp];
        matchingBz += weight * sf.Bz[idxp];
        matchingP += weight * sf.p[idxp];
      }
      
      matchingRho /= sumWeights;
      matchingU /= sumWeights;
      matchingV /= sumWeights;
      matchingW /= sumWeights;
      matchingBx /= sumWeights;
      matchingBy /= sumWeights;
      matchingBz /= sumWeights;
      matchingP /= sumWeights;
      
      m_flxPntRhos[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingRho;
      m_flxPntUs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingU;
      m_flxPntVs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingV;
      m_flxPntWs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingW;
      m_flxPntBxs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingBx;
      m_flxPntBys[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingBy;
      m_flxPntBzs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingBz;
      m_flxPntPs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingP;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::performInitialSpatialInterpolation(const std::vector<SurfaceData*>& surfaces, 
                                                                 const std::vector<FlxPntStruct>& bndFlxPnts)
{
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  const CFuint nbSurf = surfaces.size();
  cf_assert(nbSurf >= 1);
  
  RealVector tmpNode(m_dim);
  
  ClosestPointData closestPoint;
  closestPoint.surfaceIDs.resize(m_nbClosestPoints); 
  closestPoint.pointsIDs.resize(m_nbClosestPoints); 
  closestPoint.r.resize(m_nbClosestPoints);
  
  const CFuint nbBnd = bndFlxPnts.size(); 
  for(CFuint iFlx=0; iFlx<nbBnd; iFlx++) 
  {
    const FlxPntStruct& currFlxPnt = bndFlxPnts[iFlx];

    // during this preprocessing m_nbClosestPoints closest neighbors are sought
    closestPoint.reset();
      
    bool flagOut = false;
    for (CFuint is = 0; is < nbSurf && (!flagOut); ++is) 
    {
      const SurfaceData& sf = *surfaces[is];
      const CFuint nbPoints = sf.rho.size();
      for (CFuint ip = 0; ip < nbPoints; ++ip) 
      {
        sf.xyz.putRow(ip,tmpNode);

        const CFreal distance = MathFunctions::getDistance(currFlxPnt.getCentreCoordinates(), tmpNode);

        if (distance < 1.0e-8) 
        {
          // in this case we assume that the current node coincides with the mapping node
          cf_assert(ip < sf.rho.size());
          m_flxPntRhos[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.rho[ip];
          m_flxPntUs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.u[ip];
          m_flxPntVs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.v[ip];
          m_flxPntWs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.w[ip];
          m_flxPntBxs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.Bx[ip];
          m_flxPntBys[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.By[ip];
          m_flxPntBzs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.Bz[ip];
          m_flxPntPs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.p[ip];
          flagOut = true;
          break;
        }
        else 
        {
          CFint counter = -1;
          for (CFuint n = 0; n < m_nbClosestPoints; ++n) 
          {
            cf_assert(n < closestPoint.r.size());
            if (distance < closestPoint.r[n]) 
            {
              counter++;
            }
          }
            
          if (counter >= 0) 
          {
            // counter == (m_nbClosestPoints-1) corresponds to the point with min distance within the stencil
            // counter == 0                     corresponds to the point with max distance within the stencil
            for (CFuint i = 0; i < static_cast<CFuint>(counter); ++i) 
            {
              closestPoint.regressionFromTo(i+1, i);
            }
              
            cf_assert(counter < closestPoint.surfaceIDs.size());
            closestPoint.surfaceIDs[counter] = is;
            cf_assert(counter < closestPoint.pointsIDs.size());
            closestPoint.pointsIDs[counter] = ip;
            cf_assert(counter < closestPoint.r.size());
            closestPoint.r[counter] = distance;
          }
        }
      }
    }

    if (!flagOut) 
    {
      CFreal matchingRhow = 0.0;
      CFreal matchingUw = 0.0;
      CFreal matchingVw = 0.0;
      CFreal matchingWw = 0.0;
      CFreal matchingBxw = 0.0;
      CFreal matchingByw = 0.0;
      CFreal matchingBzw = 0.0;
      CFreal matchingPw = 0.0;
      CFreal sumWeights = 0.0;
      
      for (CFuint n = 0; n < m_nbClosestPoints; ++n) 
      {
        const CFuint idxs = closestPoint.surfaceIDs[n];
        cf_assert(idxs < surfaces.size());
        const SurfaceData& sf = *surfaces[idxs];
        cf_assert(closestPoint.r[n] > 0.);
        const CFreal weight = 1./closestPoint.r[n]; 
        sumWeights += weight;
        const CFuint idxp = closestPoint.pointsIDs[n];
        cf_assert(idxp < sf.rho.size());
               matchingRhow += weight*sf.rho[idxp]; 
        matchingUw += weight*sf.u[idxp]; 
        matchingVw += weight*sf.v[idxp]; 
        matchingWw += weight*sf.w[idxp]; 
        matchingBxw += weight*sf.Bx[idxp]; 
        matchingByw += weight*sf.By[idxp]; 
        matchingBzw += weight*sf.Bz[idxp]; 
        matchingPw += weight*sf.p[idxp]; 
      }
      
      matchingRhow /= sumWeights;
      matchingUw /= sumWeights;
      matchingVw /= sumWeights;
      matchingWw /= sumWeights;
      matchingBxw /= sumWeights;
      matchingByw /= sumWeights;
      matchingBzw /= sumWeights;
      matchingPw /= sumWeights;
        
      m_flxPntRhos[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingRhow;
      m_flxPntUs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingUw;
      m_flxPntVs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingVw;
      m_flxPntWs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingWw;
      m_flxPntBxs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingBxw;
      m_flxPntBys[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingByw;
      m_flxPntBzs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingBzw;  
      m_flxPntPs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingPw;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::readSurfaceData(std::vector<SurfaceData*>& surfaces, const std::string& fileName)
{ 
  const CFuint dim = (m_extractCoordZID < 0) ? PhysicalModelStack::getActive()->getDim() : 3;

  boost::filesystem::path fname;
  if (Common::StringOps::startsWith(fileName,"."))
  {
    fname = boost::filesystem::path(fileName);
  }
  else if (boost::filesystem::path(fileName).is_absolute())
  {
    // If fileName is already an absolute path, use it directly
    fname = boost::filesystem::path(fileName);
  }
  else
  {
    fname = Environment::DirPaths::getInstance().getBaseDir() / boost::filesystem::path(fileName);
  }

  // check the file type
  std::ifstream fin(fname.string().c_str());
  if(!fin) throw Common::FilesystemException (FromHere(),"Could not open file: " + fname.string());
  // The format is as follows:
  
  // NumberOfSurfaces
  // SURFACE_NAME1 NumberOfPoints
  // x y z Rho u v w Bx By Bz P
  // ...
  // SURFACE_NAME2 NumberOfPoints
  // x y z Rho u v w Bx By Bz P
  // ...
  
  CFLogInfo("BCInletHelioUnsteadyMHD::readSurfaceData() => START reading file " <<
	    fileName << "\n");
  
  CFuint nbSurf = 0;
  fin >> nbSurf;
  surfaces.resize(nbSurf);
  
  // store all the surface data
  for (CFuint is = 0; is < nbSurf; ++is) {
    std::string nameSurf = "";
    fin >> nameSurf;
    CFuint nbPoints = 0;
    fin >> nbPoints;
    
    CFLogInfo("nbSurf = " << is << "/" << nbSurf << ", NameSurf = " <<
	      nameSurf << ", nbPoints = " << nbPoints << "\n");
    SurfaceData* sf = new SurfaceData();
    sf->xyz.resize(nbPoints,dim);
    sf->rho.resize(nbPoints);
    sf->u.resize(nbPoints);
    sf->v.resize(nbPoints);
    sf->w.resize(nbPoints);
    sf->Bx.resize(nbPoints);
    sf->By.resize(nbPoints);
    sf->Bz.resize(nbPoints);
    sf->p.resize(nbPoints);
    
    for (CFuint ip = 0; ip < nbPoints; ++ip) {
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
	fin >> sf->xyz(ip,iDim);
      }
      
      if (m_xvec.size() > 0) {
	cf_assert(std::abs(m_angle) > 0.0);
	
	// rotate coordinates
	const CFreal x1 =  sf->xyz(ip,m_xvec[0])*std::cos(m_angle) + sf->xyz(ip,m_xvec[1])*std::sin(m_angle);
	const CFreal y1 = -sf->xyz(ip,m_xvec[0])*std::sin(m_angle) + sf->xyz(ip,m_xvec[1])*std::	cos(m_angle);
	sf->xyz(ip,m_xvec[0]) = x1;				   
	sf->xyz(ip,m_xvec[1]) = y1;				   
      }	
      fin >> sf->rho[ip];
      fin >> sf->u[ip];
      fin >> sf->v[ip];
      fin >> sf->w[ip];
      fin >> sf->Bx[ip];
      fin >> sf->By[ip];
      fin >> sf->Bz[ip];
      fin >> sf->p[ip];
      //CFLogInfo("ip = " << ip << " rho= " << sf->Vr[ip] << ", p= " << sf->rho[ip] << ", Vr= " << sf->p[ip] << "\n");
    }
	
    surfaces[is] = sf; //(m_extractCoordZID < 0) ?  sf : extractLineData(sf);
  }
  
  CFLogInfo("BCInletHelioUnsteadyMHD::readSurfaceData() => END reading file " <<
	    fileName << "\n"); 
}

//////////////////////////////////////////////////////////////////////////////

typename BCInletHelioUnsteadyMHD::SurfaceData* 
BCInletHelioUnsteadyMHD::extractLineData(SurfaceData* surface)
{
  // if m_extractCoordZID >= 0 a line distribution is extracted from a {x y z T} distribution for z=0
  //                           with z = x[m_extractCoordZID]
  
  CFreal minZ = MathConsts::CFrealMax();
  const CFuint nbPoints = surface->rho.size();
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    minZ = std::min(std::abs(surface->xyz(ip, m_extractCoordZID)), minZ);
  }
  
  CFuint nbPointsOnLine = 0;
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    if (MathChecks::isEqual(std::abs(surface->xyz(ip, m_extractCoordZID)), minZ)) {
      nbPointsOnLine++;
    }
  }
  
  SurfaceData* newData = new SurfaceData(); 
  newData->xyz.resize(nbPointsOnLine,2);
  newData->rho.resize(nbPointsOnLine);
  newData->u.resize(nbPointsOnLine);
  newData->v.resize(nbPointsOnLine);
  newData->w.resize(nbPointsOnLine);
  newData->Bx.resize(nbPointsOnLine);
  newData->By.resize(nbPointsOnLine);
  newData->Bz.resize(nbPointsOnLine);
  newData->p.resize(nbPointsOnLine);
  
  CFuint counter = 0;
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    if (MathChecks::isEqual(std::abs(surface->xyz(ip, m_extractCoordZID)), minZ)) {
      cout << surface->rho[ip] << endl;
      cout << surface->u[ip] << endl;
      cout << surface->v[ip] << endl;
      cout << surface->w[ip] << endl;
      cout << surface->Bx[ip] << endl;
      cout << surface->By[ip] << endl;
      cout << surface->Bz[ip] << endl;
      cout << surface->p[ip] << endl;
      newData->xyz(counter, XX) = surface->xyz(ip, m_extractCoordXYID[XX]);
      newData->xyz(counter, YY) = surface->xyz(ip, m_extractCoordXYID[YY]);
      newData->rho[counter] = surface->rho[ip];
      newData->u[counter] = surface->u[ip];
      newData->v[counter] = surface->v[ip];
      newData->w[counter] = surface->w[ip];
      newData->Bx[counter] = surface->Bx[ip];
      newData->By[counter] = surface->By[ip];
      newData->Bz[counter] = surface->Bz[ip];
      newData->p[counter] = surface->p[ip];
      counter++;
    }
  }
  
  cf_always_assert(counter == nbPointsOnLine);
  
  // remove old SurfaceData
  deletePtr(surface);
  
  //replace with new SurfaceData
  return newData;  
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::unsetup()
{
  CFAUTOTRACE;

  // Cleanup m_initialSolutionMap
  if (m_initialSolutionMap.size() > 0) 
  {
    for (CFuint i = 0; i < m_initialSolutionMap.size(); ++i) 
    {
      deletePtr(m_initialSolutionMap[i]);
    }
    m_initialSolutionMap.clear();
  }
  
  // Cleanup m_cellStatesFlxPnt
  for (CFuint iFlx = 0; iFlx < m_cellStatesFlxPnt.size(); ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt[iFlx]);
  }
  m_cellStatesFlxPnt.clear();
  
  // Cleanup time interpolation data - m_allSurfaces
  for (CFuint it = 0; it < m_allSurfaces.size(); ++it) {
    for (CFuint is = 0; is < m_allSurfaces[it].size(); ++is) {
      deletePtr(m_allSurfaces[it][is]);
    }
    m_allSurfaces[it].clear();
  }
  m_allSurfaces.clear();
  
  // Cleanup m_surfaceAtTime
  for (CFuint is = 0; is < m_surfaceAtTime.size(); ++is) {
    deletePtr(m_surfaceAtTime[is]);
  }
  m_surfaceAtTime.clear();

  // unsetup of the parent class
  BCStateComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void BCInletHelioUnsteadyMHD::updateBoundaryFluxPointValues()
{
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  if (!m_useTimeInterpolation) return;
  
  const CFuint nbSurf = m_surfaceAtTime.size();
  if (nbSurf == 0) {
    CFLog(ERROR, "No surface data available for updating boundary flux point values\n");
    return;
  }
  
  RealVector tmpNode(m_dim);
  
  ClosestPointData closestPoint;
  closestPoint.surfaceIDs.resize(m_nbClosestPoints); 
  closestPoint.pointsIDs.resize(m_nbClosestPoints); 
  closestPoint.r.resize(m_nbClosestPoints);
  
  // Rebuild the boundary flux points structure to get coordinates
  // This is necessary because we need the flux point coordinates for spatial interpolation
  std::vector<FlxPntStruct> bndFlxPnts;
  
  // Get face builder data
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  
  map< std::string , vector< vector< CFuint > > >&
      bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[m_thisTRS->getName()];
  
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1; 
  const CFuint nbTRs = m_thisTRS->getNbTRs();
  
  CFuint localFaceID = 0;
  
  // Rebuild flux points to get their coordinates
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
    for (CFuint orient = 0; orient < nbOrients; ++orient) {
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][orient+1];
      
      for (CFuint iFace = startFaceIdx; iFace < stopFaceIdx; ++iFace) {
        faceData.idx = iFace;
        GeometricEntity *const face = _faceBuilder->buildGE();
        
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx) {
          m_flxPntCoords[iFlx] = face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
          
          FlxPntStruct thisFlxPnt;
          thisFlxPnt.setCentreCoordinates(m_flxPntCoords[iFlx]);
          thisFlxPnt.setLocalFaceID(localFaceID);
          thisFlxPnt.setLocalFluxID(iFlx);
          bndFlxPnts.push_back(thisFlxPnt);
        }
        
        _faceBuilder->releaseGE();
        localFaceID++;
      }
    }
  }
  
  // Now perform spatial interpolation with current time-interpolated surface data
  const CFuint nbBnd = bndFlxPnts.size();
  for(CFuint iFlx = 0; iFlx < nbBnd; iFlx++) {
    const FlxPntStruct& currFlxPnt = bndFlxPnts[iFlx];
    closestPoint.reset();
    
    bool flagOut = false;
    for (CFuint is = 0; is < nbSurf && (!flagOut); ++is) {
      const SurfaceData& sf = *m_surfaceAtTime[is];
      const CFuint nbPoints = sf.rho.size();
      for (CFuint ip = 0; ip < nbPoints; ++ip) {
        sf.xyz.putRow(ip, tmpNode);
        const CFreal distance = MathFunctions::getDistance(currFlxPnt.getCentreCoordinates(), tmpNode);
        
        if (distance < 1.0e-8) {
          // Exact match - use this point directly
          m_flxPntRhos[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.rho[ip];
          m_flxPntUs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.u[ip];
          m_flxPntVs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.v[ip];
          m_flxPntWs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.w[ip];
          m_flxPntBxs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.Bx[ip];
          m_flxPntBys[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.By[ip];
          m_flxPntBzs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.Bz[ip];
          m_flxPntPs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = sf.p[ip];
          flagOut = true;
          break;
        } else {
          // Build closest points stencil
          CFint counter = -1;
          for (CFuint n = 0; n < m_nbClosestPoints; ++n) {
            if (distance < closestPoint.r[n]) {
              counter++;
            }
          }
          
          if (counter >= 0) {
            for (CFuint i = 0; i < static_cast<CFuint>(counter); ++i) {
              closestPoint.regressionFromTo(i+1, i);
            }
            closestPoint.surfaceIDs[counter] = is;
            closestPoint.pointsIDs[counter] = ip;
            closestPoint.r[counter] = distance;
          }
        }
      }
    }
    
    if (!flagOut) {
      // Weighted interpolation from closest points
      CFreal matchingRho = 0.0, matchingU = 0.0, matchingV = 0.0, matchingW = 0.0;
      CFreal matchingBx = 0.0, matchingBy = 0.0, matchingBz = 0.0, matchingP = 0.0;
      CFreal sumWeights = 0.0;
      
      for (CFuint n = 0; n < m_nbClosestPoints; ++n) {
        const CFuint idxs = closestPoint.surfaceIDs[n];
        const SurfaceData& sf = *m_surfaceAtTime[idxs];
        const CFreal weight = 1./closestPoint.r[n];
        sumWeights += weight;
        const CFuint idxp = closestPoint.pointsIDs[n];
        
        matchingRho += weight * sf.rho[idxp];
        matchingU += weight * sf.u[idxp];
        matchingV += weight * sf.v[idxp];
        matchingW += weight * sf.w[idxp];
        matchingBx += weight * sf.Bx[idxp];
        matchingBy += weight * sf.By[idxp];
        matchingBz += weight * sf.Bz[idxp];
        matchingP += weight * sf.p[idxp];
      }
      
      matchingRho /= sumWeights;
      matchingU /= sumWeights;
      matchingV /= sumWeights;
      matchingW /= sumWeights;
      matchingBx /= sumWeights;
      matchingBy /= sumWeights;
      matchingBz /= sumWeights;
      matchingP /= sumWeights;
      
      m_flxPntRhos[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingRho;
      m_flxPntUs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingU;
      m_flxPntVs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingV;
      m_flxPntWs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingW;
      m_flxPntBxs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingBx;
      m_flxPntBys[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingBy;
      m_flxPntBzs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingBz;
      m_flxPntPs[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPntsMax + currFlxPnt.getLocalFluxID()] = matchingP;
    }
  }
  
  CFLog(VERBOSE, "Updated " << nbBnd << " boundary flux point values with time-interpolated data\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD


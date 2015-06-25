#include <fstream>

#include "PARADE/ParadeLibrary.hh"
#include "PARADE/Parade.hh"

#include "Common/CFLog.hh"
#include "Common/DebugFunctions.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Common/BadValueException.hh"
#include "Common/Stopwatch.hh"
#include "Common/PEFunctions.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/PhysicalConsts.hh"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Parade {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ParadeLibrary,
			    RadiationLibrary,
			    ParadeModule,
			    1>
paradeLibraryProvider("Parade");

//////////////////////////////////////////////////////////////////////////////
      
void ParadeLibrary::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("EmisDom","Flag for emission dominated case.");
  options.addConfigOption< bool >("CompHeatFlux","Flag for computing the heat flux according to the tangent slab approximation.");
  options.addConfigOption< CFreal >("Tb1","Temperature at boundary 1 to solve the RTE.");
  options.addConfigOption< CFreal >("Tb2","Temperature at boundary 2 to solve the RTE.");
  options.addConfigOption< CFreal >("UndRel","Under relaxation factor.");
  options.addConfigOption< CFreal >("NDmin","Minimum number density.");
  options.addConfigOption< CFreal >("Tmin","Minimum temperature.");
  options.addConfigOption< string >("LocalDirName","Name of the local temporary directories where Parade is run.");
  options.addConfigOption< CFuint, Config::DynamicOption<> >
    ("ReuseProperties", "Reuse existing radiative data (requires the same number of processors as in the previous run).");
}
      
//////////////////////////////////////////////////////////////////////////////

ParadeLibrary::ParadeLibrary(const std::string& name) :
  RadiationLibrary(name),
  m_inFileHandle(),
  m_outFileHandle(),
  m_paradeDir(),
  m_gridFile(),
  m_tempFile(),
  m_densFile(),
  m_radFile(),
  m_library(CFNULL),
  m_spectrumSize(0),
  m_trTempID(),
  m_elTempID(),
  m_vibTempID(),
  m_wavelengths(),
  _kNu(),
  _emNu(),
  _tauNu(),
  _sourceNu(),
  m_mmasses(),
  m_avogadroOvMM(),
  m_molecularSpecies()
{
  addConfigOptionsTo(this);
  
  flag_Em = "false"; 
  setParameter("EmisDom",&flag_Em); 

  flag_HeatFlux = "false";
  setParameter("CompHeatFlux",&flag_HeatFlux); 

  _Tb1 = 0.0;
  setParameter("Tb1",&_Tb1);

  _Tb2 = 0.0;
  setParameter("Tb2",&_Tb2);
  
  _uFac = 0.0;    
  setParameter("UndRel", &_uFac);
  
  m_ndminFix = 1e+10;
  setParameter("NDmin", &m_ndminFix);
  
  m_tminFix = 300.;
  setParameter("Tmin", &m_tminFix);
  
  m_localDirName = "Parade";
  setParameter("LocalDirName", &m_localDirName);
  
  m_reuseProperties = 0;
  setParameter("ReuseProperties", &m_reuseProperties);
}
      
//////////////////////////////////////////////////////////////////////////////

ParadeLibrary::~ParadeLibrary()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::configure ( Config::ConfigArgs& args )
{
  RadiationLibrary::configure(args);
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ParadeLibrary::setup()
{
  RadiationLibrary::setup();
  
  m_inFileHandle  = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  m_outFileHandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  // create a m_localDirName-P# directory for the current process
  // cp parade executable and database inside the m_localDirName-P# directory 
  m_paradeDir = m_localDirName + "-P" + StringOps::to_str(PE::GetPE().GetRank(nsp));
  boost::filesystem::path paradeDir(m_paradeDir);
  
  boost::filesystem::path grid("grid.flo");
  m_gridFile = Environment::DirPaths::getInstance().getWorkingDir() / paradeDir / grid;
  
  boost::filesystem::path temp("temp.flo");
  m_tempFile = Environment::DirPaths::getInstance().getWorkingDir() / paradeDir / temp;
  
  boost::filesystem::path dens("dens.flo");
  m_densFile = Environment::DirPaths::getInstance().getWorkingDir() / paradeDir / dens;
  
  boost::filesystem::path rad("parade.rad");
  m_radFile = Environment::DirPaths::getInstance().getWorkingDir() / paradeDir / rad;
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  if (m_library.isNotNull()) {
    const CFuint nbSpecies = m_library->getNbSpecies();
    m_mmasses.resize(nbSpecies);
    m_library->getMolarMasses(m_mmasses);
    
    m_avogadroOvMM.resize(nbSpecies);
    m_avogadroOvMM = PhysicalConsts::Avogadro()/m_mmasses;
    
    vector<CFuint> moleculeIDs;
    m_library->setMoleculesIDs(moleculeIDs);
    cf_assert(moleculeIDs.size() > 0);
    
    m_molecularSpecies.resize(nbSpecies, false);
    for (CFuint i = 0; i < moleculeIDs.size(); ++i) {
      m_molecularSpecies[moleculeIDs[i]] = true;
    }
    
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint nbTemperatures = 1 + m_library->getNbTempVib() + m_library->getNbTe();
    // we assume here 2 (T-Tv) or 3 (T, Tv, Te) temperatures
    // translational temperature ID
    m_trTempID = nbEqs - nbTemperatures; 
    // electron temperature ID
    m_elTempID = (m_library->getNbTe() > 0) ? m_trTempID + 2 : m_trTempID + 1;
    // vibrational temperature ID
    m_vibTempID = m_trTempID + 1;
  }
  
  // if this is a parallel simulation, only ONE process at a time sets the library
  runSerial<void, ParadeLibrary, &ParadeLibrary::setLibrarySequentially>(this, nsp); 
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ParadeLibrary::setLibrarySequentially()
{ 
  // the path to the data files must be specified (there is no default)
  if (m_libPath == "") {
    CFLog(NOTICE, "ParadeLibrary::libpath NOT SET\n");
    abort();
  }
  
  if (!m_reuseProperties) {
    // create a m_localDirName-P# directory for the current process
    // cp parade executable and database inside the m_localDirName-P# directory 
    std::string command1   = "rm -fr " + m_paradeDir + " ; mkdir " + m_paradeDir;
    Common::OSystem::getInstance().executeCommand(command1);
    std::string command2   = "cp " + m_libPath + "/Binary/parade " + m_paradeDir;
    Common::OSystem::getInstance().executeCommand(command2);
    std::string command3   = "cp -R " + m_libPath + "/Data " + m_paradeDir;
    Common::OSystem::getInstance().executeCommand(command3);
  }
  
  // all the following can fail if the format of the file parade.con changes
  // initialize some private data in this class (all processors, one by one, must run this)
  ifstream fin("parade.con"); 
  string line;
  while (getline(fin,line)) {
    // look for the number of points in the spectrum
    if (m_spectrumSize == 0) {
      readValue<CFuint>(line, "npoints", m_spectrumSize);
    }
    
    // look for the minimum wavelenght if it has not been specified in the CFcase file
    if (m_wavMin < 0.) {
      readValue<CFreal>(line, "wavlo", m_wavMin);
    }
    
    // look for the maximum wavelenght if it has not been specified in the CFcase file
    if (m_wavMax < 0.) {
      readValue<CFreal>(line, "wavhi", m_wavMax);
    }
  }
  fin.close();
  
  cf_assert(m_spectrumSize > 0);
  cf_assert(m_spectrumSize < 1e8);
  cf_assert(m_wavMin > 0);
  cf_assert(m_wavMax > 0);
  
  // use the full spectrum size as stride if the stride is not specified 
  // (i.e., if the default value=1 is unchanged )
  if (m_wavStride == 1) {
    m_wavStride = m_spectrumSize;
  }
} 

//////////////////////////////////////////////////////////////////////////////
            
void ParadeLibrary::unsetup()
{
  if(isSetup()) {
    
    RadiationLibrary::unsetup();
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::runOnStagnationLine
(Common::SafePtr<std::vector<CFuint> > stagnationLineCells,
 ProxyDofIterator<CFreal>* pstates, CFreal* qrad)
{
  CFout << "ParadeLibrary::runOnStagnationLine()\n";
 
  writeStagnationLineData(stagnationLineCells, pstates);
  
  runLibrary();
 
  if (flag_Em == 0) {
    computeQradStagLine(qrad); 
  } 
  else { 
    computeQradEmDom(qrad);
  } 
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ParadeLibrary::runOnStructuredMesh
(const std::vector<std::vector<CFuint>* >& meshByLine,
 ProxyDofIterator<CFreal>* pstates, CFreal* qrad)
{
  CFout << "ParadeLibrary::runOnStructuredMesh()\n";

  writeAllData(meshByLine, pstates);
 
  runLibrary();
 
  if (flag_Em == 0) {
    computeQradMesh(meshByLine, pstates, qrad);
  } else { 
    computeQradEmDom(qrad);
  }

}

//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::writeStagnationLineData
(Common::SafePtr<std::vector<CFuint> > cellLine,
 ProxyDofIterator<CFreal>* pstates)
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // Fixes for temperature and number densities in order to avoid numerical problems 
  // when the radiation code computes the spectra
  // This needs to inverted
  const CFreal nbminFix = 1.e+10;
  const CFreal tminFix = 300.;

  // write  the mesh file
  ofstream* fout = new ofstream("grid.flo");
  *fout << "TINA" << endl;
  *fout <<  cellLine->size() << " " << "1" << endl;
  
  const CFuint nbCellsInLine = cellLine->size();
  for (CFuint i =0; i < nbCellsInLine; ++i) {
    const CFuint cellID = (*cellLine)[i];
    fout->precision(14);
    fout->setf(ios::scientific,ios::floatfield);
    
    CFreal *const node = pstates->getNode(cellID);
    
    if (dim == DIM_1D) {
      *fout << node[XX] << " " << 0.0 << " " << 0.0  << endl;
    }
    if (dim == DIM_2D) {
      *fout << node[XX] << " " << node[YY] << " " << 0.0  << endl;
    }
    if (dim == DIM_3D) {
      *fout << node << endl;
    }
  } 
  fout->close();
  
  // write the temperatures
  fout->open("temp.flo");
  const CFuint nbTemperatures = 1 + m_library->getNbTempVib() + m_library->getNbTe();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq(); 
  const CFuint tempID = nbEqs - nbTemperatures; 
  *fout <<  cellLine->size() << " "  << "1" <<  " " << nbTemperatures << endl;  

  // here it is assumed that the temperatures are the LAST variables 
  for (CFuint i =0; i < nbCellsInLine; ++i) {
    const CFuint cellID = (*cellLine)[i];
    CFreal *const currState = pstates->getState(cellID);
    for (CFuint t = 0; t < nbTemperatures; ++t) {
      fout->precision(14);
      fout->setf(ios::scientific,ios::floatfield);
      const CFreal temp = std::max(currState[tempID + t],tminFix);
      *fout << temp << " ";
//      *fout << currState[tempID + t] << " ";
    }
    *fout << endl;
  }
  fout->close();
  
  // write the number densities
  fout->open("dens.flo");
  const CFuint nbSpecies = m_library->getNbSpecies();
  *fout << cellLine->size() <<  " " << "1" << " " << nbSpecies << endl; 
 
  // here it is assumed that the species densities are the FIRST variables 
  for (CFuint i =0; i < nbCellsInLine; ++i) {
    const CFuint cellID = (*cellLine)[i];
    CFreal *const currState = pstates->getState(cellID);
    for (CFuint t = 0; t < nbSpecies; ++t) {
      fout->precision(14);
      fout->setf(ios::scientific,ios::floatfield);
      // number Density = partial density/ molar mass * Avogadro number
      const CFreal nb = std::max(currState[t]*m_avogadroOvMM[t],nbminFix);
      *fout << nb << " "; // here convert to number densities
    }
    *fout << endl;
  }
  fout->close();
  
  delete fout;  

}

//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::writeAllData(const std::vector<std::vector<CFuint>* >& meshByLine,
				 ProxyDofIterator<CFreal>* pstates)
{ 
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // Fixes for temperature and number densities in order to avoid numerical problems 
  // when the radiation code computes the spectra
  const CFreal nbminFix = 1.e+4;
  const CFreal tminFix = 300.;

  // write  the mesh file
  ofstream* fout = new ofstream("grid.flo");
  *fout << "TINA" << endl;
  
  const CFuint iMax = meshByLine[0]->size();
  const CFuint kMax = meshByLine.size();
  *fout << iMax << " " << kMax << endl;
  for (CFuint i = 0; i < iMax; ++i) {       
    for (CFuint k = 0; k < kMax; ++k) {
      vector<CFuint>& cellLine = *meshByLine[k];
      const CFuint cellID = cellLine[i];
      fout->precision(14);
      fout->setf(ios::scientific,ios::floatfield);
      CFreal *const node = pstates->getNode(cellID);
      
      if (dim == DIM_2D) {
        *fout << node[XX] << " " << node[YY] << " " << 0.0  << endl;
      }
      if (dim == DIM_3D) {
        *fout << node << endl;
      }
          
    }
  }
  fout->close();  

  // write the indices i,k and corresponding cellID's (just for a check)
  fout->open("domain.dat");
  *fout << "Cell centroid axial and radial indices, identifiers and coordinates" << endl;
  *fout << "i " << "k " << " cell " << "x " << " r" << endl; 
  for (CFuint k = 0; k < kMax; ++k) {
       vector<CFuint>& cellLine = *meshByLine[k];  
       for (CFuint i = 0; i < iMax; ++i) {
         const CFuint cellID = cellLine[i];
         CFreal *const node = pstates->getNode(cellID);
         *fout << i << " " << k << " " << cellID << " " << node[XX] << " " << node[YY] << " " << endl;                   
      }
  }
  fout->close();
  // write the temperatures
  fout->open("temp.flo");
  const CFuint nbTemperatures = 1 + m_library->getNbTempVib() + m_library->getNbTe();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq(); 
  const CFuint tempID = nbEqs - nbTemperatures; 
  *fout << iMax << " " << kMax << " " << nbTemperatures << endl;
   
  // @TODO here it is assumed that the temperatures are the LAST variables 
  for (CFuint i = 0; i < iMax; ++i) {
    for (CFuint k = 0; k < kMax; ++k) {
      vector<CFuint>& cellLine = *meshByLine[k];
      const CFuint cellID = cellLine[i];
      CFreal *const currState = pstates->getState(cellID);
      for (CFuint t = 0; t < nbTemperatures; ++t) {
	fout->precision(14);
	fout->setf(ios::scientific,ios::floatfield);
        // Fix introduced
        const CFreal temp = std::max(currState[tempID + t],tminFix);
        *fout << temp << " ";
      }
      *fout << endl;
    } 
  }
  fout->close();
  // write the number densities
  fout->open("dens.flo");
  const CFuint nbSpecies = m_library->getNbSpecies();
  *fout << iMax << " " << kMax << " " << nbSpecies << endl;
  
  // @TODO here it is assumed that the species densities are the FIRST variables 
  for (CFuint i = 0; i < iMax; ++i) {
    for (CFuint k = 0; k < kMax; ++k) {
      vector<CFuint>& cellLine = *meshByLine[k];
      const CFuint cellID = cellLine[i];
      CFreal *const currState = pstates->getState(cellID);
      for (CFuint t = 0; t < nbSpecies; ++t) { 
	fout->precision(14);
	fout->setf(ios::scientific,ios::floatfield);
	// number Density = partial density/ molar mass * Avogadro number (fix introduced)
	const CFreal nb = std::max(currState[t]*m_avogadroOvMM[t],nbminFix);
        *fout << nb << " ";  
      }
      *fout << endl;
    }
  }
  fout->close();

  delete fout; 

}
//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::computeQradEmDom(CFreal* qrad)
{
  CFLog(VERBOSE, "ParadeLibrary::computeQradEmissionDominated()");

  ifstream grid("grid.flo");

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
 
  string tmp;
  CFuint iMax = 0;
  CFuint kMax = 0;				
  CFuint nbCells = 0;
  CFreal dum = 0.;
  
  grid >> tmp; // TINA
  grid >> iMax >> kMax; // imax and kmax indices for cells in 1D and 2D structured mesh
  nbCells = iMax*kMax;

  // Check on the maximum number of slabs
  if (kMax > 1000) {
     cout << "Maximum number of slabs allowed is 1000" << endl;  
     abort();
  }
  
  //Useful constants (wavelenghts are specified in Amstrong units)
  const CFreal convA = 1.0e-10;
  const CFreal pi = 3.1415926535;

  //Reading the first slab file in order to fill the wavelenght vector
  CFuint wlPoints = 0;
  string inFile = "slab0000.dat";
  ifstream fin(inFile.c_str());
  fin >> dum >> wlPoints;

  m_wavelengths.resize(wlPoints); 
  _emNu.resize(wlPoints,1); 
 
  for (CFuint w = 0; w < wlPoints; ++w) {
       fin >>  m_wavelengths[w] >> dum >> dum; 
  }

  fin.close();
  
  // Computing the qrad loss term for each cell along each slab
  for (CFuint k = 0; k < kMax; ++k) {

      // Slab filename
      string inFile = "slab";

      if (k < 10) {
        inFile += "000";
      }
      if (k < 100 && k > 9) {
        inFile += "00";
      }
    
      if (k < 1000 && k > 99) {
        inFile += "0";
      }
    
      inFile += StringOps::to_str(k) + ".dat";    
      ifstream slab(inFile.c_str());

      // Loop over the cells contained in k-th slab
      for (CFuint i = 0; i<iMax; ++i) {

          // cell identifier
          CFuint cellID = 0; 
          if (dim == DIM_1D) {
            cellID = i;   
          }
          else {
            cellID = (nbCells - 1 - (k + i*kMax));
          }
          slab >> dum >> dum;

          // First wavelength data of the cell considered
          slab >> dum >> _emNu(0,0) >> dum;  

          // Reading the rest of the wavelength data of the cell considered
          qrad[cellID] = 0.0;
          for (CFuint w = 1; w < wlPoints; ++w) {
               slab >> dum >> _emNu(w,0) >> dum;
               const CFreal dLambda = (m_wavelengths[w] - m_wavelengths[w-1]);
               const CFreal h = 0.5*(_emNu(w,0) + _emNu(w-1,0));
               qrad[cellID] += dLambda*h;
          }

          // Needed conversion (wavelengths are given in Amstrong in spectra) 
          qrad[cellID] *= _uFac*4.0*pi*convA;
      }

      slab.close();

  } 

  CFout << "ParadeLibrary::deleting .arg .txt files for memory saving\n";
  std::string command = "rm -f *.imo *.txt *.arg slab*";
  Common::OSystem::getInstance().executeCommand(command);
  
  ofstream fout("qradEmDom.dat");
  cout << "ParadeLibrary::writing qrad source term to qradEmDom.dat\n";
  for (CFuint i = 0; i < nbCells; ++i) {
       fout << i << " " << qrad[i] << endl;
  }
  fout.close(); 
   
}
    
//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::computeQradStagLine(CFreal* qrad)
{  
  CFLog(VERBOSE, "ParadeLibrary::computeQradStagLine()");
  
  ifstream grid("grid.flo"); 

  string tmp = "";
  CFuint num = 0;				
  CFuint nbCells = 0;
  CFreal dum = 0.;
  CFuint iMax = 0;
  CFuint kMax = 0; 
  
  grid >> tmp; // TINA
  grid >> iMax >> kMax; // iax and kmax indices for cells in a 2D structured mesh
  nbCells = iMax*kMax;
  grid >> num >> nbCells; // dimension,  number of cells

  // Location of cells along stagnation line (or nozzle length)
  std::vector<CFreal> x;
  x.resize(nbCells); 
  for (CFuint i = 0; i < nbCells; ++i) {
    grid >> x[i] >> dum >> dum; 
  }
  grid.close();
  
  //Useful constants (wavelenghts are specified in Amstrong units)
  const CFreal convA = 1.0e-10;
  const CFreal pi = 3.1415926535;
  
  //Reading the first slab file in order to fill the wavelenght vector
  CFuint wlPoints = 0;
  string inFile = "slab0000.dat";
  ifstream fin(inFile.c_str());
  fin >> dum >> wlPoints;

  m_wavelengths.resize(wlPoints); 
  _kNu.resize(wlPoints,nbCells);
  _emNu.resize(wlPoints,nbCells); 
  _tauNu.resize(wlPoints,nbCells); 
  _sourceNu.resize(wlPoints,nbCells);
  _tauNu = 0.0;
  _sourceNu = 0.0;
 
  for (CFuint w = 0; w < wlPoints; ++w) {
       fin >>  m_wavelengths[w] >> dum >> dum; 
  }

  fin.close();

  // File where stagnation line data are stored  
  string slabFile = "slab0000.dat";
  ifstream stagLine(slabFile.c_str());

  for (CFuint i = 0; i < nbCells; ++i) {

       stagLine >> dum >> dum;
        
       // First wavelength data of the cell considered
       stagLine >> dum >> _emNu(0,i) >> _kNu(0,i);  
       if (i > 0) { 
         const CFreal dx = x[i] - x[i-1];
         const CFreal h = 0.5*(_kNu(0,i) + _kNu(0,i-1));
        _tauNu(0,i) = _tauNu(0,i-1) + dx*h;
       }    
         
       // Reading the rest of the wavelength data of the cell considered
       for (CFuint w = 1; w < wlPoints; ++w) {
            stagLine >> dum >> _emNu(w,i) >> _kNu(w,i);
            const CFreal dLambda = (m_wavelengths[w] - m_wavelengths[w-1]);
            const CFreal h = 0.5*(_emNu(w,i) + _emNu(w-1,i));
            qrad[i] += dLambda*h;
            if (i > 0) {
               const CFreal dx = x[i] - x[i-1]; 
               const CFreal h = 0.5*(_kNu(w,i) + _kNu(w,i-1));
               _tauNu(w,i) = _tauNu(w,i-1) + dx*h; 
            }
       }   

       // Needed conversion (wavelengths are given in Amstrong in spectra) 
       qrad[i] = 4.0*pi*convA*qrad[i];

  }
  
   // Adding the absorption contribution for tangent slab method 
   cout << "ParadeLibrary::application of tangent slab method to solve RTE\n"; 
            
   // Computation of absorption contribution 
   for (CFuint w = 0; w < wlPoints; ++w) {  
        for (CFuint i = 0; i < nbCells; ++i) {
             _sourceNu(w,i) = 0.0;
             // Left-to-right direction 
             for (CFuint j = 0; j < i; ++j) {
                  const CFreal dx = x[j+1] - x[j];
                  const CFreal tau1 = _tauNu(w,i) - _tauNu(w,j); 
                  const CFreal tau2 = _tauNu(w,i) - _tauNu(w,j+1);   
                  const CFreal e1m1 = E1(tau1);
                  const CFreal e1p1 = E1(tau2);
                  const CFreal h = 0.5*(_emNu(w,j)*e1m1 + _emNu(w,j+1)*e1p1);  
                  _sourceNu(w,i) += _kNu(w,i)*h*dx;          
             }
       
             // Right-to-left direction
             for (CFuint j = i+1; j < nbCells; ++j) {
                  const CFreal dx = x[j] - x[j-1];
                  const CFreal tau1 = _tauNu(w,j-1) - _tauNu(w,i);  
                  const CFreal tau2 = _tauNu(w,j) - _tauNu(w,i); 
                  const CFreal e1m1 = E1(tau1);
                  const CFreal e1p1 = E1(tau2); 
                  const CFreal h = 0.5*(_emNu(w,j-1)*e1m1 + _emNu(w,j)*e1p1); 
                  _sourceNu(w,i) += _kNu(w,i)*h*dx;
             }

             // Boundary conditions (side 1 and 2 of the line)
             const CFreal b1 = Planck((m_wavelengths[w]*convA), _Tb1);     // boundary 1 is the inlet  (for nozzle)
             const CFreal b2 = Planck((m_wavelengths[w]*convA), _Tb2);     // boundary 2 is the outlet (for nozzle)
             const CFreal tau1 = _tauNu(w,i) - _tauNu(w,0);
             const CFreal tau2 = _tauNu(w,nbCells-1) - _tauNu(w,i);
             const CFreal e21 = E2(tau1);
             const CFreal e22 = E2(tau2);
             _sourceNu(w,i) += _kNu(w,i)*(b1*e21 + b2*e22); 

        }
   }

   // Integration over wavelength range
   for (CFuint i = 0; i < nbCells; ++i) {
        for (CFuint w = 0; w < (wlPoints-1); ++w) { 
             const CFreal dLambda = m_wavelengths[w+1] - m_wavelengths[w];
             const CFreal h = 0.5*(_sourceNu(w+1,i) + _sourceNu(w,i));
             qrad[i] -= 2.0*pi*h*dLambda*convA;
        }
        qrad[i] = _uFac*qrad[i];
        
   } 

  CFout << "ParadeLibrary::deleting .imo .arg .txt files for memory saving\n";
  std::string command = "rm -f *.txt *.arg slab*";
  Common::OSystem::getInstance().executeCommand(command);
 
  ofstream fout("qrad.dat");
  cout << "ParadeLibrary::writing qrad source term to qrad.dat\n";
  for (CFuint i = 0; i < nbCells; ++i) {
    fout << i << " " << qrad[i] << endl;
  }
  fout.close();

  CFLog(VERBOSE, "ParadeLibrary::computeQradStagLine()");
 
}

//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::computeQradMesh
(const std::vector<std::vector<CFuint>* >& meshByLine,
 Framework::ProxyDofIterator<CFreal>* pstates,
 CFreal* qrad)
{  
  CFLog(VERBOSE, "ParadeLibrary::computeQradMesh()");
  
  ifstream grid("grid.flo");
  
  string tmp = "";
  CFreal dum = 0.; 			
  CFuint nbCells = 0;
  CFuint iMax = 0;
  CFuint kMax = 0;
  
  grid >> tmp;             // TINA
  grid >> iMax >> kMax;    // imax and kmax indices for cells in a 2D structured mesh
  nbCells = iMax*kMax;
  
  // Check on the maximum number of slabs
  if (kMax > 1000) {
     cout << "Maximum number of slabs allowed is 1000" << endl;  
     abort();
  }
 
  grid.close();
  
  //Useful constants (wavelenghts are specified in Amstrong units)
  const CFreal convA = 1.0e-10;
  const CFreal pi = 3.1415926535;
  const CFreal qradTol = 1.0e+8;

  //Reading the first slab file in order to fill the wavelenght vector
  CFuint wlPoints = 0;

  string inFile = "slab0000.dat";
  ifstream fin(inFile.c_str());

  // Reading the number of wavelength points 
  fin >> dum >> wlPoints;  

  // Allocating matrices 
  m_wavelengths.resize(wlPoints); 
  _kNu.resize(wlPoints,iMax);
  _emNu.resize(wlPoints,iMax); 
  _tauNu.resize(wlPoints,iMax);
  _sourceNu.resize(wlPoints,iMax);

  // Reading and allocating the wavelength vector
  for (CFuint w = 0; w < wlPoints; ++w) {
       fin >>  m_wavelengths[w] >> dum >> dum; 
  }

  fin.close();

  // splitting of qrad
//  ofstream qrad_fout("qrad_Em_Abs_Tot.dat");

  // Computing the qrad loss term for each cell along each slab
  for (CFuint k = 0; k < kMax; ++k) {

      cout << "ParadeLibrary::Solving the radiative transfer equation along mesh line " << k+1 << endl;

      // Slab filename
      string inFile = "slab";

      if (k < 10) {
        inFile += "000";
      }
      if (k < 100 && k > 9) {
        inFile += "00";
      }
    
      if (k < 1000 && k > 99) {
        inFile += "0";
      }
    
      inFile += StringOps::to_str(k) + ".dat";
      ifstream slab(inFile.c_str());

      // Loop over the cells contained in k-th slab
      vector<CFuint>& cellLine = *meshByLine[k];
      
      // Data inizialization
      _tauNu = 0.;
      _emNu = 0.;
      _kNu = 0.;
 
      // Reading data of the first cell of the slab
      const CFuint i = 0;

      // Cell ID
      const CFuint ID = cellLine[i]; 

      // Cell ID and number of wave lenght points (PARADE output file)
      slab >> dum >> dum;

      slab >> dum >> _emNu(0,i) >> _kNu(0,i);

      qrad[ID] = 0.0;
      for (CFuint w = 1; w < wlPoints; ++w) {
        slab >> dum >> _emNu(w,i) >> _kNu(w,i);
        const CFreal dLambda = (m_wavelengths[w] - m_wavelengths[w-1]);
        const CFreal h = 0.5*(_emNu(w,i) + _emNu(w-1,i));
        qrad[ID] += dLambda*h;
      }
      // Source term due to emission
      qrad[ID] *= 2.;

      // Reading data of the other cell cell of the slab
      for (CFuint i = 1; i<iMax; ++i) {

          // Cell identifier         
          const CFuint ID = cellLine[i]; 

          // Cell ID and number of wave lenght points (PARADE output file)
          slab >> dum >> dum;

          slab >> dum >> _emNu(0,i) >> _kNu(0,i);

          const CFuint cellID_m1 = cellLine[i-1];
          CFreal *const node1 = pstates->getNode(ID);
          CFreal *const node2 = pstates->getNode(cellID_m1);
          const CFreal ds = pow((pow((node2[XX] - node1[XX]),2) + pow((node2[YY] - node1[YY]),2)),0.5);
          const CFreal h = 0.5*(_kNu(0,i) + _kNu(0,i-1)); 
	  _tauNu(0,i) += _tauNu(0,i-1)*h*ds; 

          // Reading the rest of the wavelength data of the cell considered
          qrad[ID] = 0.0;

          for (CFuint w = 1; w < wlPoints; ++w) {
            slab >> dum >> _emNu(w,i) >> _kNu(w,i);
            const CFreal dLambda = (m_wavelengths[w] - m_wavelengths[w-1]);
            const CFreal h1 = 0.5*(_emNu(w,i) + _emNu(w-1,i));
            qrad[ID] += dLambda*h1;
               
            // Computation of plasma opticall thickness along the mesh line k 
            const CFuint cellID_m1 = cellLine[i-1];
            CFreal *const node1 = pstates->getNode(ID);
            CFreal *const node2 = pstates->getNode(cellID_m1);
            const CFreal ds = pow((pow((node2[XX] - node1[XX]),2) + pow((node2[YY] - node1[YY]),2)),0.5);
            const CFreal h2 = 0.5*(_kNu(w,i) + _kNu(w,i-1)); 
	    _tauNu(w,i) += _tauNu(w,i-1)*h2*ds;  
              
          }

          // Source term due to emission
          qrad[ID] *= 2.0;

      } 

      // Closing slab file
      slab.close();

      cout << "ParadeLibrary::Computed emission term and distribution of plasma optical thickness along mesh line " << k+1 << endl;
      
      // Loop over all mesh line cells in order to compute the asbsorption source term according to the tangent slab approximation  
      _sourceNu = 0.;
      for (CFuint i = 0; i<iMax; ++i) {

          // Cell identifier
          const CFuint ID = cellLine[i];
          for (CFuint w = 0; w < wlPoints; ++w) { 
          
              for (CFuint j = 0; j < i; ++j) {
         	const CFuint cellID = cellLine[j];
                const CFuint cellID_p1 = cellLine[j+1];
	        CFreal *const node1 = pstates->getNode(cellID);
	        CFreal *const node2 = pstates->getNode(cellID_p1);
	        const CFreal ds = pow((pow((node2[XX] - node1[XX]),2) + pow((node2[YY] - node1[YY]),2)),0.5);
	        const CFreal tau1 = _tauNu(w,i) - _tauNu(w,j); 
	        const CFreal tau2 = _tauNu(w,i) - _tauNu(w,j+1);   
	        const CFreal e1m1 = E1(tau1);
	        const CFreal e1p1 = E1(tau2);
                const CFreal h = 0.5*(_emNu(w,j)*e1m1 + _emNu(w,j+1)*e1p1);  
	        _sourceNu(w,i) += h*ds;
              }

              for (CFuint j = i+1; j < iMax; ++j) {
         	const CFuint cellID = cellLine[j];
                const CFuint cellID_m1 = cellLine[j-1];
	        CFreal *const node1 = pstates->getNode(cellID);
                CFreal *const node2 = pstates->getNode(cellID_m1);
	        const CFreal ds = pow((pow((node2[XX] - node1[XX]),2) + pow((node2[YY] - node1[YY]),2)),0.5); 	
	        const CFreal tau1 = _tauNu(w,j-1) - _tauNu(w,i);  
	        const CFreal tau2 = _tauNu(w,j) - _tauNu(w,i); 
	        const CFreal e1m1 = E1(tau1);
	        const CFreal e1p1 = E1(tau2); 
                const CFreal h = 0.5*(_emNu(w,j-1)*e1m1 + _emNu(w,j)*e1p1); 
	        _sourceNu(w,i) += h*ds;
	     }

             // Boundary conditions (side 1 and 2 of the line)
             const CFreal b1 = Planck((m_wavelengths[w]*convA), _Tb1);   // boundary 1 is the wall
             const CFreal b2 = Planck((m_wavelengths[w]*convA), _Tb2);   // boundary 2 is the shock
	     const CFreal tau1 = _tauNu(w,i) - _tauNu(w,0);
	     const CFreal tau2 = _tauNu(w,iMax-1) - _tauNu(w,i);
	     const CFreal e21 = E2(tau1);
	     const CFreal e22 = E2(tau2);
             _sourceNu(w,i) += (b1*e21 + b2*e22);

             // Multipling factor: absorption coefficient of the cell under consideration
             _sourceNu(w,i) *= _kNu(w,i);

          }

          CFreal qradEm = qrad[ID];
          CFreal qradAbs = 0.;
          // Integration over the wavelenght domain
          for (CFuint w = 1; w < wlPoints; ++w) { 
             const CFreal dLambda = m_wavelengths[w] - m_wavelengths[w-1];
             const CFreal h = 0.5*(_sourceNu(w,i) + _sourceNu(w-1,i));
 	     qrad[ID] -= h*dLambda;
             qradAbs +=  h*dLambda;
          }
          qrad[ID] *= _uFac*2.0*pi*convA;
          qradEm *= _uFac*2.0*pi*convA;
          qradAbs *= _uFac*2.0*pi*convA;

//          qrad_fout << ID << " " << qradEm << " " << qradAbs << " " << qrad[ID] << endl; 

     }

     cout << "ParadeLibrary::Computed absorption term along mesh line " << k+1 << endl << endl;

  }
 
//  qrad_fout.close();

  CFout << "ParadeLibrary::deleting .imo .arg .txt files for memory saving\n";
  std::string command = "rm -f *.txt *.arg slab*";
  Common::OSystem::getInstance().executeCommand(command);
  
  // Applying a fix to the radiative source term
  for (CFuint i = 0; i < nbCells; ++i) {
    if (qrad[i] >= qradTol) {
        qrad[i] = qradTol; 
    }
    if (qrad[i] <= (-qradTol)) {
        qrad[i] = -qradTol; 
    }
  } 
  
  ofstream fout("qrad.dat");
  cout << "ParadeLibrary::writing qrad source term to qrad.dat\n";
  for (CFuint i = 0; i < nbCells; ++i) {
    fout << i << " " << qrad[i] << endl;
  }
  fout.close();
  
  CFLog(VERBOSE, "ParadeLibrary::computeQradMesh()");
 
}

//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::computeRadWallHeatFlux
(const std::vector<std::vector<CFuint>* >& meshByLine,
 Framework::ProxyDofIterator<CFreal>* pstates)
{
  
 cout << "Not implemented yet " << endl;
 abort();
    
}

//////////////////////////////////////////////////////////////////////////////

CFreal ParadeLibrary::Planck(const CFreal& lambda, const CFreal& T)
{
  CFreal Planck = 0.0;

  CFreal c = 299792458.0;
  CFreal h = 6.626075e-34;
  CFreal kb = 1.380658e-23 ;
  CFreal eps = 1.0e-10;

  Planck = 2.0*h*pow(c,2)/pow(lambda,5)/(exp(h*c/(kb*T*lambda)) - 1.0 + eps);
  return Planck;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ParadeLibrary::E1(const CFreal& tau)
{
  CFreal e1 = 0.0;
  std::vector<CFreal> a;
  std::vector<CFreal> b;
 
  a.resize(2);
  b.resize(2);

  a[0] = 2.591;
  a[1] = 1.708;

  b[0] = -18.70;
  b[1] = -2.110;

  for (CFuint i = 0; i < 2; ++i) {
      e1 += a[i]*exp(tau*b[i]);
  } 

  return e1;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ParadeLibrary::E2(const CFreal& tau)
{
  CFreal e2 = 0.0;
  std::vector<CFreal> a;
  std::vector<CFreal> b;
 
  a.resize(2);
  b.resize(2);

  a[0] = 0.2653;
  a[1] = 0.7347;

  b[0] = -8.659;
  b[1] = -1.624;

  for (CFuint i = 0; i < 2; ++i) {
      e2 += a[i]*exp(tau*b[i]);
  } 

  return e2; 
}

//////////////////////////////////////////////////////////////////////////////

CFreal ParadeLibrary::E3(const CFreal& tau)
{
  CFreal e3 = 0.0;
  std::vector<CFreal> a;
  std::vector<CFreal> b;
 
  a.resize(2);
  b.resize(2);

  a[0] = 0.0929;
  a[1] = 0.4071;

  b[0] = -4.08;
  b[1] = -1.33;

  for (CFuint i = 0; i < 2; ++i) {
      e3 += a[i]*exp(tau*b[i]);
  } 

  return e3; 
}

//////////////////////////////////////////////////////////////////////////////

void ParadeLibrary::computeProperties(ProxyDofIterator<CFreal>* pstates,
				      RealMatrix& data, CFuint iWavRange) 
{ 
  CFLog(INFO,"ParadeLibrary::computeProperties() => START\n");
  
  if (!m_reuseProperties) {
    Stopwatch<WallTime> stp;
    
    // update the wavelength range inside parade.con
    // copy the modified files into the local Parade directories
    updateWavRange(iWavRange);
    
    // write the local grid, temperature and densities fields 
    writeLocalData(pstates);
    
    stp.start();
    
    // run concurrently Parade in each local directory, one per process 
    runLibraryInParallel();
    
    CFLog(INFO,"ParadeLibrary::runLibraryInParallel() took " << stp << "s\n");
  }
  
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  PE::GetPE().setBarrier(nsp);
  
  // read in the radiative properties from parade.rad
  readLocalRadCoeff(data); 
  
  CFLog(INFO,"ParadeLibrary::computeProperties() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ParadeLibrary::updateWavRange(CFuint iWavRange)
{
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  if (PE::GetPE().GetRank(nsp) == 0) {
    // all the following can fail if the format of the file parade.con changes
    CFreal minWav = getMinWavelength();
    CFreal maxWav = getMaxWavelength();
    if (getWavLoopSize() > 1) {
      const CFuint wavStride = getWavelengthStride();
      minWav = getMinWavelength() + iWavRange*wavStride;
      maxWav = minWav + wavStride;
    }  
    
    cout << "ParadeLibrary: computing wavelength range [" << minWav << ", " << maxWav << "]\n"; 
    
    // back up the last parade.con file
    std::string command = "cp parade.con parade.con.bkp";
    Common::OSystem::getInstance().executeCommand(command); 
    
    // read from the original parade.con file 
    ifstream fin("parade.con.bkp");
    ofstream fout("parade.con");
    
    string line;
    while (getline(fin,line)) {
      // detect and store the position in the input file where cell data will be inserted
      bool isModified = false;
      
      // update the min wavelength
      size_t pos1 = line.find("wavlo");
      if (pos1 != string::npos) {
	fout << minWav << " wavlo [A]" << endl;
	isModified = true;
      }
      
      // update the max wavelength
      size_t pos2 = line.find("wavhi");
      if (pos2 != string::npos) {
	fout << maxWav << " wavhi [A] " << endl;
	isModified = true;
      }
      
      // update the number of spectral points (= number of wavelengths considered)
      size_t pos3 = line.find("npoints");
      if (pos3 != string::npos) {
	fout << m_wavStride << " npoints" << endl;
	isModified = true;
      }
      
      if (!isModified) {
	fout << line << endl;
      }
    }
    fin.close();
    fout.close();
  }
  PE::GetPE().setBarrier(nsp);
  
  // each processor copies the newly updated parade.con to its own directory
  std::string command   = "cp parade.con " + m_paradeDir;
  Common::OSystem::getInstance().executeCommand(command);
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ParadeLibrary::writeLocalData(Framework::ProxyDofIterator<CFreal>* pstates)
{ 
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbPoints = pstates->getSize();
  
  // write the mesh file
  ofstream& foutG = m_outFileHandle->open(m_gridFile);
  foutG << "TINA" << endl;
  foutG << 1 << " " << nbPoints << endl;
  for (CFuint i =0; i < nbPoints; ++i) {
    foutG.precision(14);
    foutG.setf(ios::scientific,ios::floatfield);
    CFreal *const node = pstates->getNode(i);
    
    if (dim == DIM_1D) {
      foutG << node[XX] << " " << 0.0 << " " << 0.0  << endl;
    }
    if (dim == DIM_2D) {
      foutG << node[XX] << " " << node[YY] << " " << 0.0  << endl;
    }
    if (dim == DIM_3D) {
      foutG << node << endl;
    }
  } 
  foutG.close();
  
  // write the temperatures
  ofstream& foutT = m_outFileHandle->open(m_tempFile);
  const CFuint nbTemperatures = 1 + m_library->getNbTempVib() + m_library->getNbTe();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq(); 
  const CFuint tempID = nbEqs - nbTemperatures; 
  foutT << 1 << " " << nbPoints << " " << nbTemperatures << endl;  
  // here it is assumed that the temperatures are the LAST variables 
  for (CFuint i =0; i < nbPoints; ++i) {
    CFreal *const currState = pstates->getState(i);
    for (CFuint t = 0; t < nbTemperatures; ++t) {
      foutT.precision(14);
      foutT.setf(ios::scientific,ios::floatfield);
      const CFreal temp = std::max(currState[tempID + t], m_tminFix);
      foutT << temp << " ";
    }
    foutT << endl;
  }
  foutT.close();
  
  // write the number densities
  ofstream& foutD = m_outFileHandle->open(m_densFile);
  const CFuint nbSpecies = m_library->getNbSpecies();
  foutD << 1 << " " << nbPoints << " " << nbSpecies << endl;  
  // here it is assumed that the species densities are the FIRST variables 
  for (CFuint i =0; i < nbPoints; ++i) {
    CFreal *const currState = pstates->getState(i);
    for (CFuint t = 0; t < nbSpecies; ++t) {
      foutD.precision(14);
      foutD.setf(ios::scientific,ios::floatfield);
      // number Density = partial density/ molar mass * Avogadro number
      const CFreal nb = std::max(currState[t]*m_avogadroOvMM[t],m_ndminFix);
      foutD << nb << " ";
    }
    foutD << endl;
  }
  foutD.close();
}   
      
//////////////////////////////////////////////////////////////////////////////
  
void ParadeLibrary::readLocalRadCoeff(RealMatrix& data)
{
  // DEBUG_CONDITIONAL_PAUSE<__LINE__>();
  
  fstream& fin = m_inFileHandle->openBinary(m_radFile);
  //string nam = "file-" + StringOps::to_str(PE::GetPE().GetRank());
  //ofstream fout(nam.c_str());
  
  int one = 0;
  fin.read((char*)&one, sizeof(int));
  //fout << "one = "<< one;
  cf_assert(one == 1);
  
  int nbCells = 0;
  fin.read((char*)&nbCells, sizeof(int));
  //fout << " nbCells = " << nbCells << endl;
  cf_assert(nbCells == (int)data.nbRows());
  
  vector<int> wavptsmx(3);
  fin.read((char*)&wavptsmx[0], 3*sizeof(int));
  //fout << "wav3 =  "<< wavptsmx[0] << " " << wavptsmx[1] << " " << wavptsmx[2]<< endl;
  cf_assert(wavptsmx[0] == (int)m_wavStride);
  
  double etot = 0.;
  int wavpts = 0;
  const CFuint sizeCoeff = m_wavStride*3;
  for (int iPoint = 0 ; iPoint < nbCells; ++iPoint) {
    fin.read((char*)&etot, sizeof(double));
    //fout << "etot = " << etot << endl;
    fin.read((char*)&wavpts, sizeof(int));
    cf_assert(wavpts == (int)m_wavStride);
    //fout << "wavpt = " <<  wavpts << endl;
    // this reads [wavelength, emission, absorption] for each wavelength
    fin.read((char*)&data[iPoint*sizeCoeff], sizeCoeff*sizeof(double));
  }
  fin.close();
  //fout.close();
  
  // if emission dominated (optically thin), assign the absorption coefficient to 0
  if (flag_Em) {
    for (int iPoint = 0 ; iPoint < nbCells; ++iPoint) {
      for (int iWav = 0 ; iWav < m_wavStride; ++iWav) {
	data(iPoint, iWav*3+2) = 0.;
      }
    }
  }  
  
  // if (PE::GetPE().GetRank() == 0) {
//     ofstream fout1("inwav.txt"); 
//     fout1 << "TITLE = Original spectrum of radiative properties\n";
//     fout1 << "VARIABLES = Wavl EmCoef AbCoef \n";
//     for (CFuint i = 0; i < data.nbCols()/3; ++i) {
//       fout1 << data(0,i*3) << " " << data(0,i*3+1) << " " << data(0,i*3+2) << endl;
//     }
//     fout1.close();
//   }
  
//   // Convert spectrum to frequency space
//   const CFreal c = PhysicalConsts::LightSpeed();
//   const CFuint stride = data.nbCols()/3;
//   for (CFuint iPoint = 0 ; iPoint < nbCells; ++iPoint) {
//     for (CFuint iw = 0; iw < stride; ++iw) {
//       const CFreal lambda = data(iPoint, iw*3)*1.e-10; // conversion from [A] to [m]
//       // convert emission from wavelength spectrum to frequency spectrum
//       data(iPoint, iw*3+1) *= c/(lambda*lambda);
//       // override wavelength with frequency (nu = c/lambda)
//       data(iPoint, iw*3) = c/lambda; 
//     }
//   }
  
//   if (PE::GetPE().GetRank() == 0) {
//     ofstream fout1("infreq.txt"); 
//     fout1 << "TITLE = Original spectrum of radiative properties\n";
//     fout1 << "VARIABLES = Freq EmCoef AbCoef \n";
//     for (CFuint i = 0; i < data.nbCols()/3; ++i) {
//       fout1 << data(0,i*3) << " " << data(0,i*3+1) << " " << data(0,i*3+2) << endl;
//     }
//     fout1.close();
//   }
}

//////////////////////////////////////////////////////////////////////////////

CFuint ParadeLibrary::getWavLoopSize() const
{
  if (m_spectrumSize == m_wavStride) {
    return 1;
  }
  return (m_wavMax - m_wavMin)/m_wavStride;
}
      
//////////////////////////////////////////////////////////////////////////////
   
} // namespace Parade

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#ifndef COOLFluiD_RadiativeTransfer_RadiativeTransferFVDOM_hh
#define COOLFluiD_RadiativeTransfer_RadiativeTransferFVDOM_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Environment {
    class FileHandlerInput;
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace RadiativeTransfer {
    class RadiationPhysicsHandler;
    
//////////////////////////////////////////////////////////////////////////

/**
 * This class compute the radiative heat transfer using a Finite Volume algorithm
 *
 * @author Alejandro Alvarez (C++ version of Alan Wray's algorithm)
 * @author Andrea Lani (parallelization)
 * @author Javier Martinez
 * @author Claudio Lisco
 */
class RadiativeTransferFVDOM : public Framework::DataProcessingCom {
public:
  
  class Interpolator {
  public:
    /// constructor
    HOST_DEVICE Interpolator() {}
    
    /// destructor
    HOST_DEVICE ~Interpolator() {}
    
    /// Interpolates the values of the opacity tables
    HOST_DEVICE inline void tableInterpolate
    (const CFuint nbBins, const CFuint nbTemp,const CFuint nbPress, 
     const CFreal* Ttable, const CFreal* Ptable, const CFreal* opacities, 
     const CFreal* radSource, CFreal T, CFreal p, CFuint ib, 
     CFreal& val1, CFreal& val2);
  };
  
  /**
   * Constructor
   */
  RadiativeTransferFVDOM(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RadiativeTransferFVDOM();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);  

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();
  
  /**
   * Execute on a set of dofs
   */
  virtual void execute();
  
  /**
   * Compute the directions depending on the option selected 
   */  
  void getDirections();
  
  /**
   * Compute the advance order depending on the option selected 
   */  
  void getAdvanceOrder(const CFuint d, Framework::LocalArray<CFint>::TYPE& advanceOrder); 
  
  /**
   * Compute the advance order depending on the option selected 
   */  
  void getFieldOpacities(CFuint ib);
    
  /**
   * Compute the advance order depending on the option selected 
   */  
  void getFieldOpacities(const CFuint ib, const CFuint iCell);
  
  /**
   * Reads the binary file containing the opacities as function of the temperature, pressure
   * and wavelength  
   */
  void readOpacities();
  
  /**
   * Compute the order of the cells depending on the option selected (For the moment in the execute())
   */
  void findOrderOfAdvance();
  
 
  /**
   * Writes radial q and divQ to a file for the Sphere case 
   * (only if option RadialData is enabled)
   * 
   */
  void writeRadialData();
  
  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected: //function
  
  /// set the normal corresponding to the given face ID
  void setFaceNormal(const CFuint faceID, const CFuint elemID) 
  {
    Framework::DataHandle<CFreal> normals = socket_normals.getDataHandle();
    const CFuint startID = faceID*3;
    const CFreal factor = (static_cast<CFuint>
			   (socket_isOutward.getDataHandle()[faceID]) != elemID) ? -1. : 1.;
    for (CFuint dir = 0; dir < 3; ++dir) {
      m_normal[dir] = normals[startID+dir]*factor;
    }    
  }
  
  /// get the inner product normal*direction (the normal includes the area)
  CFreal getDirDotNA(const CFuint d, const RealVector& normal) const
  {
    return normal[XX]*m_dirs(d,XX) + normal[YY]*m_dirs(d,YY) + normal[ZZ]*m_dirs(d,ZZ);
  }	  
  
  /// compute the radiative heat flux
  /// @param ib  ID of the bin
  /// @param d   ID of the direction
  void computeQ(const CFuint ib, const CFuint d);
  
  /// diagnose problem when advance order algorithm fails
  void diagnoseProblem(const CFuint d, const CFuint m, const CFuint mLast);
  
  /// compute the dot products direction*normal for each face (sign will be adjusted on-the-fly) 
  void computeDotProdInFace(const CFuint d, 
			    Framework::LocalArray<CFreal>::TYPE& dotProdInFace);
  
  /// get the neighbor cell ID to the given face and cell
  CFuint getNeighborCellID(const CFuint faceID, const CFuint cellID) const 
  {
    // find the TRS to which the current face belong
    const Framework::TopologicalRegionSet& faceTrs = *m_mapGeoToTrs->getTrs(faceID);
    // find the local index for such a face within the corresponding TRS
    const CFuint faceIdx = m_mapGeoToTrs->getIdxInTrs(faceID);
    // first neighbor of the current face
    const CFuint sID0 = faceTrs.getStateID(faceIdx, 0);
    return (cellID == sID0) ? faceTrs.getStateID(faceIdx, 1) : sID0;
  }
    
  /// Compute radiative fluxes by looping over bins
  virtual void loopOverBins(const CFuint startBin, 
			    const CFuint endBin, 
			    const CFuint startDir,
			    const CFuint endDir);
  
  /// Compute radiative fluxes by looping over directions
  virtual void loopOverDirs(const CFuint startBin, 
			    const CFuint endBin, 
			    const CFuint startDir,
			    const CFuint endDir);
  
  /// @return the wall face ID to be used inside the qradFluxWall socket
  /// @post returns -1 if the face is not a wall face
  CFint getWallFaceID(const CFuint faceID)
  {
    if (m_isWallFace[faceID]) {
      // find the local index for such a face within the corresponding TRS
      const CFuint faceIdx = m_mapGeoToTrs->getIdxInTrs(faceID);
      const std::string wallTRSName = m_mapGeoToTrs->getTrs(faceID)->getName();
      return m_mapWallTRSToOffsetFaceID.find(wallTRSName)->second + faceIdx;
    }
    return -1;
  }


protected: //data
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of nstates (states in nodes)
  Framework::DataSocketSink<RealVector> socket_nstates;  
  
  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;  
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;  
  
  /// storage of isOutward
  Framework::DataSocketSink<CFint> socket_isOutward; 
  
  /// storage of normals
  Framework::DataSocketSink<CFreal> socket_normals; 
  
  /// storage of face centers
  Framework::DataSocketSink<CFreal> socket_faceCenters; 
  
  /// storage of the the stage of the order of advance
  Framework::DataSocketSource <CFreal> socket_CellID;
  
  /// storage of the divq 
  Framework::DataSocketSource <CFreal> socket_divq;
  
  /// storage of the qx 
  Framework::DataSocketSource <CFreal> socket_qx;
  
  /// storage of the qy
  Framework::DataSocketSource <CFreal> socket_qy;
  
  /// storage of the qz
  Framework::DataSocketSource <CFreal> socket_qz;
  
  /// storage of the radiative source table. Source[ib](it,ip)
  Framework::DataSocketSource <CFreal> socket_TempProfile; 

  /// storage of the binned opacity
  Framework::DataSocketSource <CFreal> socket_alpha_avbin;
  
  /// storage of the binned radiative source
  Framework::DataSocketSource <CFreal> socket_B_bin;
  
  /// the socket to the radiative heat flux at the wall faces
  Framework::DataSocketSource < CFreal > socket_qradFluxWall;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library; 
  
  /// pointer to the radiation library interface
  Common::SharedPtr<RadiationPhysicsHandler> m_radiation;
  
  /// map faces to corresponding TRS and index inside that TRS
  Common::SafePtr<Framework::MapGeoToTrsAndIdx> m_mapGeoToTrs;
  
  /// flag array telling whether a face is on the wall
  std::vector<bool> m_isWallFace;
  
  /// map with TRS name as key and the offset for the wall face IDs as value
  std::map<std::string, CFuint> m_mapWallTRSToOffsetFaceID;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  
  /// wall face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_wallFaceBuilder;
  
  /// temporary normal to the face
  RealVector m_normal; 
  
  /// Weights for the directions
  RealVector m_weight;
  
  /// Field source of oppacity table 
  Framework::LocalArray<CFreal>::TYPE m_fieldSource;
  
  /// Field Absorption of oppacity table used if exponential Method
  Framework::LocalArray<CFreal>::TYPE m_fieldAbsor;
  
  /// Field Absorption of oppacity table used if not Exponential Method
  Framework::LocalArray<CFreal>::TYPE m_fieldAbSrcV;
  
  /// Field Absorption of oppacity table used if not Exponential Method
  Framework::LocalArray<CFreal>::TYPE m_fieldAbV;  
  
  /// Exponent for the radiation of oppacity table
  Framework::LocalArray<CFreal>::TYPE m_In;
  
  /// Exponent for the radiation of oppacity table
  Framework::LocalArray<CFreal>::TYPE m_II;
  
  /// Opacities read from table. Stored as follows:
  /// kappa(p0,b0,T0),kappa(p0,b0,T1),...,kappa(p0,b0,Tn),
  /// kappa(p0,b1,T0),kappa(p0,b1,T1),...,kappa(p0,b1,Tn), .... 
  /// kappa(p0,bn,T0),kappa(p0,bn,T1),...,kappa(p0,bn,Tn), ....
  /// kappa(p1,b0,T0),kappa(p1,b0,T1),...,kappa(p1,b0,Tn), ...
  Framework::LocalArray<CFreal>::TYPE m_opacities;
  
  /// Radiative source reaf from table and stored the same way as the opacities
  Framework::LocalArray<CFreal>::TYPE m_radSource;
  
  /// storage of the temperatures of the opacity table. 
  Framework::LocalArray<CFreal>::TYPE m_Ttable;
  
  /// storage of the pressure of theopacity table
  Framework::LocalArray<CFreal>::TYPE m_Ptable; 
  
  /// storage of the dot products per face
  Framework::LocalArray<CFreal>::TYPE m_dotProdInFace; 
  
  /// Done status of a cell in a given direction at the end of a stage
  std::vector<bool> m_sdone;
    
  /// temporary list of cell indexes to be processed
  std::vector<CFuint> m_cdoneIdx;
  
  /// names of the TRSs of type "Wall"
  std::vector<std::string> m_wallTrsNames;
  
  /// Directions 
  RealMatrix m_dirs;  
  
  /// Create the array advance_order, which contains, for each direction (1st index), the list of cells to be advanced.  The list is divided into "stages" that consist of 
  /// a set of cells that can be advanced in parallel.These sets are terminated by a negative entry.  E.g., if there are 8 cells and advance_order(1:8,1) = (1,2,-5,3,4,-8,6,-7), 
  /// then cells (1,2,5) can be done first, and are in fact the boundary cells for direction 1,
  /// then cells (3,4,8) can be done; finally cells (6,7) can be done to complete the sweep in direction 1.
  std::vector<Framework::LocalArray<CFint>::TYPE > m_advanceOrder;   
  
  /// Radiative flux
  RealMatrix m_q;
  
  /// Radial average of q vector for a Sphere
  RealVector m_qrAv;
  
  /// Radial average of divQ for a Sphere
  RealVector m_divqAv;  
  
  /// Number of directions types
  CFuint m_nbDirTypes;   
  
  /// user defined number of directions
  CFuint m_nbDirs;
  
  /// name of the temporary local directory where Parade is run
  boost::filesystem::path m_dirName;
  
  /// name of the .dat binary table with the opacities
  std::string m_binTabName;
  
  /// File where the table is written
  std::string m_outTabName;
  
  /// bool to write the table in a file
  bool m_writeToFile;
  
  /// Use exponential method
  bool m_useExponentialMethod;
  
  /// option to print the radial q and divQ for the Sphere
  bool m_radialData;
  
  /// old algorithm just kept for comparison purposes
  bool m_oldAlgo;
  
  /// user defined number of points in the radial direction
  CFuint m_Nr;
  
  /// constant pressure
  CFreal m_constantP;
  
  /// minimum temparature
  CFreal m_Tmin;
  
  /// maximum temperature
  CFreal m_Tmax;
  
  /// temperature interval
  CFreal m_deltaT;
  
  /// ID of pressure in the state vector
  CFuint m_PID;
  
  /// ID of temperature in the state vector
  CFuint m_TID;
  
  /// input file handle
  Common::SelfRegistPtr<Environment::FileHandlerInput> m_inFileHandle;
  
  /// opacities file
  boost::filesystem::path m_binTableFile;
  
  /// start/end bin to consider
  std::pair<CFuint, CFuint> m_startEndBin;
  
  /// start/end direction to consider
  std::pair<CFuint, CFuint> m_startEndDir;
  
  /// number of bins
  CFuint m_nbBins;
  
  /// number of Temperatures
  CFuint m_nbTemp;
  
  /// number of pressures
  CFuint m_nbPress;
  
  /// total number of threads/CPUs in which the algorithm has to be split
  CFuint m_nbThreads;
  
  /// ID of the thread/CPU within the parallel algorithm
  CFuint m_threadID;
  
  /// flag telling to run a loop over bins and then over directions (or the opposite)
  bool m_loopOverBins;
  
  /// flag telling to run without solving anything, just for testing
  bool m_emptyRun;
  
  /// flag telling to read the opacity tables (default behaviour)
  bool m_readOpacityTables;
  
}; // end of class RadiativeTransferFVDOM
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::Interpolator::tableInterpolate
(const CFuint nbBins, const CFuint nbTemp, const CFuint nbPress, 
 const CFreal* Ttable, const CFreal* Ptable, const CFreal* opacities, 
 const CFreal* radSource, CFreal T, CFreal p, CFuint ib, CFreal& val1, CFreal& val2)
{
  //Find the lower bound fo the temperature and the pressure ranges
  //we assume that the temperature and pressure always fall in the bounds.
  //If they don't then the value are still interpolated from the nearest
  //two points in the temperature or pressure list
  CFuint it = nbTemp - 2;
  for (CFuint i = 1; i < (nbTemp - 2); i++){
    if(Ttable[i] > T) { it = i - 1; break;}
  }
  
  CFuint ip = nbPress - 2;
  for (CFuint i = 1; i < (nbPress - 2); i++){
    if(Ptable[i] > p) { ip = i - 1; break;}
  }
  
  //Linear interpolation for the pressure
  
  const CFuint iPiBiT           = it + ib*nbTemp + ip*nbBins*nbTemp;
  const CFuint iPplus1iBiT      = it + ib*nbTemp + (ip + 1)*nbBins*nbTemp;
  const CFuint iPiBiTplus1      = (it + 1) + ib*nbTemp + ip*nbBins*nbTemp;
  const CFuint iPplus1iBiTplus1 = (it + 1) + ib*nbTemp + (ip + 1)*nbBins*nbTemp;
  
  // Linear interpolation for the pressure
  // Interpolation of the opacities
  const CFreal bt1op = (opacities[iPplus1iBiT] - opacities[iPiBiT])*
		    (p - Ptable[ip])/(Ptable[ip + 1] - Ptable[ip]) + opacities[iPiBiT];
  
  const CFreal bt2op = (opacities[iPplus1iBiTplus1] - opacities[iPiBiTplus1])*
		    (p - Ptable[ip])/(Ptable[ip + 1] - Ptable[ip]) + opacities[iPiBiTplus1];
  
  // Interpolation of the source
  const CFreal bt1so = (radSource[iPplus1iBiT] - radSource[iPiBiT])*
		    (p - Ptable[ip])/(Ptable[ip + 1] - Ptable[ip]) + radSource[iPiBiT];
  
  const CFreal bt2so = (radSource[iPplus1iBiTplus1] - radSource[iPiBiTplus1])*
		    (p - Ptable[ip])/(Ptable[ip + 1] - Ptable[ip]) + radSource[iPiBiTplus1];    
  
  // Logarithmic interpolation for the temperature
  // Protect against log(0) and x/0 by switching to linear interpolation if either
  // bt1 or bt2 == 0.  (Note we can't allow log of negative numbers either)
  // Interpolation of the opacities   
  if(bt1op <= 0 || bt2op <= 0){
    val1 = (bt2op - bt1op)*(T - Ttable[it])/(Ttable[it + 1] - Ttable[it]) + bt1op;
//    cout <<"\nOption1 \n";
//    cout <<"T = "<< T <<"\tTi+1 = "<<Ttable[it + 1]<<"\tTi = "<<Ttable[it] <<"\n";
//    cout <<"val1 = " << val1 <<"\tbt2op ="<< bt2op <<"\tbt1op ="<< bt1op <<"\n";
  }
  else {
    val1 = std::exp((T - Ttable[it])/(Ttable[it + 1] - Ttable[it])*std::log(bt2op/bt1op))*bt1op;
//     cout <<"\nOption2 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<Ttable[it + 1]<<"\tTi = "<<Ttable[it] <<"\n";
//     cout <<"val1 = " << val1 <<"\tbt2op ="<< bt2op <<"\tbt1op ="<< bt1op <<"\n";
  }
  // Interpolation of the source
  if(bt1so <= 0 || bt2so <= 0){
    val2 = (bt2so - bt1so)*(T - Ttable[it])/(Ttable[it + 1] - Ttable[it]) + bt1so;
//     cout <<"\nOption3 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<Ttable[it + 1]<<"\tTi = "<<Ttable[it] <<"\n";
//     cout <<"val1 = " << val2 <<"\tbt2so ="<< bt2so <<"\tbt1so ="<< bt1so <<"\n";
  }
  else {
    val2 = std::exp((T - Ttable[it])/(Ttable[it + 1] - Ttable[it])*std::log(bt2so/bt1so))*bt1so;
//     cout <<"\nOption3 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<Ttable[it + 1]<<"\tTi = "<<Ttable[it] <<"\n";
//     cout <<"val2 = " << val2 <<"\tbt2so ="<< bt2so <<"\tbt1so ="<< bt1so <<"\n";
  }
  
  //cf_assert(ib == 0);
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_RadiativeTransferFVDOM_hh




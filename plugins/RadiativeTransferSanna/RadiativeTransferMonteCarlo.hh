#ifndef COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloFVMCC_hh
#define COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include <numeric>
#include <boost/random.hpp>

#include "Framework/DataProcessingData.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Config/ConfigObject.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "RadiativeTransferSanna/RandomNumberGenerator.hh"
//#include "RadiativeTransfer/MersenneTwister.hh"
#include "RadiativeTransferSanna/ParticleTracking.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Framework {
    class RadiationLibrary;
  }
  
  namespace Numerics {
    namespace FiniteVolume {
      class CellCenterFVMData;
    }
  }
  
 namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////


   
//////////////////////////////////////////////////////////////////////////////

  /**
   * This struct serves to comunicate to the entity identified with ID and type Type the radiative heat flux entering into itself
   *
   * @author Alessandro Sanna
   *
   */
struct Qout
{

  CFuint ID;
  CFint Type; // before: 1 if it is a cell, 0 if it is a face 
              // now: -1 for a cell; globalID > 0 for a face
  CFreal Q;

};

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class compute the radiativeHeatSource for each cell
   *
   * @author Alessandro Sanna
   *
   */
class RadiativeTransferMonteCarlo : public Framework::DataProcessingCom
                          
{
public:
  
  /// enumerators that define the entity type
  enum EntityType {DISAPPEARED =0, WALL_FACE =1, INTERNAL_CELL =2, PARTITION_GHOSTCELL =3, NEGLIGIBLE=4};
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
   static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
   RadiativeTransferMonteCarlo(const std::string& name);

  /**
   * Default destructor
   */
   ~RadiativeTransferMonteCarlo();

  /**
   * Configures the command.
   */
   void configure(Config::ConfigArgs& args);
  
  /**
   * PreProcess phase
   */
   void setup();

  /**
   * executeOnTrs
   */
   void executeOnTrs();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
   std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

private:


  RealMatrix m_wavReduced,m_emReduced,m_amReduced;


  /**
   * boundary knowledge
   */
   void BoundaryKnowledge();

  /**
   * MonteCarlo
   */
   void MonteCarlo();
  
  /**
   * ray tracing
   */
  CFuint rayTracing(Ray& beam, CFuint begin, EntityType& entity, const std::string& trs);
  

  inline CFreal linearInterpol(CFreal x0,CFreal y0, CFreal x1, CFreal y1, CFreal x){
    return  y0+(y1-y0)*(x-x0)/(x1-x0);
  }


  inline CFreal trapezoid(std::vector<CFreal>& x, std::vector<CFreal>& y){
    //using namespace std;
    cf_assert(x.size()==y.size());
    cf_assert(x.size()>1);
    CFreal sum=0.;
    //cout<<x[0]<<' '<<y[0]<<' ';
    for(CFuint i=0; i<x.size()-1; ++i){
       //cout<<x[i+1]<<' '<<y[i+1]<<' ' ;
       sum+= ( y[i+1]+y[i] )/2.*( x[i+1]-x[i] );

    }
    //cout<<endl;
    return sum;
  }



  /**
   * get the absorbtion coefficient of cell with cellID
   */
  CFreal getK(CFuint cellID, Ray ray)
  {

   //return 0.;

   const CFuint localID = _Global2LocalInCellIdMap.find(cellID);

   return m_amReduced(localID,m_spectralIdx) ;

   //const CFreal rayWav=ray.wavelength;

   //typeAbsMap::iterator it=absorptionCoeff2[localID]->upper_bound(rayWav);
   //ypeAbsMap::iterator it2=(--it);
   //return (it!=absorptionCoeff2[localID]->end() && rayWav<=it2->second.first)? it2->second.second : 0.;


//   CFreal absMinWav=0.,absMaxWav=0.;

//   for(CFuint wavIdx=0; rayWav>=absMaxWav-1e-2 || wavIdx<absorptionCoeff.nbCols(); wavIdx+=3){
//     absMinWav=absorptionCoeff(localID, wavIdx);
//     absMaxWav=absorptionCoeff(localID, wavIdx+1);
//     if(rayWav>= absMinWav-1e-2 && rayWav<=absMaxWav+1e-2){
//       return absorptionCoeff(localID, wavIdx+2);
//     }
//   }

//   return 0;

    // CFreal r = _CellIDmap.find(cellID);
    // CFreal absorption = 1;//(0.9 - 0.3*r);
    //if (m_radCoeffReduced.size() > 0) {
    //  // skip the first wavelength which is not averaged on the interval
    //  return m_radCoeffReduced(localID, (m_spectralIdx+1)*3+2);
    //}
    //return m_radCoeff(localID, m_spectralIdx*3+2);
  }
  
  /**
   * get the absorbtion coefficient of wall face with faceID
   */
  inline CFreal getWallK(CFuint faceID) const {return m_wallAbs;}
  
  /**
   * get the emission coefficient of cell with cellID
   */
  inline CFreal getEm(CFuint cellID)
  {
    const CFuint localID = _Global2LocalInCellIdMap.find(cellID);
    //return m_radCoeffReduced(localID, m_spectralIdx*3+1);
    return m_emReduced(localID,m_spectralIdx);
    //CFreal r = _CellIDmap.find(cellID);
    //CFreal emission = 1;//(0.9 - 0.3*r);
    // return emission;
    //if (m_radCoeffReduced.size() > 0) {
      // skip the first wavelength which is not averaged on the interval
    //  return m_radCoeffReduced(localID, (m_spectralIdx+1)*3+1);
    //}
    //return m_radCoeff(localID, m_spectralIdx*3+1);
  }
   
  /**
   * get the wavelength lambda
   */
  inline CFreal getCurrentLambda() const
  {  
    //if (m_radCoeffReduced.size() > 0) {
    //  // skip the first wavelength which is not averaged on the interval
    //  return m_radCoeffReduced(0, (m_spectralIdx+1)*3);
    //}
      return m_wavReduced(0,m_spectralIdx);
    //return m_radCoeffReduced(0, m_spectralIdx*3);
  }
  
  /**
   * get the emission coefficient of wall face with faceID
   */
  inline CFreal getWallEm(CFuint faceID) const {return m_wallEm;}
  
  /**
   * build cell IDs map
   */
   void buildCellIDradiusmap();

  /**
   * build vector of radiative heat source along a single radius in the middle of the cilinder
   */
   void getRHF();

  /**
   * build vector of radiative heat source along a radius computed as average along a circular sections in the middle of the cilinder
   */
   void getRHFcircularSections();

  /**
   * build vector of radiative heat source along a radius computed as average along a circular sections in the middle of the cilinder
   */
   void getRHFcircularSectionsBis();

  /**
   * build vector of radiative heat source along the radius of the axialsymmetric cylinder
   */
   void getRHFaxyCylinder();

  /**
   * build vector of radiative heat source along a radius computed as average along a circular sections in the middle of the cilinder
   */
   void getRHFslab();

  /**
   * build the MPI_struct of the Ray struct
   */
   void build_MPI_strut_Ray(MPI_Datatype* MPI_Ray_ptr);

  /**
   * build the MPI_struct of the Qout struct
   */
   void build_MPI_strut_Qout(MPI_Datatype* MPI_Qout_ptr);

  /**
   * compute Heat Flux
   */
   void computeHeatFlux();
  
  /// class for sorting pairs
  template <typename T1, typename T2>
  class PairSort {
  public:
    bool operator() (const std::pair<T1,T2>& p1, 
		     const std::pair<T1,T2>& p2) const
    {
      return (p1.first < p2.first) ? true : false;
    }
  };
  
  /// class for normalize an array
  void normalizeArray(CFreal* v)
  {
    int dim = _dim;
    //if(_Axi) dim = 3;
    CFreal invNorm2 = 0.0;
    if(dim == 2){
      invNorm2 = 1.0/sqrt(v[0]*v[0] + v[1]*v[1]);
      v[0] *= invNorm2;
      v[1] *= invNorm2;
    }
    if(dim == 3){
      invNorm2 = 1.0/sqrt(v[0]*v[0] + v[1]*v[1] +v[2]*v[2]);
      v[0] *= invNorm2;
      v[1] *= invNorm2;
      v[2] *= invNorm2;
    }
  }

  /// Compute cell rays
  void computeCellRays(std::vector<std::vector<Ray> >& partitionBeams);
  
  /// Compute wall faces rays
  void computeWallFaceRays(std::vector<std::vector<Ray> >& partitionBeams);
  
  /// Compute rays crossing partition faces
  void computePartitionFaceRays(std::vector<std::vector<Ray> >& partitionBeams,
				std::vector<std::vector<CFuint> >& returningPartitionBeams);
  
  /// Compute the final step of the communication
  void computeFinalCommunication(std::vector<std::vector<CFuint> >& returningPartitionBeams);
  
  /// Set the list of boundary face IDs
  void setFaceIDs(const std::vector<std::string>& names, 
		  std::vector<CFuint>& faceIDs, bool isWall);
  
  /// ray tracing for the axisymmetric case
  void axiRayTacing(Ray beam,
		    CFuint absorbingEntityID, 
		    std::vector<CFuint>& absorbingIDs,
                    std::vector<std::vector<Ray> >& partitionBeams);
  
  /// Reduce the spectrum of radiation properties
  void reduceSpectra(const RealMatrix& indata, RealMatrix& outdata);
  
  /// Reduce the spectrum of radiation properties
  void myReduceSpectra(const RealMatrix& indata, RealMatrix& outdata);

  /// Flag all cells in the overlap
  void flagOverlapCells();

  /// Build Global to Local cell ID map
  void buildCellIDmap();

  //RealVector nextAxiPoint(const RealVector old2Daxi, const Ray ray);

  void convert3Dto2DAxi(Ray ray, RealVector axiCoord);

  RealVector covert2DAxito3Dcart(Ray ray, RealVector axiCoord);

  CFreal get3Ddist(RealVector axi1, RealVector axi2, Ray ray);

  RealVector convertCylindricalToCartesian(CFreal z,CFreal r, CFreal tetha);

  //GeoEntOut myAxiRayTracing(Ray ray, CFuint currentCellID);

  RealVector get2DAxicoord(Ray ray);

  RealVector get3Dcoord(Ray ray);

  inline CFreal Plank(CFreal T,CFreal lambda){

    static CFreal c1=1.1911e-16; //  c1 = 2 h c^2 [W m ^2]
    static CFreal c2=1.439e-2 ; //   c2 = h c k_b    [m K]
    // c   : speed of light    [m/s]
    // h   : Plank constant    [J s]
    // k_b : Boltzman constant [J/k]
    return c1 / ( pow(lambda,5) * ( exp( c2/(lambda*T) ) -1 ) );
  }

  inline CFreal intPlankSimpson(CFreal T,CFreal lambda,CFreal interval){
    CFreal a=lambda-interval/2;
    CFreal b=lambda+interval/2;
    return interval/6*(Plank(T,a)+ 4*Plank(T,lambda)+Plank(T,b) );
  }
  inline CFreal intPlankSimple(CFreal T,CFreal lambda,CFreal interval){
    return Plank(T,lambda)*interval;
  }

private: 
  
  /// particle tracking algorithm 
  ParticleTracking _particleTracking;
  
  /// random number generator
  //MersenneTwister m_MTRand;

  /// Socket for the Gas Radiative Heat Source
  Framework::DataSocketSource < CFreal > socket_qrad; //GasRadiativeHeatSource
  
  /// the socket to the radiative heat flux at the wall faces
  Framework::DataSocketSource < CFreal > socket_qradFluxWall;
    
  /// handle to the face normals
  Framework::DataSocketSink< CFreal> socket_normals;

  /// storage of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// storage of the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of the nodal state's
  Framework::DataSocketSink < RealVector > socket_nstates;

  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// storage of the face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// rank of the processor containing the right state for partition faces
  Framework::DataSocketSink<CFuint> socket_rankPartitionFaces;

  /// IDs corresponding to the cell for which the normal point outward
  Framework::DataSocketSink<CFint> socket_isOutward;
    
  /// proxy pointer for the states storage
  std::auto_ptr<Framework::ProxyDofIterator<CFreal> > m_nstatesProxy;
  
  /// mapping between nodeID and stateID
  std::vector<CFuint> m_nodeIdToStateId;
  
  /// array storing absorption and emission coefficients 
  /// wavelength       = m_radCoeff(local state ID, spectral point idx*3)
  /// emission coeff   = m_radCoeff(local state ID, spectral point idx*3+1)
  /// absorption coeff = m_radCoeff(local state ID, spectral point idx*3+2)
  RealMatrix m_radCoeff;
  
  /// array storing absorption and emission coefficients 
  /// wavelength       = m_radCoeffReduced(local state ID, spectral point idx*3)
  /// emission coeff   = m_radCoeffReduced(local state ID, spectral point idx*3+1)
  /// absorption coeff = m_radCoeffReduced(local state ID, spectral point idx*3+2)
  RealMatrix m_radCoeffReduced;
  
  /// Special array for absorption coefficient in reduced spectra
  // minWav           = AbsorptionCoef(local state ID, spectral poin idx*3)
  // maxWav           = AbsorptionCoef(local state ID, spectral poin idx*3+1)
  // absorption coeff = AbsorptionCoef(local state ID, spectral poin idx*3+2)
  RealMatrix absorptionCoeff;

  typedef std::pair<CFreal,CFreal> typeMappedAbs;
  typedef std::map< CFreal,typeMappedAbs > typeAbsMap;
  std::vector<typeAbsMap*> absorptionCoeff2;

  CFreal m_deltaWav;

  /// wavelength intervals
  RealVector m_deltaWavs;
  
  /// array
  RealVector _actualPoint;
  
  /// array
  RealVector _intersectionPointOld;
   
  /// array
  RealVector _internalPoint;
  
  /// array
  RealVector _intersectionPoint;
  
  /// array
  RealVector _reflectedDirection;
  
  /// array
  RealVector _intersection;

  /// array
  RealVector _nn;
  
  /// array
  RealVector m_midPointE;
  
  /// array
  RealVector m_midPointA;
  
  /// array
  RealVector m_diffAE;
  
  /// array
  RealVector m_projARandom;
  
  ///vector that maps each cell with its radiative heat source
  std::vector<CFreal> _RHSGas;

  ///vector that maps each wall face ID with its radiative heat source
  std::vector<CFreal> _RHSWall;
  
  /// map wall local faceIDs to internal cell IDs
  Common::CFMap<CFuint, CFuint> m_mapWallFaceIDtoCellID; 
  
  /// 3D volumes for axisymmetric case
  RealVector m_axiVolumes;
  
  /// 3D surfaces for axisymmetric case
  RealVector m_wallSurfaces;
  
  /// number of dimension
  CFuint _dim;

  /// vector storing the ID of all boundary faces but the wall faces
  std::vector<CFuint> _boundaryIDs;

  /// vector storing the ID of all wall faces
  std::vector<CFuint> _wallIDs; 
    
  /// vector storing the ID of all symmetry faces
  std::vector<CFuint> _symmetryIDs; 
  
  /// vector storing the ID of the ghost cells of the partition faces 
  std::vector<CFuint> _partitionIDs;
  
  ///my process
  CFuint _myP;

  ///number processes
  CFuint _nP;

  ///MPI_COMM_WORLD
  MPI_Comm _comm;

  /// Map that stores:
  /// as key: the globalID of the state of the ghost cells of the partition faces
  /// as value: the rank of the processor in which the ghost cells are updated cells
  Common::CFMap<CFuint,CFuint> _ProcessorOfPartitionGhostCell;

  /// CFmap that maps the globalID with the getID of the cell along the inner layer partition faces
  Common::CFMap<CFuint,CFuint> _Global2LocalInCellIdMap;

  /// Global to Local cell ID map
  Common::CFMap<CFuint,CFuint> _CellIDmap;

  /// total number of wall faces (considering all wall trses)
  CFuint _nWallFaces;
  
  /// index of currently in use spectral point 
  CFuint m_spectralIdx;
  
  /// index for looping over ranges of wavelengths
  CFuint m_iWavRange;
  
  /// index indicating the cell completion
  std::vector<CFuint> m_cellCompletionIdx;
  
  /// index indicating the face completion
  std::vector<CFuint> m_faceCompletionIdx;
  
  /// map storing the ID of the emitting cells
  Common::CFMap<CFuint,CFuint> _EmittingCellIDs;

  /// vector storing the ID of the absorbing cells
  std::vector<CFint> _AbsorbingCellIDs;

  /// vector storing the rank of the processor which the absorbing cells belong to
  std::vector<CFint> _AbsorbingCellRanks;

  /// vector storing true if the absorbing element is a cell and false if it is a wall face
  std::vector<CFint> _AbsorbingCellType;
  
  /// map storing the ID of the emitting faces
  Common::CFMap<CFuint,CFuint> _EmittingWallFaceIDs;
  
  /// vector storing the ID of the absorbing cells
  std::vector<CFint> _AbsorbingWallFaceIDs;
  
  /// vector storing the rank of the processor which the absorbing cells belong to
  std::vector<CFint> _AbsorbingWallFaceRanks; 

  /// vector storing true if the absorbing element is a cell and false if it is a wall face
  std::vector<CFint> _AbsorbingWallFaceType;
  
  // create send, receive vectors
  std::vector<int> m_sendCounts;
  std::vector<int> m_recvCounts;
  std::vector<int> m_sendDispls;
  std::vector<int> m_recvDispls;
  
  /// vector storing true if the element belongs to the overlap region
  std::vector<bool> m_isOverlapCell;
  
  /// vector storing true if the wall face have absorbed something (debugging purposes)
  std::vector<bool> m_isWallFaceAbsorbing;
  
  /// table telling which processors share each overlap node
  std::vector<std::vector<CFuint> > m_overlapCellRanks;

  /// Cell IDs map
  Common::CFMap<CFuint,CFreal> _CellG2LIDmap;
  
  /// cell builder
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_cellBuilder;
  
  /// wall face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_wallFaceBuilder;
  
  /// physical properties library
  Common::SelfRegistPtr<Framework::RadiationLibrary> m_radLibrary;
  
  /// inner TRS boundary names
  std::string m_radLibraryName;
  
  /// names of the wall TRSs
  std::vector<std::string> _wallTrsNames;
  
  /// names of the boundary TRSs
  std::vector<std::string> _boundaryTrsNames;
  
  /// names of the symmetry TRSs
  std::vector<std::string> _symmetryTrsNames;
  
  /// number of rays that each element emits
  CFuint _NRAYS;
  
  /// Number of points
  //CFuint _nbPoints;
  
  /// Identifier for the testcase
  CFuint _testcaseID;
 
  /// maximum number of visited cells
  CFuint _maxVisitedCells; 
  
  /// flag tellinga to freeze the statistical pattern
  CFuint m_freezePattern; 
  
  /// For testcase 3, set the number of cell in y and x direction
  //CFuint _nCx;
  //CFuint _nCy;
  
  CFuint m_nbSteps;

  /// True if it is an axisymmetric simulation
  bool _Axi;
  
  /// Planck function at the wall
  bool m_usePlanck;
  
  /// wall emissivity
  CFreal m_wallEm; 
  
  /// wall absorption
  CFreal m_wallAbs; 
  
  /// Number of points for the virtual discretization in the normal plane for the axisymmetric case
  CFuint m_nbAxiPoints; 
  
  /// Number of points to use for reducing the spectra coming out from the radiation library
  CFuint m_reducedSpectraSize;
  
  /// Temperature ID
  CFuint m_temperatureID;
  
  /// Free stream temperature value
  CFreal m_freeStreamT;

  RandomNumberGenerator m_rand;

  CFreal m_myWav;

}; // end of class RadiativeTransferMonteCarlo

//////////////////////////////////////////////////////////////////////////////
   
  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloFVMCC_hh

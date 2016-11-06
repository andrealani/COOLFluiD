#ifndef COOLFluiD_Numerics_FiniteVolume_HartmannSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_HartmannSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"

#include "Framework/DataSocketSource.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#include "Framework/MathTypes.hh"
#include "Framework/VarSetTransformerT.hh"
#include "FiniteVolume/FluxData.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Source term for MultiFluid considering 2 fluids: plasma + neutrals
 * variables
 *
 * @author Alejandro Alvarez
 *
 */
template <class UPDATEVAR>
class HartmannSourceTerm : public ComputeSourceTermFVMCC {

public:


#ifdef CF_HAVE_CUDA

  /// nested class defining local options
  template <typename P = NOTYPE>
  class DeviceConfigOptions {
  public:
    /// constructor
    HOST_DEVICE DeviceConfigOptions() {}
    
    /// destructor
    HOST_DEVICE ~DeviceConfigOptions() {}
    
    /// initialize with another object of the same kind
    HOST_DEVICE void init(DeviceConfigOptions<P> *const in) 
    {
      isResistive = in->isResistive;
      electricalConductivity = in->electricalConductivity;
      //parameters needed
    }
    
    //parameters needed
    CFreal electricalConductivity;
    bool isResistive;
  };
  
  /// nested class defining a functor
  template <DeviceType DT, typename VS>
  class DeviceFunc {
  public:
    typedef VS MODEL;
    typedef HartmannSourceTerm BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
    
    /// Compute the source term : implementation
    HOST_DEVICE void operator()(CFreal* state, VS* model, CFreal* source); 
   

  private:
    DeviceConfigOptions<NOTYPE>* m_dco;

    
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_physicalData;
    typename MathTypes<CFreal, DT, 6>::VEC d_NonInducedEMField;
    typename MathTypes<CFreal, DT, 3>::VEC d_Btotal;
    typename MathTypes<CFreal, DT, 3>::VEC d_Etotal;
    

  };
  
  /// copy the local configuration options to the device
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {

    CFLog(NOTICE, "HartmannSourceTerm::copyConfigOptionsToDevice START1 \n \n");
    
    //CopyParameters
    CFreal electricalConductivity = getElectricalConductivity();
    bool isResistive = getIsResistive();

    CFLog(NOTICE, "HartmannSourceTerm::DeviceConfigOptions electricalConductivity = " << _electricalConductivity  << "\n");
    CFLog(NOTICE, "HartmannSourceTerm::DeviceConfigOptions isResitive = " << _isResistive  << "\n");
    CudaEnv::copyHost2Dev(&dco->electricalConductivity, &electricalConductivity, 1);
    CudaEnv::copyHost2Dev(&dco->isResistive, &isResistive, 1);


    CFLog(NOTICE, "HartmannSourceTerm::copyConfigOptionsToDevice END \n \n");
  }  
  
  /// copy the local configuration options to the device
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {

    CFreal electricalConductivity = getElectricalConductivity();
    bool isResistive = getIsResistive();

    CFLog(VERBOSE, "HartmannSourceTerm::copyConfigOptions \n");
    // consider to copy to constant memory
    dco->electricalConductivity = electricalConductivity;
    dco->isResistive = isResistive;

  } 



#endif



  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  HartmannSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~HartmannSourceTerm();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
    _sockets.template createSocketSink<RealVector>("nstates");
  }
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian);

  /**
   * Returns the DataSocket's that this command provides as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  CFreal getElectricalConductivity(){return _electricalConductivity;}
  bool getIsResistive(){return _isResistive;}
    
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;
  
  /// handle to outward normal
  Framework::DataHandle<CFint> _isOutward;
  
  /// socket for storing the divergence of magnetic field
  Framework::DataSocketSource<CFreal> socket_divB;
  
  /// socket for storing the Current
  Framework::DataSocketSource<CFreal> socket_Current;

  /// socket for storing the Bx potential
  Framework::DataSocketSource<CFreal> socket_BxPotential;
 
  /// socket for storing the By potential
  Framework::DataSocketSource<CFreal> socket_ByPotential;

  /// socket for storing the Bz potential
  Framework::DataSocketSource<CFreal> socket_BzPotential;

  /// array to store the mass fractions
  RealVector _ys;
  
  /// Euler physical data
  RealVector _physicalData;

  /// vector to store temporary result
  RealVector _temp;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;
  
  /// array of temporary values
  RealMatrix _values;

  ///Non Induced Part of the electrocmagnetic Field
  RealVector _NonInducedEMField;
  
  /// Current density vector
  RealVector _J;
  
  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients;
  
private:

  /// Electrical conductivity
  CFreal _electricalConductivity;

  /// Orszag-Tang Conductivity Flag
  bool _isResistive;

  /// Left State to compute the Orszag Conductivity
  RealVector _dataLeftState;

   /// Right State to compute the Orszag Conductivity
  RealVector _dataRightState;

  /// gradient of Bx
  RealVector _gradBx;

  /// gradient of By
  RealVector _gradBy;

 /// gradient of Bz
 RealVector _gradBz;

}; // end of class HartmannSourceTerm

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA

template <class UPDATEVAR>
template <DeviceType DT, typename VS>
HOST_DEVICE void HartmannSourceTerm<UPDATEVAR>::DeviceFunc<DT, VS>::operator()(CFreal* state, VS* model, CFreal* source) 
                                                        
{


  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
 
/* NOT IMPLEMENTED GPU 
  const vector<State*>* const states = element->getStates();
  const CFuint elemID = element->getID();
  DataHandle<CFreal> divB = socket_divB.getDataHandle();
  DataHandle<CFreal> Current = socket_Current.getDataHandle();
  DataHandle<CFreal> BxPotential = socket_BxPotential.getDataHandle();
  DataHandle<CFreal> ByPotential = socket_ByPotential.getDataHandle();
  DataHandle<CFreal> BzPotential = socket_BzPotential.getDataHandle();
*/
  
  typename VS::UPDATE_VS* updateVS = model->getUpdateVS();
  updateVS->computePhysicalData(&state[0], &d_physicalData[0]); 

  const CFuint dim = updateVS->getDim();
  const bool is2DHalf = updateVS->getIs2DHalf();
  printf("dim: %d \t is2DHalf: %d \n", dim, is2DHalf);
 // cf_assert(states->size() == 1);
  
  const CFuint firstVelocity = updateVS->getFirstVelocity(); 
 
 // State *const currState = (*states)[0]; 
 // _varSet->computePhysicalData(*(*states)[0], _physicalData);
  /*
  if(is2DHalf || PhysicalModelStack::getActive()->getDim() == 2)
  {
    _NonInducedEMField = _varSet->getModel()->getNonInducedEMField
      (currState->getCoordinates()[XX], currState->getCoordinates()[YY],0.); //The third component is null
  }
  else
  {
    _NonInducedEMField = _varSet->getModel()->getNonInducedEMField
      (currState->getCoordinates()[XX], currState->getCoordinates()[YY], currState->getCoordinates()[ZZ]); 
  }
*/

  for (CFuint i=0; i<6; i++){
    d_NonInducedEMField[i] = updateVS->getNonInducedEMField()[i]; 
  }

  d_Btotal[0] = d_physicalData[UPDATEVAR::PTERM::BX] + d_NonInducedEMField[0];
  d_Btotal[1] = d_physicalData[UPDATEVAR::PTERM::BY] + d_NonInducedEMField[1];
  d_Btotal[2] = d_physicalData[UPDATEVAR::PTERM::BZ] + d_NonInducedEMField[2];
  d_Etotal[0] = d_physicalData[UPDATEVAR::PTERM::EX] + d_NonInducedEMField[3];
  d_Etotal[1] = d_physicalData[UPDATEVAR::PTERM::EY] + d_NonInducedEMField[4];
  d_Etotal[2] = d_physicalData[UPDATEVAR::PTERM::EZ] + d_NonInducedEMField[5];

  //cout <<"NonInduced EM Field = "<< _NonInducedEMField << endl;

  //   RealVector& refData = _varSet->getModel()->getReferencePhysicalData();
  //DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  //DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  //DataHandle<CFint> isOutward =  this->socket_isOutward.getDataHandle();
  
  //const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  //const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();	// Maxwell's Eqs.+ Multifluid NavierStokes Eqs.
  //const CFuint nbEqs = eqSS.getNbEqsSS();  				
  //const CFuint iEqSS = eqSS.getEqSS();

  //Set the velocities
  const CFreal u = d_physicalData[firstVelocity];
  const CFreal v = d_physicalData[firstVelocity + 1];
  CFreal w = 0.;
  if(dim == 3 || is2DHalf == true){ w = d_physicalData[firstVelocity + 2];}
  const CFreal omega = m_dco->electricalConductivity;
  //const CFreal ovEpsilon = 1./term->getPermittivity(); OLD way of getting the permittivity
  const CFreal c_e = updateVS->getLightSpeed();
  const CFreal mu0 = updateVS->getPermeability();
  const CFreal ovEpsilon = c_e*c_e*mu0;
  CFreal Jx = 0., Jy = 0., Jz = 0.;
  
  if(dim == 2) 
  {
    Jx = omega*(d_Etotal[0] + v*d_Btotal[2]);
    Jy = omega*(d_Etotal[1] - u*d_Btotal[2]);
    Jz = omega*(d_Etotal[2] + u*d_Btotal[1] - v*d_Btotal[0]);
  }
  if(dim == 3 || is2DHalf == true)
  {
    Jx = omega*(d_Etotal[0] + v*d_Btotal[2] - w*d_Btotal[1]);
    Jy = omega*(d_Etotal[1] + w*d_Btotal[0] - u*d_Btotal[2]);
    Jz = omega*(d_Etotal[2] + u*d_Btotal[1] - v*d_Btotal[0]);  
  }
  CFreal JzCurl = 0.;

    //  cf_assert(states->size() == 1);

// NOT IMPLEMENTED
//      const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
//      const CFuint nbFacesInElem = faces.size();
//      CFreal faceAvBx = 0.0;
//      CFreal faceAvBy = 0.0;
      CFreal curlB    = 0.0;

//      for (CFuint iFace = 0; iFace < nbFacesInElem; ++iFace) {

//        const GeometricEntity *const face = element->getNeighborGeo(i);
//        const CFuint faceID = face->getID();
//        const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
//        const CFuint nbFaceNodes = face->nbNodes();
//        const CFreal ovNbFaceNodes = 1./(CFreal)nbFaceNodes;

//        const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

//        CFreal nx = normals[startID];
//        CFreal ny = normals[startID + 1];

//        if (isOutward[faceID] != static_cast<CFint>(elemID)) {
//          nx *= -1.;
//          ny *= -1.;
//        }

//        for (CFuint n = 0; n < nbFaceNodes; ++n) {
//          const CFuint nodeID = face->getNode(n)->getLocalID();
//          const RealVector& nodalState = nstates[nodeID];
//          const CFuint BxID = 0;
//          const CFuint ByID = 1;
//          // consider all 3D components here even in 2D
//          for (CFuint d = 0; d < 2; ++d) {
//            const CFreal nd = m_normal[d];
//            _gradBx[d] += nd*nodalState[BxID]*ovNbFaceNodes;
//            _gradBy[d] += nd*nodalState[BxID]*ovNbFaceNodes;
//          }
//        }


//        //const CFuint faceID = (faces)[iFace]->getID();
//        //State *const leftState = (faces)[iFace]->getState(0);
//        ////State *const rightState = (*faces)[iFace]->getState(1);
//        //const GeometricEntity *const currFace = (faces)[iFace];
//        //State *const rightState = (currFace->getState(1) != currState) ?
//        //            currFace->getState(1) : currFace->getState(0);

//        //const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

//        //CFreal nx = normals[startID];
//        //CFreal ny = normals[startID + 1];

//        //if (isOutward[faceID] != static_cast<CFint>(elemID)) {
//          //nx *= -1.;
//          //ny *= -1.;
//        //}

//        //_varSet->computePhysicalData(*leftState, _dataLeftState);
//        //_varSet->computePhysicalData(*rightState, _dataRightState);
//        //faceAvBx = 0.5*(_dataLeftState[UPDATEVAR::PTERM::BX] +
//        //        _dataRightState[UPDATEVAR::PTERM::BX]);
//        //faceAvBy = 0.5*(_dataLeftState[UPDATEVAR::PTERM::BY] +
//        //        _dataRightState[UPDATEVAR::PTERM::BY]);
//        curlB += -ny*faceAvBx + nx*faceAvBy;
//      }

// computing gradients for divB
/*
      const CFuint BxID = 0;
      const CFuint ByID = 1;
      const CFuint BzID = 2;
      const CFuint gradBxID = elemID*totalNbEqs + BxID;
      const CFuint gradByID = elemID*totalNbEqs + ByID;
      const CFuint gradBzID = elemID*totalNbEqs + BzID;
      _gradBx[XX] = this->m_ux[gradBxID];
      _gradBx[YY] = this->m_uy[gradBxID];

      _gradBy[XX] = this->m_ux[gradByID];
      _gradBy[YY] = this->m_uy[gradByID];

      if(dim == 3)
      {
        _gradBx[ZZ] = this->m_uz[gradBxID];
        _gradBy[ZZ] = this->m_uz[gradByID];

        _gradBz[XX] = this->m_ux[gradBzID];
        _gradBz[YY] = this->m_uy[gradBzID];
        _gradBz[ZZ] = this->m_uz[gradBzID];
      }
*/

/*  NOT IMPLEMENTED ON GPU
 // if (nbEqs == totalNbEqs || iEqSS == 0) {
    source[0] = 0.;			//x-Faraday's Law
    source[1] = 0.;			//y-Faraday's Law
    source[2] = 0.;			//z-Faraday's Law
    source[3] = - ovEpsilon*Jx;			//x-Ampere's Law
    source[4] = - ovEpsilon*Jy;			//y-Ampere's Law
    source[5] = - ovEpsilon*Jz;  		//z-Ampere's Law
    source[6] = 0.;			//y-Ampere's Law
    source[7] = 0.;  		//z-Ampere's Law
 // }
  
 // if (nbEqs == totalNbEqs || iEqSS == 1) {
    
    curlB = _gradBy[XX] - _gradBx[YY];
    JzCurl = curlB/mu0;
    if(is2DHalf) {
      source[8] = 0;                                    //Continuity equation
      source[9] = Jy*d_Btotal[2] -JzCurl*d_Btotal[1];			//x-momentum equation
      source[10] = JzCurl*d_Btotal[0] -Jx*d_Btotal[2];          //y-momentum equation
      source[11] = Jx*d_Btotal[1] - Jy*d_Btotal[1];          //z-momentum equation
      source[12] = Jx*d_Etotal[0] + Jy*d_Etotal[1] + JzCurl*(v*d_Btotal[0] - u*d_Btotal[1]);	//Total Energy
    }
    if(dim == 3)
    {
      Jx = (_gradBz[YY] - _gradBy[ZZ])/mu0;
      Jy = (_gradBx[ZZ] - _gradBz[XX])/mu0;
      Jz = (_gradBy[XX] - _gradBx[YY])/mu0;

      CFreal Ex = -(v*d_Btotal[2] - w*d_Btotal[1]);
      CFreal Ey = -(w*d_Btotal[0] - u*d_Btotal[2]);
      CFreal Ez = -(u*d_Btotal[1] - v*d_Btotal[0]);


      source[8]  = 0.;                                      //Continuity equation
      source[9]  = Jy*d_Btotal[2] - Jz*d_Btotal[1];                //x-momentum equation
      source[10] = Jz*d_Btotal[0] - Jx*d_Btotal[2];          	  //y-momentum equation
      source[11] = Jx*d_Btotal[1] - Jy*d_Btotal[1];          	  //z-momentum equation
      source[12] = Jx*Ex + Jy*Ey + Jz*Ez;  //Total Energy
    }
    if(dim == 2) {
      source[8] = 0;                                  //Continuity equation
      source[9] = Jy*d_Btotal[2] -JzCurl*d_Btotal[1];			//x-momentum equation
      source[10] = JzCurl*d_Btotal[0] -Jx*d_Btotal[2];		//y-momentum equation
      source[11] = Jx*d_Etotal[0] + Jy*d_Etotal[1] + JzCurl*(v*d_Btotal[0] - u*d_Btotal[1]);	//Total Energy
    }
 // }


*/

  if (m_dco->isResistive == true) {
    if(is2DHalf || dim ==3 ) {
      source[8] = 0;                                      //Continuity equation
      source[9] = Jy*d_Btotal[2] -Jz*d_Btotal[1];                 //x-momentum equation
      source[10] = Jz*d_Btotal[0] -Jx*d_Btotal[2];                //y-momentum equation
      source[11] = Jx*d_Btotal[1] - Jy*d_Btotal[1];               //z-momentum equation
      source[12] = Jx*d_Etotal[0] + Jy*d_Etotal[1] + Jz*d_Etotal[2];  //Total Energy
    }
    else {
      source[8] = 0;                                      //Continuity equation
      source[9] = Jy*d_Btotal[2] -Jz*d_Btotal[1];                 //x-momentum equation
      source[10] = Jz*d_Btotal[0] -Jx*d_Btotal[2];                //y-momentum equation
      source[11] = Jx*d_Etotal[0] + Jy*d_Etotal[1] + Jz*d_Etotal[2];  //Total Energy
    }
  }
 // source *= volumes[elemID];

/* NOT IMPLEMENTED GPU
  divB[elemID] = _gradBx[XX] + _gradBy[YY];
  Current[elemID] = curlB/mu0;
  BxPotential[elemID] = _NonInducedEMField[0];
  ByPotential[elemID] = _NonInducedEMField[1];
  BzPotential[elemID] = _NonInducedEMField[2];
*/
  //cout << "source = "<< source << endl;


}

#endif

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "HartmannSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HartmannSourceTerm_hh

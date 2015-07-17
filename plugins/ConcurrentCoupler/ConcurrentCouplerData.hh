#ifndef COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerData_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/PE.hh"

#include "Common/CFMap.hh"
#include "Common/Quartet.hh"

#include "Framework/MethodCommand.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * ConcurrentCouplerCom commands that compose @see ConcurrentCoupler
   *
   * @author Andrea Lani
   * @author Thomas Wuilbaut
   */
class ConcurrentCouplerData : public Framework::CouplerData {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Definition of a pair for defining the geometric entity index
  /// we store: 1) Name of the TRS where the face is
  ///           2) idx of Face inside the TRS
  typedef std::pair<Common::SafePtr<Framework::TopologicalRegionSet>, CFuint> GeoEntityIdx;
  
  /// Definition of a quartet for defining the geometric entity
  /// we store: 1) Name of the TRS where the face is
  ///           2) idx of Face inside the TRS
  ///           3) ShapeFunction values at the corresponding coord
  ///           4) Coordinates of the corresponding point
  typedef Common::Quartet<Common::SafePtr<Framework::TopologicalRegionSet>, CFuint, RealVector, RealVector> CoupledGeoEntity;
  
  ///vector of CoupledGeoEntity (for each otherTRS of a coupled interface)
  typedef std::vector<CoupledGeoEntity> CoupledGeoEntities;
  
  /// vector of CoupledGeoEntities (for each CoupledInterface)
  typedef std::vector<CoupledGeoEntities> CoupledInterface;
  
  /// vector of CoupledInterface's
  typedef std::vector< std::vector<CoupledInterface> > CoupledInterfaces;
  
  /**
   * Default constructor without arguments
   */
  ConcurrentCouplerData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~ConcurrentCouplerData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Sets up the FiniteElementData
   */
  void setup();

  /**
   * Does the dtaa transfer occur through files
   */
  bool isTransferFiles() {return _isTransferFiles;}
  
  /**
   * Gets the name of the socket for Coordinates (current SubSystem)
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   */
  std::vector<std::string> getThisCoupledCoordName(const std::string& interfaceName,
						   const std::string& trs);
  
  /**
   * Gets the name of the socket for Coordinates (current SubSystem)
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getThisCoupledCoordName(const std::string& interfaceName,
				      const std::string& trs,
				      const std::string& coordType);
  
  /**
   * Gets the name of the socket for Coordinates (other SubSystem)
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   */
  std::vector<std::string> getOtherCoupledCoordName(const std::string& interfaceName,
						    const std::string& trs, const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket for Coordinates (other SubSystem)
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getOtherCoupledCoordName(const std::string& interfaceName,
				       const std::string& trs,
				       const std::string& coordType, const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket for Accepted nodes
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   */
  std::vector<std::string> getThisCoupledAcceptedName(const std::string& interface,
						      const std::string& trs);
  
  /**
   * Gets the name of the socket for Accepted nodes
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param iProcessor rank of the processor from which originates the acceptance flag
   */
  std::vector<std::string> getThisCoupledAcceptedName(const std::string& interface,
						      const std::string& trs,
						      const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket for Accepted nodes
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   * @param iProcessor rank of the processor from which originates the acceptance flag
   */
  std::string getThisCoupledAcceptedName(const std::string& interface,
					 const std::string& trs,
					 const std::string& coordType,
					 const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket for Accepted nodes
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getThisCoupledAcceptedName(const std::string& interface,
					 const std::string& trs,
					 const std::string& coordType);
  
  /**
   * Gets the name of the socket for Accepted nodes
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   */
  std::vector<std::string> getOtherCoupledAcceptedName(const std::string& interface,
						       const std::string& trs, const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket for Accepted nodes
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getOtherCoupledAcceptedName(const std::string& interface,
					  const std::string& trs,
					  const std::string& coordType, const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param iProcessor rank of the processor from which originates the data
   */
  std::vector<std::string> getThisCoupledDataName(const std::string& interface,
						  const std::string& trs,
						  const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   */
  std::vector<std::string> getThisCoupledDataName(const std::string& interface,
						  const std::string& trs);
  
  
  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getThisCoupledDataName(const std::string& interface,
				     const std::string& trs,
				     const std::string& coordType);

  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getThisCoupledDataName(const std::string& interface,
				     const std::string& trs,
				     const std::string& coordType, const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   */
  std::vector<std::string> getOtherCoupledDataName(const std::string& interface,
						   const std::string& trs, const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getOtherCoupledDataName(const std::string& interface,
				      const std::string& trs,
				      const std::string& coordType, const CFuint& iProcessor);
  
  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   */
  std::vector<std::string> getThisCoupledConnectedName(const std::string& interface,
                                                const std::string& trs);

  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getThisCoupledConnectedName(const std::string& interface,
					  const std::string& trs,
					  const std::string& coordType);
  
  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   */
  std::vector<std::string> getOtherCoupledConnectedName(const std::string& interface,
							const std::string& trs);

  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType name of the type of coordinates to be used
   */
  std::string getOtherCoupledConnectedName(const std::string& interface,
					   const std::string& trs,
					   const std::string& coordType);
  
  /**
   * Gets CoupledGeoEntities for a given otherTRS of a given CoupledInterface
   * @param coupledInterfaceName name of the coupled interface
   * @param otherTRSidx idx of the TRS of the interface
   */
  CoupledGeoEntities* getCoupledInterfaces(const std::string& coupledInterfaceName,
                                           const CFuint& otherTRSidx,
                                           const CFuint& iType,
                                           const CFuint& iProc)
  {
    const CFuint idx = _coupledInterfacesMap.find(coupledInterfaceName);
    const CFuint nbProcessors = Common::PE::GetPE().GetProcessorCount(getNamespace());
    const CFuint idx2 = (nbProcessors*iType) + iProc;
    return &((_coupledInterfaces[idx])[otherTRSidx])[idx2];
  }

  /**
   * Resize the CoupledInterfaces
   * @param i new size of _coupledInterfaces
   */
  void resizeCoupledInterfaces()
  {
    _coupledInterfaces.resize(_coupledInterfacesMap.size());
  }

  /**
   * Resize the CoupledInterfaces
   */
  void resizeCoupledInterfacesTRS()
  {
    for (CFuint iInter=0;iInter < _coupledInterfaces.size();iInter++)
    {
      const CFuint nbOtherTRS = _coupledSubSystemsTRSNames[iInter].size();
      const CFuint nbProcessors = Common::PE::GetPE().GetProcessorCount(getNamespace());
      _coupledInterfaces[iInter].resize(_coupledSubSystemsTRSNames[iInter].size());

      for (CFuint iTRS=0;iTRS < nbOtherTRS;iTRS++)
      {
        const CFuint nbOtherTypes = _coupledSubSystemsCoordType[iInter].size();
        (_coupledInterfaces[iInter])[iTRS].resize(nbOtherTypes * nbProcessors);
      }
    }
  }

  CFuint getThisTransferedSize(const std::string& interfaceName)
  {
    const CFuint idx = _coupledInterfacesMap.find(interfaceName);
    return _transferedSize[idx];
  }

  /**
   * Set the names of the SubSystems to be coupled with
   */
  void setCoupledSubSystemsNames(const std::vector<std::string>& subSysNames, 
				 const std::vector<std::string>& nameSpacesNames)
  {
    _coupledSubSystemsNames = subSysNames;
    _coupledNameSpacesNames = nameSpacesNames;
    
    ///@todo this should not be here...
    _coupledSubSystemsTRSNames.resize(subSysNames.size());
    _coupledSubSystemsCoordType.resize(subSysNames.size());
    _transferedSize.resize(subSysNames.size());
    
    if(_isGeometryNonMatching.size() != subSysNames.size()) {
      _isGeometryNonMatching.resize(subSysNames.size());
    }
  }
  
  /**
   * Set the names of the TRS of the Other SubSystem
   */
  void setCoupledSubSystemTRSNames(const std::string& interfaceName, const std::vector<std::string>& trsNames)
  {
    CFuint idx = _coupledInterfacesMap.find(interfaceName);
    _coupledSubSystemsTRSNames[idx].resize(trsNames.size());
    _coupledSubSystemsTRSNames[idx] = trsNames;
  }

  /**
   * Set the type of coord passed by the Other SubSystem
   */
  void setCoupledSubSystemCoordTypes(const std::string& interfaceName, 
				     const std::vector<std::string>& coordType)
  {
    const CFuint idx = _coupledInterfacesMap.find(interfaceName);
    _coupledSubSystemsCoordType[idx].resize(coordType.size());
    _coupledSubSystemsCoordType[idx] = coordType;
  }
  
  /**
   * Set the size of the data being transfered
   */
  void setTransferedSize(const std::string& interfaceName, const CFuint dataSize)
  {
    const CFuint idx = _coupledInterfacesMap.find(interfaceName);
    _transferedSize[idx] = dataSize;
  }
  
  /**
   * Get the names of the TRS of the Other SubSystem
   * @param interfaceName name of the coupled interface
   * @return vector containing the names of the coupled subsystems
   */
  const std::vector<std::string> getCoupledSubSystemsTRSNames(const std::string& interfaceName)
  {
    CFuint idx = _coupledInterfacesMap.find(interfaceName);
    return _coupledSubSystemsTRSNames[idx];
  }

  /**
   * Get the names of the Other SubSystems
   * @param interfaceName name of the coupled interface
   * @return the name of the coupled subsystem
   */
  std::string getCoupledSubSystemName(const std::string& interfaceName)
  {
    CFuint idx = _coupledInterfacesMap.find(interfaceName);
    return _coupledSubSystemsNames[idx];
  }

  /**
   * Get the names of the coupled namespace of the Other SubSystem
   * @param interfaceName name of the coupled interface
   * @return the name of the coupled namespace
   */
  std::string getCoupledNameSpaceName(const std::string& interfaceName)
  {
    CFuint idx = _coupledInterfacesMap.find(interfaceName);
    return _coupledNameSpacesNames[idx];
  }

  /**
   * is the Interface Geometry Matching?
   */
  bool isInterfaceGeometryNonMatching(const std::string& interfaceName)
  {
    CFuint idx = _coupledInterfacesMap.find(interfaceName);
    return _isGeometryNonMatching[idx];
  }

  /**
   * Get the interface rotation to match it with the other interface
   */
  CFreal getInterfaceRotation(const std::string& interfaceName)
  {
    CFuint idx = _coupledInterfacesMap.find(interfaceName);
    return _nonMatchingGeometryRotation[idx];
  }

  /**
   * Get the interface translation vector to match it with the other interface
   */
  RealVector getInterfaceTranslation(const std::string& interfaceName)
  {
    // here we store all the vectors for the different interfaces in one large vector...
    const CFuint idx = _coupledInterfacesMap.find(interfaceName);
    const CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
    RealVector vector(dim);
    CFuint init = idx*dim;

    ///@todo Not the most efficient!!!
    for(CFuint i=0;i<dim;++i)
    {
      vector[i] = _nonMatchingGeometryVector[init+i];
    }

    return vector;
  }

  /**
   * Get the MeshMapping map
   */
  Common::SafePtr<Common::CFMap<std::string, CFuint> > getCoupledInterfacesMapPtr()
  {
    return &_coupledInterfacesMap;
  }

  /**
   * Configures the variables transfomers
   */
  void configureVariablesTransformers( Config::ConfigArgs& args );

  /**
   * Get the type of coord being transfered
   */
  std::vector<std::string> getThisCoordType(const std::string& interfaceName)
  {
    const CFuint idx = _coupledInterfacesMap.find(interfaceName);
    cf_assert(idx < _coordType.size());

    return _coordType[idx];
  }

  /**
   * Get the type of coord being transfered
   */
  std::vector<std::string> getOtherCoordType(const std::string& interfaceName)
  {
    const CFuint idx = _coupledInterfacesMap.find(interfaceName);
    cf_assert(idx < _coupledSubSystemsCoordType.size());

    return _coupledSubSystemsCoordType[idx];
  }

  /**
   * Get the threshold for non matching geometries
   */
  CFreal getNonMatchingGeometryThreshold(const std::string& interfaceName)
  {
    CFuint idx = _coupledInterfacesMap.find(interfaceName);
    if(_nonMatchingGeometryThreshold.size()!=0)
    {
      if (idx < _nonMatchingGeometryThreshold.size()) return _nonMatchingGeometryThreshold[idx];
      if (idx >= _nonMatchingGeometryThreshold.size()) return _nonMatchingGeometryThreshold[0];
    }

    // else
    return _defaultNonMatchingGeometryThreshold;
  }

  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> >
  getStdTrsGeoBuilder()
  {
    return &_stdTrsGeoBuilder;
  }

  /**
   * @return the GeometricEntity builder for faces in FVMCC
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
  getFaceTrsGeoBuilder()
  {
    return &_faceTrsGeoBuilder;
  }
  
  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "Coupler";
  }

  /**
   * Sets the interfaces in the Method data so that all commands can access it
   * All intterfaces in this method are CommandGroups
   * @param theValue is a vector with pointer to CommandGroups
   */
  void setInterfaces(const std::vector<Common::SafePtr< Framework::CommandGroup > >& theValue)
  {
    _interfaces = theValue;
  }
  
  /**
   * Gets the interfaces in the Method data so that all commands can access it
   * All intterfaces in this method are CommandGroups
   * @return a vector with pointers to CommandGroups
   */
  std::vector< Common::SafePtr< Framework::CommandGroup > > getInterfaces() const
  {
    return _interfaces;
  }

  /**
   * Sets the SpaceMethod which this Coupler uses
   */
  void setSpaceMethod(Framework::MultiMethodHandle<Framework::SpaceMethod> spaceMtd)
  {
    _spaceMethod= spaceMtd;
  }
  
  /**
   * Gets the SpaceMethod which this DataProcessing uses
   * @return pointer to the SpaceMethod
   */
  Framework::MultiMethodHandle<Framework::SpaceMethod> getSpaceMethod() const
  {
    return _spaceMethod;
  }
  
private:

  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType type of point coupled
   * @param Type type of data to be transferred
   */
  std::string thisNameBuilder(const std::string& interface,
			      const std::string& trs,
			      const std::string& coordType,
			      const std::string& type);
  
  /**
   * Gets the name of the socket to which is coupled
   * @param interface name of the coupled interface
   * @param trs name of the trs of the interface
   * @param coordType type of point coupled
   * @param Type type of data to be transferred
   */
  std::string otherNameBuilder(const std::string& interface,
			       const std::string& trs,
			       const std::string& coordType,
			       const std::string& type);
  
private:
  
  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;

  /// builder for face TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceTrsGeoBuilder;
  
  /// Coupled interfaces correspondance data
  CoupledInterfaces _coupledInterfaces;
  
  ///CFMap to get the index of the interface with its name as a key
  Common::CFMap<std::string, CFuint>  _coupledInterfacesMap;
  
  ///CFMap to get the index of a TRS with its name as a key
  Common::CFMap<std::string, CFuint>  _trsMap;

  ///The names  of the subsytems and their namespace to be coupled with
  std::vector<std::string> _coupledSubSystemsNames;
  std::vector<std::string> _coupledNameSpacesNames;

  ///The names of the TRS of the coupled subsytem
  std::vector< std::vector <std::string> > _coupledSubSystemsTRSNames;

  ///Flag for Non Matching Geometries at the interfaces
  std::vector<bool> _isGeometryNonMatching;

  ///Non Matching Geometries at the interfaces:
  ///Threshold for partially non matching TRSs
  std::vector<CFreal> _nonMatchingGeometryThreshold;

  ///Default value of threshold for partially non matching TRSs
  CFreal _defaultNonMatchingGeometryThreshold;

  ///Non Matching Geometries at the interfaces:
  ///Rotation Angle to be applied on the TRSs of the interface
  std::vector<CFreal> _nonMatchingGeometryRotation;

  ///Non Matching Geometries at the interfaces:
  ///Translation Vector to be applied on the TRSs of the interface
  std::vector<CFreal> _nonMatchingGeometryVector;

  ///Type of Coordinates (Nodal or at Quadrature Points)
  std::vector<std::string> _coordTypeStr;

  ///Type of Coord to transfer
  std::vector< std::vector<std::string> > _coordType;

  /// This assumes that all the command groups in this methos are interfaces
  std::vector<Common::SafePtr<Framework::CommandGroup> > _interfaces;

  /// handle to the space method
  Framework::MultiMethodHandle<Framework::SpaceMethod> _spaceMethod;

  ///Temporary vector
  std::vector<std::string> _tempStringVector;

  ///Vector with the type of the data being transfered from the other side (nodal, gauss...)
  std::vector< std::vector<std::string> > _coupledSubSystemsCoordType;

  ///Vector with the size of the data being transfered from the other side
  std::vector<CFuint> _transferedSize;

  ///Flag to know if we should use files to transfer the data
  bool _isTransferFiles;

}; // end of class ConcurrentCouplerData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for Coupler
typedef Framework::MethodCommand<ConcurrentCouplerData> ConcurrentCouplerCom;

/// Definition of a command provider for Coupler
typedef Framework::MethodCommand<ConcurrentCouplerData>::PROVIDER ConcurrentCouplerComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerData_hh

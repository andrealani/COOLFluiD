#ifndef COOLFluiD_UFEM_UFEMSolverData_hh
#define COOLFluiD_UFEM_UFEMSolverData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntityPool.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/VolumeIntegrator.hh"

#include "Framework/DofDataHandleIterator.hh"
#include "Framework/ProxyDofIterator.hh"

#include "UFEM/UFEM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace UFEM {

    class UFEMTerm;
    struct ElemAssembler;

    typedef std::string                                  TRSName;
    typedef CFuint                                       ElemTypeID;
    /// @todo when solution field is ready add here
    typedef std::pair < TRSName, ElemTypeID >            ElemID;
    typedef std::map  < ElemID , ElemAssembler* >        ElemAssemblerMap;

//////////////////////////////////////////////////////////////////////////////

/// This class represents data object accessed by different UFEMSolverCom's
/// @author Tiago Quintino
/// @author Pedro Maciel
class UFEM_API UFEMSolverData : public Framework::SpaceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit UFEMSolverData(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~UFEMSolverData();

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre LinearSystemSolver pointer is not constant to allow dynamic_casting
  void setLinearSystemSolver( Framework::MultiMethodHandle< Framework::LinearSystemSolver > lss )
  {
    m_lss = lss;
  }

  /// Get the linear system solver
  Framework::MultiMethodHandle< Framework::LinearSystemSolver > getLinearSystemSolver() const
  {
    cf_assert(m_lss.isNotNull());
    return m_lss;
  }

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to ConvergenceMethod is not constant to
  ///      allow dynamic_casting
  void setConvergenceMethod(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd)
  {
    m_convergenceMtd = convMtd;
  }

  /// Get the ConvergenceMethod
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> getConvergenceMethod() const
  {
    cf_assert(m_convergenceMtd.isNotNull());
    return m_convergenceMtd;
  }


  /// @return the GeometricEntity builder
  Common::SafePtr<
    Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > >
    getStdTrsGeoBuilder()
  {
    return &m_stdTrsGeoBuilder;
  }

  /// Gets the Class name
  static std::string getClassName()  { return "UFEMSolver"; }

  /// Setup the data of the method UFEMSolver
  virtual void setup();

  /// Unsetup the data of the method UFEMSolver
  virtual void unsetup();

  /// accessor to time step
  void setDt(CFreal Dt) { m_Dt=Dt; }

  /// accessor to time step
  CFreal getDt() { return m_Dt; }

  /// accessor to time step multiplicator
  CFreal getDtmult() { return m_Dtmult; }

  /// accessor to time step limit
  CFreal getDtlimit() { return m_Dtlimit; }

  /// accessor to density at element center
  CFreal getRhoElm() { return m_RhoElm; }

  /// accessor to density at element center
  void setRhoElm(CFreal rhoelm) { m_RhoElm=rhoelm; }

  /// accessor to effective viscosity at element center
  CFreal getMuElm() { return m_MuElm; }

  /// accessor to effective viscosity at element center
  void setMuElm(CFreal muelm) { m_MuElm=muelm; }

  /// accessor to 
  CFint getPrintNode() { return m_PrintNode; }

  /// accessor to 
  CFint getOExtrap() { return m_OExtrap; }

  /// accessor to elem_assembler
  ElemAssemblerMap& getelem_assemblers() { return m_elem_assemblers; }

  /// accessor wall distance calculation switch
  bool getCalcWallDistance() { return m_calcWallDistance; }

  /// accessor wall distance from file switch
  bool getReadWallDistFromFile() { return m_ReadWallDistFromFile; }

  /// accessor wall distance file name
  std::string getWallDistFileName() { return m_WallDistFileName; }
  
  
  
private:  // helper functions

  /// Configures the terms of the numerical treatment
  void setupUFEMTerms();

private:  // data

  /// Linear system solver
  Framework::MultiMethodHandle< Framework::LinearSystemSolver > m_lss;

  /// Convergence Method
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> m_convergenceMtd;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > m_stdTrsGeoBuilder;

  /// The volume integrator
  Framework::VolumeIntegrator m_volumeIntegrator;

  /// String for configuring the numerical integrator QuadratureType
  std::string m_intquadStr;

  /// String for configuring the numerical integrator Order
  std::string m_intorderStr;

  /// Time step size (at start)
  CFreal m_Dt;

  /// Time step size multiplicator
  CFreal m_Dtmult;

  /// Time step size limit
  CFreal m_Dtlimit;

  /// Density at element center
  CFreal m_RhoElm;

  /// Effective viscosity at element center
  CFreal m_MuElm;

  /// Time step size
  CFint m_PrintNode;

  /// Time step size
  CFint m_OExtrap;

  /// list of terms
  std::vector < std::string > m_terms_strs;

  /// a vector of element assemblers treatments applied per TRS, per Element Type
  ElemAssemblerMap  m_elem_assemblers;

  /// store the configure arguments for later
  Config::ConfigArgs m_stored_args;

  /// switch to calcuilate wall distance
  bool m_calcWallDistance;
  
  bool m_ReadWallDistFromFile;
  
  std::string m_WallDistFileName;

};  // end of class UFEMSolverData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a MethodCommand for UFEMSpaceMethod
typedef Framework::MethodCommand< UFEMSolverData > UFEMSolverCom;

/// Definition of a MethodStrategy for UFEMSpaceMethod
typedef Framework::MethodStrategy< UFEMSolverData > UFEMSolverStrategy;

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_UFEMSolverData_hh


#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAeroSpectralFD_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAeroSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "MathTools/FunctionParser.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/FaceTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

    namespace Framework
    {
      class FaceToCellGEBuilder;
    }
    namespace Physics
    {
      namespace NavierStokes
      {
        class EulerVarSet;
        class NavierStokesVarSet;
      }
    }
    namespace SpectralFD
    {
      class BaseBndFaceTermComputer;
      class BaseFaceTermComputer;
      class BaseVolTermComputer;
      class BCStateComputer;
      class CellToFaceGEBuilder;
      class SpectralFDMethodData;
    }

//////////////////////////////////////////////////////////////////////////////

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Wall values and aerodynamic coefficients
 * for NavierStokes2DCons with spectral difference
 *
 * @author Kris Van den Abeele
 *
 */
class NavierStokesConsComputeAeroSpectralFD : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesConsComputeAeroSpectralFD(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokesConsComputeAeroSpectralFD();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  void unsetup();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Compute the values at the wall and write them to file
   */
  void computeWall();

  /**
   * Compute the aerodynamic coefficients and write them to file
   */
  void computeAero();

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileWall();

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileAero();

  /**
   * Write the aerodynamic coefficients to file
   */
  void updateOutputFileAero();

  /// set face term data
  void setFaceAndVolumeTermData();

  /// set the local indexes of the other faces (not the current boundary face)
  /// @pre m_faces is set
  void setOtherFacesLocalIdxs();

  /// set the face neighbour states
  /// @pre setOtherFacesLocalIdxs()
  void setFaceNeighbourStates();

  /// compute the required data in the face term for the current face
  void computeFaceData();

  /// set the data related to the neighbouring cell
  /// @pre setOtherFacesLocalIdxs()
  /// @pre setFaceNeighbourStates()
  void setCellData();

  /// compute cell gradients
  void computeCellGradients();

private:

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;

  /// SpectralFD method data
  Common::SafePtr<SpectralFD::SpectralFDMethodData> m_spectralFDData;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<SpectralFD::CellToFaceGEBuilder> > m_cellBuilder;

  /// SpectralFD strategy that computes the boundary face terms
  Common::SafePtr< SpectralFD::BaseBndFaceTermComputer > m_bndFaceTermComputer;

  /// volume term computer
  Common::SafePtr< SpectralFD::BaseVolTermComputer > m_volTermComputer;

  /// face term computers
  std::vector< Common::SafePtr< SpectralFD::BaseFaceTermComputer    > > m_faceTermComputers;

  /// boundary face term computers
  std::vector< Common::SafePtr< SpectralFD::BaseBndFaceTermComputer > > m_bndFaceTermComputers;

  /// the BCStateComputer for this BC
  Common::SafePtr< SpectralFD::BCStateComputer > m_bcStateComputer;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< SpectralFD::BCStateComputer > > > m_bcStateComputers;

  /// Update variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_updateVarSet;

  /// Diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffusiveVarSet;

  /// variable for current face
  Framework::GeometricEntity* m_face;

  /// variable for current internal cell
  Framework::GeometricEntity* m_intCell;

  /// variable for (other) faces
  const std::vector< Framework::GeometricEntity* >* m_faces;

  /// variable for current face orientation
  CFuint m_orient;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;

  /// vector containing pointers to the left and right states with respect to a face
  std::vector< std::vector< std::vector< Framework::State* >* > > m_faceNghbrStates;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// internal cell gradients
  std::vector< std::vector< RealVector >* > m_cellGrads;

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// Jacobian determinants
  std::valarray<CFreal> m_solJacobDet;

  /// cell local indexes of the other faces (not the face itself)
  std::vector< CFuint > m_otherFaceLocalIdxs;

  /// pointer to booleans telling whether a face is on the boundary
  Common::SafePtr< std::vector< bool > > m_isFaceOnBoundary;

  /// pointer to neighbouring cell side vector
  Common::SafePtr< std::vector< CFuint > > m_nghbrCellSide;

  /// pointer to current cell side vector
  Common::SafePtr< std::vector< CFuint > > m_currCellSide;

  /// pointer to orientation vector
  Common::SafePtr< std::vector< CFuint > > m_faceOrients;

  /// pointer to BC index vector
  Common::SafePtr< std::vector< CFuint > > m_faceBCIdx;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// Storage for choosing when to save the wall values file
  CFuint m_saveRateWall;

  /// Storage for choosing when to save the aero coef file
  CFuint m_saveRateAero;

  /// Output File for Wall Values
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fileWall;

  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileWall;

  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileAero;

  //// Variables to store the wall values + lift/drag
  /// Storage for lift coeficient on trs
  CFreal m_lift;

  /// Storage for drag coeficient on trs
  CFreal m_drag;

  /// Storage for friction drag coeficient on trs
  CFreal m_fricDrag;

  /// force coefficient
  RealVector m_forceCoef;

  /// friction force coefficient
  RealVector m_fricForceCoef;

  /// velocity at infinity
  CFreal m_uInf;

  /// density at infinity
  CFreal m_rhoInf;

  /// pressure at infinity
  CFreal m_pInf;

  /// flow direction
  std::vector< CFreal > m_flowDir;

  /// projected surface
  CFreal m_projSurf;

  ///flag for appending iteration
  bool m_appendIter;

  ///flag for appending time
  bool m_appendTime;

  /// variable containing current gradients
  std::vector< RealVector*  > m_grads;

}; /// end of class NavierStokesConsComputeAeroSpectralFD

//////////////////////////////////////////////////////////////////////////////

    } /// namespace AeroCoef

} /// namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif /// COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAeroSpectralFD_hh

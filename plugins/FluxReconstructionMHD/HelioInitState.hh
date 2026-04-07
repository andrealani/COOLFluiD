#ifndef COOLFluiD_FluxReconstructionMHD_HelioInitState_hh
#define COOLFluiD_FluxReconstructionMHD_HelioInitState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Initializes the heliospheric domain from a boundary surface data file
 * using Parker solar wind radial extrapolation.
 *
 * Reads the same .dat file format as BCInletHelioUnsteadyMHD, interpolates
 * angularly to each solution point, and applies radial decay:
 *   rho ~ r^-2, Br ~ r^-2, p ~ r^(-2*gamma), Parker spiral for Bphi.
 *
 * @author Rayan Dhib
 */
class HelioInitState : public FluxReconstructionSolverCom {
public:

  static void defineConfigOptions(Config::OptionList& options);

  explicit HelioInitState(const std::string& name);

  ~HelioInitState();

  virtual void setup();

  virtual void unsetup();

  virtual void configure(Config::ConfigArgs& args);

  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  void executeOnTrs();

  /// Local nested class for grouping surface data
  class SurfaceData {
  public:
    RealMatrix xyz;
    RealVector rho;
    RealVector u;
    RealVector v;
    RealVector w;
    RealVector Bx;
    RealVector By;
    RealVector Bz;
    RealVector p;
  };

  /// Local nested class for closest point search
  class ClosestPointData {
  public:
    std::vector<CFint> surfaceIDs;
    std::vector<CFint> pointsIDs;
    RealVector r;

    void reset()
    {
      surfaceIDs.assign(surfaceIDs.size(), -1);
      pointsIDs.assign(pointsIDs.size(), -1);
      r = MathTools::MathConsts::CFrealMax();
    }

    void regressionFromTo(CFuint start, CFuint end)
    {
      cf_assert(end < surfaceIDs.size());
      cf_assert(start < surfaceIDs.size());
      surfaceIDs[end] = surfaceIDs[start];
      cf_assert(end < pointsIDs.size());
      cf_assert(start < pointsIDs.size());
      pointsIDs[end] = pointsIDs[start];
      cf_assert(end < r.size());
      cf_assert(start < r.size());
      r[end] = r[start];
    }
  };

  /// Read boundary surface data from .dat file
  void readSurfaceData(std::vector<SurfaceData*>& surfaces,
                       const std::string& fileName);

  /// Interpolate surface data at a given angular position using k-nearest neighbors
  void interpolateSurface(const RealVector& queryCoords,
                          CFreal& rho, CFreal& u, CFreal& v, CFreal& w,
                          CFreal& Bx, CFreal& By, CFreal& Bz, CFreal& p);

  /// Apply Parker solar wind radial decay from r_bc to r
  void applyParkerDecay(const RealVector& coords, CFreal r, CFreal theta,
                        CFreal& rho, CFreal& u, CFreal& v, CFreal& w,
                        CFreal& Bx, CFreal& By, CFreal& Bz, CFreal& p);

protected: // data

  /// path to boundary .dat file
  std::string m_fileNameTw;

  /// number of closest points for interpolation
  CFuint m_nbClosestPoints;

  /// inner boundary radius (0 = auto-detect from file)
  CFreal m_innerRadius;

  /// radial velocity scaling exponent (default 0 = constant Vr)
  CFreal m_velExponent;

  /// include Parker spiral winding for Bphi
  bool m_parkerSpiral;

  /// solar rotation rate (non-dimensional)
  CFreal m_omegaSun;

  /// input variable set name
  std::string m_inputVarStr;

  /// loaded surface data
  std::vector<SurfaceData*> m_surfaces;

  /// detected/configured inner boundary radius
  CFreal m_rBoundary;

  /// number of equations
  CFuint m_nbrEqs;

  /// dimensionality
  CFuint m_dim;

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> m_varSet;

  /// Transformer from input to update Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVar;

  /// states in initialization points
  std::vector< Framework::State* > m_initPntsStates;

  /// input state
  Framework::State* m_inputState;

  /// coordinates of initialization point
  RealVector m_initPntCoords;

  /// temporary node for distance computation
  RealVector m_tmpNode;

}; // class HelioInitState

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMHD_HelioInitState_hh

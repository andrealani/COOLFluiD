#ifndef COOLFluiD_Numerics_FluctSplit_CarbuncleFix2_hh
#define COOLFluiD_Numerics_FluctSplit_CarbuncleFix2_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeJacobianFix.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Common/CFMultiMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  

    
    namespace FluctSplit {
      class FluctuationSplitData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes a jacobian fix to cure entropy violation or carbuncle
 * for compressible flows
 *
 * @author Andrea Lani
 */
class CarbuncleFix2 : public ComputeJacobianFix {
public:
  
  /**
   * Default constructor without arguments
   */
  CarbuncleFix2(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~CarbuncleFix2();
  
  /**
   * Set up
   */
  virtual void setup();
  
  /**
   * Compute the jacobian fix
   */
  virtual void computeFix(const InwardNormalsData& normalsData,
                          RealVector& delta);
private:
  
  /**
   * Update the correction
   */
  void updateDelta(CFuint iState);

  /**
   * Compute the correction
   */
  void computeCellDelta(Framework::GeometricEntity *const cell,
			CFint iStateCurr);

private:
  
  /// normals data for the current cell
  InwardNormalsData* m_normals;
  
  /// correction values for each state
  RealVector* m_delta;
  
  /// max pressure gradient*volume*dim
  RealVector m_pGradTimesVolDim;
  
  /// normal shock unit normal
  RealVector m_nShock;

  /// physical data associeted to each state 
  std::vector<RealVector> m_pData;
  
  /// u*n - a eigenvalue for each state
  std::vector<RealVector> m_lambda;
  
  /// average cell speed
  RealVector m_avSpeed;

  /// array storing the local IDs in the 2D cell
  /// corresponding to the two nodes defining the shock face
  std::vector<CFuint> m_shockFaceNodes;
    
  /// cell builder
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> m_cellBuilder;

  /// node to cell connecivity
  Common::CFMultiMap<CFuint,CFuint> m_mapNode2CellID;
  
}; // end of class CarbuncleFix2
      
//////////////////////////////////////////////////////////////////////////////

    }  // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_CarbuncleFix2_hh

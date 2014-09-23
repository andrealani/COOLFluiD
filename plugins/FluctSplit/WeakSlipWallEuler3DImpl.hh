#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler3DImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler3DImpl_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for Euler3D
 *
 * @author Andrea Lani
 *
 */
class WeakSlipWallEuler3DImpl : public FluctuationSplitCom {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSlipWallEuler3DImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallEuler3DImpl();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

 protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Compute the normal flux and the corresponding jacobian
   */
  void computeNormalFluxAndJacob(const Framework::State& state,
				 const RealVector& normal,
				 RealVector& flux,
				 RealMatrix& fluxJacob);
  
  
  /**
   * Compute the face normal
   */
  void setFaceNormal(const CFuint faceID,
		     RealVector& normal);
  
  /**
   * @return the accumulator corresponding to the current face
   * @pre the id of the accumulator is assumed to be 0 for triangles
   *      and 1 for quadrilaters
   */
  Framework::BlockAccumulator* getAccumulator
  (const CFuint nbStatesInFace) const
  {
    return _accumulators[nbStatesInFace-3];
  }

 private:

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;
  
  /// Euler var set
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> _varSet;
  
  /// distribution coefficient
  CFreal                         _alpha;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im0;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in0;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im1;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in1;

  /// list of accumulators (you can have one for triangular faces
  /// and one for quadrilateral)
  std::vector<Framework::BlockAccumulator*> _accumulators;

  /// temporary storage for the fluxes
  std::vector<RealVector>     _fluxes;

  /// temporary storage for the fluxes
  std::valarray<RealMatrix>     _fluxJacobs;
  
  /// physical data
  RealVector _physicalData;
  
  /// temporary flux jacobian matrix
  RealMatrix _tJacob;
  
}; // end of class WeakSlipWallEuler3DImpl
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler3DImpl_hh

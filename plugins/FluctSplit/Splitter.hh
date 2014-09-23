#ifndef COOLFluiD_Numerics_FluctSplit_Splitter_hh
#define COOLFluiD_Numerics_FluctSplit_Splitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodStrategy.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace FluctSplit {

      class FluctuationSplitData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a generic RDS splitter
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API Splitter : public Framework::MethodStrategy<FluctuationSplitData> {
public:

  typedef Framework::BaseMethodStrategyProvider  <FluctuationSplitData,Splitter> PROVIDER;

  /// Constructor
  Splitter(const std::string& name);

  /// Virtual destructor
  virtual ~Splitter();

  /// @returns the class name as a string
  static std::string getClassName() { return "Splitter"; }

  /// Configures the object
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the object
  virtual void setup();

  /// Unsets up the object
  virtual void unsetup();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Compute the K, K+, K-
  virtual void computeK(const std::vector<Framework::State*>& states,
                        const InwardNormalsData* const normalsData) = 0;

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual)
  {
    throw Common::NotImplementedException (FromHere(),"Splitter::distribute()");
  }

  /// Distribute the residual
  virtual void distributePart(std::vector<RealVector>& residual)
  {
    throw Common::NotImplementedException (FromHere(),"Splitter::distributePart()");
  }

  /// Compute all the contributions to the Picard jacobian
  virtual void computePicardJacob(std::vector<RealMatrix*>& jacob)
  {
    throw Common::NotImplementedException (FromHere(),"Splitter::computePicardJacob()");
  }

  /// Compute part of the contributions to the Picard jacobian
  virtual void computePicardJacobPart(std::vector<RealMatrix*>& jacob)
  {
    throw Common::NotImplementedException (FromHere(),"Splitter::computePicardJacobPart()");
  }

private:

  /// Sets the correct block limits for a Scalar Splitter
  /// Called by the RDS_Splitter constructor.
  /// @see _nbEquations
  /// @see _firstVarID
  /// @see _lastVarID
  virtual void setBlockData() = 0;

protected: // member function

  /// Set the inward normal
  void setAdimensionalNormal(const CFuint iState)
  {
    using namespace Framework;
    cf_assert(m_normals != CFNULL);

    for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
      _adimNormal[iDim] = m_normals->getNodalNormComp(iState, iDim);
    }
    _adimNormal *= 1. / m_normals->getAreaNode(iState);
  }

protected: // data

  /// nb dimensions of the problem
  CFuint                        m_dim;

  /// size of the block in which to operate
  CFuint                        _nbEquations;

  /// first index of the block in which to operate
  CFuint                        _firstVarID;

  /// last index plus one of the block in which to operate
  CFuint                        _lastVarID;

  /// Number of states in the cell
  CFuint                        _nbStatesInCell;

  /// size of block separator
  CFuint _blockSeparator;

  /// Normals of the cell being currently treated.
  /// Should be set when computing the K
  const InwardNormalsData* m_normals;

  /// adimensionalized normal vector
  RealVector               _adimNormal;

  /// The sockets for the flags telling if states are on the boundary
  Framework::DataSocketSink<bool> socket_isBState;

  /// The sockets to use in this strategy for the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// The sockets for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // end of class Splitter

//////////////////////////////////////////////////////////////////////////////

    }  // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_Splitter_hh

#ifndef COOLFluiD_Numerics_FluctSplit_ComputeSourceTermFSM_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeSourceTermFSM_hh

//////////////////////////////////////////////////////////////////////////////



#include "Framework/ComputeSourceTerm.hh"
#include "Framework/Storage.hh"

#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /// This class offers and abstract interface to computes the source
  /// term for RD schemes.
  /// @author Andrea Lani
  /// @author Tiago Quintino
class FluctSplit_API ComputeSourceTermFSM : public Framework::ComputeSourceTerm<FluctuationSplitData> {
public:

  typedef Framework::BaseMethodStrategyProvider
  <FluctuationSplitData,ComputeSourceTermFSM> PROVIDER;

  /// Constructor
  ComputeSourceTermFSM(const std::string& name) :
    Framework::ComputeSourceTerm<FluctuationSplitData>(name),
    socket_volumes("volumes"),
    _jacobMatrix()
  {
  }
  
  /// Default destructor
  virtual ~ComputeSourceTermFSM()
  {
  }
  
  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup()
  {
    Framework::ComputeSourceTerm<FluctuationSplitData>::setup();
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    _jacobMatrix.resize(nbEqs, nbEqs);
  }

  /// Compute the source term
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian)
  {
    throw Common::NotImplementedException (FromHere(),"ComputeSourceTermFSM::computeSource()");
  }
  
  /// Compute the source term
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
				RealVector& source,
				const FluctSplit::InwardNormalsData& normalsData) = 0;

  /// Compute the source term
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
				std::vector<RealVector>& source,
				const FluctSplit::InwardNormalsData& normalsData)
  {
    throw Common::NotImplementedException (FromHere(),"ComputeSourceTermFSM::computeSource()");
  }
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "FluctSplit::ComputeSourceTermFSM";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
      Framework::ComputeSourceTerm<FluctuationSplitData>::needsSockets();
    result.push_back(&socket_volumes);

    return result;
  }

protected:

  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// jacobian matrix
  RealMatrix _jacobMatrix;
  
}; // end of class ComputeSourceTermFSM

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeSourceTermFSM_hh

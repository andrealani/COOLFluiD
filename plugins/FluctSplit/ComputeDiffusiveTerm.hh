#ifndef COOLFluiD_Numerics_FluctSplit_ComputeDiffusiveTerm_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeDiffusiveTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "MathTools/RealVector.hh"
#include "Environment/ConcreteProvider.hh"

#include "Common/NullableObject.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework
  {
    class DiffusiveVarSet;
    class GeometricEntity;
    class ConvectiveVarSet;
  }



    namespace FluctSplit {

      class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents an object computing a diffusive term
/// with Fluctuation Split method
/// @author Andrea Lani
class FluctSplit_API ComputeDiffusiveTerm : public Framework::MethodStrategy<FluctuationSplitData> {

public:

  typedef Framework::BaseMethodStrategyProvider<FluctuationSplitData,ComputeDiffusiveTerm> PROVIDER;
 
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Constructor
  ComputeDiffusiveTerm(const std::string& name);

  /// Default destructor
  virtual ~ComputeDiffusiveTerm();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args );

  /// Set up private data to prepare the simulation
  virtual void setup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
    
    result.push_back(&socket_updateCoeff);
    result.push_back(&socket_normals);
    result.push_back(&socket_volumes);

    return result;
  }

  /// Set mesh data
  void setMeshData();

  /// Set the diffusive variable set
  virtual void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar) {}
  
  /// Set the update variable set
  virtual void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet> updateVar) {}
  
  /// Compute the flux in the current face
  virtual void computeDiffusiveTerm(Framework::GeometricEntity *const geo,
				    std::vector<RealVector>& result,
				    const bool updateCoeff) {}
  
  /// Compute the cell gradient and average state and put them
  /// into @see DistributionData
  virtual void computeCellGradientsAndAverageState
  (Framework::GeometricEntity *const geo, const RealVector& pdata)
  {
    throw Common::NotImplementedException
      (FromHere(), "ComputeDiffusiveTerm::computeCellGradientsAndAverageState()");
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeDiffusiveTerm";
  }
 
  virtual void computePicardDiffJacob(Framework::GeometricEntity *const cell,std::vector<RealMatrix*>& _diffjacob) 
  {
    throw Common::NotImplementedException(FromHere(),"ComputeDiffusiveTerm::computePicardDiffJacob()");
  }
  
protected:
  
  /// Tells if this term can be added to another derived one
  virtual bool addToDerivedTerm() {return false;}
  
protected:

  /// socket for storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;
  
  /// socket for the inward normals storage
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// socket for the volume storage
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// socket with theta clending coefficient
  Framework::DataSocketSource<CFreal> socket_Pe_cell;

  /// Export the cell Peclet number
  bool _store_Pe_cell;
  
}; // end of class ComputeDiffusiveTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeDiffusiveTerm_hh

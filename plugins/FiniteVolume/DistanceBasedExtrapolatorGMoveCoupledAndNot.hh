#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveCoupledAndNot_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveCoupledAndNot_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/DistanceBasedExtrapolatorGMove.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object to be used
 * in combination with BCs that moves the ghost states (nodes) like
 * @see NoSlipWallAdiabaticTurb2D
 *
 * @author Thomas Wuilbaut
 *
 */
class DistanceBasedExtrapolatorGMoveCoupledAndNot : public DistanceBasedExtrapolatorGMove {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMoveCoupledAndNot(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveCoupledAndNot();

  /**
   * Configuration of the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Extrapolate the solution in all mesh nodes
   */
  virtual void extrapolateInAllNodes();

  /**
   * Extrapolate the solution in the given nodes
   */
  virtual void extrapolateInNodes(const std::vector<Framework::Node*>& nodes);

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
      DistanceBasedExtrapolatorGMove::needsSockets();

    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result2 = _sockets.getAllSinkSockets();
    for(CFuint i=0;i<result2.size();++i)
    {
      result.push_back(result2[i]);
    }

    return result;
  }

protected:

  virtual void checkTRSList(const CFuint count)
  {
  }

private:

  void setIndex();

private:

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  // indices of values taken from the coupler
  std::vector<CFuint> _wCoupledValuesIdx;

  //Names of the interfaces to which belong the TRS
  std::vector<std::string> _interfaceNames;
  std::vector<std::string> _coupledTRSNames;

  // Maps to get back the index of the node in the TRS list from its LocalID
  std::vector<Common::CFMap<CFuint, CFuint>  > _trsNodeIDMap;

  //temporary vector holding the coupled values
  RealVector _coupledValues;

  //index of the data related to index of the node for each of the TRSs
  std::vector<std::vector<CFint> > _coupledDataID;

  //Flag to know if the setup of the index list has been performed
  bool _isSetIndex;

  //Values to use when rejected node
  std::vector<CFreal> _defaultCoupledValues;

  CFuint _defaultIterations;

  CFuint _nbNonCoupledValues;
}; // end of class DistanceBasedExtrapolatorGMoveCoupledAndNot

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveCoupledAndNot_hh

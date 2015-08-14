#ifndef COOLFluiD_PhysicalModelDummy_CouplingModelDummySendToRecv_hh
#define COOLFluiD_PhysicalModelDummy_CouplingModelDummySendToRecv_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace PhysicalModelDummy {

    class CouplingModelDummy;
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a default transformer of variables from send to recv
 * when @see CouplingModelDummy is used
 *
 * @author Andrea Lani
 *
 */
class CouplingModelDummySendToRecv : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  CouplingModelDummySendToRecv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~CouplingModelDummySendToRecv();
   
  /// Transform an array of variables into another one of potentially different length
  void transform(const RealVector& state, RealVector& result);
  
  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
private:
  
  Common::SafePtr<CouplingModelDummy> m_model;
  
}; // end of class CouplingModelDummySendToRecv

//////////////////////////////////////////////////////////////////////////////

    } // namespace PhysicalModelDummy

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_PhysicalModelDummy_CouplingModelDummySendToRecv_hh

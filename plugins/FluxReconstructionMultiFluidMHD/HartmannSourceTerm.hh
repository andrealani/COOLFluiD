#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_HartmannSourceTerm_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_HartmannSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/StdSourceTerm.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/StdSourceTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Physics {
    namespace MultiFluidMHD {
      class DiffMFMHD2DVarSet;
    }

    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Source term for MultiFluid considering 2 fluids: plasma + neutrals
 * variables
 */
template <class UPDATEVAR>
class HartmannSourceTerm : public StdSourceTerm {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  HartmannSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~HartmannSourceTerm();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * get data required for source term computation
   */
  void getSourceTermData();

  /**
   * add the source term
   */
  void addSourceTerm(RealVector& resUpdates);

  CFreal getElectricalConductivity(){return _electricalConductivity;}
  bool getIsResistive(){return _isResistive;}

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // data

  /// variable for physical data of sol
  RealVector m_solPhysData;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

  /// dimensionality
  CFuint m_dim;

  /// physical model (in conservative variables)
  Common::SafePtr<UPDATEVAR > m_varSet;

  /// corresponding diffusive variable set
  Common::SafePtr<Framework::DiffusiveVarSet> m_diffVarSet;

  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;
  
  /// Euler physical data
  RealVector _physicalData;

  /// Euler physical data
  RealVector _sockets;

  /// handle to outward normal
  Framework::DataHandle<CFint> _isOutward;

  /// the gradients in the neighbouring cell
  std::vector< std::vector< RealVector >* > m_cellGrads;

  ///Dummy vector for the gradients
  std::vector<RealVector*> m_dummyGradients;
  
  /// socket for storing the divergence of magnetic field
  Framework::DataSocketSource<CFreal> socket_divB;
  
  /// socket for storing the Current
  Framework::DataSocketSource<CFreal> socket_Current;

  /// socket for storing the Bx potential
  Framework::DataSocketSource<CFreal> socket_BxPotential;
 
  /// socket for storing the By potential
  Framework::DataSocketSource<CFreal> socket_ByPotential;

  /// socket for storing the Bz potential
  Framework::DataSocketSource<CFreal> socket_BzPotential;

  /// array to store the mass fractions
  RealVector _ys;

  /// vector to store temporary result
  RealVector _temp;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;
  
  /// array of temporary values
  RealMatrix _values;

  ///Non Induced Part of the electromagnetic Field
  RealVector _NonInducedEMField;
  
  /// Current density vector
  RealVector _J;
  
  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients;
 
  /// Source term
  RealVector m_source;
  
private:

  /// Electrical conductivity
  CFreal _electricalConductivity;

  /// Orszag-Tang Conductivity Flag
  bool _isResistive;

  /// Left State to compute the Orszag Conductivity
  RealVector _dataLeftState;

   /// Right State to compute the Orszag Conductivity
  RealVector _dataRightState;

  /// gradient of Bx
  std::vector<CFreal> _gradBx;

  /// gradient of By
  std::vector<CFreal> _gradBy;

 /// gradient of Bz
 std::vector<CFreal> _gradBz;

}; // end of class HartmannSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "HartmannSourceTerm.ci"

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_HartmannSourceTerm_hh


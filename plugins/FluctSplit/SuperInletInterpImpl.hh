#ifndef COOLFluiD_Numerics_FluctSplit_SuperInletInterpImpl_hh
#define COOLFluiD_Numerics_FluctSplit_SuperInletInterpImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {
    template <class KEY, class VALUE> class LookUpTable;
  }

  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a supersonic inlet command
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API SuperInletInterpImpl : public SuperInlet{
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  SuperInletInterpImpl(const std::string& name);

  /// Default destructor
  ~SuperInletInterpImpl();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Set up the member data
  virtual void setup();

protected:

  /// Execute on the current TRS
  void executeOnTrs();

  /// fill the lookup table
  void fillTable();

protected: // data

  /// the socket to the data handle of the boundary state neighbors
  Framework::DataSocketSink<std::valarray<Framework::State*> >
  socket_bStatesNeighbors;

  /// the socket to the data handle of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// the socket to the data handle of arrays of flags specifying if the
  /// time jacobian contribution of certain variables in boundary states
  /// have to be discarded
  Framework::DataSocketSink<std::vector<bool> > socket_discardTimeJacob;


  /// RealVector holding the BC variables
  RealVector m_variables;

  /// factor that multiplies the coefficient for the diagonal entry
  /// in the linear system matrix
  CFreal _diagCoeffFactor;

   /// input data file name
  std::string m_infile;
  
  /// look up table for u(y)
  std::vector< Common::LookUpTable< CFreal, CFreal>* > m_lookupState;

}; // end of class SuperInletInterpImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SuperInletInterpImpl_hh

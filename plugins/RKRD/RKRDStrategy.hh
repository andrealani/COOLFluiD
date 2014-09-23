#ifndef COOLFluiD_RKRD_RKRDStrategy_hh
#define COOLFluiD_RKRD_RKRDStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "MathTools/MatrixInverter.hh"

#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Mario Ricchiuto
/// @author Tiago Quintino
class RKRDStrategy : public FluctSplit::FluctuationSplitStrategy {

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  RKRDStrategy(const std::string& name);

  /// Destructor.
  virtual ~RKRDStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args );

  /// Set up private data and data
  virtual void setup();

  /// Unsetup the private data and data of the aggregated classes
  /// in this strategy after the  processing phase
  virtual void unsetup();

  /// function called at each prepare phase before the computation of the RHS
  virtual void prepare ();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Compute the fluctuation
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation(std::vector<RealVector>& residual);

protected:

  /// computes the physical flux
  void flux ( RealVector& flux, const RealVector& state, const RealVector& normal );

  /// splits the jacobian into positive and negative eigen vectors projected along normal
  void splitJacobian( RealMatrix& jacobPlus,
                      const RealVector& state,
                      const RealVector& normal );
private: // data

  /// socket for states of k Runge-Kutta stages
  Framework::DataSocketSink< RealMatrix > socket_kstates;
  /// datahandle to socket
  Framework::DataHandle< RealMatrix > dh_kstates;

  /// cell volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  /// datahandle to cell volumes
  Framework::DataHandle< CFreal > dh_volumes;

//  /// The socket to use in this strategy for the update coefficients
//  Framework::DataSocketSink<CFreal> socket_updateCoeff;
//  /// datahandle to socket
//  Framework::DataHandle< CFreal > dh_update_coeff;

  /// lumping factor parameter used in switching
  /// from full lumping to partial lumping ( either one or zero )
  bool lump;

  /// current k-step
  CFuint kstep;

  /// RK order
  CFuint order;

  /// time step
  CFreal dt;

  /// temporary storage of the flux
  RealVector vflux;

  /// matrix of alpha parameters ( got from RKRD method )
  RealMatrix alpha;

  /// matrix of beta parameters ( got from RKRD method )
  RealMatrix beta;

  /// vector for delta u_s
  std::vector<RealVector> du;

  /// normals
  std::vector<RealVector> normals;
  std::vector<RealVector> adim_normals;

  /// current kstae used on integration
  RealVector tkstate;

  RealMatrix sum_kplus;

  RealMatrix inv_k;

  RealVector u_tmp;

  std::vector<RealMatrix> kplus;

  RealMatrix right_eigenv;
  RealMatrix left_eigenv;

  RealVector eValues;
  RealVector eValuesP;

  MathTools::MatrixInverter*       inverter;

//  /// quadrature points
//  std::vector<RealVector> m_qd_pdata;

}; // class RKRDStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RKRD_RKRDStrategy_hh

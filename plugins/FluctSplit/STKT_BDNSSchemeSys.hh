#ifndef COOLFluiD_Numerics_FluctSplit_STKT_BDNSSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_BDNSSchemeSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/STKT_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTime B scheme for RDS space discretization
 * We consider it as the Space-Time LDA + theta * dissipation used for N scheme.
 * This approach is developped
 * In the thesis of Mario Ricchiuto (see pp 139)
 * @author Nadege Villedieu
 *
 */
class STKT_BDNSSchemeSys : public STKT_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STKT_BDNSSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STKT_BDNSSchemeSys();

  /**
   * Set up
   */
  void setup();

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);


  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  // static void defineConfigOptions(Config::OptionList& options);

 /// Returns the DataSockets that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
private:

/// temporary data for computation ofthe sum of the positive part of the "space-time"upwind coeff
  RealMatrix m_sumKplus;

  RealMatrix _invK;

  RealVector m_uTemp;

  std::vector<RealVector> m_diss;

  RealVector m_sumKplusU;

  std::vector<RealVector> m_phiN;

  RealVector m_phitot;

  RealVector m_absphitot;

  RealVector m_phiT;

  CFreal m_r0;

  CFreal m_beta;

  /// the socket stores the data of the damping coefficient
    Framework::DataSocketSink<CFreal> socket_dampingCoeff;
}; // end of class STKT_BDNSSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_BDNSSchemeSys_hh

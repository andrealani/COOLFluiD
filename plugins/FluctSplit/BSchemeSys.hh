#ifndef COOLFluiD_Numerics_FluctSplit_BSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_BSchemeSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterSys.hh"
#include "FluctSplit/BSchemeBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the N scheme for RDS space discretization
   *
   * @author Andrea Lani
   */
class BSchemeSys : public BSchemeBase<RDS_SplitterSys> {
public:

  /// Constructor
  explicit BSchemeSys(const std::string& name);

  /// Destructor
  virtual ~BSchemeSys();

  /// Setup this object with data depending on the mesh
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  ///  Distribute part of the residual
  virtual void distributePart(std::vector<RealVector>& residual);
  
  /// Compute all the contributions for the Picard jacobian
  virtual void computePicardJacob(std::vector<RealMatrix*>& jacob);
  
protected: // functions

  /// Compute the blending coefficients
  virtual void computeBlendingCoeff();

  /**
   * Add an isotropic dissipative term.
   */
  virtual void addExtraDissipation();

protected:

  RealVector                 _sumKminU;

  RealMatrix                 _sumKmin;

  RealMatrix                 _sumKplus;

  RealMatrix                 _k;

  RealMatrix                 _invKmin;
  
  RealMatrix                 _invKplus;
  
  RealVector                 _uInflow;

  RealVector                 _phi;

  RealVector                 _phiLDA;

  RealVector                 _temp;

  RealMatrix                 _betaLDA;
  
  RealMatrix _tempMat;
  
}; // end of class BSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BSchemeSys_hh

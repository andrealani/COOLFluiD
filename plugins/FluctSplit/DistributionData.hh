#ifndef COOLFluiD_Numerics_FluctSplit_DistributionData_hh
#define COOLFluiD_Numerics_FluctSplit_DistributionData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



      namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class holds some data needed for the fluctuation distribution
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API DistributionData : public Common::NonCopyable<DistributionData>
{
public: // functions

  /// Constructor.
  DistributionData() : Common::NonCopyable<DistributionData>()
  {
    states = CFNULL;
    tStates = CFNULL;
    cell = CFNULL;
    currBetaMat = CFNULL;
    cellID = 0;
    sourceTermID = 0;
    computeBetas = false;
    isPerturb = false;
    iVar = 0;
    sourceComputeGradients = false;
    needDiss = false;
}

  /// Destructor.
  ~DistributionData()
  {
    for (CFuint i = 0; i< gradients.size(); ++i) {
      deletePtr(gradients[i]);
    }
  }

  /// Setup the data in this object
  void setup()
  {
    using namespace COOLFluiD::Framework;

    nbeqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint maxStates = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

    phi.resize(nbeqs);
    phiS.resize(nbeqs);

    nodalST.resize(maxStates);
    for (CFuint i = 0; i < nodalST.size(); ++i)
    {
      nodalST[i].resize(nbeqs);
    }

    betaMats.resize(1);
    betaMats[0].resize(maxStates);
    for (CFuint i = 0; i <  betaMats[0].size(); ++i)
    {
      betaMats[0][i].resize(nbeqs, nbeqs);
    }

    // by default the current beta matrix points to the first entry in betaMats
    currBetaMat = &betaMats[0];
    cf_assert(currBetaMat.isNotNull());

    // array of gradients
    gradients.resize(nbeqs);
    for (CFuint i = 0; i< nbeqs; ++i) {
      gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
    }

    avState.resize(nbeqs);
  }

public: // data

  /// locally store the number of equations
  CFuint nbeqs;

  /// all states in current cell
  std::vector<Framework::State*>* states;

  /// transformed states in current cell
  std::vector<Framework::State*>* tStates;

  /// This is the states of the sub-cell high order discretization
  std::vector<Framework::State*>* subStates;

  /// Pointer to the current cell being processed.
  Framework::GeometricEntity* cell;

  /// current beta matrices
  Common::SafePtr<std::vector<RealMatrix> > currBetaMat;

  /// index of the current cell
  CFuint cellID;

  /// ID of the source term
  CFuint sourceTermID;

  /// flag telling if to compute the distribution coefficients
  bool computeBetas;

  /// flag telling if the residual is being perturbed
  bool isPerturb;

  /// ID of the variable currently perturbed
  CFuint iVar;

  /// flag telling if the source term computes the cell gradients
  /// (instead of the diffusive term computer)
  bool sourceComputeGradients;

  /// flag telling if we need to compute adissipation for the Space-TimeU
  bool needDiss;

  /// pointer to the temporary storage of the fluctuation term
  RealVector phi;

  /// temporary storage of the source term
  RealVector phiS;

  /// nodal source terms in the current cell ( this may not be needed )
  std::vector<RealVector> nodalST;

  /// beta matrices
  typedef std::vector<std::vector<RealMatrix> > BetaMatrices;
  BetaMatrices betaMats;

  /// array of gradients
  std::vector<RealVector*> gradients;
  
  /// average state in the cell
  RealVector avState;
  
  /// inward normals data ( not used yet )
  FluctSplit::InwardNormalsData* inward_normals;

  /// Pointer to the past residuals, contributed by the past states.
  /// RealVector is sized number of equations times number of states.
  RealVector*  past_residuals;

  /// Pointer to the intermediate residuals, contributed by the inter states.
  /// RealVector is sized number of equations times number of states.
  RealVector*  inter_residuals;

  /// Pointer to the past residuals, contributed by the past states,
  /// computed with 1st order scheme it is used with Space-Time blending
  /// RealVector is sized number of equations times number of states.
  RealVector* past_residuals_order1;

  /// Pointer to the time component of the temporary storage of the fluctuation term
  /// (\int_T \farc{\partial u}{\partial t}
  RealVector phi_time;

  /// Physical time that we are considering. This can be the time of the present
  /// or the time if the past. This is used for the unsteady source term
  CFreal time;

  /// This boolean is used to tell the time splitter that we are doing HO
  /// Then the area of the volume should be devided by 4
  bool isHO;

	//sub iterations
  CFint subiter;


  /// Booleanused in unsteady to see if we were initializing one ore two level
  bool isdbInit;


  /// Boolean used in high order unsteady to know if we are in the first iteration done 
  /// with P1 discretization
  bool isfirstP1;


  /// Matrix of the normals of the sub-tertra in a P2 tetra
  RealMatrix FacesMat_tet;

}; // class DistributionData

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_DistributionData_hh

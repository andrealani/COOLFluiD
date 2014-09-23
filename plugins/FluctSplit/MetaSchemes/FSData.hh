#ifndef COOLFluiD_FluctSplit_FSData_hh
#define COOLFluiD_FluctSplit_FSData_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
// #include "Framework/GeometricEntity.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DataHandle.hh"
#include "Framework/VarSetMatrixTransformer.hh"
#include "FluctSplit/DistributionData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// Storage of subelement data
/// @author Tiago Quintino
template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
struct SubElem
{
  SubElem() : residual        (PHYSICS::NEQS),
              nodal_normals   (ELEMGEO::NBSUBSTATES),
              states          (ELEMGEO::NBSUBSTATES)
  {
    nodal_areas.resize(ELEMGEO::NBSUBSTATES);
    for (CFuint i = 0; i < ELEMGEO::NBSUBSTATES; ++i)
      nodal_normals[i].resize(PHYSICS::DIM);
  }

  /// subelement residual
  RealVector residual;
  /// subelement nodal unit normals
  std::vector<RealVector> nodal_normals;
  /// subelement nodal areas
  std::vector<CFreal> nodal_areas;
  /// subelement states
  std::vector<Framework::State*> states;
};

/// Storage of quadrature point data
/// @author Tiago Quintino
template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
struct QdPt
{
  QdPt() : state      (RealVector(PHYSICS::NEQS)),
           node       (RealVector(PHYSICS::DIM),true),
           physdata   (PHYSICS::PHYSDATA),
           extra_vars (PHYSICS::XVARS),
           flux       (PHYSICS::NEQS, PHYSICS::DIM)
  {
    state.setSpaceCoordinates(&node);
    CFout << "node " << &node << " " << node << "\n";
  }

  /// quadrature point state variables
  Framework::State state;
  /// quadrature point nodal coordinates
  Framework::Node node;
  /// quadrature point extra variables to be linearized
  RealVector physdata;
  /// quadrature point extra variables to be linearized
  RealVector extra_vars;
  /// quadrature point fluxes
  RealMatrix flux;
};

/// class to share data between parts of the algorithm
/// @author Tiago Quintino
template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
struct FSData
{

  enum { NBNODES  = ELEMGEO::NBNODES  };
  enum { NBSTATES = ELEMGEO::NBSTATES };
  enum { NEQS     = PHYSICS::NEQS };
  enum { DIM      = PHYSICS::DIM };

  // information about subelement
  enum { NBSUBELEM   = ELEMGEO::NBSUBELEM };
  enum { NBSUBSTATES = ELEMGEO::NBSUBSTATES };

  typedef ELEMGEO ElemGeo_type;
  typedef PHYSICS Physics_type;
  typedef NUMERICS Neumerics_type;

  typedef SubElem<ELEMGEO,PHYSICS,NUMERICS> SubElem_type;
  typedef QdPt<ELEMGEO,PHYSICS,NUMERICS> QdPt_type;

  typedef typename NUMERICS::Splitter_type::template Splitter<ELEMGEO,PHYSICS> Splitter_type;
  typedef typename NUMERICS::Integration_type::template Integrator<ELEMGEO,PHYSICS> Integrator_type;

  /// Constructor
  FSData() : sub_elem (ELEMGEO::NBSUBELEM),
             qdpts    (Integrator_type::NBQDPT),
             nodes    (ELEMGEO::NBNODES),
             states   (ELEMGEO::NBSTATES),
             kplus    (ELEMGEO::NBSUBSTATES),
             kmin     (ELEMGEO::NBSUBSTATES),
             evalues  (ELEMGEO::NBSUBSTATES),
             updateCoeff     (CFNULL),
             method_data     (CFNULL)

  {
//     CFout << "+++++  sub_elem size " << sub_elem.size() << "\n";
//     CFout << "+++++  qdpts size " << qdpts.size() << "\n";
//     for (CFuint i = 0; i < Integrator_type::NBQDPT; ++i)
//     {
//       CFout << "satte [" << i << "] " << & qdpts[i].state << " " << qdpts[i].state << "\n";
//       CFout << "node [" << i << "] " << & qdpts[i].node << " " << qdpts[i].node << "\n";
//     }
    for (CFuint i = 0; i < ELEMGEO::NBSUBSTATES; ++i)
    {
      kplus[i].resize(PHYSICS::NEQS,PHYSICS::NEQS);
      kmin[i].resize(PHYSICS::NEQS,PHYSICS::NEQS);
      evalues[i].resize(PHYSICS::NEQS);
    }
  }

  /// storage of subelement structures
  std::vector<SubElem_type> sub_elem;
  /// storage of quadrature points
  std::vector<QdPt_type> qdpts;
  /// element nodes
  std::vector< Framework::Node* > nodes;
  /// element states
  std::vector< Framework::State* > states;
  /// transformed states in current subcell
  std::vector<Framework::State*>* t_states;
  /// temporary storage of kplus upwind coefficients
  std::vector<RealMatrix> kplus;
  /// temporary storage of kmin upwind coefficients
  std::vector<RealMatrix> kmin;
  /// temporary storage of eigen values
  std::vector<RealVector> evalues;
  /// input residual to change and return, NBSTATES x NEQS
  std::vector<RealVector> * return_residual;
  /// object that splits the element residual, distributing it by the states
  Splitter_type splitter;
  /// distribute varset
  Common::SafePtr<Framework::ConvectiveVarSet>   distribute_varset;
  /// linearzation varset
  Common::SafePtr<Framework::ConvectiveVarSet>   linear_varset;
  /// jacobian linearizer
  Common::SafePtr<Framework::JacobianLinearizer> jacob_linearizer;
  /// vector transformation from update to linear variables
  Common::SafePtr<Framework::VarSetTransformer>  update_to_linear;
  /// matrix transformation from linear to distribution variables
  Common::SafePtr<Framework::VarSetMatrixTransformer>  linear_to_distrib;
  /// datahandle to the update coefficients
  Framework::DataHandle< CFreal > updateCoeff;
  /// transformed linear states in current cell
  std::vector<Framework::State*> * linear_states;
  /// access to the Fluctsplit method data
  Common::SafePtr<FluctSplit::FluctuationSplitData> method_data;
  /// access to the distribution data
  Common::SafePtr<FluctSplit::DistributionData> dist_data;
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_FSData_hh

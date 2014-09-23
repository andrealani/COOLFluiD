#ifndef COOLFluiD_FluctSplit_FSMHO_hh
#define COOLFluiD_FluctSplit_FSMHO_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/static_assert.hpp>

#include "Common/Meta/Loop.hh"
#include "Common/CFLog.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/MetaSchemes/FSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a template strategy to compute the element
/// contribution of the Residual Distribution method
///
/// @todo Possible improvements to the performance of this class:
///       - use Lobatto quadrature in contour integration
///         to save interpolation at the nodes
///       - remove PlaceSubStates
///       - numerical vectors and matrices with static size
///       - storage of normals
///       - accumulating loops may be implemented with boost::accumulator
///       - non isoparametric elements P1P2
///       - avoid copy of states and nodes pointers (?)
///
/// @author Tiago Quintino
///
template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
class FSMHO : public COOLFluiD::FluctSplit::FluctuationSplitStrategy
{

  public: // typedefs

    enum { NBNODES  = ELEMGEO::NBNODES  };
    enum { NBSTATES = ELEMGEO::NBSTATES };
    enum { NEQS     = PHYSICS::NEQS };
    enum { DIM      = PHYSICS::DIM };

    typedef ELEMGEO  ElemGeo_type;
    typedef PHYSICS  Physics_type;
    typedef NUMERICS Numerics_type;

    typedef FSData < ELEMGEO, PHYSICS, NUMERICS> FSData_type;

    typedef typename NUMERICS::Splitter_type::template Splitter<ELEMGEO,PHYSICS> Splitter_type;
    typedef typename NUMERICS::Integration_type::template Integrator<ELEMGEO,PHYSICS> Integrator_type;

  public: // functions

  /// Constructor
  FSMHO (const std::string& name);

  /// Destructor
  virtual ~FSMHO();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Compute the fluctuation
  /// @pre residual size is number of states
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation (std::vector<RealVector>& residual);

  static std::string getSchemeName()
  {
    return ELEMGEO::GeoSF_type::getShapeName() +
            + "P" + Common::StringOps::to_str((int) ELEMGEO::GEOORDER)
            + "P" + Common::StringOps::to_str((int) ELEMGEO::SOLORDER)
            + "_" + Physics_type::PhysicalModel_type::getClassName()
            + "_" + Splitter_type::getClassName()
            + "_" + Numerics_type::getClassName() + Integrator_type::getClassName();
  };

  /// @return a std::string with details of the scheme
  static std::string details ()
  {
    std::ostringstream oss;
    oss << "FSMHO ELEMENT [" << getSchemeName() << "]\n";
    oss << "SHAPE    : [" << ELEMGEO::SHAPE     << "]\n";
    oss << "DIM      : [" << ELEMGEO::DIM       << "]\n";
    oss << "NBNODES  : [" << ELEMGEO::NBNODES   << "]\n";
    oss << "NBSTATES : [" << ELEMGEO::NBSTATES  << "]\n";
    oss << "GEOORDER : [" << ELEMGEO::GEOORDER  << "]\n";
    oss << "SOLORDER : [" << ELEMGEO::SOLORDER  << "]\n";
    oss << "NBSUBELEM: [" << ELEMGEO::NBSUBELEM << "]\n";
    return oss.str();
  }

  template < unsigned int SELEM >
  struct PlaceSubStates
  {
    template < unsigned int SSTATE, typename FSDATA >
    static inline void exec ( FSDATA& data )
    {
      data.sub_elem[SELEM].states[SSTATE] = data.states[ ELEMGEO::SolSF_type::template SubElem<SELEM>::template SNode<SSTATE>::ID ] ;
    }
  };

  struct DistributeSubElemResidual
  {
    template < unsigned int SELEM, typename FSDATA >
    static inline void exec ( FSDATA& data )
    {
      // distributes the subelement residual
      NUMERICS::template distribute_residual <SELEM, FSDATA> ( data );
    }
  };

  /// Prepares the subelements by placing the states and getting the normals
  /// @todo In the future the computation of the normals
  ///       should be made a Trait class that either computes
  ///       on the fly or gets them from the storage
  struct PrepareSubElem
  {
    template < unsigned int SELEM, typename FSDATA >
    static inline void exec ( FSDATA& data )
    {
      // place the states in the subelements
      Common::Loop1<PlaceSubStates<SELEM>, ELEMGEO::NBSUBSTATES>::run(data);
      // computes the nodal normals
      ELEMGEO::SolSF_type::template SubElem<SELEM>::compute_nodal_normals (data);
      // compute the consistent states
      // transform the update states to linear varset
      data.linear_states = data.update_to_linear->transform(&data.sub_elem[SELEM].states);
      // computes an average state in which jacobians will be linearized
      data.jacob_linearizer->linearize( *(data.linear_states) );
      // transformation from linear to consistent variables with the distribution varset
      data.t_states = data.linear_to_distrib->transformFromRef(data.linear_states);
    }
  };

private: // helper functions

  /// Prints each subelem normals
  void print_subelem_normals ()
  {
    for (CFuint i = 0; i < ELEMGEO::NBSUBELEM; ++i)
    {
      for (CFuint j = 0; j < ELEMGEO::NBSUBSTATES; ++j)
        CFout << "[" << m_fsdata.sub_elem[i].nodal_normals[j] << "]\n";
    }
  }

public: // data

    /// data to share between the algorithms
    FSData_type m_fsdata;
    /// handle to the update coefficient
    Framework::DataSocketSink< CFreal>  socket_updateCoeff;

}; // end of class FSMHO

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
FSMHO<ELEMGEO,PHYSICS,NUMERICS>::FSMHO (const std::string& name)
  : FluctuationSplitStrategy(name),
    socket_updateCoeff("updateCoeff")
{
  // some sanity checks
  BOOST_STATIC_ASSERT ((int) ELEMGEO::DIM               == (int) PHYSICS::DIM);
  BOOST_STATIC_ASSERT ((int) ELEMGEO::GeoSF_type::DIM   == (int) ELEMGEO::SolSF_type::DIM);
  BOOST_STATIC_ASSERT ((int) ELEMGEO::GeoSF_type::SHAPE == (int) ELEMGEO::SolSF_type::SHAPE);
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
FSMHO<ELEMGEO,PHYSICS,NUMERICS>::~FSMHO ()
{
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
void FSMHO<ELEMGEO,PHYSICS,NUMERICS>::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitStrategy::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
void FSMHO<ELEMGEO,PHYSICS,NUMERICS>::setup ()
{
  FluctuationSplitStrategy::setup();

  FluctuationSplitData& mdata = getMethodData();

  m_fsdata.method_data = & mdata;
  m_fsdata.dist_data   = & mdata.getDistributionData();

  m_fsdata.linear_varset     = mdata.getLinearVar();
  m_fsdata.distribute_varset = mdata.getDistribVar();
  m_fsdata.jacob_linearizer  = mdata.getLinearizer();
  m_fsdata.update_to_linear  = mdata.getUpdateToLinearVecTrans();
  m_fsdata.linear_to_distrib = mdata.getLinearToDistribMatTrans();

  m_fsdata.updateCoeff = socket_updateCoeff.getDataHandle();
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
void FSMHO<ELEMGEO,PHYSICS,NUMERICS>::unsetup ()
{
  m_fsdata.method_data = CFNULL;

  m_fsdata.distribute_varset = CFNULL;
  m_fsdata.jacob_linearizer  = CFNULL;
  m_fsdata.update_to_linear  = CFNULL;
  m_fsdata.linear_to_distrib = CFNULL;

  FluctuationSplitStrategy::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
FSMHO<ELEMGEO,PHYSICS,NUMERICS>::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO, typename PHYSICS, typename NUMERICS >
void FSMHO<ELEMGEO,PHYSICS,NUMERICS>::computeFluctuation (std::vector<RealVector>& residual)
{
  using namespace COOLFluiD::Framework;

  // check preconditions
  cf_assert (residual.size() == ELEMGEO::NBSTATES);

//   CFout << getSchemeName() << CFendl;
// CF_DEBUG_POINT;

  const DistributionData& ddata = *m_fsdata.dist_data;

/// @todo not using external nodes because we only do isoparemetric elements now
//   const std::vector<Node*>& nodes = * ddata.cell->getNodes();
//   cf_assert (nodes.size()  == ELEMGEO::NBNODES);
  const std::vector<State*>& states = * ddata.states;
  cf_assert (states.size() == ELEMGEO::NBSTATES);

  // put the cell states in the FSData
  for (CFuint i = 0; i < NBSTATES; ++i) { m_fsdata.states[i] = states[i]; }
//   for (CFuint i = 0; i < NBSTATES; ++i) {  CFout << "state " << *states[i] << "\n";  }
  // put the cell nodes in the FSData
  for (CFuint i = 0; i < NBNODES; ++i) { m_fsdata.nodes[i] = states[i]->getNodePtr(); }
//   for (CFuint i = 0; i < NBNODES; ++i) {  CFout << "node " << states[i]->getCoordinates() << "\n";  }
  // reset the residual to zero and cache it in the FSdata
  m_fsdata.return_residual = &residual;
  for (CFuint i = 0; i < NBSTATES; ++i) { residual[i] = 0.0; }

  // prepares the subelements by computin normals, volumes ...
  Common::Loop1<PrepareSubElem, ELEMGEO::NBSUBELEM>::run(m_fsdata);

  // computes the subelem fluctuation
  NUMERICS::compute_residual ( m_fsdata );
//   CF_DEBUG_EXIT;

  // distribute subelem fluctuation by all sub elements
  Common::Loop1<DistributeSubElemResidual, ELEMGEO::NBSUBELEM>::run(m_fsdata);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_FSMHO_hh

#ifndef COOLFluiD_FluctSplit_Nscheme_hh
#define COOLFluiD_FluctSplit_Nscheme_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/MetaSchemes/BaseSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// Template Nscheme class that is used as Trait class for the FSMHO strategy
/// @author Tiago Quintino
class NschemeT
{
  public:
  template < typename ELEMGEO , typename PHYSICS >
  class Splitter : public BaseSplitter < ELEMGEO, PHYSICS >
  {
    public: // typedefs

    enum { NBSTATES = ELEMGEO::NBSTATES };
    enum { NEQS     = PHYSICS::NEQS };

    public: // functions

    /// Constructor
    Splitter();
    /// Destructor
    ~Splitter();
    /// returns the class name
    static std::string getClassName () { return "Nscheme"; };
    /// distributes the subelement residual
    template < unsigned int SELEM, typename FSDATA > inline void distribute_intres ( FSDATA& data );
    /// distributes the subelement residual
    template < unsigned int SELEM, typename FSDATA > inline void distribute_linres ( FSDATA& data );

    private: // data

        /// sum of k ( plus | minus ) times u
        RealVector sk_u;
        /// sum of k ( plus | minus )
        RealMatrix sk;
        /// inverse of sum of k
        RealMatrix inv_k;
        /// upwind inlet state
        RealVector u_in;
        /// matrix inverter
        MathTools::MatrixInverterT<NEQS> inverter;

    private: // helper classes

      /// Adds the subelement nodal residual to the master element residual
      template < unsigned int SELEM >
      struct NResidual
      {
        template < unsigned int SSTATE >
        static inline void exec ( std::vector<RealVector>& res, const std::vector<RealMatrix>& kplus, const std::vector<Framework::State*>& states, const RealVector& u_in )
        {
          res [ ELEMGEO::SolSF_type::template SubElem<SELEM>::template SNode<SSTATE>::ID ] += kplus[SSTATE] * ( *states[SSTATE] - u_in );
        }
      };

  }; // end class Splitter

}; // end of class NschemeT

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO , typename PHYSICS >
NschemeT::Splitter<ELEMGEO,PHYSICS>::Splitter()
 : sk_u(NEQS), sk(NEQS,NEQS), inv_k(NEQS,NEQS), u_in(NEQS)
{
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO , typename PHYSICS >
NschemeT::Splitter<ELEMGEO,PHYSICS>::~Splitter()
{
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO , typename PHYSICS >
template < unsigned int SELEM, typename FSDATA >
inline void NschemeT::Splitter<ELEMGEO,PHYSICS>::distribute_intres ( FSDATA& data )
{
  using namespace Framework;

  PRINT_DEBUG ( CFout << "NschemeT::distribute_intres() selem [" << SELEM << "] \n" << CFendl; )

  const std::vector<RealMatrix>& kplus =  data.kplus;
  const RealVector&       sub_residual =  data.sub_elem[SELEM].residual;
  const std::vector<State*>&    states = *data.t_states;
  std::vector<RealVector>&    residual = *data.return_residual;

  sk   = kplus[0];
  sk_u = kplus[0] * *states[0];
  for ( CFuint i = 1; i < ELEMGEO::NBSUBSTATES; ++i )
  {
    sk   += kplus[i];
    sk_u += kplus[i] * *states[i];
  }

  inverter.invert(sk,inv_k);
  // compute inlet state
  u_in = inv_k * (sk_u - sub_residual);
  // res = kplus ( u - u_in)
  Common::Loop4<NResidual<SELEM>, ELEMGEO::NBSUBSTATES>::run(residual,kplus,states,u_in);
}

//////////////////////////////////////////////////////////////////////////////

template < typename ELEMGEO , typename PHYSICS >
template < unsigned int SELEM, typename FSDATA >
inline void NschemeT::Splitter<ELEMGEO,PHYSICS>::distribute_linres ( FSDATA& data )
{
  using namespace Framework;

// CFout << "NschemeT::distribute_linres() selem [" << SELEM << "] \n" << CFendl;

  const std::vector<RealMatrix>& kmin   =  data.kmin;
  const std::vector<RealMatrix>& kplus  =  data.kplus;
  const std::vector<State*>&     states = *data.t_states;
  std::vector<RealVector>&     residual = *data.return_residual;

  // compute sum kmin times U
  sk   = kmin[0];
  sk_u = kmin[0] * *states[0];
  for ( CFuint i = 1; i < ELEMGEO::NBSUBSTATES; ++i )
  {
    sk   += kmin[i];
    sk_u += kmin[i] * *states[i];
  }

  // invert sum kmin U
  inverter.invert(sk,inv_k);
  // compute inlet state
  u_in = inv_k * sk_u;
  // res = kplus ( u - u_in)
  Common::Loop4<NResidual<SELEM>, ELEMGEO::NBSUBSTATES>::run(residual,kplus,states,u_in);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_Nscheme_hh


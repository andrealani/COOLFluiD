#include <boost/progress.hpp>
#include <vector>

#include "Common/SwapEmpty.hh"

#include "AutoTemplateLoader/BaseRunLib.hh"

using namespace COOLFluiD;

template < typename SCHEME >
class RunNscheme : public BaseRunLib
{
  public:

    typedef typename SCHEME::Matrix_t Matrix_t;
    typedef typename SCHEME::Vector_t Vector_t;

    enum { NSTATES = SCHEME::E_NSTATES };
    enum { NEQS    = SCHEME::E_NEQS };

  RunNscheme ()
  {
    cell_kplus.resize(NSTATES);
    cell_states.resize(NSTATES);
    state_residuals.resize(NSTATES);
  }

  protected: // helper functions

  void compute()
  {    
    CFLogInfo ( "\nRunning benchmark\n" );
    for( CFuint e = 0; e < nbelems; ++e)
    {
      for( CFuint s = 0; s < NSTATES; ++s)
      {
        CFuint idx = e*NSTATES+s;
        cell_kplus  [s] = &gkplus[idx];
        cell_states [s] = &gstates[idx];
        state_residuals[s] = &gsresiduals[idx];
      }
      m_scheme.distribute(cell_kplus, cell_states, gresidual[e], state_residuals);
    }
  }

  void init()
  {
    CFLogInfo ( "\nAllocating benchmark data ..." );
    boost::progress_display progress (4);
    gkplus.resize(nbelems*NSTATES);      ++progress;
    gstates.resize(nbelems*NSTATES);     ++progress;
    gresidual.resize(nbelems);           ++progress;
    gsresiduals.resize(nbelems*NSTATES); ++progress;

    CFLogInfo ( "\nInitializing benchmark data ..." );
    boost::progress_display progress_init (nbelems);
    for( CFuint i = 0; i < nbelems; ++i)
    {
      ++progress_init;
      for( CFuint s = 0; s < NSTATES; ++s)
      {
        CFuint idx = i*NSTATES+s;
        SCHEME::resize(gkplus[idx],NEQS,NEQS);
        SCHEME::init(gkplus[idx]);
//         gkplus[idx].resize(NEQS,NEQS);
//         for( CFuint e = 0; e < NEQS; ++e) gkplus[idx](e,e) = 1.0; // identity matrix
        SCHEME::resize(gstates[idx],NEQS);
        SCHEME::init(gstates[idx]);
        SCHEME::resize(gsresiduals[idx],NEQS);
        SCHEME::init(gsresiduals[idx]);
      }
      SCHEME::resize(gresidual[i],NEQS);
    }
  }

  void finalize()
  {
    using namespace COOLFluiD::Common;
    SwapEmpty (gkplus);
    SwapEmpty (gstates);
    SwapEmpty (gresidual);
    SwapEmpty (gsresiduals);
  }

  protected:

    void do_init ();

  private: // data

    /// global data
    std::vector<Matrix_t> gkplus;
    std::vector<Vector_t> gstates;
    std::vector<Vector_t> gresidual;
    std::vector<Vector_t> gsresiduals;
    
    /// cell wise data
    std::vector<Matrix_t*> cell_kplus;
    std::vector<Vector_t*> cell_states;
    std::vector<Vector_t*> state_residuals;

    /// the scheme
    SCHEME m_scheme;
};




#ifndef COOLFluiD_Post_Process_LSF_hh
#define COOLFluiD_Post_Process_LSF_hh

#include "PostProcess.hh"



using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;


namespace COOLFluiD {

namespace RadiativeTransfer {

/////////////////////////////////////////////////////////////////////////////
//
// Calculates the radiative heat flux inside, discarding the boundary values
// This post-process routine assumes regular geometry
//
/////////////////////////////////////////////////////////////////////////////

class PostProcessLSFexp : public PostProcess
{
public:

    PostProcessLSFexp(const string &name);
    void configure ( Config::ConfigArgs& args );

    /// Default destructor
    ~PostProcessLSFexp();
    static void defineConfigOptions(Config::OptionList& options);

    void runPostProcess(vector<CFreal> dataVector);
    private:
    bool m_mirrorPoints;
    vector<string> m_regectionTRS;
    vector<CFreal> m_direction;
    vector<CFreal> m_regectionDists;
    CFuint m_nbCoeffs;
    CFuint m_nbPointsPlot;
};



}
}

#endif

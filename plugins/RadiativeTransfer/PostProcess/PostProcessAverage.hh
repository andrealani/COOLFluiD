#ifndef COOLFluiD_Post_Process_Average_hh
#define COOLFluiD_Post_Process_Average_hh

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

class PostProcessAverage : public PostProcess
{
public:

    PostProcessAverage(const string &name);
    void configure ( Config::ConfigArgs& args );

    /// Default destructor
    ~PostProcessAverage();
    static void defineConfigOptions(Config::OptionList& options);

    void runPostProcess(DataHandle<CFreal> dataVector);
    private:
    vector<string> m_regectionTRS;
    vector<CFreal> m_direction;
    vector<CFreal> m_regectionDists;
    CFuint m_nbIntervals;
    bool m_outNormals;
    bool m_mirrorOnOrigin;
    bool m_isRevolutionVolume3d;
};



}
}

#endif

#ifndef COOLFluiD_Post_Process_Null_hh
#define COOLFluiD_Post_Process_Null_hh

#include "PostProcess.hh"

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;


namespace COOLFluiD {

namespace RadiativeTransfer {

/////////////////////////////////////////////////////////////////////////////
//
// Null PostProcess Routine
//
/////////////////////////////////////////////////////////////////////////////

class PostProcessNull : public PostProcess
{
public:

    PostProcessNull(const string &name):PostProcess(name){;}
    void configure ( Config::ConfigArgs& args ){;}

    /// Default destructor
    ~PostProcessNull(){;}
    static void defineConfigOptions(Config::OptionList& options){;}

    void runPostProcess(DataHandle<CFreal> dataVector){;}
};

}
}

#endif

#ifndef COOLFluiD_Post_Process_LSF_hh
#define COOLFluiD_Post_Process_LSF_hh

#include "PostProcess.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

/////////////////////////////////////////////////////////////////////////////
//
// Calculates the radiative heat flux inside, discarding the boundary values
// This post-process routine assumes regular geometry
//
/////////////////////////////////////////////////////////////////////////////

class PostProcessLSF : public PostProcess
{
public:

  PostProcessLSF(const std::string &name);
  void configure ( Config::ConfigArgs& args );
  
  /// Default destructor
  ~PostProcessLSF();
  static void defineConfigOptions(Config::OptionList& options);
  
  void runPostProcess(Framework::DataHandle<CFreal> dataVector);
private:
  bool m_isRevolutionVolume3d;
  bool m_mirrorOnOrigin;
  std::vector<std::string> m_regectionTRS;
  std::vector<CFreal> m_direction;
  std::vector<CFreal> m_regectionDists;
  CFuint m_nbCoeffs;
  CFuint m_nbPointsPlot;
  bool m_outNormals;
  CFreal m_scalingFactor;
  
};

  
}
}

#endif

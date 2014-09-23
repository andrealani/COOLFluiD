#ifndef COOLFluiD_Muffin_SystemMITReMPLaS_hh
#define COOLFluiD_Muffin_SystemMITReMPLaS_hh

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"

#include "PLaS/PLaSTracking.hh"
#include "Muffin/SystemPLaS.hh"

namespace COOLFluiD {
  namespace Framework { class Node;          }
  namespace PLaS      { class PLaSTracking;  }
  namespace Muffin    { struct GasOnSurface; }
}

namespace COOLFluiD {
  namespace Muffin {


/// Description of bubble in flow medium
struct Bubble {
  double x,y,z;  // original position
  double u,v,w;  // initial velocity
  double d;      // diameter
  double t;      // temperature
};


/// Class for using PLaS (dispersed two-phase flow solver) by gas-evolution
/// driven bubbles
class SystemMITReMPLaS : public SystemPLaS {

 public:  // non-virtual functions

  /// MITReMPLaS system constructor
  SystemMITReMPLaS(const std::string& name);

  /// MITReMPLaS system destructor
  ~SystemMITReMPLaS() {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);


 public:  // virtual functions implementations

  /// Set private data before processing phase
  void setup();

  /// Iterate over the system
  void execute();


 private:  // bubble-related functions

  /// Generate a random integer between min and max
  inline int generateRInteger(int min=0, int max=1) const
  {
    double r = (double) max - (double) min + 1.;
    return min + (int)(r*rand()/(RAND_MAX+1.));
  }

  /// Generate a random double between min and max
  inline double generateRDouble(double min=0., double max=1.) const
  {
    int r = generateRInteger(0,100000);
    return (double)(r/100000.) * (max-min) + min;
  }

  /// Polar Box-Muller transformation for random Gaussian distribution
  ///
  /// (c) Copyright 1994, Everett F. Carter Jr.
  /// Permission is granted by the author to use this software for any
  /// application provided this copyright notice is preserved.
  /// input variables: mean m, standard deviation s
  double generateRGaussian(double m, double s) const;

  /// Generate new bubble
  Bubble generateNewBubble(
    Framework::Node* p1, Framework::Node* p2, Framework::Node* p3 );


 private:  // functions

  /// Write surface gas tracking/bubbles information
  void writeAttachedBubbles(
    const std::string& fname, const unsigned i, const double t );


 private:  // sockets

  /// Socket to access surface gas tracking (per boundary TRS, per boundary element)
  Framework::DataSocketSink< std::vector< GasOnSurface > > s_gasonsurf;


 public:  // sockets

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r(SystemPLaS::needsSockets());
    r.push_back(&s_gasonsurf);
    return r;
  }

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > r(SystemPLaS::providesSockets());
    r.push_back(&s_mn_voidfraction);
    return r;
  }


 public:  // data

  /// Gas evolution surface gas fraction maximum (configurable)
  double m_ge_fmax;

  /// Bubbles distribution type (as string) (configurable)
  std::string m_bubbles_type_str;

  /// Bubbles distribution mean (configurable)
  double m_bubbles_mu;

  /// Bubbles distribution standard deviation (configurable)
  double m_bubbles_sig;

  /// Bubbles distribution minimum diameter (configurable)
  double m_bubbles_min;

  /// Bubbles distribution type
  int m_bubbles_type;

  /// Generated bubbles temperature
  double m_bubbles_t;

  /// Bubbles on surface tracking (per TRS, per b. element)
  std::vector< std::vector< Bubble > > m_vbubbles;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_SystemMITReMPLaS_hh


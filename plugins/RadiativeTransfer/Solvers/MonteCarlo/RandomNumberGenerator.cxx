#include <boost/random.hpp>
#include "RandomNumberGenerator.hh"
#include "Common/COOLFluiD.hh"
namespace COOLFluiD {

namespace RadiativeTransfer {
  using namespace std;

  CFreal RandomNumberGenerator::uniformRand(const CFreal i0, const CFreal i1){
    boost::uniform_real<CFreal> uniformDist(i0,i1);
    boost::variate_generator<typeGenerator&, boost::uniform_real<CFreal> >
             uniform(m_generator, uniformDist);
    return uniform();
  }

  void RandomNumberGenerator::seed(CFuint seedNumber){
    m_generator.seed(seedNumber);
  }
}
}

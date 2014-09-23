#ifndef RANDOMNUMBERGENERATOR_HH
#define RANDOMNUMBERGENERATOR_HH

#include "MathTools/MathFunctions.hh"
#include <boost/random.hpp>
#include "Common/COOLFluiD.hh"
#include <vector>

/*  Wrapper class for the Boost Random library
*   This works for the 1.42 version, for the new version,
*   the interface is broken
*   Authors: Pedro Santos
*/
namespace COOLFluiD {

namespace RadiativeTransfer {

typedef boost::mt19937 typeGenerator; //Marsenne Twister generator

class RandomNumberGenerator{

public:

  template<typename Tout>
  void sphereDirections(CFuint dim, Tout &directions);

  template<typename Tin, typename Tout>
  void hemiDirections(CFuint dim , Tin faceNormals, Tout &directions);

  CFreal uniformRand(const CFreal i0=0., const CFreal i1=1.);

  void seed(CFuint seedNumber);

private:

  typeGenerator m_generator;
};

template<class Tout>
void RandomNumberGenerator::sphereDirections(CFuint dim, Tout &directions){
  boost::uniform_on_sphere<CFreal, std::vector<CFreal> > uniformOnSphere(dim);
  boost::variate_generator<typeGenerator&, boost::uniform_on_sphere<CFreal , std::vector<CFreal> > >
       randSphere(m_generator, uniformOnSphere);
  std::vector<CFreal> temp=randSphere();
  for (CFuint i=0; i<dim;++i){
    directions[i]=temp[i];
  }
}

template<typename Tin, typename Tout>
void RandomNumberGenerator::hemiDirections(CFuint dim , Tin faceNormals, Tout &directions){
  //generate spherical directions;
  sphereDirections(dim, directions);
  //if the direction is in the wrong half of the sphere
  CFreal mult=(MathTools::MathFunctions::innerProd(directions,faceNormals)<0.)?-1:1;
  //reverse the directions
  for(CFuint i=0; i<dim;++i){
    directions[i]*=mult;
  }
}


}
}

#endif // RANDOMNUMBERGENERATOR_HH

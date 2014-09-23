#include "Acoustic_ArtDiffStrategy_HO.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/InwardNormalsData.hh"

////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Acoustic_ArtDiffStrategy_HO,
                       FluctuationSplitData,
                       ArtificialDiffusionStrategy,
                       FluctSplitModule>
                       acousticArtHODiffStrategyProvider("AcousticHO");

////////////////////////////////////////////////////////////////////////////

void Acoustic_ArtDiffStrategy_HO::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Artnu","Set up the cutoff frequency of the filter.");
}

////////////////////////////////////////////////////////////////////////////

Acoustic_ArtDiffStrategy_HO::Acoustic_ArtDiffStrategy_HO(const std::string& name) :
  ArtificialDiffusionStrategy(name)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_nu = 1e-5;
  setParameter("Artnu",&m_nu);
}

////////////////////////////////////////////////////////////////////////////

Acoustic_ArtDiffStrategy_HO::~Acoustic_ArtDiffStrategy_HO()
{
}

////////////////////////////////////////////////////////////////////////////


void Acoustic_ArtDiffStrategy_HO::setup()
{
  CFAUTOTRACE;

 // first call parent method
  ArtificialDiffusionStrategy::setup();

  m_dim = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_stateNormal.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
      m_stateNormal[iState].resize(m_dim);
  }

}

////////////////////////////////////////////////////////////////////////////

/**
 * 
 * @param residual 
 */
void Acoustic_ArtDiffStrategy_HO::addArtificialDiff(std::vector<RealVector>& residual)
{

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbNodes = residual.size();
  const CFuint sbNodes = nbNodes/2.;


  DistributionData& ddata = getMethodData().getDistributionData();
  DataHandle<bool> isBState = socket_isBState.getDataHandle();
  vector<State*>& states = *ddata.states;
  const CFreal volume = ddata.cell->computeVolume();

//   RealVector c;
//   c.resize(sbNodes);
  CFreal Rs;
  RealVector Damping;
  Damping.resize(nbEqs);
  RealVector SmoothRes;
  SmoothRes.resize(nbEqs);
  RealVector r;
  r.resize(sbNodes);

  CFuint subnodes[3];

  const CFuint cell_ID = getMethodData().getDistributionData().cellID;

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[cell_ID]);

  Rs = m_nu/volume;


// Triangle 1 : nodes 0-3-5

  subnodes[0] = 0;
  subnodes[1] = 3;
  subnodes[2] = 5;

  for (CFuint iState = 0; iState < sbNodes; ++iState) {
    residual[iState] = 0.;
  }

   for (CFuint node = 0; node < sbNodes; ++node) {
     const Node& currnode = states[node]->getCoordinates();
     if(!isBState[states[node]->getLocalID()]) {

       for (CFuint other = 0; other < sbNodes; ++other) {
         r[other] =0.;
       }
 
       for (CFuint other = 0; other < sbNodes; ++other) {
          const Node& othernode = states[other]->getCoordinates();
          for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
            r[other] += (currnode[iDim]-othernode[iDim])*(currnode[iDim]-othernode[iDim]);
         }
 	 r[other] = -sqrt(r[other]);
       }

       CFreal sum = 0.;
       for (CFuint other = 0; other < sbNodes; ++other)
         sum += r[other];
       r[node] = -sum;

       for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          Damping[iEq] = 0.;
          SmoothRes[iEq] = 0.;

          for (CFuint other = 0; other < sbNodes; ++other) {
            SmoothRes[iEq] += r[other]*((*states[other])[iEq]);
          }
          Damping[iEq] = Rs*SmoothRes[iEq];
          residual[node][iEq] = Damping[iEq];
        }
      }
      else {
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          residual[node][iEq] = 0.0;
      }
    }


// Triangle 2 : nodes 3-1-4

  subnodes[0] = 3;
  subnodes[1] = 1;
  subnodes[2] = 4;

  for (CFuint iState = 0; iState < sbNodes; ++iState) {
    residual[iState] = 0.;
  }

   for (CFuint node = 0; node < sbNodes; ++node) {
     const Node& currnode = states[node]->getCoordinates();
     if(!isBState[states[node]->getLocalID()]) {

       for (CFuint other = 0; other < sbNodes; ++other) {
         r[other] =0.;
       }
 
       for (CFuint other = 0; other < sbNodes; ++other) {
          const Node& othernode = states[other]->getCoordinates();
          for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
            r[other] += (currnode[iDim]-othernode[iDim])*(currnode[iDim]-othernode[iDim]);
         }
 	 r[other] = -sqrt(r[other]);
       }

       CFreal sum = 0.;
       for (CFuint other = 0; other < sbNodes; ++other)
         sum += r[other];
       r[node] = -sum;

       for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          Damping[iEq] = 0.;
          SmoothRes[iEq] = 0.;

          for (CFuint other = 0; other < sbNodes; ++other) {
            SmoothRes[iEq] += r[other]*((*states[other])[iEq]);
          }
          Damping[iEq] = Rs*SmoothRes[iEq];
          residual[node][iEq] = Damping[iEq];
        }
      }
      else {
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          residual[node][iEq] = 0.0;
      }
    }

// Triangle 3 : nodes 5-4-2

  subnodes[0] = 5;
  subnodes[1] = 4;
  subnodes[2] = 2;

  for (CFuint iState = 0; iState < sbNodes; ++iState) {
    residual[iState] = 0.;
  }

   for (CFuint node = 0; node < sbNodes; ++node) {
     const Node& currnode = states[node]->getCoordinates();
     if(!isBState[states[node]->getLocalID()]) {

       for (CFuint other = 0; other < sbNodes; ++other) {
         r[other] =0.;
       }
 
       for (CFuint other = 0; other < sbNodes; ++other) {
          const Node& othernode = states[other]->getCoordinates();
          for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
            r[other] += (currnode[iDim]-othernode[iDim])*(currnode[iDim]-othernode[iDim]);
         }
 	 r[other] = -sqrt(r[other]);
       }

       CFreal sum = 0.;
       for (CFuint other = 0; other < sbNodes; ++other)
         sum += r[other];
       r[node] = -sum;

       for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          Damping[iEq] = 0.;
          SmoothRes[iEq] = 0.;

          for (CFuint other = 0; other < sbNodes; ++other) {
            SmoothRes[iEq] += r[other]*((*states[other])[iEq]);
          }
          Damping[iEq] = Rs*SmoothRes[iEq];
          residual[node][iEq] = Damping[iEq];
        }
      }
      else {
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          residual[node][iEq] = 0.0;
      }
    }

// Triangle 4 : nodes 4-5-3

  subnodes[0] = 4;
  subnodes[1] = 5;
  subnodes[2] = 3;

  for (CFuint iState = 0; iState < sbNodes; ++iState) {
    residual[iState] = 0.;
  }

   for (CFuint node = 0; node < sbNodes; ++node) {
     const Node& currnode = states[node]->getCoordinates();
     if(!isBState[states[node]->getLocalID()]) {

       for (CFuint other = 0; other < sbNodes; ++other) {
         r[other] =0.;
       }
 
       for (CFuint other = 0; other < sbNodes; ++other) {
          const Node& othernode = states[other]->getCoordinates();
          for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
            r[other] += (currnode[iDim]-othernode[iDim])*(currnode[iDim]-othernode[iDim]);
         }
 	 r[other] = -sqrt(r[other]);
       }

       CFreal sum = 0.;
       for (CFuint other = 0; other < sbNodes; ++other)
         sum += r[other];
       r[node] = -sum;

       for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          Damping[iEq] = 0.;
          SmoothRes[iEq] = 0.;

          for (CFuint other = 0; other < sbNodes; ++other) {
            SmoothRes[iEq] += r[other]*((*states[other])[iEq]);
          }
          Damping[iEq] = Rs*SmoothRes[iEq];
          residual[node][iEq] = Damping[iEq];
        }
      }
      else {
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          residual[node][iEq] = 0.0;
      }
    }

  }

////////////////////////////////////////////////////////////////////////////

  }// End namespace FluctSplit

}// End namespace COOLFluiD

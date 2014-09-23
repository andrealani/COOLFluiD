
#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "NavierStokes/NavierStokes.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/SpaceMethodData.hh"

#include "GradComputer.hh"
#include "FluctSplit/MaximumGradientMagnitude.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace GradComputer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MaximumGradientMagnitude,
                      DataProcessingData,
		                  GradComputerFSModule>
aMaximumGradientMagnitudeProvider("MaximumGradientMagnitude");

//////////////////////////////////////////////////////////////////////////////

void MaximumGradientMagnitude::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

MaximumGradientMagnitude::MaximumGradientMagnitude( const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  m_sockets()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

MaximumGradientMagnitude::~MaximumGradientMagnitude()
{
}

//////////////////////////////////////////////////////////////////////////////

void MaximumGradientMagnitude::setup()
{
  DataProcessingCom::setup(); // first call setup of parent class

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> grads = m_sockets.getSocketSource<CFreal>("MaximumGradientMagnitude")->getDataHandle();
  grads.resize(PhysicalModelStack::getActive()->getNbEq()*states.size());
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
MaximumGradientMagnitude::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MaximumGradientMagnitude::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();
  return result;
}


//////////////////////////////////////////////////////////////////////////////

void MaximumGradientMagnitude::execute()
{
  CFAUTOTRACE;
  CFout << "Computing Maximum Gradient Magnitudes.\n";

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> grads = m_sockets.getSocketSource<CFreal>("MaximumGradientMagnitude")->getDataHandle();
  CFreal nbstates = states.size();
  if (grads.size()!=nbstates) grads.resize(nbEqs*nbstates);
  grads=0.;
  std::vector<CFreal> sumvol(nbstates,0);

  // builder for standard TRS GeometricEntity's that will be used
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setup();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.trs = getCurrentTRS();
  CFuint nbCells = getCurrentTRS()->getLocalNbGeoEnts();

  // cell loop to contribute
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    geoData.idx = iCell;
    GeometricEntity & currCell = *geoBuilder.buildGE();
    CFreal vol=currCell.computeVolume();
    RealVector centroid(nbDim);
    centroid=0.;
    for(int iState=0; iState<(const int)currCell.nbStates(); ++iState)
      for(int iDim=0; iDim<(const int)(nbDim); ++iDim)
        centroid[iDim]+=currCell.getState(iState)->getCoordinates()[iDim];
    centroid/=(CFreal)currCell.nbStates();
    //centroid=currCell.computeCentroid();
    vector<RealVector> fcoords(1,centroid);

    const std::vector<RealMatrix> cellGradients= currCell.computeSolutionShapeFunctionGradients(fcoords);
    vector<Framework::State*>& cellStates=  (*(currCell.getStates()));

    // compute gradients
    std::vector< RealVector > gradvec(0);
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq){
      gradvec.push_back(RealVector(0.,nbDim));
      for (CFuint Istate = 0; Istate < currCell.nbStates() ; ++Istate)
        for (CFuint iDim = 0; iDim < nbDim ; ++iDim)
          gradvec[iEq][iDim] += (cellGradients[0])(Istate,iDim)*(*(cellStates[Istate]))[iEq];
    }
    RealVector gradmag(0.,nbEqs);
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      gradmag[iEq]=gradvec[iEq].norm2();

    // sum from elements to nodes
    for (CFuint Istate = 0; Istate < currCell.nbStates() ; ++Istate)
    {
      const CFuint LID=cellStates[Istate]->getLocalID();
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        grads[LID*nbEqs+iEq]+=gradmag[iEq]*vol;
      sumvol[LID]+=vol;
    }

    geoBuilder.releaseGE();
  }

  // divide by sum of cell volumes around the state
  for (CFuint Istate = 0; Istate < getCurrentTRS()->getNbStatesInTrs() ; ++Istate) {
    const CFreal invvol=1./sumvol[Istate];
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      grads[Istate*nbEqs+iEq]*=invvol;
  }

}

//////////////////////////////////////////////////////////////////////////////

void MaximumGradientMagnitude::unsetup()
{
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void MaximumGradientMagnitude::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );
  m_sockets.createSocketSource<CFreal>("MaximumGradientMagnitude");
}

//////////////////////////////////////////////////////////////////////////////

  } //

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




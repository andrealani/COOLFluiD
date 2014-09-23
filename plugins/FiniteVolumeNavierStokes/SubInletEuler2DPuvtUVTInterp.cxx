#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubInletEuler2DPuvtUVTInterp.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/FilesystemException.hh"
#include "NavierStokes/Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletEuler2DPuvtUVTInterp, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletEuler2DPuvtUVTInterpFVMCCProvider("SubInletEuler2DPuvtUVTInterpFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DPuvtUVTInterp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Datafile","input file with u(y),v(y),T(y)");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler2DPuvtUVTInterp::SubInletEuler2DPuvtUVTInterp
(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _lookupTableU(),
  _lookupTableV(),
  _lookupTableT()
{
  addConfigOptionsTo(this);

  _datafile = "";
  setParameter("Datafile",&_datafile);
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler2DPuvtUVTInterp::~SubInletEuler2DPuvtUVTInterp()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DPuvtUVTInterp::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DPuvtUVTInterp::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const CFreal yCoord = 0.5*(ghostState->getCoordinates()[YY] +
			     innerState->getCoordinates()[YY]);

  CFreal uIn = _lookupTableU.get(yCoord);
  CFreal vIn = _lookupTableV.get(yCoord);
  CFreal TIn = _lookupTableT.get(yCoord);

  uIn /= _varSet->getModel()->getVelRef();
  vIn /= _varSet->getModel()->getVelRef();
  TIn /= _varSet->getModel()->getTempRef();

  (*ghostState)[0] = (*innerState)[0];
  (*ghostState)[1] = 2.*uIn - (*innerState)[1];
  (*ghostState)[2] = 2.*vIn - (*innerState)[2];
  (*ghostState)[3] = 2.*TIn - (*innerState)[3];
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DPuvtUVTInterp::setup()
{
  FVMCC_BC::setup();

  // open the given file with the following format
  // nbPoints
  // y, u(y), v(y), T(y)

  ifstream file(_datafile.c_str());

  if (!file) {
    throw Common::FilesystemException
      (FromHere(), "SubInletEuler2DPuvtUVTInterp::setup() => trying to open file " + _datafile);
  }

  CFuint nbPoints = 0;
  file >> nbPoints;

  _lookupTableU.reserve(nbPoints);
  _lookupTableV.reserve(nbPoints);
  _lookupTableT.reserve(nbPoints);

  for (CFuint i = 0; i < nbPoints; ++i) {
    CFreal y = 0.0;
    CFreal u = 0.0;
    CFreal v = 0.0;
    CFreal T = 0.0;

    file >> y >> u >> v >> T;

    _lookupTableU.insert(y,u);
    _lookupTableV.insert(y,v);
    _lookupTableT.insert(y,T);
  }

  _lookupTableU.sortKeys();
  _lookupTableV.sortKeys();
  _lookupTableT.sortKeys();

  // output the table to file
  // ofstream* fileU =  new ofstream("Uy.dat");
  //   _lookupTableU.saveToStream(*fileU);
  //   delete fileU;

  //   ofstream* fileV =  new ofstream("Vy.dat");
  //   _lookupTableV.saveToStream(*fileV);
  //   delete fileV;

  //   ofstream* fileT =  new ofstream("Ty.dat");
  //   _lookupTableT.saveToStream(*fileT);
  //   delete fileT;

  _varSet = getMethodData().getUpdateVar().d_castTo<Physics::NavierStokes::Euler2DVarSet>();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

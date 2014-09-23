#include "StructMech/StructMech.hh"
#include "StructMech2DDisp.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMech2DDisp, ConvectiveVarSet, StructMechModule, 1>
structMech2DDispProvider("StructMech2DDisp");

//////////////////////////////////////////////////////////////////////////////

StructMech2DDisp::StructMech2DDisp(Common::SafePtr<Framework::BaseTerm> term) :
  ConvectiveVarSet(term),
  _model(CFNULL)
{
  vector<std::string> names(2);
  names[0] = "u";
  names[1] = "v";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

StructMech2DDisp::~StructMech2DDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMech2DDisp::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
result.resize(state.size());
extra.resize(6);

// Output the states
result[0] = state[0];
result[1] = state[1];

CFuint stateID = state.getLocalID();

CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
if(iter > 0){

  const std::string namespaceName = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  std::string dataHandleName = namespaceName + "_stress";
  bool exists = Framework::MeshDataStack::getActive()->getDataStorage()->checkData(dataHandleName);

  if(exists){
    DataHandle< RealVector> stress = MeshDataStack::getActive()->getDataStorage()->
      getData<RealVector>(dataHandleName);

    RealVector& stresses = stress[stateID];
    //Output the stresses
    extra[0] = stresses[0];
    extra[1] = stresses[1];
    extra[2] = stresses[2];
  }

  dataHandleName = namespaceName + "_strain";
  exists = Framework::MeshDataStack::getActive()->getDataStorage()->checkData(dataHandleName);

  if(exists){
    DataHandle< RealVector> strain = MeshDataStack::getActive()->getDataStorage()->
      getData<RealVector>(dataHandleName);

    RealVector& strains = strain[stateID];
    //Output the stresses
    extra[3] = strains[0];
    extra[4] = strains[1];
    extra[5] = strains[2];
  }
}
else{
  extra = 0.;
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> StructMech2DDisp::getExtraVarNames() const
{
  vector<std::string> names(6);
  names[0] = "sigmaxx";
  names[1] = "sigmayy";
  names[2] = "sigmaxy";
  names[3] = "eps_xx";
  names[4] = "eps_yy";
  names[5] = "eps_xy";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

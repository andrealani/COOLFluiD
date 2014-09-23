#include "MeshRigidMove/MeshRigidMove.hh"

#include "RigidMoveData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModelImpl.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<RigidMoveData>, RigidMoveData, MeshRigidMoveModule> nullRigidMoveComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void RigidMoveData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("ExpansionRatioY","Expansion Ratio in Y axis");
   options.addConfigOption< CFreal >("VY","Speed of translation along Y coord");
   options.addConfigOption< CFreal >("OX","Center of Rotation X coord");
   options.addConfigOption< CFreal >("OY","Center of Rotation Y coord");
   options.addConfigOption< CFreal >("VOY","Speed of translation of Y coord of the center of rotation");
   options.addConfigOption< CFreal >("VOX","Speed of translation of X coord of the center of rotation");
   options.addConfigOption< CFreal >("ExpansionRatioX","Expansion Ratio in X axis");
   options.addConfigOption< CFreal >("RotationSpeed","Rotation Speed");
   options.addConfigOption< CFreal >("VX","Speed of translation along X coord");
   //  options.addConfigOption< CFreal >("VZ","Speed of translation along Z coord");
}

//////////////////////////////////////////////////////////////////////////////

RigidMoveData::RigidMoveData(Common::SafePtr<Framework::Method> owner)
  : MeshAdapterData(owner),
//     _rotationCenter(PhysicalModelStack::getActive()->getDim()),
//     _rotationCenterSpeed(PhysicalModelStack::getActive()->getDim()),
//     _translationSpeed(PhysicalModelStack::getActive()->getDim()),
//     _expansionRatio(PhysicalModelStack::getActive()->getDim()),
    _rotationCenter(DIM_2D),
    _rotationCenterSpeed(DIM_2D),
    _translationSpeed(DIM_2D),
    _expansionRatio(DIM_2D),
    _rotationSpeed()
{
   addConfigOptionsTo(this);
    _expansionRatio[0] = 1.;
   setParameter("ExpansionRatioX",&_expansionRatio[0]);

    _expansionRatio[1] = 1.;
   setParameter("ExpansionRatioY",&_expansionRatio[1]);

    _rotationSpeed = 0.;
   setParameter("RotationSpeed",&_rotationSpeed);

    _rotationCenter[0] = 0.;
   setParameter("OX",&_rotationCenter[0]);

    _rotationCenterSpeed[0] = 0.;
   setParameter("VOX",&_rotationCenterSpeed[0]);

    _translationSpeed[0] = 0.;
   setParameter("VX",&_translationSpeed[0]);

//    if (PhysicalModelStack::getActive()->getDim()>1){
   _rotationCenter[1] = 0.;
   setParameter("OY",&_rotationCenter[1]);

   _rotationCenterSpeed[1] = 0.;
   setParameter("VOY",&_rotationCenterSpeed[1]);
   
   _translationSpeed[1] = 0.;
   setParameter("VY",&_translationSpeed[1]);
   
   //      }

   /// AL: this is HORRIBLE!   all RealVector should be vector<CFreal> and there should a "V" option with size == DIM
   
//     if (PhysicalModelStack::getActive()->getDim()>2){
//       _translationSpeed[2] = 0.;
//    setParameter("VZ",&_translationSpeed[2]);
//       }

}

//////////////////////////////////////////////////////////////////////////////

RigidMoveData::~RigidMoveData()
{
}

//////////////////////////////////////////////////////////////////////////////

void RigidMoveData::configure ( Config::ConfigArgs& args )
{
  MeshAdapterData::configure(args);
}


//////////////////////////////////////////////////////////////////////////////

void RigidMoveData::setup()
{
  MeshAdapterData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void RigidMoveData::unsetup()
{
  MeshAdapterData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////


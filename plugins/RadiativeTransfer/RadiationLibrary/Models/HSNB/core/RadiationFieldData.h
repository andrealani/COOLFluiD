#ifndef COOLFluiD_RadiativeTransfer_RADIATIONFIELDDATA_H
#define COOLFluiD_RadiativeTransfer_RADIATIONFIELDDATA_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Constants.h"
#include "Common/StringOps.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Common/CFMap.hh"
#include "Environment/CFEnv.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

class RadiationFieldData
{
public:
    ///Empty constructor
    RadiationFieldData();

    ///For init with Mesh States
    RadiationFieldData(const CFuint pressureID, const CFuint trID, const CFuint tvID, const CFuint emIndex, const Common::SafePtr<RealVector> stateVector, const CFuint &nbSpecies, const Common::SafePtr<RealVector> avogadroOvMM=NULL);

    //  NOTE: TODO
    CFreal* X();
    CFreal X(CFuint i);
    CFreal p();
    CFreal Tr();
    CFreal Tv();

    CFreal *getNumberDensities();

    ~RadiationFieldData();
    CFreal N(CFuint i) const;
    
private:
    RealVector m_numberDensity;
    Common::SafePtr<RealVector> m_stateVector;


    ///Position of the total pressure in the state vector
    CFuint m_pressureID;

    ///Position of the translational temp. in the state vector
    CFuint m_trID;

    ///Position of the vibrational temp. in the state vector
    CFuint m_tvID;

    CFuint m_nbSpecies;

    bool m_convertPartialDensity;

};

}
}

#endif // RADIATIONFIELD_H

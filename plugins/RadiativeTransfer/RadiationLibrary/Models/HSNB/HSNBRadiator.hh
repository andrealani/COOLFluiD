#ifndef COOLFluiD_RadiativeTransfer_HSNBRadiator_hh
#define COOLFluiD_RadiativeTransfer_HSNBRadiator_hh

#include <vector>
#include <string>


#include "MathTools/RealMatrix.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"
#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SnbDiatomicSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SnbAtomicSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SnbContinuumSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/AtomicLines.h"
//#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/HSNBSharedPtr.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ProbDistFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/StringUtils.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/PhotonPath.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SpeciesLoadData.h"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBPhotonTrace.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBLocalParameterSet.hh"

#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBPhotonData.hh"

#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Common/OSystem.hh"
#include "boost/filesystem.hpp"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Common/StringOps.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Common/CFMap.hh"
#include "Common/Stopwatch.hh"
#include "Common/SafePtr.hh"
#include "Environment/FileHandlerInput.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"



//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Environment {
class FileHandlerInput;
}

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////


struct BBCE {
    std::vector<double> TTimesLamda;
    std::vector<double> cumulative;
};

struct TestStruct {
    CFreal* pressure;
    CFreal* temperature;
};

//!
//! \brief Define computation of emission / absorption as introduced
//! in JB's code
//!
class HSNBRadiator : public Radiator
{
public:
    
    HSNBRadiator(const std::string& name);
    ~HSNBRadiator();

    //! Compute the sum of emission over all mechanisms associated
    //! with the current cell
    //! BT: corresponds to sum(mech(i).row())
    virtual CFreal getSpectraLoopPower();

    /// Associates a photon with a emission mechanisms and a subset of the mechanisms emissive power
    void generateCellPhotonData(CFuint stateID, CFreal& energyFraction, CFuint &mechanismIndex, CFint &mechType, CFreal &lambda);

    /// Computes the emission of all mechanisms for the current state / cell and associates a subset
    /// of the total number of rays for this cell over all mechanisms according to a probability distribution
    /// (the higher the emission by a mechanisms the higher the probability that rays will be associated with
    /// that mechanism
    void setState(CFuint nbRays);

    /// Computes the energy absorbed in a cell given a photon and it's trajectory / trace
    /// in which all necessary HSNB parameters have been stored during the raytracing process
    /// First the transmissivity for the path until the curCell is computed. If transmissivity
    /// is lower than a certain tolerance
    CFreal computeAbsorbedEnergy(HSNBDataContainer &traceSet, CFint curCell);

    void addStateParams(HSNBDataContainer& photonTraceSet, CFuint cellID, CFreal raydistance);

    void setAbsorptionTolerance(const CFreal tol);

    void initTrace(HSNBPhotonTrace &trace, CFuint mechID, CFreal startDistance);

    /// Computes transmission as the negative exponential of the optical thickness.
    /// The optical thickness is computed using the HSNBRadiator using an HSNB model
    /// for radiative properties
    CFreal computeTransmission(HSNBDataContainer &traceSet);

    CFreal getAbsorption(CFreal lambda, RealVector &s_o);

    CFreal getAbsorption(HSNBDataContainer &traceSet);

    CFreal tau(HSNBDataContainer &traceSet);

    void absorb(HSNBPhotonTrace& trace, HSNBMechanismType mechType);

    void setup();

    static std::string getClassName() { return "HSNBRadiator"; }
    
    /**
     * Defines the Config Option's of this class
     * @param options a OptionList where to add the Option's
     */
    static void defineConfigOptions(Config::OptionList& options);

    /**
     * Configures this configurable object.
     */
    virtual void configure ( Config::ConfigArgs& args );

    /**
     * Unsetups the data of the library
     */
    virtual void unsetup();

    void setupSpectra(CFreal wavMin, CFreal wavMax);

    CFreal getEmission( CFreal lambda, RealVector &s_o );

    void computeEmissionCPD();

    //lamda generation in generatePhotonData
    void getRandomEmission(CFreal &lambda, RealVector &s_o );

    void getRandomDirection(RealVector &s_o);

    void getRandomBBSigma(CFreal &sig, CFreal temperature);

    void addNewCellAbsorption(CFreal& ku, CFreal& param2, CFreal& param3, CFreal &wavenumber, CFreal &distance, CFuint& mechanismID);

    /**
     * Returns the emissive power as summed up over all cells living in this radiator
     */
    CFreal getTotalEmissivePower() const;

    CFreal getCurrentCellEmissivePower() const;

    CFuint getNbNonThickDiatomics() const;

    CFuint getNbThickDiatomics() const;

    CFuint getNbAtoms() const;

    CFuint getNbContinua() const;
private:

    std::string m_libPath;

    CFreal m_tolerance;

    NonUniformPath m_nupTreatment;

    boost::filesystem::path m_hsnbDir;

    std::string m_localDirName;

    CFreal m_totalPower;

    Common::SafePtr<HSNBPhotonData> m_curPhoton;
    Common::SafePtr<HSNBPhotonTrace> m_curTrace;

    /// iterator for the state vector
    Common::SafePtr<Framework::DofDataHandleIterator<CFreal, Framework::State, Framework::GLOBAL> > m_pstates;

    void setUpRadiationField();

    void setUpMechanisms();

    void load_BBCE_data();

    void determineBandRange();


//    std::vector<>
//    /// thermodynamic library BT: In our case mutation
//    Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
    /// indicates at which pos. of the state vector a species composition (X,rho_i) is specified
    std::vector<CFuint> m_speciesCompID;

    /// Thermodynamic library (convert partial pressures into mole frac. if needed)
    Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

    /// array with molar masses
    RealVector m_mmasses;

    std::vector<CFuint> m_nBRaysPerSys;

    CFreal testTemperature;

    /// array with Avogadro number/molar masses <> 1/m_molec
    /// used for the computation of mole fractions
    RealVector m_avogadroOvMM;

    /// namespace within which HSNBRadiator is run in parallel
    std::string m_namespace;

    /// Define how species composition is defined
    /// "partialDensity", "moleFractions"
    std::string m_compositionType;

    ///Vector containing the states id of the TRS;
    std::vector<CFuint> m_statesID;

    /// Reuse existing radiative data (requires the same number of processors as in the previous run).
    bool m_reuseProperties;

    ///Minimum number density
    CFreal m_ndminFix = 1e+10;

    CFuint m_atomicSelfAbsorption;

    CFuint m_tvI;
    CFuint m_trI;
    CFuint m_pI;
    
    CFuint m_nbCells;
    CFuint m_nbFaces;
    CFuint m_nbSpecies;

    std::size_t m_bmin;
    std::size_t m_bmax;
    std::size_t m_nbands;

    CFuint m_nbMechanisms;

    //Specific to the current state, has to be reset if the state is changed
    CFuint m_curMechanismID;
    CFuint m_curMechanismRayCount;

    ProbDistFunc lineEmission_pdf;
    std::vector<CFreal> m_emissionByMechanisms;

    CFreal m_totalEmissionCurrentCell;
    CFuint m_totalNbRaysCurrentCell;


    boost::filesystem::path m_processFile;
    boost::filesystem::path m_blackBodyFile;
    boost::filesystem::path m_HSNBPath;

    /// input file handle
    Common::SelfRegistPtr<Environment::FileHandlerInput> m_inFileHandle;

    CFuint m_nbRaysPerCell;
    CFuint m_nbRaysPerFace;



    //Keep indices to conveniently access vector elements
    std::vector<CFuint> m_thickDiatomicsIndices;
    std::vector<CFuint> m_thinDiatomicsIndices;

    std::map<CFuint,CFuint> thickDiatomicIndices;
    std::map<CFuint,CFuint> thinDiatomicIndices;

    std::vector<SnbDiatomicSystem> m_diatomics;
    std::vector<SnbContinuumSystem> m_continua;
    std::vector<AtomicLines>  m_atoms;

    CFuint m_nbDiatomics;
    CFuint m_nbContinua;
    CFuint m_nbAtoms;

    CFuint m_nbThickDiatomics;
    CFuint m_nbNonThickDiatomics;

    std::vector<CFreal> m_line_emis;

    /// BlackBody Dataset
    BBCE m_BBCE;

    CFuint m_rank;
    CFuint m_nbProc;

    std::string m_nuPathTreatment;

    ThermoData m_thermoData;

    bool m_convertPartialDensity=false;

    //Keep memory allocated for local variables:
    CFuint m_tempLocalID;
    Common::SafePtr<HSNBThickParameterSet> m_tempThickParamSet;
    CFreal m_tempWavenumber;
    CFreal m_tempAbsorbedEnergy;
    CFreal m_tempTransp;
    CFuint m_tempMechIndex;
    CFuint m_tempNewStateMechIndex;
    CFreal m_tempEnergyFraction;
    CFreal m_rayDistance;

    CFuint m_nbSteps;
    CFreal m_avgStepDistance;

    //Needed for the computation of transmission
    CFreal m_tempShortenedStartDistance;
    CFuint m_tempCurMechIndex;
    CFreal m_tempTau;
    CFreal m_debugTau;
    CFint m_tempEmittingMechanism;
    CFreal m_tempTauk;
    CFreal m_tempds;
    CFreal m_Taukds;

    Common::SafePtr<HSNBNonThickParameterSet> m_tempNonThickParamSet;
    Common::SafePtr<HSNBAtomicParameterSet> m_tempAtomicParamSet;

    //Assumes that m_avogadroOvMM is initialised and that the
    //species defined in the convertVector are the same / in the same order
    //as in the m_library
    void moleFractionsFromPartialDensity(std::vector<CFreal>& convertVector);

    CFreal updateLineEmission(bool printDebug=0);

    bool mechanismIsDiatomic(CFuint mechanismID);
    bool mechanismIsContinous(CFuint mechanismID);

    //!
    //! \brief Load the participating species and emission processes
    //!
    //! Data is stored in the "processes" file in the /data directory
    void loadProcessData();

    //DEBUGGING ONLY
    CFreal emissivePowerDebug(CFuint stateID);
};

//////////////////////////////////////////////////////////////////////////////

}
}

#endif // HSNBRADIATOR_H

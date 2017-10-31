#ifndef COOLFluiD_RadiativeTransfer_HSNBTRACE_hh
#define COOLFluiD_RadiativeTransfer_HSNBTRACE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/PhotonData.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBPhotonData.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBLocalParameterSet.hh"
#include "Framework/PhysicalModel.hh"
#include <vector>
#include "Common/SafePtr.hh"
#include <boost/shared_ptr.hpp>


//using namespace boost;

namespace COOLFluiD {

    namespace RadiativeTransfer {



//////////////////////////////////////////////////////////////////////////////

    const CFuint TraceStaticRealCount=4;
    const CFuint TraceDynamicRealCount=0;
    const CFreal diffFactor=0.001;

///Create a trace of photon absorption parameters for all cells crossed and all mechanisms
struct HSNBPhotonTrace {
public:


    HSNBPhotonTrace();
    ~HSNBPhotonTrace();

    //No need since we got rid of nasty ptrs
    //HSNBPhotonTrace(const HSNBPhotonTrace& obj);

    void setup(CFuint nbContinua, CFuint nbAtoms, CFuint nbSpecies);
    void addTrackingState();

    void addDiatomicThickParams(CFuint mechID, CFreal ku, CFreal betaD, CFreal betaL);
    void addDiatomicNonThickParams(CFuint mechID, CFreal ku);
    void addAtomicParams(CFuint mechID,CFreal optThick);
    void addContinuumParams(CFuint mechID, CFreal ku);

    void addThickDiatomicState(bool isEmittingMechanism);
    void addThinDiatomicState(bool isEmittingMechanism);

    void clear();
    bool isUninitialised() const;

    void print(bool printAll=false);

    virtual std::string name() const {
        return "HSNBPhotonTrace";
    }

    CFuint getBufferRealCount();
    CFuint byteSize();



/////////////////////////////
public:

    std::vector<CFuint> cellID;

    CFuint m_nbSpecies;
    CFuint m_nbDiatomics;
    CFuint m_nbContinua;
    CFuint m_nbAtoms;

    //Static header variables
    CFuint m_mechLocalID;
    CFuint m_nbCellsCrossed;
    CFreal m_kappa0;
    CFreal m_distance0;

    //Use Polymorphism for diatomic systems
    std::vector<HSNBThickParameterSet> m_thickDiatomics;
    std::vector<HSNBNonThickParameterSet> m_thinDiatomics;

    std::vector<HSNBNonThickParameterSet> m_continua;
    std::vector<HSNBAtomicParameterSet> m_atoms;




};

typedef boost::shared_ptr<HSNBPhotonTrace> TracePtr;

//Associate Photon traces with a photon header.
struct HSNBDataContainer {

    void reset() {
        //std::cout << "RESET \n";
        //trace=TracePtr(new HSNBPhotonTrace());
        //CFLog(INFO, "USE COUNT " << trace.use_count() << "\n");
        trace.reset(new HSNBPhotonTrace());
        headerData=NULL;
    }



    void print(bool printAll=false) {
        if (headerData!=NULL) {
            if (printAll) {
                std::cout << "HSNBDataContainer::print => headerData->mechanismIndex= " << headerData->mechanismIndex << " \n";
                std::cout << "HSNBDataContainer::print => headerData->mechType= " << headerData->mechType<< " \n";
                std::cout << "HSNBDataContainer::print => headerData->nbCrossedCells= " << headerData->nbCrossedCells<< " \n";
                std::cout << "HSNBDataContainer::print => headerData->energyFraction= " <<  std::setprecision(16)<< headerData->energyFraction<< " \n";
                std::cout << "HSNBDataContainer::print => headerData->energyResiduum= " <<  std::setprecision(16)<< headerData->energyResiduum<< " \n";
                std::cout << "HSNBDataContainer::print => transm*energy= " <<  std::setprecision(16) << headerData->energyFraction*headerData->transm<< " \n";
                std::cout << "HSNBDataContainer::print => headerData->prevProcessId= " << headerData->prevProcessId<< " \n";
                std::cout << "HSNBDataContainer::print => headerData->fatherProcessId= " << headerData->fatherProcessId<< " \n";
                std::cout << "HSNBDataContainer::print => headerData->wavelength= " <<  std::setprecision(16) << headerData->wavelength<< " \n";
                std::cout << "HSNBDataContainer::print => headerData->wavenumber= " <<  std::setprecision(16) << (1/headerData->wavelength)<< " \n";
            }
            else {
                CFLog(INFO, "HSNBDataContainer::print => headerData->mechanismIndex= " << headerData->mechanismIndex << " \n");
                CFLog(INFO, "HSNBDataContainer::print => headerData->mechType= " << headerData->mechType<< " \n");
                CFLog(INFO, "HSNBDataContainer::print => headerData->nbCrossedCells= " << headerData->nbCrossedCells<< " \n");
                CFLog(INFO, "HSNBDataContainer::print => headerData->energyFraction= " <<  std::setprecision(16)<< headerData->energyFraction<< " \n");
                CFLog(INFO, "HSNBDataContainer::print => headerData->energyResiduum= " <<  std::setprecision(16)<< headerData->energyResiduum<< " \n");
                CFLog(INFO, "HSNBDataContainer::print => transm*energy= " <<  std::setprecision(16) << headerData->energyFraction*headerData->transm<< " \n");
                CFLog(INFO, "HSNBDataContainer::print => headerData->prevProcessId= " << headerData->prevProcessId<< " \n");
                CFLog(INFO, "HSNBDataContainer::print => headerData->wavelength= " <<  std::setprecision(16) << headerData->wavelength<< " \n");
                CFLog(INFO, "HSNBDataContainer::print => headerData->wavenumber= " <<  std::setprecision(16) << (1/headerData->wavelength)<< " \n");
            }
        }
        if (trace.use_count()!=0) {
            trace->print(printAll);
        }
    }

    boost::shared_ptr<HSNBPhotonTrace> trace;
    Common::SafePtr<HSNBPhotonData> headerData;
};

    }
}

#endif // HSNBLOCALPARAMETERSET_H

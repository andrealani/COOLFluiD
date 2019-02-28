#ifndef COOLFluiD_RadiativeTransfer_HSNBLOCALPARAMETERSET_hh
#define COOLFluiD_RadiativeTransfer_HSNBLOCALPARAMETERSET_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "HSNBPhotonData.hh"
#include "Framework/PhysicalModel.hh"
#include <vector>
#include "Common/SafePtr.hh"


namespace COOLFluiD {

    namespace RadiativeTransfer {


//////////////////////////////////////////////////////////////////////////////


    ///
    /// \brief Simple public container to store parameters necessary in order to compute absorption
    ///
struct HSNBLocalParameterSet {

    virtual std::string name() const {
        return "HSNBLocalParameterSet";
    }

    virtual ~HSNBLocalParameterSet() {

    }

    //For debugging
    virtual void print(bool printAll=false) const {
        CFLog(INFO, "HSNBLocalParameter::print=> BaseClass \n");
    }

//    HSNBLocalParameterSet(const HSNBLocalParameterSet& obj) {
//        nbCells=obj.nbCells;
//    }


    HSNBLocalParameterSet(){
        nbCells=0;
    }

    void addState() {
        nbCells++;
    }

    bool isEmpty() const {
        return (nbCells==0);
    }

    CFuint backIndex() const {
        return (nbCells-1);
    }

    virtual CFuint getRealCount() const {
        return 1;
    }


    CFuint nbCells=0;
};

struct HSNBAtomicParameterSet: HSNBLocalParameterSet {

    virtual std::string name() const {
        return "HSNBAtomicParameterSet";
    }


    ~HSNBAtomicParameterSet() {

    }


    //For debugging
    virtual void print(bool printAll=false) const {
        if (printAll) {
            std::cout <<"HSNBAtomicParameterSet::print=> Parameters: \n";
            std::cout <<"HSNBAtomicParameterSet::print=> optThick: " << optThick << " \n";
            std::cout <<"\n";
        }
        else {
            CFLog(INFO, "HSNBAtomicParameterSet::print=> Parameters: \n");
            CFLog(INFO, "HSNBAtomicParameterSet::print=> optThick: " << optThick << " \n");
            CFLog(INFO, "\n");
        }
    }

    HSNBAtomicParameterSet(): HSNBLocalParameterSet(){
        optThick=0.0;
    }

    void addState(CFreal newOptThick) {
        optThick=newOptThick;
        HSNBLocalParameterSet::addState();
    }

    virtual CFuint getRealCount() const {
        return 1;
    }

    CFreal optThick;

};

struct HSNBNonThickParameterSet: HSNBLocalParameterSet {

    virtual std::string name() const {
        return "HSNBNonThickParameterSet";
    }



    ~HSNBNonThickParameterSet() {

    }

    virtual void print(bool printAll=false) const {
        if (printAll) {
            std::cout <<"HSNBNonThickParameterSet::print=> Parameters: kappa=" << kappa << " \n";
            std::cout <<"\n";
        }
        else {
            CFLog(INFO, "HSNBNonThickParameterSet::print=> Parameters: kappa=" << kappa << " \n");
            CFLog(INFO, "\n");
        }
    }

//    HSNBNonThickParameterSet(const HSNBNonThickParameterSet& obj): HSNBLocalParameterSet(obj) {
//        kappa=obj.kappa;
//    }

    HSNBNonThickParameterSet(): HSNBLocalParameterSet(){
        kappa=0.0;
    }

    void addState(CFreal newKappa) {
        kappa=newKappa;
        HSNBLocalParameterSet::addState();
    }

    virtual CFuint getRealCount() const {
        return 1;
    }

    //FOR NONTHICK SYSTEMS WE ONLY NEED THE PREVIOUS VALUE, NOT THE ENTIRE ARRAY
    CFreal kappa=0.0;
};

struct HSNBThickParameterSet: HSNBLocalParameterSet {


    ~HSNBThickParameterSet() {

    }

//    HSNBThickParameterSet(const HSNBThickParameterSet& obj): HSNBLocalParameterSet(obj) {
//        kappa=obj.kappa;
//        betaD=obj.betaD;
//        betaL=obj.betaL;
//    }


    HSNBThickParameterSet():HSNBLocalParameterSet()
    {
        startDistance=0.0;
    }

    virtual void print(bool printAll=false) const {
        if (printAll) {
            std::cout <<"HSNBThickParameterSet::print=> Parameters: \n";
            for (int i=0; i<nbCells; i++) {
                std::cout <<"HSNBThickParameterSet::print=> Cell[" << i << "],kappa="
                      << kappa[i] <<", betaD=" << betaD[i] << ", betaL=" << betaL[i] <<"\n";
            }
            std::cout <<"\n";
        }
        else {
            CFLog(INFO, "HSNBThickParameterSet::print=> Parameters: \n");
            for (int i=0; i<nbCells; i++) {
                CFLog(INFO, "HSNBThickParameterSet::print=> Cell[" << i << "],kappa="
                      << kappa[i] <<", betaD=" << betaD[i] << ", betaL=" << betaL[i] <<"\n");
            }
            CFLog(INFO, "\n");
        }
    }

    virtual std::string name() const {
        return "HSNBThickParameterSet";
    }


    void addState(CFreal newKappa, CFreal newBetaD,
                  CFreal newBetaL) {

        kappa.push_back(newKappa);
        betaD.push_back(newBetaD);
        betaL.push_back(newBetaL);

        HSNBLocalParameterSet::addState();
    }

    virtual CFuint getRealCount() const {
        return nbCells*3;
    }

    CFreal startDistance;
    std::vector<CFreal> kappa;
    //NEEDED ONLY FOR OPTICALLY THICK SYSTEMS
    std::vector<CFreal> betaD;
    std::vector<CFreal> betaL;
};

struct HSNBCO2ParameterSet: HSNBLocalParameterSet {


    ~HSNBCO2ParameterSet() {

    }


    void reset() {
        kappa=0.0;
        nbCells=0;
    }


    HSNBCO2ParameterSet():HSNBLocalParameterSet()
    {
        kappa=0.0;
    }

    virtual void print(bool printAll=false) const {
        if (printAll) {
            std::cout <<"HSNBCO2ParameterSet::print=> Parameters: kappa=" << kappa << " \n";
            std::cout <<"\n";
        }
        else {
            CFLog(INFO, "HSNBCO2ParameterSet::print=> Parameters: kappa=" << kappa << " \n");
            CFLog(INFO, "\n");
        }
    }

    virtual std::string name() const {
        return "HSNBCO2ParameterSet";
    }


    void addState(CFreal newKappa) {
        kappa=newKappa;
        HSNBLocalParameterSet::addState();
    }

    virtual CFuint getRealCount() const {

        return 1;
    }

    //Store previous value (sum over the optical depth along previous trajectory)
    CFreal kappa;

};


    }
}

#endif // HSNBLOCALPARAMETERSET_H

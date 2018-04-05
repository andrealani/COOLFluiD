#include "HSNBPhotonTrace.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {



HSNBPhotonTrace::HSNBPhotonTrace()
{
m_nbAtoms=0;
m_nbDiatomics=0;
m_nbCellsCrossed=0;
m_nbContinua=0;
m_nbSpecies=0;
m_kappa0=0;
m_mechLocalID=0;

}

HSNBPhotonTrace::~HSNBPhotonTrace()
{

    //std::cout << "Trace Destructor called... \n ";

//    for (int i=0; i<m_diatomics.size(); i++) {
//        if (m_diatomics[i]!=NULL) {
//            delete m_diatomics[i];
//        }
//    }

    //std::cout <<" ....and finished \n";
}

void HSNBPhotonTrace::setup(CFuint nbContinua, CFuint nbAtoms, CFuint nbSpecies)
{
        m_nbCellsCrossed=0;
        m_nbAtoms=nbAtoms;

        m_nbContinua=nbContinua;
        m_nbSpecies=nbSpecies;
        m_nbDiatomics=0;

        m_thickDiatomics.clear();
        m_thinDiatomics.clear();
        m_continua.clear();
        m_atoms.clear();

        m_continua.resize(m_nbContinua);
        m_atoms.resize(m_nbAtoms);
        //For diatomics we have to distinguish between thick and thin diatomics
        //Thus diatomics are added seperately using addThickDiatomicState / addThinDiatomicState
}

void HSNBPhotonTrace::addTrackingState()
{
    m_nbCellsCrossed++;
}

void HSNBPhotonTrace::addDiatomicThickParams(CFuint mechID, CFreal ku, CFreal betaD, CFreal betaL)
{

    m_thickDiatomics[mechID].addState(ku,betaD,betaL);
}

void HSNBPhotonTrace::addDiatomicNonThickParams(CFuint mechID, CFreal ku)
{
    m_thinDiatomics[mechID].addState(ku);
}

void HSNBPhotonTrace::addAtomicParams(CFuint mechID, CFreal optThick)
{
    m_atoms[mechID].addState(optThick);
}

void HSNBPhotonTrace::addContinuumParams(CFuint mechID, CFreal ku)
{

    m_continua[mechID].addState(ku);
}



void HSNBPhotonTrace::addThickDiatomicState(bool isEmittingMechanism)
{

    //If the current mechanism is the emitting diatomic system we keep the
    //ID to be able to compute the correction for thick diatomics in the computation
    //of transmissivity
    if (isEmittingMechanism) {
        m_mechLocalID=m_thickDiatomics.size();
    }
    m_nbDiatomics++;
    m_thickDiatomics.push_back(HSNBThickParameterSet());

}


void HSNBPhotonTrace::addThinDiatomicState(bool isEmittingMechanism)
{

    //If the current mechanism is the emitting diatomic system we keep the
    //ID to be able to compute the correction for thick diatomics in the computation
    //of transmissivity
    if (isEmittingMechanism) {
        m_mechLocalID=m_thinDiatomics.size();
    }

    m_nbDiatomics++;


    m_thinDiatomics.push_back(HSNBNonThickParameterSet());
    //std::cout << "create thin state : m_thinDiatomics.size()=" <<  m_thinDiatomics.size()
    //          << ", trying to access..." ;
    //std::cout << m_thinDiatomics.back().name()  << " \n";
}




void HSNBPhotonTrace::clear()
{

    cellID.clear();

    m_nbDiatomics=0;
    m_nbCellsCrossed=0;

    m_thinDiatomics.clear();
    m_thickDiatomics.clear();
    m_continua.clear();
    m_atoms.clear();

    m_continua.resize(m_nbContinua);
    m_atoms.resize(m_nbAtoms);

}

bool HSNBPhotonTrace::isUninitialised() const
{
    return ((m_thinDiatomics.size()+m_thickDiatomics.size())==0);
}

void HSNBPhotonTrace::print(bool printAll)
{
    if (printAll) {
        std::cout << "HSNBPhotonTrace::print()=> Print Trace: \n";
        std::cout << "HSNBPhotonTrace::print()=> Total Bytesize:"<< byteSize() << " // " << getBufferRealCount()<<  " REAL values \n";
        std::cout << "HSNBPhotonTrace::print()=> nbSpecies" << m_nbSpecies << " \n";
        std::cout << "HSNBPhotonTrace::print()=> nbDiatomics" << m_nbDiatomics << " \n";
        std::cout << "HSNBPhotonTrace::print()=> nbContinua" << m_nbContinua << " \n";
        std::cout << "HSNBPhotonTrace::print()=> nbAtoms" << m_nbAtoms << " \n";
        std::cout << "HSNBPhotonTrace::print()=> starting Distance" << m_distance0 << " \n";
        std::cout << "HSNBPhotonTrace::print()=> nbCellsCrossed" << m_nbCellsCrossed << " \n";

        std::cout << "HSNBPhotonTrace::print()=> cellIDs: \n";
        for (int i=0; i<cellID.size(); i++) {
            std::cout << "HSNBPhotonTrace::print()=>cellID["<<i<< "]: "<< cellID[i] << " \n";
        }


        std::cout << "HSNBPhotonTrace::print()=> Thin Diatomics: \n";
        for (int i=0; i<m_thinDiatomics.size(); i++) {
            m_thinDiatomics[i].print(true);
        }

        std::cout << "HSNBPhotonTrace::print()=> Thick Diatomics: \n";
        for (int i=0; i<m_thickDiatomics.size(); i++) {
            m_thickDiatomics[i].print(true);
        }

        std::cout << "HSNBPhotonTrace::print()=> Continua: \n";
        for (int i=0; i<m_continua.size(); i++) {
            m_continua[i].print(true);
        }

        std::cout << "HSNBPhotonTrace::print()=> Atoms: \n";
        for (int i=0; i<m_atoms.size(); i++) {
            m_atoms[i].print(true);
        }
    }
    else {

        CFLog(INFO, "HSNBPhotonTrace::print()=> Print Trace: \n");
        CFLog(INFO, "HSNBPhotonTrace::print()=> Total Bytesize:"<< byteSize() << " // " << getBufferRealCount()<<  " REAL values \n");
        CFLog(INFO, "HSNBPhotonTrace::print()=> nbSpecies" << m_nbSpecies << " \n");
        CFLog(INFO, "HSNBPhotonTrace::print()=> nbDiatomics" << m_nbDiatomics << " \n");
        CFLog(INFO, "HSNBPhotonTrace::print()=> nbContinua" << m_nbContinua << " \n");
        CFLog(INFO, "HSNBPhotonTrace::print()=> nbAtoms" << m_nbAtoms << " \n");
        CFLog(INFO, "HSNBPhotonTrace::print()=> nbCellsCrossed" << m_nbCellsCrossed << " \n");

        CFLog(INFO, "HSNBPhotonTrace::print()=> cellIDs: \n");
        for (int i=0; i<cellID.size(); i++) {
            CFLog(INFO, "HSNBPhotonTrace::print()=>cellID["<<i<< "]: "<< cellID[i] << " \n");
        }


        CFLog(INFO, "HSNBPhotonTrace::print()=> Thin Diatomics: \n");
        for (int i=0; i<m_thinDiatomics.size(); i++) {
            m_thinDiatomics[i].print();
        }

        CFLog(INFO, "HSNBPhotonTrace::print()=> Thick Diatomics: \n");
        for (int i=0; i<m_thickDiatomics.size(); i++) {
            m_thickDiatomics[i].print();
        }

        CFLog(INFO, "HSNBPhotonTrace::print()=> Continua: \n");
        for (int i=0; i<m_continua.size(); i++) {
            m_continua[i].print();
        }

        CFLog(INFO, "HSNBPhotonTrace::print()=> Atoms: \n");
        for (int i=0; i<m_atoms.size(); i++) {
            m_atoms[i].print();
        }

    }

    cf_assert(m_atoms.size()==m_nbAtoms);
    cf_assert(m_continua.size()==m_nbContinua);


    cf_assert((m_thinDiatomics.size()+m_thickDiatomics.size())==m_nbDiatomics);
}


CFuint HSNBPhotonTrace::getBufferRealCount()
{
    CFuint count=0;

    //Real Variables in diatomic traces
    if (m_thinDiatomics.empty()==false) {
        count+=m_thinDiatomics[0].getRealCount()*m_thinDiatomics.size();
    }

    //Real Variables in diatomic traces
    if (m_thickDiatomics.empty()==false) {
        count+=m_thickDiatomics[0].getRealCount()*m_thickDiatomics.size();
    }

    //Real variables in continous traces
    if (m_continua.empty()==false) {
        count+=m_continua[0].getRealCount()*m_continua.size();
    }

    //Real Variables in atomic traces
    if (m_atoms.empty()==false) {
        count+=m_atoms[0].getRealCount()*m_atoms.size();
    }

    //Static header variables
    count +=TraceStaticRealCount;

    return count;
}

CFuint HSNBPhotonTrace::byteSize()
{
    return getBufferRealCount()*sizeof(CFreal);
}






    }
}

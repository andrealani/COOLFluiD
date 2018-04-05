#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////



typedef std::numeric_limits< CFreal > realVal;

ThermoData::ThermoData(): m_p(0), m_tr(0), m_tv(0), m_currentStateID(0) {

}


void ThermoData::setup(const CFuint pressureID, const CFuint trID, const CFuint tvID, RealVector* avogadroOvMM, bool convertPartialPressure)
{
    m_pressureID=pressureID;
    m_trID=trID;
    m_tvID=tvID;
    m_convertPartialPressure=convertPartialPressure;
    m_avogadroOvMM=avogadroOvMM;
}

CFuint ThermoData::getLocalCellID(CFuint cellID)
{
    stateMapIt=stateIDMap.find(cellID);

    try {
        if (isValidLocalID(stateMapIt)) {
            return stateMapIt->second;
        }
        else {
            throw COOLFluiD::Common::BadValueException (FromHere(),"ThermoData::setState(): The stateID provided does not refer to any state held by this instance of ThermoData");
        }
    }
    catch (Common::Exception& e) {
      CFout << e.what() << "\n";
      throw;
    }
}

double ThermoData::N(int i) const
{
    if (m_externalState.active==false) {
        return (i >= 0 ? m_states[m_currentStateID].N(i) : 0.0);
    }
    else {
        return (i >= 0 ? m_externalState.m_N[i] : 0.0);
    }

}

CFreal* ThermoData::N()
{
    if (m_externalState.active==false) {
        return m_states[m_currentStateID].getNumberDensities();
    }
    else {
        return m_externalState.m_N;
    }
}

CFuint ThermoData::nCells() const
{
    return m_states.size();
}

void ThermoData::setState(CFuint stateIndex)
{
    m_externalState.active=false;
    stateMapIt=stateIDMap.find(stateIndex);

    //cf_assert(stateMapIt->second==stateIndex);
    //CFLog(VERBOSE, "ThermoData::setState => stateIndex=" << stateIndex << "\n");

    try {
        if (isValidLocalID(stateMapIt)) {
            m_tr  = m_states[stateMapIt->second].Tr();
            m_tv  = m_states[stateMapIt->second].Tv();
            m_p   = m_states[stateMapIt->second].p();

            m_currentStateID=stateMapIt->second;
            //CFLog(VERBOSE, "ThermoData::setState => m_currentStateID=" << m_currentStateID << "\n");
        }
        else {
            throw COOLFluiD::Common::BadValueException (FromHere(),"ThermoData::setState(): The stateID provided does not refer to any state held by this instance of ThermoData");
        }
    }
    catch (Common::Exception& e) {
      CFout << e.what() << "\n";
      throw;
    }

}

void ThermoData::setExternalState(CFreal newTr, CFreal newTv, CFreal newP, CFreal *newN)
{
    m_tr=newTr;
    m_tv=newTv;
    m_p=newP;
    m_externalState.activateState(newN);
}

void ThermoData::deactivateExternalState()
{
   m_externalState.active=false;
}

void ThermoData::addState(RealVector *stateVector, CFuint localStateID)
{
    stateIDMap.insert(std::make_pair(localStateID, m_states.size()));
    m_states.push_back(RadiationFieldData(m_pressureID,m_trID,m_tvID,emIndex,stateVector,m_nbSpecies, m_avogadroOvMM));
}

void ThermoData::printState(CFuint stateID)
{
   setState(stateID);
   CFLog(INFO, "ThermoData::printState => Printing state " << stateID << " with local ID " << m_currentStateID << " \n");
   CFLog(INFO, "ThermoData::printState => Tr=" <<  std::setprecision(16) << m_tr << "; \n");
   CFLog(INFO, "ThermoData::printState => Tv=" <<  std::setprecision(16) << m_tv << "; \n");
   CFLog(INFO, "ThermoData::printState => P=" << std::setprecision(16) <<  m_p << "; \n\n");

   for (int i=0; i<nSpecies();i++) {
       CFLog(INFO, "ThermoData::printState => Species " << i << ": " <<  speciesName(i) << " \n");
   }

   CFLog(INFO, " \n");

   for (int i=0; i<nSpecies();i++) {
       CFLog(INFO, "ThermoData::printState => N[" << i << "]=" << std::setprecision(16) << N(i) << "; \n");
   }


}

bool ThermoData::isValidLocalID(std::map<CFuint, CFuint>::iterator stateMapIt) const
{
    return (stateMapIt!= stateIDMap.end());
}

}
}

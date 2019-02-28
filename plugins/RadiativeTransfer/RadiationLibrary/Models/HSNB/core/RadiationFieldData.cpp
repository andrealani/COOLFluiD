#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/RadiationFieldData.h"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////


RadiationFieldData::RadiationFieldData()
{

}

RadiationFieldData::RadiationFieldData(const CFuint pressureID, const CFuint trID, const CFuint tvID, const CFuint emIndex, const Common::SafePtr<RealVector> stateVector, const CFuint& nbSpecies, const Common::SafePtr<RealVector> avogadroOvMM): m_nbSpecies(nbSpecies), m_pressureID(pressureID),
    m_trID(trID), m_tvID(tvID)
{
    m_numberDensity.resize(nbSpecies);
//    m_numberDensity = new RealVector(nbSpecies);
    m_stateVector = stateVector;

    CFreal p,Tr,Tv;

//    std::cout << "RadiationFieldData::RadiationFieldData => m_stateVector->size=" << m_stateVector->size() << std::endl;
//    std::cout << "RadiationFieldData::RadiationFieldData => pressureID=" << pressureID << " trID="  <<trID << " tvID=" << tvID <<std::endl;



    p= (*m_stateVector)[pressureID];
    Tr=(*m_stateVector)[trID];
    Tv=(*m_stateVector)[tvID];


    CFLog(DEBUG_MAX, "RadiationFieldData::RadiationFieldData => New state. P: " <<p << ", Tr " << Tr << ", Tv " << Tv << "\n");
//    std::cout<< "RadiationFieldData::RadiationFieldData => New state. P: " <<p << ", Tr " << Tr << ", Tv " << Tv << std::endl;

    for (CFuint si=0; si<m_nbSpecies; si++) {
        //Need to convert all part. densities to nD first
//       std::cout<< "RadiationFieldData::RadiationFieldData => m_stateVector).at(" << si << "):"<<si <<" nbSpecies=" << nbSpecies << std::endl;

       m_numberDensity[si]=((*m_stateVector).at(si)*(*avogadroOvMM).at(si));

//       if (m_numberDensity[si]<0.0) {
//           CFLog(INFO, "RadiationFieldData::RadiationFieldData => m_stateVector[" << si << "]=" << (*m_stateVector).at(si) << "\n");
//           CFLog(INFO, "RadiationFieldData::RadiationFieldData => avogadroOvMM[" << si << "]=" << (*avogadroOvMM).at(si) << "\n");
//           CFLog(INFO, "RadiationFieldData::RadiationFieldData => m_numberDensity[" << si << "]=" << m_numberDensity[si]<< "\n");
//       }

//       std::cout << "NUMBER DENSITY m_numberDensity[" << si << "]=" << m_numberDensity[si] << " Avogadro=" << (*avogadroOvMM).at(si) << " STATE=" <<(*m_stateVector).at(si) << std::endl;
//       cf_assert(m_numberDensity[si]>=0.0);

    }


//    std::cout<< "\n";




//    for (int i=0; i<m_nbSpecies; i++) {
//       std::cout << "AvogadroOvMM("<< i << ")=" << (*avogadroOvMM).at(i) <<"\n";
//    }
//     std::cout<< "\n";

//    if (convertPartialDensity) {
//         //Convert to molar Fraction first
//        X_e=(*m_stateVector).at(emIndex)*(*avogadroOvMM).at(emIndex);
//        n = p / (KB * (Tr + X_e*(Tv - Tr)));

//        for (CFuint si=0; si<m_nbSpecies; si++) {
//            //Need to convert all part. densities to MF first
//            moleFraction_i=((*m_stateVector).at(si)*(*avogadroOvMM).at(si));
////            std::cout << "moleFraction_i " << (*m_stateVector).at(si) << " // "<< moleFraction_i << std::endl;
//             m_numberDensity[si]=moleFraction_i*n;
//        }

//    }
//    else {
//        X_e=(*m_stateVector).at(emIndex);
//        n = p / (KB * (Tr + X_e*(Tv - Tr)));

//        for (CFuint si=0; si<m_nbSpecies; si++) {
//            //No need to convert
//            moleFraction_i=((*m_stateVector).at(si));
//            m_numberDensity[si]=moleFraction_i*n;
//        }
//    }

}

CFreal *RadiationFieldData::X()
{
    return NULL;
}

CFreal RadiationFieldData::X(CFuint i)
{
    return 1.0;
}

CFreal RadiationFieldData::p()
{
    return (*m_stateVector).at(m_pressureID);
}

CFreal RadiationFieldData::Tr()
{
    return (*m_stateVector).at(m_trID);
}

CFreal RadiationFieldData::Tv()
{
    return (*m_stateVector).at(m_tvID);
}

CFreal* RadiationFieldData::getNumberDensities()
{
    return m_numberDensity.ptr();
}


RadiationFieldData::~RadiationFieldData()
{

}

CFreal RadiationFieldData::N(CFuint i) const
{
    return m_numberDensity[i];
}
    
}
}

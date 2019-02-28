#ifndef COOLFluiD_RadiativeTransfer_SPECIESLOADDATA_H
#define COOLFluiD_RadiativeTransfer_SPECIESLOADDATA_H

#include <string>
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/StringUtils.h"
#include "Environment/CFEnv.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

struct SpeciesLoadData {
    SpeciesLoadData(std::string basePath, std::string loadDataset) {

        std::vector<std::string> tokens;
        String::tokenize(loadDataset,tokens,"/");
        //TODO ASSERT tokens.size==2
        this->speciesName=tokens[0];
        this->systemName=tokens[1];
        this->baseDirectory=basePath;

    }

    std::string baseDirectory;
    std::string speciesName;
    std::string systemName;


    void print(){
        std::cout << "Loaded species data "<< baseDirectory <<"/" << speciesName << systemName << " \n" << std::endl;
    }
};

}
}
#endif // SPECIESLOADDATA_H

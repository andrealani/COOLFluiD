#ifndef COOLFluiD_RadiativeTransfer_SPECIESDATA_H
#define COOLFluiD_RadiativeTransfer_SPECIESDATA_H

#include <string>
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

struct SpeciesData {
    SpeciesData(const std::string& n, double c, double m)
        : name(n), charge(c), mw(m)
    { }

    std::string name;
    double      charge;
    double      mw;
};

}
}

#endif // SPECIESDATA_H

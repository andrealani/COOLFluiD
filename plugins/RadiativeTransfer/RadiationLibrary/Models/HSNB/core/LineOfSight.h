#ifndef COOLFluiD_RadiativeTransfer_LINE_OF_SIGHT_H
#define COOLFluiD_RadiativeTransfer_LINE_OF_SIGHT_H

#include <vector>
#include <cassert>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/FieldData.h"

//enum Geometry
//{
//    CARTESIAN,
//    SPHERICAL
//};


/**
 * @class LineOfSight 
 * @brief 
 */
class LineOfSight
{
public:
    
    LineOfSight(const FieldData& field,
                const int& index, const double& mu);
    
    LineOfSight(const FieldData& field, int i1, int i2);

    LineOfSight() { }

    ~LineOfSight() {}

    /**
      * Returns the number of cells crossed by the path
      */
    int nPath() const {
        return m_path.size();
    }


    /**
      * Returns the length associated to the path index
      */
    double length(int index) const {
        return m_path[index].first;
    }


    /**
      * Returns the cell number associated to the path index
      */
    int cellNum(int index) const {
        return m_path[index].second;
    }

    /**
      * Returns the wall index associated to the path
      */
    int wallNum() const {
        return m_wall;
    }

    // Adds a new cell to this distance
    void addCell(int i, double d) {
        assert(d > 0);
        m_path.push_back(std::make_pair(d,i));
    }


private:

    void buildCartesianPath(const int& index, const double& mu);
    void buildSphericalPath(const int& index, const double& mu);
    
private:
    
    std::vector<std::pair<double, int> > m_path;
    int m_wall;
    
    FieldData m_field;
};

#endif // LINE_OF_SIGHT_H

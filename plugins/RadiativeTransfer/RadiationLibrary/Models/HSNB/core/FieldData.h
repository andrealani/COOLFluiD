#ifndef COOLFluiD_RadiativeTransfer_FIELD_DATA_H
#define COOLFluiD_RadiativeTransfer_FIELD_DATA_H

#include <vector>
#include <iostream>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Constants.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

enum Geometry
{
  CARTESIAN,
  SPHERICAL
};

/**
 * @class CellData
 * @brief Stores information regarding a cell of the computational domain.
 */
class CellData
{
public:
    
    /**
     * Creates a new CellData from the pressure, rotational and vibrational
     * temperatures and mixture composition in species mole fractions.
     */
    CellData(double P, double Tr, double Tv, const double composition [],
	     COOLFluiD::RadiativeTransfer::ThermoData& thermo, int nbSpecies=0)
        : m_p(P), m_tr(Tr), m_tv(Tv), 
          m_mf(nbSpecies),
          m_mn(nbSpecies)
    {
        double sum = 0.0;
        for (int i = 0; i < nbSpecies; ++i) {
            m_mf[i] = std::max(1.0e-20, composition[i]);
            sum += m_mf[i];
        }
        for (int i = 0; i < nbSpecies; ++i)
            m_mf[i] /= sum;


        int ie = thermo.speciesIndex("e-");
        double n = P / (KB * (Tr + X(ie)*(Tv - Tr)));

        for (int i = 0; i < nbSpecies; ++i)
            m_mn[i] = m_mf[i]*n;
    }
    

    /**
     * Returns the pressure at this cell.
     */
    double P() const { return m_p; }
    
    /**
     * Returns the rotational temperature at this cell.
     */
    double Tr() const { return m_tr; }
    
    /**
     * Returns the vibrational temperature at this cell.
     */
    double Tv() const { return m_tv; }
    
    /**
     * Returns the pointer to the species mole fraction array at this cell.
     */
    const double* const X() const { return m_mf.data(); }
    double X(int i) const { return (i >= 0 ? m_mf[i] : 0.0); }

    /**
     * Returns the pointer to the species number density array at this cell.
     */
    const double* const N() const { return m_mn.data(); }
    double N(int i) const { return (i >= 0 ? m_mn[i] : 0.0); }

private:
    
    double  m_p;
    double  m_tr;
    double  m_tv;
    std::vector<double> m_mf;
    std::vector<double> m_mn;
};

/**
 * @class FieldData
 * @brief Stores information regarding mixture state at all the cells of the computational domain.
 */
class FieldData
{
public:
    
    /**
     * Empty constructor.
     */
    FieldData() {}
    
    /**
     * Loads the field data from a file.
     */
    FieldData(const std::string& file);
    
    /**
     * Destructor.
     */
    ~FieldData() {}
    
    /**
     * Returns the number of points.
     */
    int nPoints() const { return m_x.size(); }

    /**
     * Returns the number of cells.
     */
    int nCells() const { return m_cells.size(); }

    Geometry geom() const { return m_geom; }

    void setGeometry(const std::string& geom) {
        if (geom == "CARTESIAN") { m_geom = CARTESIAN; }
        else if (geom == "SPHERICAL") { m_geom = SPHERICAL; }
        else { std::cout << "ERROR GEOMETRY INITIALIZATION" << std::endl; }
    }
    
    /**
     * Adds a new point to the domain.
     */
    void addPoint(const double& loc) {
        m_x.push_back(loc);
    }

    /**
     * Adds a new cell to the domain.
     */
    void addCell(const CellData& loccell) {
        m_cells.push_back(loccell);
    }

    /**
     * Adds a boundary condition.
     */
    void addBoundaryCondition(double epsilon, double twall) {
        m_bc.push_back(std::make_pair(epsilon, twall));
    }

    /**
     * Returns a reference to the CellData at the given index.
     */
    const CellData& cell(int index) const {
        return m_cells[index];
    }
    
    /**
     * Returns the location of the point.
     */
    double loc(int index) const {
        return m_x[index];
    }

    /**
      * Returns the emissivity of the wall point.
      */
    double epsilon(int index) const {
        return m_bc[index].first;
    }

    /**
      * Returns the temperature of the wall point.
      */
    double twall(int index) const {
        return m_bc[index].second;
    }
    
    friend std::ostream& operator << (
        std::ostream& os, const FieldData& field);
    friend std::istream& operator >> (
        std::istream& is, FieldData& field);

    
private:

    Geometry m_geom;    
    std::vector<std::pair< double, double> > m_bc;
    std::vector<CellData> m_cells;
    std::vector<double> m_x;
};

/**
 * Output operator for a FieldData object. Prints a small header and the data
 * for each cell in the domain.
 */
std::ostream& operator << (std::ostream& os, const FieldData& field);

/**
 * Input operator for a FieldData object.  Reads in the data in from columns.
 */
std::istream& operator >> (std::istream& is, FieldData& field);

#endif // FIELD_DATA_H

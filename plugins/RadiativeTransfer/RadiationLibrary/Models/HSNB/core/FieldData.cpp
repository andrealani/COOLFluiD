#include <cstdlib>
#include <iomanip>
#include <fstream>

#include "FieldData.h"
#include "StringUtils.h"

FieldData::FieldData(const std::string& file)
{
    std::ifstream is(file.c_str());
    is >> (*this);
    is.close();
}

std::ostream& operator << (std::ostream& os, const FieldData& field)
{
    using std::setw;
    const GlobalThermoData& thermo = GlobalThermoData::getInstance();

    // Print the number of points and cells
    os << "Number of points - number of cells \n";
    os << setw(10) << field.m_x.size() 
       << setw(10) << field.m_cells.size() << "\n";
    
    // Print the number of points and cells
    os << "Boundary conditions, eps1, Tw1(K), eps2, Tw2(K) \n";
    os << setw(10) << field.m_bc[0].first 
       << setw(10) << field.m_bc[0].second
       << setw(10) << field.m_bc[1].first 
       << setw(10) << field.m_bc[1].second << "\n";

    // Print the header
    os << setw(7)  << "Loc (cm)";
    os << setw(10) << "P (Pa)";
    os << setw(10) << "Tr (K)";
    os << setw(10) << "Tv (K)";
    
    for (int i = 0; i < thermo.nSpecies(); ++i)
        os << setw(12) << thermo.speciesName(i);
    
    os << "\n";
    
    // Print each cell and the last point
    for (int i = 0; i < field.m_cells.size(); ++i) {
        os << setw(7)  << field.m_x[i];
        os << setw(10) << field.m_cells[i].P();
        os << setw(10) << field.m_cells[i].Tr();
        os << setw(10) << field.m_cells[i].Tv();
        
        for (int j = 0; j < thermo.nSpecies(); ++j)
            os << setw(12) << field.m_cells[i].X()[j];
        
        os << "\n";
    }
    os << setw(7)  << field.m_x[field.m_x.size()-1];
        
    return os;
}

std::istream& operator >> (std::istream& is, FieldData& field)
{
    std::string line;
    std::vector<std::string> tokens;

    // Read in the number of points
    getline(is, line); 
    getline(is, line); 
    String::tokenize(String::trim(line, " \t"), tokens, " \t");
    const int np = atoi(tokens[0].c_str());

    // Read in the 1D geometry type
    getline(is, line); 
    getline(is, line); 
    String::tokenize(String::trim(line, " \t"), tokens, " \t");
    field.setGeometry(tokens[0]);

    // Read in the boundary conditions
    double epsilon, twall;
    getline(is, line);
    getline(is, line);
    String::tokenize(String::trim(line, " \t"), tokens, " \t");
    epsilon = atof(tokens[0].c_str());
    twall = atof(tokens[1].c_str());
    field.addBoundaryCondition(epsilon, twall);
    epsilon = atof(tokens[2].c_str());
    twall = atof(tokens[3].c_str());
    field.addBoundaryCondition(epsilon, twall);

    // Read in the species data from the table header
    const GlobalThermoData& thermo = GlobalThermoData::getInstance();
    getline(is, line); 
    String::tokenize(String::trim(line, " \t"), tokens, " \t");
    
    for (int i = 4; i < tokens.size(); ++i)
        GlobalThermoData::getInstance().addSpecies(tokens[i]);

    // Read in each of the cells
    const int ns = thermo.nSpecies();
    double loc, tr, tv, p, x[ns];
    
    getline(is, line);
    for (int j=0; j < np-1; j++) {
        tokens.clear();
        String::tokenize(String::trim(line, " \t"), tokens, " \t");
        
        loc = atof(tokens[0].c_str());
        tr  = atof(tokens[1].c_str());
        tv  = atof(tokens[2].c_str());
        p   = atof(tokens[3].c_str());
        
        for (int i = 0; i < thermo.nSpecies(); ++i) {
            x[i] = atof(tokens[i+4].c_str());
        }
        
        field.addPoint(loc);
        field.addCell(CellData(p, tr, tv, x));
        getline(is, line);
    }
    tokens.clear();
    String::tokenize(String::trim(line, " \t"), tokens, " \t");
    loc = atof(tokens[0].c_str());
    field.addPoint(loc);
    
    return is;
}



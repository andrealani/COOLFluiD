#ifndef COOLFluiD_RadiativeTransfer_LORENTZ_TOPBASE_H
#define COOLFluiD_RadiativeTransfer_LORENTZ_TOPBASE_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LorentzModel.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

using namespace COOLFluiD::RadiativeTransfer;

/**
 * Implements the Lorentz model using TOPBASE data.
 */
class LorentzTopbase : public LorentzModel
{
public:
    /// Constructor
    LorentzTopbase(
        const std::string& species, const std::vector<int>& topbase_index) :
        m_topbase_id(topbase_index)
    {
        // Compute the maximum id
        m_ntopbase = 0;
        for (int i = 0; i < m_topbase_id.size(); ++i)
            m_ntopbase = std::max(m_ntopbase, m_topbase_id[i]);
        m_ntopbase++;

        // Open the Lorentz line-width database for this atom
        std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
        std::string database = datadir + "/lbl_data/widths/atoms/" + species + ".lor";
        std::ifstream file(database.c_str());

        if (!file.is_open()) {
            std::cout << "Could not open database file " << database << "!" << std::endl;
            std::exit(1);
        }

        // First read the type of lorentz data
        int lorentz_type;
        file >> lorentz_type;
        m_lorentz_type = static_cast<LorentzType>(lorentz_type);

        // Allocate storage for Lorentz parameters
        switch (m_lorentz_type) {
        case NEUTRAL:
            m_nglor = 11; break;
        case SINGLE_ION:
            m_nglor = 24; break;
        default:
            std::cout << "Lorentz type not implemented!" << std::endl;
            exit(1);
        }
        mp_lorentz_data = new double [m_ntopbase*m_nglor];

        // Read in the lorentz data
        double* pg;
        size_t i; file >> i;
        while (i <= m_ntopbase && !file.eof()) {
            pg = mp_lorentz_data + (i-1)*m_nglor;
            for (int k = 0; k < m_nglor; ++k)
                file >> pg[k];
            file >> i;
        }

        // Close Lorentz file
        file.close();
    }

    /// Destructor
    virtual ~LorentzTopbase() {}

    void update(const ThermoData& thermo)
    {
        // Get the species indices
        static int species_idx [9];
        static bool first_call = true;
        if (first_call) {
            species_idx[0] = thermo.speciesIndex("O");
            species_idx[1] = thermo.speciesIndex("N");
            species_idx[2] = thermo.speciesIndex("O2");
            species_idx[3] = thermo.speciesIndex("N2");
            species_idx[4] = thermo.speciesIndex("NO");
            species_idx[5] = thermo.speciesIndex("C");
            species_idx[6] = thermo.speciesIndex("CO");
            species_idx[7] = thermo.speciesIndex("CO2");
            species_idx[8] = thermo.speciesIndex("e-");
            first_call = false;
        }

        double Th = thermo.Th();
        double Te = thermo.Te();
        double Tcal = std::min(std::max(Th/1000.0, 2.0), 25.0);
        m_omega = Tcal - (int) Tcal;
        m_ioni = (int) Tcal - 2;
        if ((int) Tcal == 25) {
            m_omega = 1.0;
            m_ioni  = 22;
        }

        switch (m_lorentz_type) {
        case NEUTRAL:
            m_lorentz_fac[8] = 0.5*std::pow(Th/2000.0, 0.3);
            for (int i = 0; i < 8; ++i)
                m_lorentz_fac[i] = m_lorentz_fac[8]*thermo.N(species_idx[i]);
            m_lorentz_fac[8] = 0.0; // ignoring this factor
            m_lorentz_fac[10] = 0.5e-22*thermo.N(species_idx[8]);
            m_lorentz_fac[9] = m_lorentz_fac[10]*std::pow(Te/2000.0, 1.0/6.0);
            break;
        case SINGLE_ION:
            m_lorentz_fac[0] = 1.0e-22*thermo.N(species_idx[8]);
            break;
        default:
            std::cout << "this lorentz type is not supported yet!" << std::endl;
            //exit(1);
        }
    }

    double hwhm(int i)
    {
        double gaml = 0.0;
        double *pg = mp_lorentz_data + m_nglor*m_topbase_id[i];

        switch (m_lorentz_type) {
        case NEUTRAL:
            for (int k = 0; k < 11; ++k)
                gaml += std::abs(pg[k])*m_lorentz_fac[k];
            break;
        case SINGLE_ION:
            gaml = m_lorentz_fac[0]*(pg[m_ioni]*(1.0-m_omega)+pg[m_ioni+1]*m_omega);
            break;
        default:
            std::cout << "This Lorentz type is not handled here." << std::endl;
        }

        return gaml;
    }

    void setHWHM(std::vector<LineData>& lines)
    {
        double *pg;

        switch (m_lorentz_type) {
        case NEUTRAL:
            for (int i = 0; i < lines.size(); ++i) {
                pg = mp_lorentz_data + m_nglor*m_topbase_id[i];
                LineData& line = lines[i];
                line.gaml = 0.0;
                for (int k = 0; k < 11; ++k)
                    line.gaml += std::abs(pg[k])*m_lorentz_fac[k];
            }
            break;
        case SINGLE_ION:
            for (int i = 0; i < lines.size(); ++i) {
                pg = mp_lorentz_data + m_nglor*m_topbase_id[i];
                lines[i].gaml =
                    m_lorentz_fac[0]*(pg[m_ioni]*(1.0-m_omega)+pg[m_ioni+1]*m_omega);
            }
            break;
        default:
            std::cout << "This Lorentz type is not handled here." << std::endl;
        }
    }

    int topbaseID(int i) const {
        return m_topbase_id[i];
    }

private:

    enum LorentzType
    {
        NEUTRAL,
        SINGLE_ION,
        DOUBLE_ION
    } m_lorentz_type;

    double m_lorentz_fac [11];

    int m_nglor;
    int m_ntopbase;
    double* mp_lorentz_data;

    double m_omega;
    int m_ioni;

    std::vector<int> m_topbase_id;
};

#endif // LORENTZ_TOPBASE_H

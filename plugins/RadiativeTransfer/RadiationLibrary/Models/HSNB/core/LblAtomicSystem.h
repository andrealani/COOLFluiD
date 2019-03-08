#ifndef COOLFluiD_RadiativeTransfer_LBL_ATOMIC_SYSTEM_H
#define COOLFluiD_RadiativeTransfer_LBL_ATOMIC_SYSTEM_H

#include <string>
#include <vector>
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/AtomicPartFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/QssAtoms.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

class LblSpectralGrid;

enum SpectralDataType {
    ETA,
    KAPPA,
    ETA_KAPPA
};

struct AtomicLineData {
    AtomicLineData(double gl, double gd, double sc, double e)
        : gaml(gl), gamd(gd), sigc(sc), eta(e) {};
    double gaml;
    double gamd;
    double sigc;
    double eta;
};

/**
 * Represents an atomic LBL radiative system.
 */
class LblAtomicSystem 
{
public:

    /**
     * Constructor taking the name of the atom to load.
     */
    LblAtomicSystem(const COOLFluiD::RadiativeTransfer::ThermoData& thermo, const std::string& name, const std::string datadir, const bool& qss = false);
    
    /**
     * Copy constructor.
     */
    LblAtomicSystem(const LblAtomicSystem& sys)
        : m_name(sys.m_name),
          m_datadir(sys.m_datadir),
          m_index(sys.m_index),
          m_nlines(sys.m_nlines),
          mp_line_data(
              sys.mp_line_data == NULL ? NULL : new LineData [sys.m_nlines]),
          m_lorentz_type(sys.m_lorentz_type),
          m_nglor(sys.m_nglor),
          m_ntopbase(sys.m_ntopbase),
          mp_lorentz_data(
              sys.mp_lorentz_data == NULL ? NULL :
              new double [sys.m_ntopbase*sys.m_nglor]),
          m_part_func(sys.m_part_func),
          m_test_qss(sys.m_test_qss)
    {
        std::copy(
            sys.mp_line_data, sys.mp_line_data+m_nlines, mp_line_data);
        std::copy(
            sys.mp_lorentz_data, sys.mp_lorentz_data+m_ntopbase*m_nglor,
            mp_lorentz_data);
    }

    /**
     * Destructor.
     */
    ~LblAtomicSystem();
    
    /**
     * Assignment operator.
     */
    LblAtomicSystem& operator= (LblAtomicSystem sys)
    {
        swap(*this, sys);
        return *this;
    }

    /**
     * Returns the name of this atom.
     */
    const std::string& name() const { return m_name; }

    /**
     * Computes emmission and absorption coefficient spectra for the atomic LBL
     * system.
     */
    void spectra(COOLFluiD::RadiativeTransfer::ThermoData &, const LblSpectralGrid&, double* const p_spectra,
        const SpectralDataType type);
    
    /**
     * Computes emmission and absorption coefficient spectra for the atomic LBL
     * system. This version ensure that equilibrium is retrieved.
     */
    void spectraEq
      (COOLFluiD::RadiativeTransfer::ThermoData&,
       const LblSpectralGrid&, double* const p_spectra);
    
    /**
     * Generates a list of atomic line data
     */
    void getLineData
      (COOLFluiD::RadiativeTransfer::ThermoData &thermo,
       std::vector<AtomicLineData>& line_data);
    
 private:

    /**
     * Load the line data from the NIST files.
     */
    void loadLineData();

    /**
     * Loads the data necessary for computing Lorentz half-widths.
     */
    void loadLorentzData();

    void loadQssInterface();

    /**
     * Computes the Lorentz HWHM values for all lines and stores in the array.
     */
    void computeLorentzWidths(const COOLFluiD::RadiativeTransfer::ThermoData&, double* const p_gamma);
    
    friend void swap(LblAtomicSystem& s1, LblAtomicSystem& s2);

private:

    std::string m_name;
    size_t      m_index;

    /**
     * Necessary data to generate spectra for a single atomic line
     */
    struct LineData {
        size_t topbase; // line number in TOPBASE database
        double sigma;   // [cm-1] wavenumber of the line
        double Eu;      // [cm-1] upper energy level of the line
        double lngf;    // [-] natural log of g_l*f_lu
        size_t Gqssu;   // Group index associated with the upper level (Johnston's QSS model) 
        size_t Gqssl;   // Group index associated with the lower level (Johnston's QSS model) 
    };
    
    std::size_t m_nlines;
    LineData*   mp_line_data;
    
    enum LorentzType
    {
        NEUTRAL,
        SINGLE_ION,
        DOUBLE_ION,
        HYDROGEN
    } m_lorentz_type;

    std::size_t m_nglor;
    std::size_t m_ntopbase;
    double* mp_lorentz_data;

    AtomicPartFunc m_part_func;

    bool        m_test_qss;
    std::string m_datadir;
};

void swap(LblAtomicSystem& s1, LblAtomicSystem& s2);

#endif // LBL_ATOMIC_SYSTEM_H

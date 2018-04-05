#ifndef COOLFluiD_RadiativeTransfer_LINE_DATA_H
#define COOLFluiD_RadiativeTransfer_LINE_DATA_H

/// Necessary line data
struct LineData
{
    /// Line center in per cm
    double sigc;
    /// Log of product of lower state degeneracy and oscillator strength
    double lnglflu;
    /// Upper state energy level in per cm
    double Eu;
    /// Half-width at half-maximum for Lorentz profile in per cm
    double gaml;
    /// Half-width at half-maximum for Doppler profile in per cm
    double gamd;
    /// Lower state non-Boltzmann population correction factor
    double nonbolt_l;
    /// Upper state non-Boltzmann population correction factor
    double nonbolt_u;
};

#endif // LINE_DATA_H

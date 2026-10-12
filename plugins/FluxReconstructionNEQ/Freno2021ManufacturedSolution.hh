#ifndef COOLFluiD_FluxReconstructionMethod_Freno2021ManufacturedSolution_hh
#define COOLFluiD_FluxReconstructionMethod_Freno2021ManufacturedSolution_hh

//////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "Common/BadValueException.hh"
#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Manufactured solution of B. A. Freno, B. R. Carnes, V. G. Weirs, J. Comput. Phys. 425 (2021) 109752, Eq. (19),
 * for the inviscid two-temperature five-species air equations in 2D. With s = sin, c = cos, e the amplitude
 * (option Epsilon), ub the velocity scale (vb = ub as in the paper), Tb and Tvb the temperature scales:
 *   rho_N2 = rhob_N2 [1 - e s(5/4 pi x) (s(pi y) + c(pi y))]
 *   rho_O2 = rhob_O2 [1 + e s(3/4 pi x) (s(pi y) + c(pi y))]
 *   rho_NO = rhob_NO [1 + e s(pi x) s(pi y)]
 *   rho_N  = rhob_N  [1 + e s(pi x) c(1/4 pi y)]
 *   rho_O  = rhob_O  [1 + e s(pi x) (s(pi y) + c(1/4 pi y))]
 *   u      = ub      [1 + e s(1/4 pi x) (s(pi y) + c(pi y))]
 *   v      = ub      [  - e s(5/4 pi x) s(pi y)]
 *   T      = Tb      [1 + e s(5/4 pi x) (s(pi y) + c(pi y))]
 *   Tv     = Tvb     [1 + e s(3/4 pi x) (s(5/4 pi y) + c(3/4 pi y))]
 * x and y in metres. The species are given by name in the order of the chemistry library, so the shape of each
 * partial density follows its name, not its position. The primitive vector is
 * (rho_0, ..., rho_{ns-1}, u, v, T, Tv), ns being the number of species, in SI units.
 *
 * @author Rayan Dhib
 */
class Freno2021ManufacturedSolution {
public:

  /// Constructor
  Freno2021ManufacturedSolution() :
    m_names(), m_rhoBar(), m_uBar(0.), m_TBar(0.), m_TvBar(0.), m_eps(0.)
  {
  }

  /// molar mass [kg/mol] of a species of the model, used to check the names against the library
  static CFreal molarMass(const std::string& name)
  {
    if (name == "N2") return 28.0134e-3;
    if (name == "O2") return 31.9988e-3;
    if (name == "NO") return 30.0061e-3;
    if (name == "N")  return 14.0067e-3;
    if (name == "O")  return 15.9994e-3;
    throw Common::BadValueException (FromHere(), "Freno2021ManufacturedSolution: unknown species " + name +
                                     " (the solution is defined for N2, O2, NO, N, O)");
    return 0.;
  }

  /**
   * Set the parameters
   * @param names  species names in the order of the chemistry library
   * @param rhoBar bar partial densities [kg/m3], same order as names
   * @param uBar   velocity scale ub [m/s], also used for v
   * @param TBar   temperature scale Tb [K]
   * @param TvBar  vibrational temperature scale Tvb [K]
   * @param eps    amplitude e of the perturbations
   */
  void setup(const std::vector<std::string>& names, const std::vector<CFreal>& rhoBar,
             const CFreal uBar, const CFreal TBar, const CFreal TvBar, const CFreal eps)
  {
    if (names.size() != 5 || rhoBar.size() != names.size())
    {
      throw Common::BadValueException (FromHere(), "Freno2021ManufacturedSolution: give 5 species names and 5 bar densities");
    }
    m_names = names;
    m_rhoBar = rhoBar;
    m_uBar = uBar;
    m_TBar = TBar;
    m_TvBar = TvBar;
    m_eps = eps;
    for (CFuint i = 0; i < names.size(); ++i)
    {
      molarMass(names[i]); // throws for an unknown name
      for (CFuint j = 0; j < i; ++j)
      {
        if (names[j] == names[i])
        {
          throw Common::BadValueException (FromHere(), "Freno2021ManufacturedSolution: species " + names[i] + " given twice");
        }
      }
    }
  }

  /**
   * Check the species names against the molar masses of the library (relative difference below 1 percent; the
   * closest pair of the model, NO and O2, differs by 6.6 percent)
   * @param mm molar masses [kg/mol] in the order of the library
   */
  void checkMolarMasses(const RealVector& mm) const
  {
    if (mm.size() != m_names.size())
    {
      throw Common::BadValueException (FromHere(), "Freno2021ManufacturedSolution: the library has a different number of species");
    }
    for (CFuint i = 0; i < m_names.size(); ++i)
    {
      const CFreal ref = molarMass(m_names[i]);
      if (std::abs(mm[i] - ref) > 0.01*ref)
      {
        std::ostringstream msg;
        msg << "Freno2021ManufacturedSolution: species " << i << " is named " << m_names[i]
            << " but the library molar mass is " << mm[i] << " kg/mol (expected " << ref << ")";
        throw Common::BadValueException (FromHere(), msg.str());
      }
    }
  }

  /// number of species
  CFuint nbSpecies() const {return m_names.size();}

  /// names of the primitive fields: species names, then u, v, T, Tv
  std::vector<std::string> fieldNames() const
  {
    std::vector<std::string> f = m_names;
    f.push_back("u"); f.push_back("v"); f.push_back("T"); f.push_back("Tv");
    return f;
  }

  /// scale of primitive field i (bar density, ub for u and v, Tb, Tvb), used to normalise the errors
  CFreal barValue(const CFuint i) const
  {
    const CFuint ns = m_names.size();
    if (i < ns) return m_rhoBar[i];
    if (i < ns + 2) return m_uBar;
    if (i == ns + 2) return m_TBar;
    return m_TvBar;
  }

  /**
   * Evaluate the primitive vector at a point
   * @param x, y coordinates [m]
   * @param prim output (rho_0, ..., rho_{ns-1}, u, v, T, Tv), size ns + 4
   */
  void evaluate(const CFreal x, const CFreal y, RealVector& prim) const
  {
    const CFreal pi = 3.141592653589793;
    const CFreal spx  = std::sin(pi*x);
    const CFreal spy  = std::sin(pi*y);
    const CFreal cpy  = std::cos(pi*y);
    const CFreal s54x = std::sin(1.25*pi*x);
    const CFreal s34x = std::sin(0.75*pi*x);
    const CFreal s14x = std::sin(0.25*pi*x);
    const CFreal c14y = std::cos(0.25*pi*y);
    const CFreal s54y = std::sin(1.25*pi*y);
    const CFreal c34y = std::cos(0.75*pi*y);

    const CFuint ns = m_names.size();
    for (CFuint i = 0; i < ns; ++i)
    {
      const std::string& n = m_names[i];
      CFreal g = 0.;
      if      (n == "N2") g = -s54x*(spy + cpy);
      else if (n == "O2") g =  s34x*(spy + cpy);
      else if (n == "NO") g =  spx*spy;
      else if (n == "N")  g =  spx*c14y;
      else if (n == "O")  g =  spx*(spy + c14y);
      prim[i] = m_rhoBar[i]*(1. + m_eps*g);
    }
    prim[ns]     = m_uBar*(1. + m_eps*s14x*(spy + cpy));
    prim[ns + 1] = -m_uBar*m_eps*s54x*spy;
    prim[ns + 2] = m_TBar*(1. + m_eps*s54x*(spy + cpy));
    prim[ns + 3] = m_TvBar*(1. + m_eps*s34x*(s54y + c34y));
  }

private:

  /// species names in the library order; the name selects the shape of the partial density
  std::vector<std::string> m_names;

  /// bar partial densities [kg/m3], library order
  std::vector<CFreal> m_rhoBar;

  /// velocity scale ub [m/s] (vb = ub)
  CFreal m_uBar;

  /// temperature scale Tb [K]
  CFreal m_TBar;

  /// vibrational temperature scale Tvb [K]
  CFreal m_TvBar;

  /// amplitude e of the perturbations
  CFreal m_eps;

}; // class Freno2021ManufacturedSolution

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_Freno2021ManufacturedSolution_hh

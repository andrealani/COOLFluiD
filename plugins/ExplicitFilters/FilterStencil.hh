#ifndef COOLFluiD_Numerics_ExplicitFilters_FilterStencil_hh
#define COOLFluiD_Numerics_ExplicitFilters_FilterStencil_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MathConsts.hh"
#include <deque>


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  	namespace ExplicitFilters {
      
//////////////////////////////////////////////////////////////////////////////

class FilterStencil
{
  public:
    FilterStencil() : m_mustCompute(true), m_distanceToBoundary(MathTools::MathConsts::CFrealInf()) {}
    FilterStencil(const FilterStencil& S) : m_cellWidth(S.m_cellWidth), 
                                       m_mustCompute(S.m_mustCompute), 
                                       m_radius(S.m_radius), 
                                       m_distanceToBoundary(S.m_distanceToBoundary),
                                       m_stencilElements(S.m_stencilElements),
                                       m_isGhost(S.m_isGhost){}
    CFuint                                  getID() const { return m_stencilElements[0]; }
    bool                                    mustCompute()const { return m_mustCompute; }
    void                                    clear() {m_stencilElements.resize(0); m_stencilElements.clear(); m_isGhost.resize(0); m_isGhost.clear(); }
    void                                    setCompute(const bool mustCompute) { m_mustCompute = mustCompute;}
    Common::SafePtr<std::vector<CFuint> >   getElements() { return &m_stencilElements; }
    CFuint                                  getElement(const CFuint& idx) const { return m_stencilElements[idx];}
    void                                    addElement(const CFuint& elementID, const bool isGhost = false) {m_stencilElements.push_back(elementID); m_isGhost.push_back(isGhost);}
    CFreal                                  getRadius() const { return m_radius; }
    void                                    setRadius(const CFreal& radius) { m_radius = radius;}
    void                                    enlargeRadiusWithFactor(const CFreal& enlargementFactor);
    void                                    setDistanceToBoundary(const CFreal& distance) { m_distanceToBoundary = std::min(distance,m_distanceToBoundary); }
    CFreal                                  getDistanceToBoundary() const { return m_distanceToBoundary; }
    CFreal                                  getCellWidth() const { return m_cellWidth; }
    void                                    setCellWidth(const CFreal& cellWidth) { m_cellWidth = cellWidth;}
    CFuint                                  getNbElements() const { return m_stencilElements.size(); }
    bool                                    isGhost(const CFuint& idx) const{return m_isGhost[idx]; }
    bool                                    isNotIncluded(const CFuint& elementID, const bool isGhost = false) const;
    void                                    addElement(const Common::SafePtr<FilterStencil>& stencil, const CFuint& idx) { m_stencilElements.push_back(stencil->getElement(idx)); m_isGhost.push_back(stencil->isGhost(idx));  }
    void                                    addElement(const FilterStencil& stencil, const CFuint& idx) { m_stencilElements.push_back(stencil.getElement(idx)); m_isGhost.push_back(stencil.isGhost(idx));  }

  protected:
    CFreal              m_cellWidth;
    bool                m_mustCompute;
    CFreal              m_radius;
    CFreal              m_distanceToBoundary;
    std::vector<CFuint> m_stencilElements;
    std::deque<bool>    m_isGhost;
};
  
//////////////////////////////////////////////////////////////////////////////
  
inline bool FilterStencil::isNotIncluded(const CFuint& elementID, const bool isGhost) const
{
  CFuint nbElems = getNbElements();
  for (CFuint i=0; i<nbElems; i++) {
    if (elementID == m_stencilElements[i] && isGhost == m_isGhost[i]) {
      return false;
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

inline void FilterStencil::enlargeRadiusWithFactor(const CFreal& enlargementFactor)
{
  m_radius *= enlargementFactor;
}

//////////////////////////////////////////////////////////////////////////////

		} // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_FilterStencil_hh

#ifndef COOLFluiD_client_ShutdownConfirmationPanel_h
#define COOLFluiD_client_ShutdownConfirmationPanel_h

#include "ClientServer/client/CloseConfirmationPanel.hh"

class QComboBox;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    
    struct CloseConfirmationInfos;
    
    /////////////////////////////////////////////////////////////////////////////
    
    class ShutdownConfirmationPanel : public CloseConfirmationPanel
    {
      
    public:
      
      ShutdownConfirmationPanel(QDialog * parent = NULL);
      
      ~ShutdownConfirmationPanel();
      
      virtual bool isAccepted() const;
      
      virtual void getData(CloseConfirmationInfos & infos) const;
      
      virtual void setData(const CloseConfirmationInfos & infos);
      
    protected:
      
      virtual void hideComponents(bool hide);
      
    private:
      
      QComboBox * m_comboBox;
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_ShutdownConfirmationPanel_h

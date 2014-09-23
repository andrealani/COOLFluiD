#ifndef COOLFluiD_client_CommitConfirmationPanel_h
#define COOLFluiD_client_CommitConfirmationPanel_h

#include "ClientServer/client/CloseConfirmationPanel.hh"

class QComboBox;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    
    struct CloseConfirmationInfos;
    /////////////////////////////////////////////////////////////////////////////
    
    class CommitConfirmationPanel : public CloseConfirmationPanel
    {
      Q_OBJECT 
      
    public:
      
      CommitConfirmationPanel(QDialog * parent = NULL);
      
      ~CommitConfirmationPanel();
      
      virtual bool isAccepted() const;
      
      virtual void getData(CloseConfirmationInfos & infos) const;
      
      virtual void setData(const CloseConfirmationInfos & infos);
      
    protected:
      
      virtual void hideComponents(bool hide);
      
      private slots:
      
      void showDetails();
      
    private:
      
      QComboBox * comboBox;
      
      QString details;
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_CommitConfirmationPanel_h

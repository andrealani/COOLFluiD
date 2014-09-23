#ifndef COOLFluiD_client_AboutCFDialog_h
#define COOLFluiD_client_AboutCFDialo.hh

#include <QDialog>
#include <QList>

class QFormLayout;
class QLabel;
class QPushButton;
class QVBoxLayout;
class QWidget;

namespace COOLFluiD
{
  namespace client
  {
    class AboutCFDialog : public QDialog
    {
      struct CFInfo
      {
      public:
        QLabel * labName;
        QLabel * labValue;
        
        CFInfo(const QString & name, const QString & value, QFormLayout * parent);
        
        ~CFInfo();
      };
      
    public:
      
      AboutCFDialog(QWidget * parent = NULL);
      
      ~AboutCFDialog();
      
    private: // data
      
      QVBoxLayout * m_mainLayout;
      
      QPushButton * m_btOK;
      
      QFormLayout * m_infoLayout;
      
      QList<CFInfo *> m_infoList;
      
    };
  }
}

#endif // COOLFluiD_client_AboutCFDialo.hh

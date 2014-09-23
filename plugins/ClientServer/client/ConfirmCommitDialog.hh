#ifndef COOLFluiD_client_ConfirmCommitDialog_h
#define COOLFluiD_client_ConfirmCommitDialog_h

#include <QDialog>
#include <QHash>
#include <QDialogButtonBox>

class QPushButton;
class QLabel;
class QVBoxLayout;
class QTableView;
class QString;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD 
{
  namespace client 
  {
    class CommitDetails;
    
    ////////////////////////////////////////////////////////////////////////////
    
    class ConfirmCommitDialog : public QDialog
    {
      
    public:
      
      enum CommitConfirmation 
      {
        COMMIT,
        
        DONT_COMMIT,
        
        YES,
        
        NO,
        
        CANCEL
      };
      
      Q_OBJECT
      
    public:
      
      ConfirmCommitDialog(QWidget * parent = NULL);
      
      ~ConfirmCommitDialog();
      
      ConfirmCommitDialog::CommitConfirmation show(CommitDetails & commitDetails);
      
      private slots:
      
      void buttonClicked();
      
    private: // methods
      
      void createButton(const QString & text, CommitConfirmation commConf, 
                        QDialogButtonBox::ButtonRole role);
      
    private: // data
      QDialogButtonBox * m_buttonBox;
      
      QVBoxLayout * m_mainLayout;
      
      QTableView * m_detailsView;
      
      QLabel * m_labText;
      
      CommitConfirmation m_answer;
      
      QHash<CommitConfirmation, QPushButton *> m_buttons; 
      
    }; // class ConfirmCommitDialog
    
    ////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

////////////////////////////////////////////////////////////////////////////// 

#endif // COOLFluiD_client_ConfirmCommitDialog_h
#ifndef COOLFluiD_client_CommitDetailsDialog_h
#define COOLFluiD_client_CommitDetailsDialog_h

#include <QDialog>

class QPushButton;
class QDialogButtonBox;
class QVBoxLayout;
class QTableView;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD 
{
  namespace client
  {
    class CommitDetails;
    
    //////////////////////////////////////////////////////////////////////////////
    
    class CommitDetailsDialog : public QDialog
    {
      Q_OBJECT
      
    public:
      
      CommitDetailsDialog(QWidget * parent = NULL);
      
      ~CommitDetailsDialog();
      
      void setCommitDetails(CommitDetails * details);
      
      void show(CommitDetails & details);
      
      public slots:
      
      void show();
      
    private:
      
      QPushButton * m_btOk;
      
      QTableView * m_view;
      
      QDialogButtonBox * m_buttonBox;
      
      QVBoxLayout * m_mainLayout;
      
      CommitDetails * m_commitDetails;
      
    }; // class CommitDetailsDialog
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_CommitDetailsDialog_h

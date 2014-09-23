#ifndef COOLFluiD_client_StatusPanel_h
#define COOLFluiD_client_StatusPanel_h

#include <QTreeView>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    class StatusModel;
    
    //////////////////////////////////////////////////////////////////////////////
    
    class StatusPanel : public QTreeView
    {
      Q_OBJECT
      
    public:
      
      StatusPanel(StatusModel * model, QWidget * parent = NULL);
      
      ~StatusPanel();
      
      private slots:
      
      void subSystemAdded(const QModelIndex & index);
      
    private:
      
      StatusModel * m_model;
      
      
    };
    
    //////////////////////////////////////////////////////////////////////////////
    
  }
}

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_StatusPanel_h
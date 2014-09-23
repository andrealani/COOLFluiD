#ifndef COOLFluiD_client_CommitDetails_h
#define COOLFluiD_client_CommitDetails_h

#include <QHash>
#include <QAbstractItemModel>

#include "ClientServer/client/OptionTypes.hh"

class QStringList;
class QString;

namespace COOLFluiD
{
  namespace client
  {
    class CommitDetailsItem;
    
    class CommitDetails : public QAbstractItemModel
    {
      
    public:
      
      CommitDetails(QObject * parent = NULL, const QString & nodePath = QString());
      
      QVariant data(const QModelIndex &index, int role) const;
      
      QVariant headerData(int section, Qt::Orientation orientation,
                          int role = Qt::DisplayRole) const;
      
      QModelIndex index(int row, int column,
                        const QModelIndex &parent = QModelIndex()) const;
      
      QModelIndex parent(const QModelIndex &index) const;
      
      int rowCount(const QModelIndex &parent = QModelIndex()) const;
      
      int columnCount(const QModelIndex &parent = QModelIndex()) const;
      
      void setOption(const QString & optionName, const QString & oldValue, 
                     const QString & currentValue);
      
      void setNewOption(const QString & optionName, const QString & value,
                        TOptionTypes type = NO_TYPE);
      
      bool setOptionNewValue(const QString & optionName, const QString & value);
      
      bool setOptionOldValue(const QString & optionName, const QString & value);
      
      bool setNewOptionValue(const QString & optionName, const QString & value);
      
      bool removeOption(const QString & optionName);
      
      bool removeNewOption(const QString & optionName);
      
      bool getOption(const QString & optionName, QString & oldValue, 
                     QString & newValue) const;
      
      bool getNewOption(const QString & optionName, QString & value) const;
      
      TOptionTypes getNewOptionType(const QString & optionName) const;
      
      QString getNodePath() const;
      
      void setNodePath(const QString & nodePath);
      
      bool contains(const QString & optionName, bool * isNewOption = NULL) const;
      
      bool isEmpty() const;
      
      bool hasOptions() const;
      
      bool hasNewOptions() const;
      
      // clears node path
      void clear();
      
      void clearOptions();
      
      void clearNewOptions();
      
      int getOptionCount() const;
      
      int getNewOptionCount() const;
      
      QString toString() const;
      
    private:
      
      QString m_nodePath;
      
      QStringList m_options;
      
      QStringList m_newOptions;
      
      QHash<QString, QString> m_optionsOldValues;
      
      QHash<QString, QString> m_optionsNewValues;
      
      QHash<QString, QString> m_newOptionsValues;
      
      QHash<QString, TOptionTypes> m_newOptionsTypes;
      
      QList<CommitDetailsItem *> m_items;
      
    };
  }
}

#endif // COOLFluiD_client_CommitDetails_h
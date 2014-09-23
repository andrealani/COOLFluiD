// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * HierarchyMaintainer.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_HIERARCHYMAINTAINER_HH
#define CF_LOGCPP_HIERARCHYMAINTAINER_HH

#include <logcpp/Portability.hh>
#include <string>
#include <map>
#include <vector>
#include <logcpp/Category.hh>
#include <logcpp/Threading.hh>

namespace logcpp {

    /**
     * HierarchyMaintainer is an internal logcpp class. It is responsible
     * for maintaining the hierarchy of Categories. Applications should
     * not have to use this class directly.
     **/
    class HierarchyMaintainer {
        friend class Log4cppCleanup;

        public:
        typedef std::map<std::string, Category*> CategoryMap;
  
        static HierarchyMaintainer& getDefaultMaintainer();

        HierarchyMaintainer();
        virtual ~HierarchyMaintainer();
        virtual Category* getExistingInstance(const std::string& name);
        virtual Category& getInstance(const std::string& name);
        virtual std::vector<Category*>* getCurrentCategories() const;
        virtual void shutdown();
        virtual void deleteAllCategories();

        protected:
        virtual Category* _getExistingInstance(const std::string& name);
        virtual Category& _getInstance(const std::string& name);
        CategoryMap _categoryMap;
        mutable threading::Mutex _categoryMutex;

        private:
        static HierarchyMaintainer* _defaultMaintainer;
    };        
}

#endif // CF_LOGCPP_HIERARCHYMAINTAINER_HH

/*
 * HierarchyMaintainer.cpp
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#include "PortabilityImpl.hh"

#ifdef LOGCPP_HAVE_IO_H
#    include <io.h>
#endif
#ifdef CF_HAVE_UNISTD_H
#    include <unistd.h>
#endif

#include <cstdio>
#include <logcpp/HierarchyMaintainer.hh>
#include <logcpp/FileAppender.hh>

namespace logcpp {

    HierarchyMaintainer* HierarchyMaintainer::_defaultMaintainer = NULL;

    HierarchyMaintainer& HierarchyMaintainer::getDefaultMaintainer() {
        if (!_defaultMaintainer)
            _defaultMaintainer = new HierarchyMaintainer();

        return *_defaultMaintainer;
    }

    HierarchyMaintainer::HierarchyMaintainer() {
    }

    HierarchyMaintainer::~HierarchyMaintainer() {
        shutdown();
        deleteAllCategories();
    }

    Category* HierarchyMaintainer::getExistingInstance(const std::string& name) {
        threading::ScopedLock lock(_categoryMutex);
        return _getExistingInstance(name);
    }

    Category* HierarchyMaintainer::_getExistingInstance(const std::string& name) {
	Category* result = NULL;

        CategoryMap::iterator i = _categoryMap.find(name);
        if (_categoryMap.end() != i) {
            result = (*i).second;
        }

	return result;
    }

    Category& HierarchyMaintainer::getInstance(const std::string& name) {
        threading::ScopedLock lock(_categoryMutex);
        return _getInstance(name);
    }

    /* assume lock is held */
    Category& HierarchyMaintainer::_getInstance(const std::string& name) {
        Category* result;
        result = _getExistingInstance(name);
        
        if (NULL == result) {            
            if (name == "") {
                result = new Category(name, NULL, Priority::INFO);
                result->addAppender(new FileAppender("_", ::dup(fileno(stderr))));
            } else {
                std::string parentName;
                size_t dotIndex = name.find_last_of('.');
                if (name.length() <= dotIndex) {
                    parentName = "";
                } else {
                    parentName = name.substr(0, dotIndex);
                }
                Category& parent = _getInstance(parentName);
                result = new Category(name, &parent, Priority::NOTSET);
            }	  
            _categoryMap[name] = result; 
        }
        return *result;
    }

    std::vector<Category*>* HierarchyMaintainer::getCurrentCategories() const {
        std::vector<Category*>* categories = new std::vector<Category*>;

        threading::ScopedLock lock(_categoryMutex);
        {
            for(CategoryMap::const_iterator i = _categoryMap.begin(); i != _categoryMap.end(); i++) {
                categories->push_back((*i).second);
            }
        }

        return categories;
    }

    void HierarchyMaintainer::shutdown() {
        threading::ScopedLock lock(_categoryMutex);
        {
            for(CategoryMap::const_iterator i = _categoryMap.begin(); i != _categoryMap.end(); i++) {
                ((*i).second)->removeAllAppenders();
            }
        }
    }

    void HierarchyMaintainer::deleteAllCategories() {
        threading::ScopedLock lock(_categoryMutex);
        {
            for(CategoryMap::const_iterator i = _categoryMap.begin(); i != _categoryMap.end(); i++) {
                delete ((*i).second);
            }
        }
    }

}

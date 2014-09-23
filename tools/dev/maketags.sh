#!/bin/bash
ctags -R --C++-types=+pxl --excmd=pattern --langmap=C++:+.ci.cpp
exit 0

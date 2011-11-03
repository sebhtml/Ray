/*

Ray
Copyright (C)  2010  Sébastien Boisvert

http://DeNovoAssembler.SourceForge.Net/

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You have received a copy of the GNU General Public License
along with this program (COPYING).
see <http://www.gnu.org/licenses/>

*/

#ifndef __MALLOC_TYPES
#define __MALLOC_TYPES

extern const char* MALLOC_TYPES[];

#define MACRO_LIST_ITEM(x) x,

enum{
#include <memory/malloc_type_macros.h>
MALLOC_TYPES_H_DUMMY
};

#undef MACRO_LIST_ITEM

#endif

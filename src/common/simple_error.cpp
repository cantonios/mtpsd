/*

This file is part of MTPSD.

MTPSD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MTPSD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MTPSD.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "simple_error.h"
#include <string.h>

ERR::ERR(){
    this->msg[0]='\0';
}

ERR::ERR(const char *errmsg){
    strncpy(msg,errmsg,30);
}
void ERR::getmsg(char *errmsg){
    strncpy(errmsg,this->msg,30);
}
const char *ERR::getmsg(){
    return msg;
}

LAPACK_ERROR::LAPACK_ERROR(const char *errmsg) : ERR(errmsg){};
LAPACK_ERROR::LAPACK_ERROR() : ERR(){};
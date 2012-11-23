/*

Copyright (c) 2010 C. Antonio Sanchez

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


#ifndef MTPSD_SIMPLE_ERROR
#define MTPSD_SIMPLE_ERROR

//Error classes for LAPACK, and general
class ERR{
    public:    
        ERR();
        ERR(const char *msg);           //sets error message
        void getmsg(char *errmsg);      //copies error message to errmsg
        const char *getmsg();           //returns pointer to error message (for use with printf)
    protected:
        char msg[30];
};

//LAPACK_ERROR inherits everything from ERR
//    Defined so it can be caught separately
class LAPACK_ERROR : public ERR{
    public:
        LAPACK_ERROR();
        LAPACK_ERROR(const char *errmsg);
};

#endif
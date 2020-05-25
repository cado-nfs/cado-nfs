//%LICENSE////////////////////////////////////////////////////////////////
//
// Licensed to The Open Group (TOG) under one or more contributor license
// agreements.  Refer to the OpenPegasusNOTICE.txt file distributed with
// this work for additional information regarding copyright ownership.
// Each contributor licenses this file to you under the OpenPegasus Open
// Source License; you may not use this file except in compliance with the
// License.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//////////////////////////////////////////////////////////////////////////
/*
//
// The above license is compatible with the LGPL 2.1+ under which CADO-NFS
// is distributed; it is essentially the Expat license:
// http://directory.fsf.org/wiki/License:Expat
//
//%/////////////////////////////////////////////////////////////////////////////
*/

#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <string.h>
#include "portability.h"
#include "macros.h"

#define PEGASUS_TEST_ASSERT(X) ASSERT_ALWAYS(X)
#define Strlcpy strlcpy

int main()
{
    {
        char buffer[1024];
        size_t n;
        strcpy(buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

        n = Strlcpy(buffer, "abc", sizeof(buffer));
        PEGASUS_TEST_ASSERT(n == 3);
        PEGASUS_TEST_ASSERT(strcmp(buffer, "abc") == 0);
    }

    {
        char buffer[1024];
        size_t n;
        strcpy(buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

        n = Strlcpy(buffer, "abcdefg", 4);
        PEGASUS_TEST_ASSERT(n == 7);
        PEGASUS_TEST_ASSERT(strcmp(buffer, "abc") == 0);
    }

    {
        char buffer[1024];
        size_t n;
        strcpy(buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

        n = Strlcpy(buffer, "", 4);
        PEGASUS_TEST_ASSERT(n == 0);
        PEGASUS_TEST_ASSERT(strcmp(buffer, "") == 0);
    }

    {
        char buffer[1024];
        size_t n;
        strcpy(buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

        n = Strlcpy(buffer, "", 1);
        PEGASUS_TEST_ASSERT(n == 0);
        PEGASUS_TEST_ASSERT(strcmp(buffer, "") == 0);
    }

    {
        char buffer[1024];
        size_t n;
        strcpy(buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

        n = Strlcpy(buffer, "1234567890", 1);
        PEGASUS_TEST_ASSERT(n == 10);
        PEGASUS_TEST_ASSERT(strcmp(buffer, "") == 0);
    }

    {
        char buffer[1024];
        size_t n;
        strcpy(buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

        n = Strlcpy(buffer, "1234567890", 2);
        PEGASUS_TEST_ASSERT(n == 10);
        PEGASUS_TEST_ASSERT(strcmp(buffer, "1") == 0);
    }

    {
        char buffer[1024];
        size_t n;
        strcpy(buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

        n = Strlcpy(buffer, "1234567890", 0);
        PEGASUS_TEST_ASSERT(n == 10);
        PEGASUS_TEST_ASSERT(strcmp(
            buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx") == 0);
    }

    printf("+++++ passed all tests\n");

    return 0;
}

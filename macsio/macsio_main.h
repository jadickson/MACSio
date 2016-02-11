#ifndef MACSIO_MAIN_H
#define MACSIO_MAIN_H
/*
Copyright (c) 2015, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Mark C. Miller

LLNL-CODE-676051. All rights reserved.

This file is part of MACSio

Please also read the LICENSE file at the top of the source code directory or
folder hierarchy.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License (as published by the Free Software
Foundation) version 2, dated June 1991.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <json-cwx/json_object.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_MPI
extern MPI_Comm MACSIO_MAIN_Comm;
#else
extern int MACSIO_MAIN_Comm;
#endif
extern int MACSIO_MAIN_Size;
extern int MACSIO_MAIN_Rank;

extern int MACSIO_MAIN_GetRankOwningPart(json_object *main_obj, int chunkId);

#ifdef __cplusplus
}
#endif

#endif /* MACSIO_MAIN_H */

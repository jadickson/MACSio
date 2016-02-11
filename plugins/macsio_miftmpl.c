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

#include <json-cwx/json.h>

#include <macsio_clargs.h>
#include <macsio_iface.h>
#include <macsio_log.h>
#include <macsio_main.h>
#include <macsio_mif.h>
#include <macsio_utils.h>

#include <stdio.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*!
\defgroup plugins Plugins
@{
*/

#warning ADD PWRITE OPTION
#warning ADD MPI_FILE_WRITE OPTION
#warning ADD AGGREGATION OPTION

/*!
\addtogroup MIF Template
\brief A simple MIF Plugin Template

Writing a new plugin for MACSio involves a mininum of two new files in the plugins directory.
One is the C or C++ source code for the plugin. The other is a \c .make file that includes
various plugin-specific variable definitions and logic to decode library dependencies that
effect.

We'll describe the structure and content of the \c .make using a generic plugin identified
by the moniker \c plgn. In any plugin's \c .make file, all the variable names should be
prepended with the plugin's name except for these three (examples below).

\code
PLUGIN_OBJECTS += $(PLGN_SOURCES:.c=.o)
PLUGIN_LDFLAGS += $(PLGN_LDFLAGS)
PLUGIN_LIST += plgn
\endcode

The first variable is the \c BUILD_ORDER variable. This variable is used to sort
the order in which plugin object files as well as their dependent libraries appear
on the link line when linking the MACSio main executable. The numerical value
assigned to the \c BUILD_ORDER variable is a floating point number. Smaller numbers
appear earlier on the link line. If plugin A depends on libraries that are also 
used by plugin B, then plugin A should be assigned a \c BUILD_ORDER value that is
larger than plugin B.

\code
PLGN_BUILD_ORDER = 1.7
\endcode

The next variables in the \c .make file define the version number, tarfile and URL
for zero or more dependent libraries used by the plugin.

\code
PLGN_VERSION = 1.8.11
PLGN_FILE = plgn-$(PLGN_VERSION).tar.gz
PLGN_URL = http://www.greatplugins.org/ftp/PLGN/releases/$(PLGN_FILE)
\endcode

In addition, for each third party library a plugin uses, it is assumed there is
an associated \em home variable that specifies a top-level directory underneath
which lives \c include and \c lib directories for the library's header files and
library files respectively. If you have a package that does not conform to this
standard installation structure, the work-around is to use symlinks or explicit
copies to create some \em proxy home directory for the library that is
structured in the way MACSio's Makefiles need it. So, for our generic plugin,
there is a \c PLGN_HOME variable that defines where the library is installed.
Note that this variable is then typically defined in the \c config-site file.

The next section of the makefile is the conditional logic necessary to decide
if the plugin can be built based on whether its depndent library(s), if any,
are defined.

\code
ifneq ($(PLGN_HOME),)

PLGN_LDFLAGS = -L$(PLGN_HOME)/lib -lplgn
PLGN_CFLAGS = -I$(PLGN_HOME)/include

PLGN_SOURCES = macsio_plgn.c

PLGN_LDFLAGS += -lz -lm

PLUGIN_OBJECTS += $(PLGN_SOURCES:.c=.o)
PLUGIN_LDFLAGS += $(PLGN_LDFLAGS)
PLUGIN_LIST += plgn 

endif
\endcode

The next section of the \c .make file indicates how to make the plugin
object file.

\code
macsio_plgn.o: ../plugins/macsio_plgn.c
	$(CXX) -c $(PLGN_CFLAGS) $(MACSIO_CFLAGS) $(CFLAGS) ../plugins/macsio_plgn.c
\endcode

The final section of the \c .make file defines some pre-defined targets which having
to do with obtaining the third party library(s) associated with the plugin.

\code
$(PLGN_FILE):
	$(DLCMD) $(PLGN_FILE) $(PLGN_URL)

list-tpls-plgn:
	@echo "$(PLGN_FILE) ($(PLGN_URL))"

download-tpls-plgn: $(PLGN_FILE)
\endcode

The \em miftmpl plugin is intended to serve as a template for how to create a basic MIF-mode plugin
for MACSio.  This template code does indeed actually function correctly as a MACSio plugin.
It does so by writing MACSio's internal JSON objects repesenting each mesh part as ascii
strings to the individual files.

Each processor in a MIF group serializes each JSON object representing a mesh part to an
ascii string. Then, each of these strings is appended to the end of the file. For each
such string, the plugin maintains knowledge of the mesh part's ID, the filename it was
written to and the offset within the file.

The filenames, offsets and mesh part IDs are then written out as a separate JSON object
to a root or master file. Currently, there is no plugin in VisIt to read these files
and display them. But, this example code does help to outline the basic work to write a
MIF plugin.

In practice, this plugin could have simply written the entire JSON object from each
processor to its MIF group's file. However, in doing that, the resulting file would
not "know" such things as how many mesh parts there are or where a given
mesh part is located in the file set. So, we wind up writing JSON objects for each
part individually so that we can keep track of where they all are in the fileset.

Some of the aspects of this plugin code exist here only to serve as an example in
writing a MIF plugin and are non-essential to the proper operation of this plugin.

MACSio uses a static load approach to its plugins. The MACSio main executable must
be linked with all the plugins it expects to use.

In writing any MACSio plugin (MIF or SIF) be sure to declare *all* of your plugin's
symbols (functions, local variables, etc.) as static. Each plugin is being linked into
MACSio's main and any symbols that are not static file scope will wind up appearing
in and therefore being vulnerable too global namespace collisions. The plugin's
main interface methods to MACSio are handled via registration of a set of function
pointers.
@{
*/

static char const *iface_name = "miftmpl"; /**< Name of the interface this plugin uses */
static char const *iface_ext = "json";     /**< Default file extension for files generated by this plugin */
static int json_as_html = 0;               /**< Use HTML output instead of raw ascii */
static int my_opt_one;                     /**< Example of a static scope, plugin-specific variable to be set in
                                                process_args to control plugin behavior */
static int my_opt_two;                     /**< Another example variable to control plugin behavior */
static char *my_opt_three_string;          /**< Another example variable to control plugin behavior */
static float my_opt_three_float;           /**< Another example variable to control plugin behavior */

/*!
\brief Process command-line arguments specific to this plugin

Uses MACSIO_CLARGS_ProcessCmdline() to do its work.

This example plugin is implemented to route command line arguments to memory locations
(e.g. static variables) here in the plugin.  Alternatively, a plugin can choose to
route the results of MACSIO_CLARGS_ProcessCmdline() to a JSON object. MACSio's main is
implemented that way.
*/
static int process_args(
    int argi,      /**< [in] Argument index of first argument that is specific to this plugin */
    int argc,      /**< [in] argc as passed into main */
    char *argv[]   /**< [in] argv as passed into main */
)
{
    /* Can use MACSIO_CLARGS_TOJSON here instead in which case pass the pointer to
       a json_object* as first arg and eliminate all the pointers to specific
       variables. The args will be returned as a json-c object. */
    const MACSIO_CLARGS_ArgvFlags_t argFlags = {MACSIO_CLARGS_WARN, MACSIO_CLARGS_TOMEM};

    MACSIO_CLARGS_ProcessCmdline(0, argFlags, argi, argc, argv,
        "--json_as_html", "",
            "Write files as HTML instead of raw ascii [false]",
            &json_as_html,
        "--my_opt_one", "",
            "Help message for my_opt_one which has no arguments. If present, local\n"
            "var my_opt_one will be assigned a value of 1 and a value of zero otherwise.",
            &my_opt_one,
        "--my_opt_two %d", MACSIO_CLARGS_NODEFAULT,
            "Help message for my_opt_two which has a single integer argument",
            &my_opt_two,
        "--my_opt_three %s %f", MACSIO_CLARGS_NODEFAULT,
            "Help message for my_opt_three which has a string argument and a float argument",
            &my_opt_three_string, &my_opt_three_float,
           MACSIO_CLARGS_END_OF_ARGS);

    return 0;
}

/*!
\brief CreateFile MIF Callback

This implments the MACSIO_MIF_CreateFile callback needed for a MIF mode plugin.

\return A void pointer to the plugin-specific file handle
*/
static void *CreateMyFile( 
    const char *fname,     /**< [in] Name of the MIF file to create */
    const char *nsname,    /**< [in] Name of the namespace within the file for caller should use. */
    void *userData         /**< [in] Optional plugin-specific user-defined data */
)
{
    FILE *file = fopen(fname, "w");
    return (void *) file;
}

/*!
\brief OpenFile MIF Callback

This implments the MACSIO_MIF_OpenFile callback needed for a MIF mode plugin.

\return A void pointer to the plugin-specific file handle
*/
static void *OpenMyFile(
    const char *fname,            /**< [in] Name of the MIF file to open */
    const char *nsname,           /**< [in] Name of the namespace within the file caller should use */
    MACSIO_MIF_ioFlags_t ioFlags, /**< [in] Various flags indicating behavior/options */
    void *userData                /**< [in] Optional plugin-specific user-defined data */
)
{
    FILE *file = fopen(fname, "a+");
    return (void *) file;
}

/*!
\brief CloseFile MIF Callback

This implments the MACSIO_CloseFile callback needed for a MIF mode plugin.
*/
static void CloseMyFile(
    void *file,      /**< [in] A void pointer to the plugin specific file handle */
    void *userData   /**< [in] Optional plugin specific user-defined data */
)
{
    fclose((FILE*) file);
}

/*!
\brief Write a single mesh part to a MIF file

All this method does is serialize the JSON object for the given mesh
part to an ASCII string and then appends/writes that string at the
end of the current file.

After serializing the object to an ASCII string and writing it to the
file, the memory for the ASCII string is released by json_object_free_printbuf().

\return A tiny JSON object holding the name of the file, the offset at
which the JSON object for this part was written in the file and the part's ID.
*/
static json_object *write_mesh_part(
    FILE *myFile,          /**< [in] The file handle being used in a MIF dump */
    char const *fileName,  /**< [in] Name of the MIF file */
    json_object *part_obj  /**< [in] The json object representing this mesh part */
)
{
    json_object *part_info = json_object_new_object();

#warning SOMEHOW SHOULD INCLUDE OFFSETS TO EACH VARIABLE
    /* Write the json mesh part object as an ascii string */
    fprintf(myFile, "%s\n", json_object_to_json_string_ext(part_obj, JSON_C_TO_STRING_PRETTY));
    json_object_free_printbuf(part_obj);

    /* Form the return 'value' holding the information on where to find this part */
    json_object_object_add(part_info, "partid",
#warning CHANGE NAME OF KEY IN JSON TO PartID
        json_object_new_int(json_object_path_get_int(part_obj, "Mesh/ChunkID")));
    json_object_object_add(part_info, "file",
        json_object_new_string(fileName));
    json_object_object_add(part_info, "offset",
        json_object_new_double((double) ftello(myFile)));

    return part_info;
}

/*!
\brief Main MIF dump implementation for this plugin

This is the function MACSio main calls to do the actual dump of data with this plugin.

It uses \ref MACSIO_MIF twice; once for the main dump and a second time to create the
root (or master) file. However, in the second use, the file count is set to 1. That
means that the root file is effectively written using serial (e.g. non-parallel) I/O.

It is a useful exercise to ask how we might improve the implementation here to avoid
writing the root file using serial I/O.
*/
static void main_dump(
    int argi,               /**< [in] Command-line argument index at which first plugin-specific arg appears */
    int argc,               /**< [in] argc from main */
    char **argv,            /**< [in] argv from main */
    json_object *main_obj,  /**< [in] The main json object representing all data to be dumped */
    int dumpn,              /**< [in] The number/index of this dump. Each dump in a sequence gets a unique,
                                      monotone increasing index starting from 0 */
    double dumpt            /**< [in] The time to be associated with this dump (like a simulation's time) */
)
{
    int i, rank, numFiles;
    char fileName[256];
    FILE *myFile;
    MACSIO_MIF_ioFlags_t ioFlags = {MACSIO_MIF_WRITE, JsonGetInt(main_obj, "clargs/exercise_scr")&0x1};
    MACSIO_MIF_baton_t *bat;
    json_object *parts;
    json_object *part_infos = json_object_new_array();

    /* process cl args */
    process_args(argi, argc, argv);

    /* ensure we're in MIF mode and determine the file count */
#warning SIMPLIFY THIS LOGIC USING NEW JSON INTERFACE
    json_object *parfmode_obj = json_object_path_get_array(main_obj, "clargs/parallel_file_mode");
    if (parfmode_obj)
    {
        json_object *modestr = json_object_array_get_idx(parfmode_obj, 0);
        json_object *filecnt = json_object_array_get_idx(parfmode_obj, 1);
        if (!strcmp(json_object_get_string(modestr), "SIF"))
        {
            MACSIO_LOG_MSG(Die, ("miftmpl plugin cannot currently handle SIF mode"));
        }
        else
        {
            numFiles = json_object_get_int(filecnt);
        }
    }
    else
    {
        char const * modestr = json_object_path_get_string(main_obj, "clargs/parallel_file_mode");
        if (!strcmp(modestr, "SIF"))
        {
            MACSIO_LOG_MSG(Die, ("miftmpl plugin cannot currently handle SIF mode"));
        }
        else if (!strcmp(modestr, "MIFMAX"))
            numFiles = json_object_path_get_int(main_obj, "parallel/mpi_size");
        else if (!strcmp(modestr, "MIFAUTO"))
        {
            /* Call MACSio utility to determine optimal file count */
        }
    }

    bat = MACSIO_MIF_Init(numFiles, ioFlags, MACSIO_MAIN_Comm, 3,
        CreateMyFile, OpenMyFile, CloseMyFile, 0);

    rank = json_object_path_get_int(main_obj, "parallel/mpi_rank");

    /* Construct name for the silo file */
    sprintf(fileName, "%s_json_%05d_%03d.%s",
        json_object_path_get_string(main_obj, "clargs/filebase"),
        MACSIO_MIF_RankOfGroup(bat, rank),
        dumpn,
        json_object_path_get_string(main_obj, "clargs/fileext"));

    myFile = (FILE *) MACSIO_MIF_WaitForBaton(bat, fileName, 0);

    parts = json_object_path_get_array(main_obj, "problem/parts");
    for (int i = 0; i < json_object_array_length(parts); i++)
    {
        json_object *this_part = json_object_array_get_idx(parts, i);
        json_object_array_add(part_infos, write_mesh_part(myFile, fileName, this_part));
    }

    /* Hand off the baton to the next processor. This winds up closing
     * the file so that the next processor that opens it can be assured
     * of getting a consistent and up to date view of the file's contents. */
    MACSIO_MIF_HandOffBaton(bat, myFile);

    /* We're done using MACSIO_MIF for these files, so finish it off */
    MACSIO_MIF_Finish(bat);

    /* Use MACSIO_MIF a second time to manage writing of the master/root
       file contents. This winds up being serial I/O but also means we
       never collect all info on all parts to any single processor. */
#warning THERE IS A BETTER WAY TO DO USING LOOP OF NON-BLOCKING RECIEVES
    bat = MACSIO_MIF_Init(1, ioFlags, MACSIO_MAIN_Comm, 5,
        CreateMyFile, OpenMyFile, CloseMyFile, 0);

    /* Construct name for the silo file */
    sprintf(fileName, "%s_json_root_%03d.%s",
        json_object_path_get_string(main_obj, "clargs/filebase"),
        dumpn,
        json_object_path_get_string(main_obj, "clargs/fileext"));

    /* Wait for MACSIO_MIF to give this processor exclusive access */
    myFile = (FILE *) MACSIO_MIF_WaitForBaton(bat, fileName, 0);

#warning FIX THE STRING THAT WE PRODUCE HERE SO ITS A SINGLE JSON ARRAY OBJECT
    /* This processor's work on the file is just to write its part_infos */
    fprintf(myFile, "%s\n", json_object_to_json_string_ext(part_infos, JSON_C_TO_STRING_PRETTY));

    MACSIO_MIF_HandOffBaton(bat, myFile);

    MACSIO_MIF_Finish(bat);

    /* decriment ref-count (and free) part_infos */
    json_object_put(part_infos);
}

/*!
\brief Method to register this plugin with MACSio main

Due to its use to initialize a file-scope, static const variable, this
function winds up being called at load time (e.g. before main is even called).

Its purpose is to add key information about this plugin to MACSio's global
interface table.
*/
static int register_this_interface()
{
    MACSIO_IFACE_Handle_t iface;

    if (strlen(iface_name) >= MACSIO_IFACE_MAX_NAME)
        MACSIO_LOG_MSG(Die, ("Interface name \"%s\" too long", iface_name));

    /* Populate information about this plugin */
    strcpy(iface.name, iface_name);
    strcpy(iface.ext, iface_ext);
    iface.dumpFunc = main_dump;
    iface.processArgsFunc = process_args;

    /* Register this plugin */
    if (!MACSIO_IFACE_Register(&iface))
        MACSIO_LOG_MSG(Die, ("Failed to register interface \"%s\"", iface_name));

    return 0;
}

/*!
\brief Dummy initializer to trigger register_this_interface by the loader

This one statement is the only statement requiring compilation by
a C++ compiler. That is because it involves initialization and non
constant expressions (a function call in this case). This function
call is guaranteed to occur during *initialization* (that is before
even 'main' is called) and so will have the effect of populating the
iface_map array merely by virtue of the fact that this code is linked
with a main.
*/
static int const dummy = register_this_interface();

/*!@}*/

/*!@}*/


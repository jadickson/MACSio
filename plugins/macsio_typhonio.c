#include <macsio_clargs.h>
#include <macsio_iface.h>
#include <macsio_log.h>
#include <macsio_main.h>
#include <macsio_mif.h>
#include <macsio_utils.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <hdf5.h>
#include <typhonio.h>


/* Disable debugging messages */

/*!
\addtogroup plugins
@{
*/

/*!
\addtogroup TyphonIO
@{
*/

	static char const *iface_name = "typhonio";
	static char const *iface_ext = "h5";
	static char *filename;

	char  errstr[TIO_STRLEN];
	TIO_t errnum;

	static int show_errors = 0;

	#define TIO_Call(r,s) if((errnum=r) != TIO_SUCCESS) {TIO_Get_Error(errnum, errstr); printf("%s\n%s\n",s,errstr); exit(EXIT_FAILURE);}



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
       	MACSIO_CLARGS_END_OF_ARGS);

       return 0;
   }

   static char *getDate()
   {
   	time_t     now;
   	char       date[TIO_STRLEN];
   	struct tm  *ts;

   	time(&now);
   	ts = localtime(&now);
   	strftime(date, sizeof(date), "%a %d-%m-%Y %H:%M", ts);
   	return date;
   }


/*!
\brief CreateFile MIF Callback

This implments the MACSIO_MIF_CreateFile callback needed for a MIF mode plugin.

\return A void pointer to the plugin-specific file handle
*/

   static void *CreateTyphonIOFile(
    const char *fname,     /**< [in] Name of the MIF file to create */
    const char *nsname,    /**< [in] Name of the namespace within the file for caller should use. */
    void *userData         /**< [in] Optional plugin-specific user-defined data */
   	)
   {
    TIO_File_t *retval = 0;
   	TIO_File_t file_id;
   	char *date = getDate();
   	TIO_Call( TIO_Create(fname, &file_id, TIO_ACC_REPLACE, "MACSio", 
   		"0.9", date, fname, MACSIO_MAIN_Comm, MPI_INFO_NULL, MACSIO_MAIN_Rank),
   	"File Creation Failed\n");
    if (file_id >= 0)
    {
        retval = (TIO_File_t *) malloc(sizeof(TIO_File_t));
        *retval = file_id;
    }

   	return (void *) retval;
   }

   static void *OpenTyphonIOFile(
   	const char *fname,
   	const char *nsname,
   	MACSIO_MIF_ioFlags_t ioFlags, 
   	void *userData
   	)
   {
    TIO_File_t *retval = 0;
   	TIO_File_t file_id;

   	char *date = getDate();
   	TIO_Call( TIO_Open(fname, &file_id, TIO_ACC_READWRITE, "MACSio", 
   		"0.9", date, (char*)fname, MACSIO_MAIN_Comm, MPI_INFO_NULL, MACSIO_MAIN_Rank),
   	"File Open Failed\n");
    if (file_id >= 0)
    {
        retval = (TIO_File_t *) malloc(sizeof(TIO_File_t));
        *retval = file_id;
    }
   	return (void *) retval;
   }

   static void CloseTyphonIOFile(
   	void *file,
   	void *userData
   	)
   {
   	TIO_Call( TIO_Close(*(TIO_File_t*)file), 
   		"File Close Failed\n");
   }

   static void write_mesh_part(
   	TIO_File_t file_id,
   	json_object *part_obj
   	)
   {
		TIO_Object_t state_id, variable_id;
		char *state_name = "state0";
		char var_namr[TIO_STRLEN];
		json_object *vars_array = json_object_path_get_array(part_obj, "Vars");

		TIO_Call( TIO_Create_State(file_id, state_name, &state_id, 1, (TIO_Time_t)0.0, "us"),
			"State Create Failed\n");

		for (int i = 0; i < json_object_array_length(vars_array); i++)
		{
			int j;
			TIO_Size_t var_dims[3];
			TIO_Object_t fspace_id, ds_id, var_id;
			json_object *var_obj = json_object_array_get_idx(vars_array, i);
			json_object *data_obj = json_object_path_get_extarr(var_obj, "data");
			char const *varname = json_object_path_get_string(var_obj, "name");
			int ndims = json_object_extarr_ndims(data_obj);
			void const *buf = json_object_extarr_data(data_obj);

			TIO_Dims_t ndims_tio = (TIO_Dims_t)ndims;

			TIO_Data_t dtype_id = json_object_extarr_type(data_obj)==json_extarr_type_flt64? 
			TIO_DOUBLE:TIO_INT;

			for (j = 0; j < ndims; j++)
				var_dims[j] = json_object_extarr_dim(data_obj, j);

			TIO_Call( TIO_Create_Variable(file_id, state_id, varname, &var_id, TIO_DOUBLE, ndims_tio, var_dims, NULL),
				"Create variable failed\n");
			TIO_Call( TIO_Write_Variable(file_id, var_id, TIO_DOUBLE, buf),
				"Write variable failed\n");

			TIO_Call( TIO_Close_Variable(file_id, var_id),
				"Close variable failed\n");
        }

   }

   typedef struct _user_data {
   	TIO_t groupId;
   } user_data_t;

   static void main_dump_mif(json_object *main_obj, int numFiles, int dumpn, double dumpt)
   {
		int size, rank;
		TIO_t *tioFile_ptr;
		TIO_File_t tioFile;
		TIO_Object_t tioGroup;
		char fileName[256];
		int i, len;
		int *theData;
		user_data_t userData;
		MACSIO_MIF_ioFlags_t ioFlags = {MACSIO_MIF_WRITE,
			JsonGetInt(main_obj, "clargs/exercise_scr")&0x1};

#warning SET FILE AND DATASET PROPERTIES
#warning DIFFERENT MPI TAGS FOR DIFFERENT PLUGINS AND CONTEXTS
			MACSIO_MIF_baton_t *bat = MACSIO_MIF_Init(numFiles, ioFlags, MACSIO_MAIN_Comm, 3,
				CreateTyphonIOFile, OpenTyphonIOFile, CloseTyphonIOFile, &userData);

			rank = json_object_path_get_int(main_obj, "parallel/mpi_rank");
			size = json_object_path_get_int(main_obj, "parallel/mpi_size");

    /* Construct name for the silo file */
			sprintf(fileName, "%s_typhonio_%05d_%03d.%s",
				json_object_path_get_string(main_obj, "clargs/filebase"),
				MACSIO_MIF_RankOfGroup(bat, rank),
				dumpn,
				json_object_path_get_string(main_obj, "clargs/fileext"));

			tioFile_ptr = (TIO_t *) MACSIO_MIF_WaitForBaton(bat, fileName, 0);
			tioFile = *tioFile_ptr;
			tioGroup = userData.groupId;

			json_object *parts = json_object_path_get_array(main_obj, "problem/parts");

			for (int i = 0; i < json_object_array_length(parts); i++)
			{
				char domain_dir[256];
				json_object *this_part = json_object_array_get_idx(parts, i);
				TIO_t domain_group_id;

				snprintf(domain_dir, sizeof(domain_dir), "domain_%07d",
					json_object_path_get_int(this_part, "Mesh/ChunkID"));

				//domain_group_id = H5Gcreate1(h5File, domain_dir, 0);

				write_mesh_part(tioFile, this_part);

				//H5Gclose(domain_group_id);
			}

    /* Hand off the baton to the next processor. This winds up closing
     * the file so that the next processor that opens it can be assured
     * of getting a consistent and up to date view of the file's contents. */
     MACSIO_MIF_HandOffBaton(bat, tioFile_ptr);

    /* We're done using MACSIO_MIF, so finish it off */
     MACSIO_MIF_Finish(bat);

   }

   static void main_dump_sif(json_object *main_obj, int dumpn, double dumpt)
   {
// #ifdef HAVE_MPI
//  	int ndims;
//  	int i, v, p;
//  	char const *mesh_type = json_object_path_get_string(main_obj, "clargs/part_type");
//  	char fileName[256];
//  	int use_part_count;

//  	hid_t h5file_id;
//  	hid_t fapl_id = make_fapl();
//  	hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
//  	hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
//  	hid_t null_space_id = H5Screate(H5S_NULL);
//  	hid_t fspace_nodal_id, fspace_zonal_id;
//  	hsize_t global_log_dims_nodal[3];
//  	hsize_t global_log_dims_zonal[3];

//  	MPI_Info mpiInfo = MPI_INFO_NULL;

// #warning WE ARE DOING SIF SLIGHTLY WRONG, DUPLICATING SHARED NODES
// #warning INCLUDE ARGS FOR ISTORE AND K_SYM
// #warning INCLUDE ARG PROCESS FOR HINTS
// #warning FAPL PROPS: ALIGNMENT 

// #warning FOR MIF, NEED A FILEROOT ARGUMENT OR CHANGE TO FILEFMT ARGUMENT
//     /* Construct name for the HDF5 file */
//  	sprintf(fileName, "%s_typhonio_%03d.%s",
//  		json_object_path_get_string(main_obj, "clargs/filebase"),
//  		dumpn,
//  		json_object_path_get_string(main_obj, "clargs/fileext"));

//  	h5file_id = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);

//     /* Create an HDF5 Dataspace for the global whole of mesh and var objects in the file. */
//  	ndims = json_object_path_get_int(main_obj, "clargs/part_dim");
//  	json_object *global_log_dims_array =
//  	json_object_path_get_array(main_obj, "problem/global/LogDims");
//  	json_object *global_parts_log_dims_array =
//  	json_object_path_get_array(main_obj, "problem/global/PartsLogDims");
//     /* Note that global zonal array is smaller in each dimension by one *ON*EACH*BLOCK*
//        in the associated dimension. */
//  	for (i = 0; i < ndims; i++)
//  	{
//  		int parts_log_dims_val = JsonGetInt(global_parts_log_dims_array, "", i);
//  		global_log_dims_nodal[ndims-1-i] = (hsize_t) JsonGetInt(global_log_dims_array, "", i);
//  		global_log_dims_zonal[ndims-1-i] = global_log_dims_nodal[ndims-1-i] -
//  		JsonGetInt(global_parts_log_dims_array, "", i);
//  	}
//  	fspace_nodal_id = H5Screate_simple(ndims, global_log_dims_nodal, 0);
//  	fspace_zonal_id = H5Screate_simple(ndims, global_log_dims_zonal, 0);

//     /* Get the list of vars on the first part as a guide to loop over vars */
//  	json_object *part_array = json_object_path_get_array(main_obj, "problem/parts");
//  	json_object *first_part_obj = json_object_array_get_idx(part_array, 0);
//  	json_object *first_part_vars_array = json_object_path_get_array(first_part_obj, "Vars");

//     /* Dataset transfer property list used in all H5Dwrite calls */
// #if H5_HAVE_PARALLEL
//  	if (no_collective)
//  		H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
//  	else
//  		H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
// #endif


//     /* Loop over vars and then over parts */
//     /* currently assumes all vars exist on all ranks. but not all parts */
//     for (v = -1; v < json_object_array_length(first_part_vars_array); v++) /* -1 start is for Mesh */
//  	{

// #warning SKIPPING MESH
//         if (v == -1) continue; /* All ranks skip mesh (coords) for now */

//         /* Inspect the first part's var object for name, datatype, etc. */
//  		json_object *var_obj = json_object_array_get_idx(first_part_vars_array, v);
//  		char const *varName = json_object_path_get_string(var_obj, "name");
//  		char *centering = strdup(json_object_path_get_string(var_obj, "centering"));
//  		json_object *dataobj = json_object_path_get_extarr(var_obj, "data");
// #warning JUST ASSUMING TWO TYPES NOW. CHANGE TO A FUNCTION
//  		hid_t dtype_id = json_object_extarr_type(dataobj)==json_extarr_type_flt64? 
//  		H5T_NATIVE_DOUBLE:H5T_NATIVE_INT;
//  		hid_t fspace_id = H5Scopy(strcmp(centering, "zone") ? fspace_nodal_id : fspace_zonal_id);
//  		hid_t dcpl_id = make_dcpl(compression_alg_str, compression_params_str, fspace_id, dtype_id);

//         /* Create the file dataset (using old-style H5Dcreate API here) */
// #warning USING DEFAULT DCPL: LATER ADD COMPRESSION, ETC.

//  		hid_t ds_id = H5Dcreate1(h5file_id, varName, dtype_id, fspace_id, dcpl_id); 
//  		H5Sclose(fspace_id);
//  		H5Pclose(dcpl_id);

//         /* Loop to make write calls for this var for each part on this rank */
// #warning USE NEW MULTI-DATASET API WHEN AVAILABLE TO AGLOMERATE ALL PARTS INTO ONE CALL
//  		use_part_count = (int) ceil(json_object_path_get_double(main_obj, "clargs/avg_num_parts"));
//  		for (p = 0; p < use_part_count; p++)
//  		{
//  			json_object *part_obj = json_object_array_get_idx(part_array, p);
//  			json_object *var_obj = 0;
//  			hid_t mspace_id = H5Scopy(null_space_id);
//  			void const *buf = 0;

//  			fspace_id = H5Scopy(null_space_id);

//             /* this rank actually has something to contribute to the H5Dwrite call */
//  			if (part_obj)
//  			{
//  				int i;
//  				hsize_t starts[3], counts[3];
//  				json_object *vars_array = json_object_path_get_array(part_obj, "Vars");
//  				json_object *mesh_obj = json_object_path_get_object(part_obj, "Mesh");
//  				json_object *var_obj = json_object_array_get_idx(vars_array, v);
//  				json_object *extarr_obj = json_object_path_get_extarr(var_obj, "data");
//  				json_object *global_log_origin_array =
//  				json_object_path_get_array(part_obj, "GlobalLogOrigin");
//  				json_object *global_log_indices_array =
//  				json_object_path_get_array(part_obj, "GlobalLogIndices");
//  				json_object *mesh_dims_array = json_object_path_get_array(mesh_obj, "LogDims");
//  				for (i = 0; i < ndims; i++)
//  				{
//  					starts[ndims-1-i] =
//  					json_object_get_int(json_object_array_get_idx(global_log_origin_array,i));
//  					counts[ndims-1-i] =
//  					json_object_get_int(json_object_array_get_idx(mesh_dims_array,i));
//  					if (!strcmp(centering, "zone"))
//  					{
//  						counts[ndims-1-i]--;
//  						starts[ndims-1-i] -=
//  						json_object_get_int(json_object_array_get_idx(global_log_indices_array,i));
//  					}
//  				}

//                 /* set selection of filespace */
//  				fspace_id = H5Dget_space(ds_id);
//  				H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, starts, 0, counts, 0);

//                 /* set dataspace of data in memory */
//  				mspace_id = H5Screate_simple(ndims, counts, 0);
//  				buf = json_object_extarr_data(extarr_obj);
//  			}

//  			H5Dwrite(ds_id, dtype_id, mspace_id, fspace_id, dxpl_id, buf);
//  			H5Sclose(fspace_id);
//  			H5Sclose(mspace_id);

//  		}

//  		H5Dclose(ds_id);
//  		free(centering);
//  	}

//  	H5Sclose(fspace_nodal_id);
//  	H5Sclose(fspace_zonal_id);
//  	H5Sclose(null_space_id);
//  	H5Pclose(dxpl_id);
//  	H5Pclose(fapl_id);
//  	H5Fclose(h5file_id);

// #endif
   }

   static void main_dump(int argi, int argc, char **argv, json_object *main_obj,
   	int dumpn, double dumpt)
   {
 	int rank, size, numFiles;

#warning SET ERROR MODE OF HDF5 LIBRARY

    /* Without this barrier, I get strange behavior with Silo's MACSIO_MIF interface */
 	mpi_errno = MPI_Barrier(MACSIO_MAIN_Comm);

    /* process cl args */
 	process_args(argi, argc, argv);

 	rank = json_object_path_get_int(main_obj, "parallel/mpi_rank");
 	size = json_object_path_get_int(main_obj, "parallel/mpi_size");

#warning MOVE TO A FUNCTION
    /* ensure we're in MIF mode and determine the file count */
 	json_object *parfmode_obj = json_object_path_get_array(main_obj, "clargs/parallel_file_mode");
 	if (parfmode_obj)
 	{
 		json_object *modestr = json_object_array_get_idx(parfmode_obj, 0);
 		json_object *filecnt = json_object_array_get_idx(parfmode_obj, 1);
#warning ERRORS NEED TO GO TO LOG FILES AND ERROR BEHAVIOR NEEDS TO BE HONORED
 		if (!strcmp(json_object_get_string(modestr), "SIF"))
 		{
 			main_dump_sif(main_obj, dumpn, dumpt);
 		}
 		else
 		{
 			numFiles = json_object_get_int(filecnt);
 			main_dump_mif(main_obj, numFiles, dumpn, dumpt);
 		}
 	}
 	else
 	{
 		char const * modestr = json_object_path_get_string(main_obj, "clargs/parallel_file_mode");
 		if (!strcmp(modestr, "SIF"))
 		{
 			float avg_num_parts = json_object_path_get_double(main_obj, "clargs/avg_num_parts");
 			if (avg_num_parts == (float ((int) avg_num_parts)))
 				main_dump_sif(main_obj, dumpn, dumpt);
 			else
 			{
#warning CURRENTLY, SIF CAN WORK ONLY ON WHOLE PART COUNTS
 				MACSIO_LOG_MSG(Die, ("HDF5 plugin cannot currently handle SIF mode where "
 					"there are different numbers of parts on each MPI rank. "
 					"Set --avg_num_parts to an integral value." ));
 			}
 		}
 		else if (!strcmp(modestr, "MIFMAX"))
 			numFiles = json_object_path_get_int(main_obj, "parallel/mpi_size");
 		else if (!strcmp(modestr, "MIFAUTO"))
 		{
            /* Call utility to determine optimal file count */
#warning ADD UTILIT TO DETERMINE OPTIMAL FILE COUNT
 		}
 		main_dump_mif(main_obj, numFiles, dumpn, dumpt);
 	}
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


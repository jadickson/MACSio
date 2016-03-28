#include <math.h>
#include <time.h>

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

static int no_collective = 0;

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
		"--no_collective", "",
			"Use independent, not collective I/O calls in SIF mode.",
			&no_collective,
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
	char *ret = (char*) malloc(20*sizeof(char));
	strcpy(ret, date);
	return ret;
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
	//MPI_Comm *groupComm = (MPI_Comm*)userData;
	char *date = getDate();
	TIO_Call( TIO_Create(fname, &file_id, TIO_ACC_REPLACE, "MACSio",
	                     "0.9", date, (char*)fname, MPI_COMM_SELF, MPI_INFO_NULL, MACSIO_MAIN_Rank),
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
	                   "0.9", date, (char*)fname, MPI_COMM_SELF, MPI_INFO_NULL, MACSIO_MAIN_Rank),
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

static void write_quad_mesh_part(
	TIO_File_t file_id,
	TIO_Object_t state_id,
	json_object *part_obj,
	TIO_Mesh_t tio_mesh_type
)
{
	TIO_Object_t mesh_id;
	json_object *coordobj;
	void const *coords[3];
	int ndims = JsonGetInt(part_obj, "Mesh/GeomDim");
	int dims[3] = {1,1,1};
	int dimsz[3] = {1,1,1};

	dims[0] = JsonGetInt(part_obj, "Mesh/LogDims", 0);

	if (tio_mesh_type == TIO_MESH_QUAD_COLINEAR)
	{
		coordobj = JsonGetObj(part_obj, "Mesh/Coords/XAxisCoords"); // Rect Mesh
	}
	else
	{
		coordobj = JsonGetObj(part_obj, "Mesh/Coords/XCoords");	// Curv Mesh
		/* For non-colinear mesh, TyphonIO misses out the end cell in each dimension so increment by one */
		//dims[0]++;
	}

	coords[0] = json_object_extarr_data(coordobj);


	if (ndims > 1)
	{
		dims[1] = JsonGetInt(part_obj, "Mesh/LogDims", 1);

		if (tio_mesh_type == TIO_MESH_QUAD_COLINEAR)
		{
			coordobj = JsonGetObj(part_obj, "Mesh/Coords/YAxisCoords");
		} 
		else
		{
			coordobj = JsonGetObj(part_obj, "Mesh/Coords/YCoords");
			/* When creating a non-colinear mesh, TyphonIO takes the dimension as an index and missed out the end cell in each dimension so it needs to be incremented */
			//dims[1]++;
		}
		coords[1] = json_object_extarr_data(coordobj);
	}
	if (ndims > 2)
	{
		dims[2] = JsonGetInt(part_obj, "Mesh/LogDims", 2);

		if (tio_mesh_type == TIO_MESH_QUAD_COLINEAR)
		{
			coordobj = JsonGetObj(part_obj, "Mesh/Coords/ZAxisCoords");
		} 
		else
		{
			coordobj = JsonGetObj(part_obj, "Mesh/Coords/ZCoords");
			//dims[2]++;
		}
		coords[2] = json_object_extarr_data(coordobj);
	}

	TIO_Call( TIO_Create_Mesh(file_id, state_id, "mesh", &mesh_id, tio_mesh_type, 
							TIO_COORD_CARTESIAN, TIO_FALSE, "mesh_group", (TIO_Size_t)1,
							TIO_DATATYPE_NULL, TIO_DOUBLE, (TIO_Dims_t)ndims,
							(TIO_Size_t)dims[0], (TIO_Size_t)dims[1], (TIO_Size_t)dims[2],
							TIO_NULL, (TIO_Size_t)1,
							NULL, NULL, NULL,
			                NULL, NULL, NULL),
						"Create Mesh Failed\n");

	if (tio_mesh_type == TIO_MESH_QUAD_COLINEAR){
		TIO_Call( TIO_Set_Quad_Chunk(file_id, mesh_id, (TIO_Size_t)0, (TIO_Dims_t)ndims,
							0, dims[0], 0, dims[1], 0, dims[2],
							0, 0),
			"Set Quad Mesh Chunk Failed");
		TIO_Call( TIO_Write_QuadMesh_All(file_id, mesh_id, TIO_DOUBLE, coords[0], coords[1], coords[2]),
					"Write Mesh Coords failed\n");
	} 
	else 
	{
		TIO_Call( TIO_Set_Quad_Chunk(file_id, mesh_id, (TIO_Size_t)0, (TIO_Dims_t)ndims,
							0, dims[0]+1, 0, dims[1]+1, 0, dims[2]+1,
							0, 0),
			"Set Quad Mesh Chunk Failed");
		TIO_Call( TIO_Write_QuadMesh_Chunk(file_id, mesh_id, 0, TIO_XFER_INDEPENDENT, 
											TIO_DOUBLE, coords[0], coords[1], coords[2]),
					"Write Non-Colinear Mesh Coords failed\n");
	}

	TIO_Object_t variable_id;
	json_object *vars_array = json_object_path_get_array(part_obj, "Vars");

	for (int i = 0; i < json_object_array_length(vars_array); i++)
	{
		int j;
		TIO_Size_t var_dims[3];
		TIO_Object_t var_id;
		json_object *var_obj = json_object_array_get_idx(vars_array, i);
		json_object *data_obj = json_object_path_get_extarr(var_obj, "data");
		char const *varname = json_object_path_get_string(var_obj, "name");
		int ndims = json_object_extarr_ndims(data_obj);
		void const *buf = json_object_extarr_data(data_obj);

		TIO_Dims_t ndims_tio = (TIO_Dims_t)ndims;

		TIO_Data_t dtype_id = json_object_extarr_type(data_obj) == json_extarr_type_flt64 ?
		                      TIO_DOUBLE : TIO_INT;

		for (j = 0; j < ndims; j++)
			var_dims[j] = json_object_extarr_dim(data_obj, j);

		TIO_Call( TIO_Create_Variable(file_id, mesh_id, varname, &var_id, dtype_id, ndims_tio, var_dims, NULL),
		          "Create variable failed\n");
		TIO_Call( TIO_Write_Variable(file_id, var_id, dtype_id, buf),
		          "Write variable failed\n");

		TIO_Call( TIO_Close_Variable(file_id, var_id),
		          "Close variable failed\n");
	}
	TIO_Call( TIO_Close_Mesh(file_id, mesh_id),
		"Close Mesh failed\n");

}

static void write_ucdzoo_mesh_part(
	TIO_File_t file_id,
	TIO_Object_t state_id,
	json_object *part_obj,
	char const *topo_name
)
{
	TIO_Object_t mesh_id;
	json_object *coordobj, *topoobj;
	char const *coordnames[] = {"X", "Y", "Z"};
	void const *coords[3];
	int ndims = JsonGetInt(part_obj, "Mesh/GeomDim");
	int nnodes = 1, nzones = 1;
	int dims[3] = {1,1,1};
	int dimsz[3] = {1,1,1};

	coordobj = JsonGetObj(part_obj, "Mesh/Coords/XCoords");
    coords[0] = json_object_extarr_data(coordobj);
    dims[0] = JsonGetInt(part_obj, "Mesh/LogDims", 0);
    dimsz[0] = dims[0]-1;
    nnodes *= dims[0];
    nzones *= dimsz[0];
    if (ndims > 1)
    {
        coordobj = JsonGetObj(part_obj, "Mesh/Coords/YCoords");
        coords[1] = json_object_extarr_data(coordobj);
        dims[1] = JsonGetInt(part_obj, "Mesh/LogDims", 1);
        dimsz[1] = dims[1]-1;
        nnodes *= dims[1];
        nzones *= dimsz[1];
    }
    if (ndims > 2)
    {
        coordobj = JsonGetObj(part_obj, "Mesh/Coords/ZCoords");
        coords[2] = json_object_extarr_data(coordobj);
        dims[2] = JsonGetInt(part_obj, "Mesh/LogDims", 2);
        dimsz[2] = dims[2]-1;
        nnodes *= dims[2];
        nzones *= dimsz[2];
    }

    if (ndims == 1 || !strcmp(topo_name, "ucdzoo"))
    /* UCDZOO */
    {
        json_object *topoobj = JsonGetObj(part_obj, "Mesh/Topology");
        json_object *nlobj = JsonGetObj(topoobj, "Nodelist");
        void const *nodelist = (void const*) json_object_extarr_data(nlobj);
        int lnodelist = json_object_extarr_nvals(nlobj);
        TIO_Shape_t shapetype; 
        int shapesize;
        int shapecnt = nzones;
        int ncells = nzones;

        if (!strcmp(JsonGetStr(topoobj, "ElemType"), "Beam2"))
        {
            shapesize = 2;
            shapetype = TIO_SHAPE_BAR2;
        }
        else if (!strcmp(JsonGetStr(topoobj, "ElemType"), "Quad4"))
        {
            shapesize = 4;
            shapetype = TIO_SHAPE_QUAD4;
        }
        else if (!strcmp(JsonGetStr(topoobj, "ElemType"), "Hex8"))
        {
            shapesize = 8;
            shapetype = TIO_SHAPE_HEX8;
        }

 		TIO_Size_t nshapes = 1;
 		TIO_Size_t nconnectivity = 0;//ncells*shapesize;

 		 /* For unstructured: n1=nnodes, n2=ncells, n3=nshapes, n4=nconnectivity */

        TIO_Call( TIO_Create_Mesh(file_id, state_id, "mesh", &mesh_id, TIO_MESH_UNSTRUCT,
        						TIO_COORD_CARTESIAN, TIO_FALSE, "mesh_group", (TIO_Size_t)1,
        						TIO_INT, TIO_DOUBLE, (TIO_Dims_t)ndims, 
        						(TIO_Size_t)nnodes, (TIO_Size_t)ncells, (TIO_Size_t)nshapes,
								TIO_NULL, (TIO_Size_t)1,
								NULL, NULL, NULL,
				                NULL, NULL, NULL),
							"Create Mesh Failed\n");

        TIO_Call( TIO_Set_Unstr_Chunk(file_id, mesh_id, (TIO_Size_t)0, (TIO_Dims_t)ndims, (TIO_Size_t)nnodes,
        						ncells, nshapes, nconnectivity, 0, 0, 0, 0, 0, 0),
        					"Set UCDZOO Mesh Chunk Failed\n");

        TIO_Call( TIO_Write_UnstrMesh_Chunk(file_id, mesh_id, (TIO_Size_t)0, TIO_XFER_INDEPENDENT,
        								TIO_INT, TIO_DOUBLE, nodelist, nodelist, &shapetype, 
        								&ncells, (const void*)NULL, coords[0], coords[1], coords[2]),
        					"Write unstructured Mesh Failed\n");
    }
    else if (!strcmp(topo_name, "arbitrary"))
    /* ARBITRARY */
    {
 		json_object *topoobj = JsonGetObj(part_obj, "Mesh/Topology");
        json_object *nlobj = JsonGetObj(topoobj, "Nodelist");
        void const *nodelist = (void const*) json_object_extarr_data(nlobj);
        json_object *ncobj = JsonGetObj(topoobj, "NodeCounts");
        int const *nodecnt = (int const *) json_object_extarr_data(ncobj);
        int nfaces = json_object_extarr_nvals(ncobj);
        json_object *flobj = JsonGetObj(topoobj, "Facelist");
        int const *facelist = (int const *) json_object_extarr_data(flobj);
        int lfacelist = json_object_extarr_nvals(flobj);
        json_object *fcobj = JsonGetObj(topoobj, "FaceCounts");
        int const *facecnt = (int const *) json_object_extarr_data(fcobj);

        int ncells = nzones;
        TIO_Size_t nshapes = MACSIO_MAIN_Size;
        TIO_Size_t nconnectivity = 0;//ncells*shapesize;
        TIO_Shape_t *shapetype = (TIO_Shape_t*)malloc(nshapes * sizeof(TIO_Shape_t));
        for (int s = 0; s<nshapes; s++){
        	shapetype[s] = (TIO_Shape_t)4; // Abribtrary polygon with 4 nodes
        }  

    	TIO_Call( TIO_Create_Mesh(file_id, state_id, "mesh", &mesh_id, TIO_MESH_UNSTRUCT,
    							TIO_COORD_CARTESIAN, TIO_FALSE, "mesh_group", (TIO_Size_t)1,
    							TIO_INT, TIO_DOUBLE, (TIO_Dims_t)ndims,
    							(TIO_Size_t)nnodes, (TIO_Size_t)ncells, (TIO_Size_t)nshapes,
    							TIO_NULL, (TIO_Size_t)1,
    							NULL, NULL, NULL,
    							NULL, NULL, NULL),
    						"Create Arbitrary Unstructured Mesh Failed\n");

    	TIO_Call( TIO_Set_Unstr_Chunk(file_id, mesh_id, (TIO_Size_t)0, (TIO_Dims_t)ndims, (TIO_Size_t)nnodes,
    							ncells, nshapes, nconnectivity, 0, 0, 0, 0, 0, 0),
    						"Set Arbitrary Mesh Chunk Failed");

    	 TIO_Call( TIO_Write_UnstrMesh_Chunk(file_id, mesh_id, (TIO_Size_t)0, TIO_XFER_INDEPENDENT,
        								TIO_INT, TIO_DOUBLE, nodelist, nodelist, shapetype, 
        								&ncells, (const void*)NULL, coords[0], coords[1], coords[2]),
        					"Write unstructured Mesh Failed\n");
    	 free(shapetype);
    }	

	TIO_Object_t variable_id;
	json_object *vars_array = json_object_path_get_array(part_obj, "Vars");

	for (int i = 0; i < json_object_array_length(vars_array); i++)
	{
		int j;
		TIO_Size_t var_dims[3];
		TIO_Object_t var_id;
		json_object *var_obj = json_object_array_get_idx(vars_array, i);
		json_object *data_obj = json_object_path_get_extarr(var_obj, "data");
		char const *varname = json_object_path_get_string(var_obj, "name");
		int ndims = json_object_extarr_ndims(data_obj);
		void const *buf = json_object_extarr_data(data_obj);

		TIO_Data_t dtype_id = json_object_extarr_type(data_obj) == json_extarr_type_flt64 ?
		                      TIO_DOUBLE : TIO_INT;

		for (j = 0; j < ndims; j++)
			var_dims[j] = json_object_extarr_dim(data_obj, j);

		TIO_Call( TIO_Create_Variable(file_id, mesh_id, varname, &var_id, dtype_id, (TIO_Dims_t)ndims, var_dims, NULL),
		          "Create variable failed\n");
		TIO_Call( TIO_Write_Variable(file_id, var_id, dtype_id, buf),
		          "Write variable failed\n");

		TIO_Call( TIO_Close_Variable(file_id, var_id),
		          "Close variable failed\n");
	}	    


    TIO_Call( TIO_Close_Mesh(file_id, mesh_id),
				"Close Mesh Failed\n");
	
}

static void write_mesh_part(
    TIO_File_t file_id,
    TIO_Object_t state_id,
    json_object *part_obj
)
{
    if (!strcmp(JsonGetStr(part_obj, "Mesh/MeshType"), "rectilinear"))
        write_quad_mesh_part(file_id, state_id, part_obj, TIO_MESH_QUAD_COLINEAR);
    else if (!strcmp(JsonGetStr(part_obj, "Mesh/MeshType"), "curvilinear"))
        write_quad_mesh_part(file_id, state_id, part_obj, TIO_MESH_QUAD_NONCOLINEAR);
    else if (!strcmp(JsonGetStr(part_obj, "Mesh/MeshType"), "ucdzoo"))
        write_ucdzoo_mesh_part(file_id, state_id, part_obj, "ucdzoo");
    else if (!strcmp(JsonGetStr(part_obj, "Mesh/MeshType"), "arbitrary"))
        write_ucdzoo_mesh_part(file_id, state_id, part_obj, "arbitrary");
}

typedef struct _user_data {
	TIO_t groupId;
	MPI_Comm groupComm;
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
	                                JsonGetInt(main_obj, "clargs/exercise_scr") & 0x1
	                               };

#warning SET FILE AND DATASET PROPERTIES
#warning DIFFERENT MPI TAGS FOR DIFFERENT PLUGINS AND CONTEXTS
	MACSIO_MIF_baton_t *bat = MACSIO_MIF_Init(numFiles, ioFlags, MACSIO_MAIN_Comm, 3,
	                          CreateTyphonIOFile, OpenTyphonIOFile, CloseTyphonIOFile, &userData);

	rank = json_object_path_get_int(main_obj, "parallel/mpi_rank");
	size = json_object_path_get_int(main_obj, "parallel/mpi_size");

	/* Construct name for the typhonio file */
	sprintf(fileName, "%s_typhonio_%05d_%03d.%s",
	        json_object_path_get_string(main_obj, "clargs/filebase"),
	        MACSIO_MIF_RankOfGroup(bat, rank),
	        dumpn,
	        "h5");//json_object_path_get_string(main_obj, "clargs/fileext"));

	tioFile_ptr = (TIO_t *) MACSIO_MIF_WaitForBaton(bat, fileName, 0);
	tioFile = *tioFile_ptr;
	tioGroup = userData.groupId;

	json_object *parts = json_object_path_get_array(main_obj, "problem/parts");

	for (int i = 0; i < json_object_array_length(parts); i++)
	{
		char domain_dir[256];
		json_object *this_part = json_object_array_get_idx(parts, i);
		TIO_Object_t domain_group_id;

		snprintf(domain_dir, sizeof(domain_dir), "domain_%07d",
		         json_object_path_get_int(this_part, "Mesh/ChunkID"));

		TIO_Call( TIO_Create_State(tioFile, domain_dir, &domain_group_id, 1, (TIO_Time_t)0.0, "us"),
		          "State Create Failed\n");

		write_mesh_part(tioFile, domain_group_id, this_part);

		TIO_Call( TIO_Close_State(tioFile, domain_group_id),
		          "State Close Failed\n");
	}

	/* Hand off the baton to the next processor. This winds up closing
	 * the file so that the next processor that opens it can be assured
	 * of getting a consistent and up to date view of the file's contents. */
	MACSIO_MIF_HandOffBaton(bat, tioFile_ptr);

	/* We're done using MACSIO_MIF, so finish it off */
	MACSIO_MIF_Finish(bat);

}

static void write_quad_mesh_whole(
	TIO_File_t file_id, 
	TIO_Object_t state_id,
	json_object *main_obj,
	TIO_Mesh_t mesh_type)
{
	TIO_Object_t mesh_id;
	TIO_Object_t object_id;
	int ndims;
	int i, v, p;
	int use_part_count;
	const TIO_Xfer_t TIO_XFER = no_collective ? TIO_XFER_INDEPENDENT : TIO_XFER_COLLECTIVE;

	TIO_Size_t dims[3];
	TIO_Size_t dimsz[3];

	ndims = json_object_path_get_int(main_obj, "clargs/part_dim");
	json_object *global_log_dims_array = json_object_path_get_array(main_obj, "problem/global/LogDims");
	json_object *global_parts_log_dims_array = json_object_path_get_array(main_obj, "problem/global/PartsLogDims");
	

	for (i = 0; i < ndims; i++)
	{
		dims[i] = (TIO_Size_t) JsonGetInt(global_log_dims_array, "", i);
		dimsz[i] = dims[i] - JsonGetInt(global_parts_log_dims_array, "", i);
	}

	/* Get the list of vars on the first part as a guide to loop over vars */
	json_object *part_array = json_object_path_get_array(main_obj, "problem/parts");
	json_object *first_part_obj = json_object_array_get_idx(part_array, 0);
	json_object *first_part_vars_array = json_object_path_get_array(first_part_obj, "Vars");


	/* Loop over vars and then over parts */
	/* currently assumes all vars exist on all ranks. but not all parts */
	for (v = -1; v < json_object_array_length(first_part_vars_array); v++) /* -1 start is for Mesh */
	{
		/* Inspect the first part's var object for name, datatype, etc. */
		json_object *var_obj = json_object_array_get_idx(first_part_vars_array, v);

		char *centering = (char*)malloc(sizeof(char)*5);
		TIO_Data_t dtype_id;

		//MESH
		TIO_Dims_t ndims_tio = (TIO_Dims_t)ndims;

		if (v == -1) { 
			TIO_Call( TIO_Create_Mesh(file_id, state_id, "mesh", &mesh_id, mesh_type,
			                          TIO_COORD_CARTESIAN, TIO_FALSE, "mesh_group", (TIO_Size_t)1,
			                          TIO_DATATYPE_NULL, TIO_DOUBLE, ndims_tio,
			                          dims[0], dims[1], dims[2],
			                          TIO_GHOSTS_NONE, (TIO_Size_t)MACSIO_MAIN_Size,
			                          NULL, NULL, NULL,
			                          NULL, NULL, NULL),

			          "Mesh Create Failed\n");
		} else{
			char const *varName = json_object_path_get_string(var_obj, "name");
			centering = strdup(json_object_path_get_string(var_obj, "centering"));
			json_object *dataobj = json_object_path_get_extarr(var_obj, "data");
			dtype_id = json_object_extarr_type(dataobj) == json_extarr_type_flt64 ? TIO_DOUBLE : TIO_INT;
			TIO_Centre_t tio_centering = strcmp(centering, "zone") ? TIO_CENTRE_NODE : TIO_CENTRE_CELL;

			TIO_Call( TIO_Create_Quant(file_id, mesh_id, varName, &object_id, dtype_id, tio_centering,
										TIO_GHOSTS_NONE, TIO_FALSE, "qunits"),
					"Quant Create Failed\n");
		}

		/* Loop to make write calls for this var for each part on this rank */
		use_part_count = (int) ceil(json_object_path_get_double(main_obj, "clargs/avg_num_parts"));
		for (p = 0; p < use_part_count; p++)
		{
			json_object *part_obj = json_object_array_get_idx(part_array, p);
			json_object *var_obj = 0;

			void const *buf = 0;
			void const *x_coord = 0;
			void const *y_coord = 0;
			void const *z_coord = 0;

			void *x_coord_root = 0;
			void *y_coord_root = 0;
			void *z_coord_root = 0;

			if (part_obj)
			{
				if (v == -1) {
					//Mesh
					TIO_Size_t starts[3], counts[3];
					json_object *mesh_obj = json_object_path_get_object(part_obj, "Mesh");
					json_object *global_log_origin_array = json_object_path_get_array(part_obj, "GlobalLogOrigin");
					json_object *global_log_indices_array = json_object_path_get_array(part_obj, "GlobalLogIndices");
					json_object *mesh_dims_array = json_object_path_get_array(mesh_obj, "LogDims");
					int local_mesh_dims[3];

					
					for (i = 0; i < ndims; i++)
					{
						local_mesh_dims[i] = json_object_get_int(json_object_array_get_idx(mesh_dims_array, i));
						starts[i] = json_object_get_int(json_object_array_get_idx(global_log_origin_array, i));
						counts[i] = json_object_get_int(json_object_array_get_idx(mesh_dims_array, i)) - 1;
					}

					TIO_Size_t local_chunk_indices[6] = {0,0,0,0,0,0};	/* local_chunk_indices [il, ih, jl, jh, kl, kh] */

					// Sets the upper and lower index in each dimension for the current ranks chunk
					local_chunk_indices[0] = starts[0];
					local_chunk_indices[1] = starts[0] + counts[0];

					if (ndims > 1) {
						local_chunk_indices[2] = starts[1];
						local_chunk_indices[3] = starts[1] + counts[1];
						if (ndims > 2) {
							local_chunk_indices[4] = starts[2];
							local_chunk_indices[5] = starts[2] + counts[2];
						}
					}
									
					TIO_Size_t chunk_indices[MACSIO_MAIN_Size][6];
					MPI_Allgather(local_chunk_indices, 6, MPI_DOUBLE, chunk_indices, 6, MPI_DOUBLE, MACSIO_MAIN_Comm);

					for (int k=0; k<MACSIO_MAIN_Size; k++){
						TIO_Call( TIO_Set_Quad_Chunk(file_id, mesh_id, k, ndims_tio,
												chunk_indices[k][0], chunk_indices[k][1],
												chunk_indices[k][2], chunk_indices[k][3],
												chunk_indices[k][4], chunk_indices[k][5],
	                            				 (TIO_Size_t)0, (TIO_Size_t)0),
	         				 "Set Quad Mesh Chunk Failed\n");
					}

					json_object *coords = json_object_path_get_object(mesh_obj, "Coords");

					int coord_array_size[] = {0, 0, 0};

					int block_size = 1;
					switch (ndims) {
						case 3: block_size *= json_object_get_int(json_object_array_get_idx(mesh_dims_array, 2));
						case 2: block_size *= json_object_get_int(json_object_array_get_idx(mesh_dims_array, 1));
						case 1: block_size *= json_object_get_int(json_object_array_get_idx(mesh_dims_array, 0));
					}
					
					if (mesh_type == TIO_MESH_QUAD_COLINEAR){
						x_coord = json_object_extarr_data(json_object_path_get_extarr(coords, "XAxisCoords"));
						y_coord = json_object_extarr_data(json_object_path_get_extarr(coords, "YAxisCoords"));
						z_coord = json_object_extarr_data(json_object_path_get_extarr(coords, "ZAxisCoords"));	
						coord_array_size[0] = dims[0];
						coord_array_size[1] = dims[1];
						coord_array_size[2] = dims[2];	

						int parts = json_object_path_get_int(main_obj, "problem/global/TotalParts");;

						if (MACSIO_MAIN_Rank == 0){
							x_coord_root = malloc(coord_array_size[0] * sizeof(double));
							if (ndims > 1) y_coord_root = malloc(coord_array_size[1] * sizeof(double));
							if (ndims > 2) z_coord_root = malloc(coord_array_size[2] * sizeof(double));
						}

						json_object *bounds = json_object_path_get_array(mesh_obj, "Bounds");

						int color = (JsonGetInt(bounds,"",1) == 0 && JsonGetInt(bounds,"",2) ==0) ? 1: MPI_UNDEFINED;
						MPI_Comm comm;
						MPI_Comm_split(MACSIO_MAIN_Comm, color, MACSIO_MAIN_Rank, &comm);
						if (color == 1){
						MPI_Gather(x_coord, local_mesh_dims[0], MPI_DOUBLE, x_coord_root, local_mesh_dims[0], MPI_DOUBLE, 0, comm);
						}
						if (ndims > 1){
							color = (JsonGetInt(bounds, "", 0)==0 && JsonGetInt(bounds,"",2)==0) ? 1: MPI_UNDEFINED;
							MPI_Comm_split(MACSIO_MAIN_Comm, color, MACSIO_MAIN_Rank, &comm);
							if (color == 1){
							MPI_Gather(y_coord, local_mesh_dims[1], MPI_DOUBLE, y_coord_root, local_mesh_dims[1], MPI_DOUBLE, 0, comm);
							}					
						}
						if (ndims > 2){
							color = (JsonGetInt(bounds, "", 0)==0 && JsonGetInt(bounds,"",1)==0) ? 1: MPI_UNDEFINED;
							MPI_Comm_split(MACSIO_MAIN_Comm, color, MACSIO_MAIN_Rank, &comm);
							if (color == 1)
							MPI_Gather(z_coord, local_mesh_dims[2], MPI_DOUBLE, z_coord_root, local_mesh_dims[2], MPI_DOUBLE, 0, comm);

						}
					} else {
						x_coord = json_object_extarr_data(json_object_path_get_extarr(coords, "XCoords"));
						y_coord = json_object_extarr_data(json_object_path_get_extarr(coords, "YCoords"));
						z_coord = json_object_extarr_data(json_object_path_get_extarr(coords, "ZCoords"));
					}

				} else {
					//Variable 
					json_object *vars_array = json_object_path_get_array(part_obj, "Vars");
					json_object *var_obj = json_object_array_get_idx(vars_array, v);
					json_object *extarr_obj = json_object_path_get_extarr(var_obj, "data");

					buf = json_object_extarr_data(extarr_obj);
				}
			}

			if (v == -1) { 

				if (mesh_type == TIO_MESH_QUAD_COLINEAR){
					if (MACSIO_MAIN_Rank == 0){					
						TIO_Call( TIO_Write_QuadMesh_All(file_id, mesh_id, TIO_DOUBLE, x_coord_root, y_coord_root, z_coord_root),
							"Write Quad Mesh All Failed\n");
					}
				} else {
					TIO_Call( TIO_Write_QuadMesh_Chunk(file_id, mesh_id, MACSIO_MAIN_Rank, TIO_XFER,
														TIO_DOUBLE, x_coord, y_coord, z_coord),
								"Write Non-Colinear Mesh Coords failed\n");
				}
	    
			} else {				
				TIO_Call( TIO_Write_QuadQuant_Chunk(file_id, object_id, MACSIO_MAIN_Rank, 
												TIO_XFER, dtype_id, buf, (void*)TIO_NULL),
					"Write Quad Quant Chunk Failed\n");

				TIO_Call( TIO_Close_Quant(file_id, object_id),
					"Close Quant Failed\n");
			}

			free(x_coord_root);
			free(y_coord_root);
			free(z_coord_root);

		}

		free(centering);
	}

	TIO_Call( TIO_Close_Mesh(file_id, mesh_id),
          "Close Mesh Failed\n");

}

static void write_ucd_mesh_whole(
	TIO_File_t file_id, 
	TIO_Object_t state_id,
	json_object *main_obj)
{

}

static void main_dump_sif(json_object *main_obj, int dumpn, double dumpt)
{
#ifdef HAVE_MPI
	char const *mesh_type = json_object_path_get_string(main_obj, "clargs/part_type");
	char fileName[256];
	TIO_Object_t state_id, variable_id;
	char *state_name = "state0";

	TIO_File_t tiofile_id;

	MPI_Info mpiInfo = MPI_INFO_NULL;

	/* Construct name for the HDF5 file */
	sprintf(fileName, "%s_typhonio_%03d.%s",
	        json_object_path_get_string(main_obj, "clargs/filebase"),
	        dumpn,
	        "h5"); //json_object_path_get_string(main_obj, "clargs/fileext"));

	char *date = (char*)getDate();

	TIO_Call( TIO_Create(fileName, &tiofile_id, TIO_ACC_REPLACE, "MACSio",
	                     "0.9", date, fileName, MACSIO_MAIN_Comm, MPI_INFO_NULL, MACSIO_MAIN_Rank),
	          "File Creation Failed\n");

	TIO_Call( TIO_Create_State(tiofile_id, state_name, &state_id, 1, (TIO_Time_t)0.0, "us"),
	          "State Create Failed\n");

	json_object *part_array = json_object_path_get_array(main_obj, "problem/parts");
	json_object *part_obj = json_object_array_get_idx(part_array, 0);
	json_object *mesh_obj = json_object_path_get_object(part_obj, "Mesh");

	if (!strcmp(JsonGetStr(mesh_obj, "MeshType"), "rectilinear"))
		write_quad_mesh_whole(tiofile_id, state_id, main_obj, TIO_MESH_QUAD_COLINEAR);
    else if (!strcmp(JsonGetStr(mesh_obj, "MeshType"), "curvilinear"))
        write_quad_mesh_whole(tiofile_id, state_id, main_obj, TIO_MESH_QUAD_NONCOLINEAR);
    else if (!strcmp(JsonGetStr(mesh_obj, "MeshType"), "ucdzoo"))
        write_ucd_mesh_whole(tiofile_id, state_id, main_obj);
    else if (!strcmp(JsonGetStr(mesh_obj, "MeshType"), "arbitrary"))
        write_ucd_mesh_whole(tiofile_id, state_id, main_obj);

	TIO_Call( TIO_Close_State(tiofile_id, state_id),
			"State Close Failed\n");

	TIO_Call( TIO_Close(tiofile_id),
	          "Close File failed\n");


#endif
}

static void main_dump(int argi, int argc, char **argv, json_object * main_obj, int dumpn, double dumpt)
{
	int rank, size, numFiles;
	/* Without this barrier, I get strange behavior with Silo's MACSIO_MIF interface */
	mpi_errno = MPI_Barrier(MACSIO_MAIN_Comm);

	/* process cl args */
	process_args(argi, argc, argv);

	rank = json_object_path_get_int(main_obj, "parallel/mpi_rank");
	size = json_object_path_get_int(main_obj, "parallel/mpi_size");

	/* ensure we're in MIF mode and determine the file count */
	json_object *parfmode_obj = json_object_path_get_array(main_obj, "clargs/parallel_file_mode");
	if (parfmode_obj) {
		json_object *modestr = json_object_array_get_idx(parfmode_obj, 0);
		json_object *filecnt = json_object_array_get_idx(parfmode_obj, 1);

		if (!strcmp(json_object_get_string(modestr), "SIF")) {
			main_dump_sif(main_obj, dumpn, dumpt);
		}
		else {
			numFiles = json_object_get_int(filecnt);
			main_dump_mif(main_obj, numFiles, dumpn, dumpt);
		}
	}
	else {
		char const * modestr = json_object_path_get_string(main_obj, "clargs/parallel_file_mode");
		if (!strcmp(modestr, "SIF")) {
			float avg_num_parts = json_object_path_get_double(main_obj, "clargs/avg_num_parts");
			if (avg_num_parts == (float ((int) avg_num_parts)))
				main_dump_sif(main_obj, dumpn, dumpt);
			else {
				// CURRENTLY, SIF CAN WORK ONLY ON WHOLE PART COUNTS
				MACSIO_LOG_MSG(Die, ("HDF5 plugin cannot currently handle SIF mode where "
				                     "there are different numbers of parts on each MPI rank. "
				                     "Set --avg_num_parts to an integral value." ));
			}
		}
		else if (!strcmp(modestr, "MIFMAX"))
			numFiles = json_object_path_get_int(main_obj, "parallel/mpi_size");
		else if (!strcmp(modestr, "MIFAUTO")) {
			/* Call utility to determine optimal file count */
			// ADD UTILIT TO DETERMINE OPTIMAL FILE COUNT
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


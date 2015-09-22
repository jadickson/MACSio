#include <macsio_clargs.h>
#include <macsio_iface.h>
#include <macsio_log.h>
#include <macsio_main.h>
#include <macsio_mif.h>
#include <macsio_utils.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <typhonio.h>
#include <hdf5.h>

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

static const char *filename;

static int process_args(
		int argi,
		int argc, 
		char *argv[]
)
{
	const MACSIO_CLARGS_ArgvFlags_t argFlags = {MACSIO_CLARGS_WARN, MACSIO_CLARGS_TOMEM};

	MACSIO_CLARGS_ProcessCmdline(0, argFlags, argi, argc, argv, 

			MACSIO_CLARGS_END_OF_ARGS);

}

char errstr[TIO_STRLEN];
TIO_t errnum;

#define TIO_Call(rtn,emsg) if ((errnum=rtn) != TIO_SUCCESS){\
	TIO_Get_Error(errnum, errstr);\
	std::cerr << emsg << "\n" << errstr << "\n" << std::endl;\
	exit(EXIT_FAILURE);}

static void *CreateTyphonIOFile(
		const char *fname, 
		const char *nsname, 
		void *userData
		)
{
	TIO_File_t fileID;

	TIO_Call( TIO_Create(fname, &fileID, TIO_ACC_REPLACE, codeName, 
				version, date, fname, MACSIO_MAIN_Comm, MPI_INFO_NULL, rank),
			"File Creation Failed\n");

	return (void *) fileID;
}

static void *OpenTyphonIOFile(
		const char *fname,
		const char *nsname,
		MACSIO_MIF_ioFlags_t ioFlags, 
		void *userData
		)
{
	TIO_File_t fileID;

	TIO_Call( TIO_Open(fname, &fileID, TIO_ACC_READWRITE, codeName, 
				version, date, fname, MACSIO_MAIN_Comm, MPI_INFO_NULL, rank),
			"File Open Failed\n");

	return (void *) fileID;
}

static void CloseTyphonIOFile(
		void *file,
		void *userData
		)
{
	TIO_Call( TIO_Close(file), 
			"File Close Failed\n");
}



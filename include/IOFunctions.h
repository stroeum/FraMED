/*
 *  IOfunctions.h
 *  Created by Jeremy Riousset on 10/26/07.
 */

#ifndef IOFUNCTIONS_H
#define IOFUNCTIONS_H
#include "Input.h"
#include <sstream>

class IO
{
public:
    static void print(FILE*,                            string);

    static ListDouble read(char *);
	static void write(ListLink&,						char *);
	static void write(double,							char *);
	static void write(int,int, int,						char *);
	static void write(double,double, double,			char *);
	static void write(double, double, double, double,	char *);
	static void write(ListDouble&,						char *);

    static void read(int&,                              char *);
    static void read(double&,                           char *);
    static void read(SizeGrid&,                         char *);
    static void read(ResGrid&,                          char *);
    static void read(ListLink&,                         char *);
    static void read(ListDouble&,                       char *);
    static void read(CMatrix1D &,                       char *);
    static void read(CMatrix2D &,                       char *);
    static SizeGrid read(CMatrix3D &,                   char *);

    /* SAM function.  Identical to 'write(ListDouble&, char*)', except
	 * writes a list of ints.
	 */
	static void write(ListInt&,							char *);
	static void write(ListVector&,						char *);
	static void write(ListCMatrix1D&,					char *);
	static void write(CMatrix1D&,						char *);
	static void write(CMatrix2D&,						char *);
	static void write(CMatrix3D&,						char *);

    
	/* SAM functions */

	/* Used to change the path name to a new location.
	 * All IO functions write to the current subfolder.
	 * The old subfolder (if there is one) is stored in the
	 * 'path_name_old' static variable.
	 */
	static void setPathName(char * input_path_name);
	/* Return the current path name. */
	static void getPathName(char * ret_path_name);

	/* Used to append the '/' to 'input_path_name'
	 * before appending 'fname'.
	 */
	static void formFileName(char * fname);
	/* Similar to 'fopen()', except appends the file name to the path name. */
	static FILE * openFile(char * fname, string mode);

	/* Deletes the previous subfolder.  WARNING: uses system-specific
	 * command line commands.
	 */
	static void removeOldDirectory();
	/* Copies the old subfolder into the new directory. */
	static void copyFromOldDirectory();
	/* Calls 'copyFromOldDirectory()', then 'removeOldDirectory()'. */
	static void moveFromOldDirectory();

	static char * path_name; /* Specifies current subfolder. */
	static char * path_name_old; /* Specifies PREVIOUS subfolder. */
	static char * file_name; /* Specifies current file name. */
	static FILE * fp; /* Pointer to the current file. */
};

#endif

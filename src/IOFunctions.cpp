/*
 *  IOFunctions.cpp
 *  Created by Jeremy Riousset on 10/26/07.
 */

#include "IOFunctions.h"

/**************************************************************************************/
/* Read a file with doubles, store in a list of double								  */
/**************************************************************************************/
void IO::read(int& nn, char * fname)
{
    fp = openFile(fname, "r");
    double _d(0);
    if(fp) fscanf(fp,"%lf\n", &_d);
    nn=fabs((int)_d);
    fclose(fp);
}

void IO::read(double& _d, char * fname)
{
    fp = openFile(fname, "r");
    if(fp) fscanf(fp,"%lf\n", &_d);
    fclose(fp);
}

ListDouble IO::read(char * fname)
{
	ListDouble LL;
	double xx;
	ifstream inFile;

	inFile.open(fname);
	if(!inFile)
		cerr<<"cannot open file";
	//nrerror("cannot open file");

	while(inFile>>xx)
		LL.push_back(xx);
	inFile.close();
	return LL;
}

void IO::read(ListDouble& LList, char * ffile)
{
    char *  line        = NULL;
    size_t  line_length = 0;
    ssize_t nread;
    double  value;
    
    fp = openFile(ffile, "r");
    if(fp)
        while ((nread = getline(&line, &line_length, fp)) != -1) {
            sscanf(line,"%lf\n", &value );
            LList.push_back(value);
        }
    free(line);
    fclose(fp);
}

void IO::read(ListLink& TTree, char * ffile)
{
    char *  line        = NULL;
    size_t  line_length = 0;
    ssize_t nread;
    Link    LL;

    fp = openFile(ffile, "r");
    if(fp)
         while ((nread = getline(&line, &line_length, fp)) != -1) {
            sscanf(line,"%5d %5d %5d %5d %5d %5d %lf %lf %lf %lf\n", &(LL.start.i), &(LL.start.j), &(LL.start.k), &(LL.end.i), &(LL.end.j), &(LL.end.k), &(LL.l), &(LL.efield), &(LL.deltaV), &(LL.proba) );
            TTree.push_back(LL);
        }
    free(line);
    fclose(fp);
}


void IO::read(SizeGrid& _N, char * ffile)
{
    fp = openFile(ffile, "r");
    double N[3];
    if(fp) fscanf(fp,"%lf\n%lf\n%lf\n", &N[0],  &N[1], &N[2] );
    fclose(fp);
    _N.init((int)N[0],(int)N[1],(int)N[2]);
}

void IO::read(ResGrid& _d, char * ffile)
{
    fp = openFile(ffile, "r");
    double d[3];
    if(fp) fscanf(fp,"%lf\n%lf\n%lf\n", &d[0],  &d[1], &d[2] );
    fclose(fp);
    _d.init(d[0],d[1],d[2]);
}

SizeGrid IO::read(CMatrix3D & M, char * ffile)
{

    // Read Matrix Size //
    fp = openFile(ffile, "r");
    char	tmp_c, tmp_cc;
    string	tmp_s;
    double tmp_d;
    SizeGrid N;
    
    N.x = 0;	N.y = 0;	N.z = 0;
    while((tmp_c=fgetc(fp)) != EOF)
	{
		if(tmp_c == ' ')
			N.y++;
		if(tmp_c == '\n')
		{
			N.x++;
			tmp_cc = fgetc(fp);
			if(tmp_cc == '\n' && tmp_cc !=EOF)
				N.z++;
		}
	}
	N.x  = N.x/N.z;
	N.y /= N.z*N.x;
    cout<<N.x<<" "<<N.y<<" "<<N.z<<endl;
    
    rewind(fp);
    
    int ii(0),jj(0),kk(0);

    while((tmp_c=fgetc(fp)) != EOF)
	{
		tmp_s += tmp_c;
		if(tmp_c == ' ')
		{
			tmp_d = atof(tmp_s.c_str());
			M[ii][jj][kk] = tmp_d;
			jj++;
			tmp_s = "";
		}
		if(tmp_c == '\n')
		{
			ii++;
			jj=0;
			tmp_cc = fgetc(fp);
			if(tmp_cc == '\n' && tmp_cc !=EOF)
			{
				ii=0;
				jj=0;
				kk++;
			}
		}
	}
    fclose(fp);
    return N;
}

void IO::read(CMatrix2D & M, char * ffile)
{
    
    // Read Matrix Size //
    fp = openFile(ffile, "r");
    char	tmp_c;
    string	tmp_s;
    int m(0),n(0);
    
    while((tmp_c=fgetc(fp)) != EOF)
    {
        if(tmp_c == ' ' || tmp_c == '\n')
            n++;
        if(tmp_c == '\n')
            m++;
    }
    n /= m;
    
    M.init(m,n);
    
    rewind(fp);
    
    
    int ii(0),jj(0);
    
    char *  line        = NULL;
    size_t  line_length = 0;
    ssize_t nread;
    fp = openFile(ffile, "r");
    int nn(0);
    if(fp){
        ii = -1;
        while ((nread = getline(&line, &line_length, fp)) != -1) {
            if(nread > 1) {
                ii++;
            }
            istringstream iss(line);
            string word;
            
            while(iss >> word) {
                nn++;
                M[ii][jj] = atof(word.c_str());
                jj++;
            }
            jj=0;
        }
    }
    free(line);
    fclose(fp);
}

void IO::read(CMatrix1D & M, char * ffile)
{
    // Read Matrix Size //
    fp = openFile(ffile, "r");
    if(fp == 0) {
        perror("fopen");
        exit(1);
    }
    char	tmp_c;
    string	tmp_s;
    int m(0);
    
    while((tmp_c=fgetc(fp)) != EOF)
    {
        if(tmp_c == '\n')
            m++;
    }
    M.init(m);
    
    rewind(fp);
  
    int ii(0);
    
    char *  line        = NULL;
    size_t  line_length = 0;
    ssize_t nread;
    fp = openFile(ffile, "r");
    if(fp){
        ii = -1;
        while ((nread = getline(&line, &line_length, fp)) != -1) {
            //cout << nread<<endl;
            if(nread > 1) {
                ii++;
            }
            istringstream iss(line);
            string word;
            
            while(iss >> word) {
                M[ii] = atof(word.c_str());
            }
        }
    }
    free(line);
    fclose(fp);
}

/**************************************************************************************/

/* SAM functions: The structure of the 'IO::write' functions has been changed.
 * Now, a subfolder must be specified (see 'setPathName'), and the specified file name
 * (i.e. 'fname) will be appended to that folder name.  Also, there is one file pointer
 * ('fp') which is shared between all functions.  Every time the file name is 'set', the
 * the specified file is opened and assigned to 'fp'.
 */

/**************************************************************************************/
/* write a double, 3 ints, 3 doubles, 4 doubles,									  */
/* a list of doubles, vectors, 1-D matrices											  */
/**************************************************************************************/
void IO::write(char * _s, char * fname)
{
	fp = openFile(fname, "w");

	if(fp) fprintf(fp,"%s\n", _s);
	fclose(fp);
}

void IO::write(double _d, char * fname)
{
	fp = openFile(fname, "w");

	if(fp) fprintf(fp,"%f\n", _d);
	fclose(fp);
}

void IO::write(int nn, int mm, int pp,	char * fname)
{
	fp = openFile(fname, "w");

	if(fp) fprintf(fp,"%d\n%d\n%d\n", nn, mm, pp);
	fclose(fp);
}

void IO::write(double aa, double bb, double cc,	char * fname)
{
	fp = openFile(fname, "w");

	if(fp) fprintf(fp,"%f\n%f\n%f\n", aa, bb, cc);
	fclose(fp);
}

void IO::write(double aa, double bb, double cc, double _d,	char * fname)
{
	fp = openFile(fname, "w");

	if(fp) fprintf(fp,"%f\n%f\n%f\n%f\n", aa, bb, cc, _d);
	fclose(fp);
}

void IO::write(ListDouble& LL, char * fname)
{
	fp = openFile(fname, "w");
	ListDouble::iterator it;

	if(fp)
		for (it=LL.begin() ; it!=LL.end() ; it++)
			fprintf(fp,"%f\n", *it);
	fclose(fp);
}

/* SAM function: The following function was added.  Almost identical to the
 * 'IO::write(ListDouble& LL, char * fname)' file above it.
 */
void IO::write(ListInt& LL, char * fname)
{
	fp = openFile(fname, "w");
	ListInt::iterator it;

	if(fp)
		for(it=LL.begin() ; it!=LL.end() ; it++)
			fprintf(fp,"%d\n", *it);
	fclose(fp);
}

void IO::write(ListVector& LL, char * fname)
{
	fp = openFile(fname, "w");
	ListVector::iterator it;

	if(fp)
		for (it=LL.begin() ; it!=LL.end() ; it++)
			fprintf(fp,"%f %f %f\n", it->x,it->y,it->z);
	fclose(fp);
}

void IO::write(ListCMatrix1D& LL, char * fname)
{
	fp = openFile(fname, "w");
	ListCMatrix1D::iterator it;

	if(fp)
		for (it=LL.begin() ; it!=LL.end() ; it++)
		{
			for (int ii=0 ; ii<it->getNbElem() ; ii++)
				fprintf(fp,"%f ",it->getElem(ii));
			fprintf(fp,"\n");
		};
	fclose(fp);
}

void IO::write(ListLink& LL, char * fname)
{
	fp = openFile(fname, "w");
	ListLink::iterator	it;

	if(fp)
		for (it=LL.begin() ; it!=LL.end() ; it++)
			fprintf(fp,"%5d %5d %5d %5d %5d %5d %+12.2f %+12.2f %+12.2f %5.4f\n", it->start.i, it->start.j, it->start.k, it->end.i, it->end.j, it->end.k, it->l, it->efield, it->deltaV, it->proba );
	fclose(fp);
}

void IO::write(CMatrix1D& MM, char * fname)
{
	formFileName(fname);
	MM.fwrite(file_name);
}

void IO::write(CMatrix2D& MM, char * fname)
{
	formFileName(fname);
	MM.fwrite(file_name);
}

void IO::write(CMatrix3D& MM, char * fname)
{
	formFileName(fname);
	MM.fwrite(file_name);
}

/* Change the current subfolder.  Also, store the old path name */
void IO::setPathName(char input_path_name[])
{
	/* If the path name isn't NULL, then we need to store it. */
	if(path_name != NULL)
	{
		if(path_name_old != NULL)
			free(path_name_old);
		path_name_old = (char *) malloc(strlen(path_name) + 1);
		strcpy(path_name_old, path_name);
		free(path_name);
	}
	path_name = (char *) malloc(strlen(input_path_name) + 1);
	strcpy(path_name, input_path_name);
	file_name = (char *) malloc(strlen(input_path_name) + 35);
}

/* Return the current path name. */
void IO::getPathName(char ret_path_name[])
{
	if(path_name != NULL)
		strcpy(ret_path_name, path_name);
	else
		strcpy(ret_path_name, "");
}

void IO::formFileName(char * fname)
{
	/* First allocate space for the file name.  Remember, we need spots
	 * for the '/' and the NULL character at the end. */
	file_name = (char *) malloc(strlen(path_name) + strlen(fname) + 2);
	strcpy(file_name, path_name);
	strcat(file_name, "/");
	strcat(file_name, fname);
   
}
/* Similar to 'fopen()', except appends the file name to the path name. */
FILE * IO::openFile(char * fname, string mode)
{
	formFileName(fname);
	fp = fopen(file_name, mode.c_str());
	free(file_name);
	file_name = NULL;

	return fp;
}

void IO::removeOldDirectory()
{
	if(path_name_old != NULL)
	{
		/* To save time, quotation marks in the string are replaced with backslashes. */
		char * command = (char *) malloc(	strlen(	"test -d // && echo /ee: Directory") +
											strlen( "cannot be removed./ || echo /++: Directory") +
											strlen( "successfully removed!/") +
											strlen(path_name_old) + 1);

		printf("..: Removing interim results directory: %s.\n", path_name_old);
		/* WARNING: 'rm', 'test', and 'echo' command line programs used. */

		/* First, remove the old directory. */
		strcpy(command, "rm -r \"");
		strcat(command, path_name_old);
		strcat(command, "\"");
		system(command);

		/* Then check to see if it still exists.  If it does, there is a problem. */
		strcpy(command, "test -d \"");
		strcat(command, path_name_old);
		strcat(command, "\" && echo \"ee: Directory cannot be removed.\" || ");
		strcat(command, "echo \"++: Directory successfully removed!\"");
		system(command);


		free(path_name_old);
		free(command);
		path_name_old = NULL;
	}
}

void IO::copyFromOldDirectory()
{
	if(path_name_old != NULL)
	{
		char * command = (char *) malloc(	strlen("cp -r // //") + strlen(path_name_old) +
											strlen(path_name) + 1);
		char * new_dir_name;
		char * name_ptr;

		/* WARNING: 'cp', 'test', and 'echo' command line programs used. */

		/* Copy the subfolder over to the new directory. */
		strcpy(command, "cp -r \"");
		strcat(command, path_name_old);
		strcat(command, "\" \"");
		strcat(command, path_name);
		strcat(command, "\"");
		system(command);

		/* Form the new path name.  When the subfolder was copied, a NEW subfolder was created.
		 *
		 * Graphically: C:/A copied to C:/B would yield C:/B/A
		 */

		/* strrchr finds the last occurence of a given character and returns the substring bounded
		 * by it.
		 */
		name_ptr = strrchr(path_name_old, '/');
		new_dir_name = (char *) malloc(strlen(path_name) + strlen(name_ptr) + 1);
		strcpy(new_dir_name, path_name);
		strcat(new_dir_name, name_ptr);

		free(command);
		command = (char *) malloc(	strlen(	"test -d // && echo /++: Results successfully") +
									strlen( "copied!/ || echo /ee: Copy unsuccessful./") +
									strlen(new_dir_name) + 1);

		printf("..: Copying results to final directory: %s.\n", new_dir_name);

		/* Test to see if the new directory exists.  If it doesn't, there is a problem. */
		strcpy(command, "test -d \"");
		strcat(command, new_dir_name);
		strcat(command, "\" && echo \"++: Results successfully copied!\"");
		strcat(command, " || echo \"ee: Copy unsuccessful.\"");
		system(command);

		free(command);
		free(new_dir_name);
	}
}

void IO::moveFromOldDirectory()
{
	/* Pretty self-explanatory. */
	copyFromOldDirectory();
	removeOldDirectory();
}

void IO::print(FILE* file, string msg)
{
    printf(		  "%s", msg.c_str());
    fprintf(file, "%s", msg.c_str());
}


/* Basic static variables.  For more information, see 'IOFunctions.h'. */
char * IO::path_name = NULL;
char * IO::file_name = NULL;
char * IO::path_name_old = NULL;
FILE * IO::fp = NULL;
/**************************************************************************************/

/**************************************************************************************/

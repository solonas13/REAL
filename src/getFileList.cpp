/**
    REAL: An efficient REad ALigner for next generation sequencing reads.
    Copyright (C) 2010 German Tischler, Solon Pissis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#if defined(HAVE_CONFIG_H)
#include "real_config.hpp"
#endif

#include "getFileList.hpp"

#if defined(HAVE_WINDOWS_H)
#include <windows.h>
#endif

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_DIRENT_H)
#include <dirent.h>
#endif

#if defined(HAVE_SYS_STAT_H)
#include <sys/stat.h>
#endif

#include <stdexcept>

#include "stringFunctions.hpp"

#if defined(HAVE_WINDOWS_H)
bool isFile(std::string const & filename)
{
	WIN32_FIND_DATA file_data;	
	HANDLE h = FindFirstFile(filename.c_str(),&file_data);
	
	if ( h == INVALID_HANDLE_VALUE )
		return false;
	
	unsigned int filefound = 0, dirfound = 0;
	
	do
	{
		if (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) dirfound++;
		else filefound++;
	}
	while ( FindNextFile(h,&file_data) != 0 );

	return (dirfound == 0) && (filefound == 1);
}

bool isDirectory(std::string const & filename)
{
	WIN32_FIND_DATA file_data;	
	HANDLE h = FindFirstFile(filename.c_str(),&file_data);
	
	if ( h == INVALID_HANDLE_VALUE )
		return false;
	
	unsigned int filefound = 0, dirfound = 0;
	
	do
	{
		if (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) dirfound++;
		else filefound++;
	}
	while ( FindNextFile(h,&file_data) != 0 );

	return (dirfound == 1) && (filefound == 0);	
}

void enumerateFilesInDirectory(std::string const & dirname, std::vector<std::string> & filenames, std::string const & suffix)
{
	WIN32_FIND_DATA file_data;	
	HANDLE h = FindFirstFile((dirname + "\\*").c_str(),&file_data);
	
	if ( h != INVALID_HANDLE_VALUE )
		do
		{
			if (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
			{
				std::string subdirname = file_data.cFileName;
				if ( subdirname != "." && subdirname != ".." )
					enumerateFilesInDirectory(dirname + "\\" + subdirname, filenames, suffix);
			}
			else
			{
			        std::string const & compfilename = (dirname + "\\" + file_data.cFileName);
			        if ( stringFunctions::endsOn(compfilename,suffix) )
        				filenames.push_back(compfilename);
			}                                                                        
		}
		while ( FindNextFile(h,&file_data) != 0 );
}
#else
bool isFile(std::string const & filename)
{
	struct stat mystat;
	if ( stat(filename.c_str(),&mystat) == 0 )
	{
		if ( S_ISREG(mystat.st_mode) )
			return true;
		else
			return false;
	}
	else
		return false;
}
bool isDirectory(std::string const & filename)
{
	struct stat mystat;
	if ( stat(filename.c_str(),&mystat) == 0 )
	{
		if ( S_ISDIR(mystat.st_mode) )
			return true;
		else
			return false;
	}
	else
		return false;
}
void enumerateFilesInDirectory(std::string const & dirname, std::vector<std::string> & files, std::string const & suffix)
{
	DIR * mydir = opendir(dirname.c_str());

	if ( ! mydir )
		throw std::runtime_error("Could not open file/directory");

	struct dirent * direntry = 0;

	while ( (direntry = readdir(mydir)) )
	{
		std::string const filename = dirname+"/"+direntry->d_name;
		struct stat mystat;
		if ( stat(filename.c_str(),&mystat) == 0 )
		{
			if ( S_ISREG(mystat.st_mode) )
			{
			        if ( stringFunctions::endsOn(filename,suffix) )
        				files.push_back(filename);
                        }
			else if ( S_ISDIR(mystat.st_mode) && 
				std::string(direntry->d_name) != "." &&
				std::string(direntry->d_name) != ".."
			)
				enumerateFilesInDirectory(filename+"/",files,suffix);
		}
	}

	closedir(mydir);
}
#endif

void getFileList(std::string const & dirname, std::vector<std::string> & files, std::string const & suffix)
{
	if ( isFile(dirname) && stringFunctions::endsOn(dirname,suffix) )
		files.push_back(dirname);
	else if ( isDirectory(dirname) )
		enumerateFilesInDirectory(dirname,files,suffix);
}

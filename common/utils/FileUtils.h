#ifndef FILE_UTILS_H_
#define FILE_UTILS_H_

#include <iostream>
#include <fstream>
#include <string>
#include "sys/fcntl.h"
#include "sys/mman.h"
#include <stdio.h>

#include <eerno.h>
extern int errno;

using namespace std;

bool FileExists(string &fileName) {

  FILE *fp = fopen(fileName.c_str(),"r");
  if( fp ) {
    // exists
    fclose(fp);
    return true;
  } 
  else {
    return false;
  }
}

inline void CriticalOpenRead(string &fileName, ifstream &file, std::ios::openmode mode=std::ios::in) {
	file.open(fileName.c_str(), mode | std::ios::in);
	if (!file.good()) {
		cerr << "Could not open file:"  << fileName << endl;
		exit(1);
	}
}

inline int OpenRead(string &fileName, ifstream &file, std::ios::openmode mode=std::ios::in) {
	file.open(fileName.c_str(), mode | std::ios::in);
	return file.good();
}


inline void CriticalOpenWrite(string &fileName, ofstream &file, std::ios::openmode mode=std::ios::out) {
	file.open(fileName.c_str(), mode | std::ios::out);
	if (!file.good()) { 
		cerr << "Could not open file: " << fileName << endl;
		exit(1);
	}
	
}

inline int OpenWrite(string &fileName, ofstream &file, std::ios::openmode mode=std::ios::out) {
	file.open(fileName.c_str(), mode | std::ios::out);
	return file.good();
}


int CountLinesInFile(string fileName) {
  
  char* filePtr;
  long fileSize;
  int fileDes;
  fileDes = open(fileName.c_str(), O_RDONLY);  
  fileSize = lseek(fileDes, 0, SEEK_END);
  lseek(fileDes, 0, SEEK_SET);
  filePtr = (char*) mmap(0, fileSize, PROT_READ, MAP_PRIVATE, fileDes, 0);
  if (filePtr==MAP_FAILED) {
    cout << "mmap memory assignment in CountLinesInFile failed for filesize " << fileSize << " for file "<<fileName <<endl;
    cout << "Error number " << errno << "for file name" << fileName << endl;
    exit(1);
  }
  long pos;
  int numLines = 0;
  for (pos = 0; pos < fileSize; pos++, filePtr++) {
    if (*filePtr == '\n') {
      numLines++;
    }
  }
  return numLines;
}


#endif

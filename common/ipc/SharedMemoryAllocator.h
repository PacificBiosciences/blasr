#ifndef SHARED_MEMORY_ALLOCATOR_H_
#define SHARED_MEMORY_ALLOCATOR_H_

#include <iostream>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

using namespace std;
template<typename T_Data>
int AllocateMappedShare(string &handle, int dataLength, T_Data *&dataPtr, int &shmId) {
	cout << "opening shm" << endl;
	shmId = shm_open(handle.c_str(), O_CREAT| O_RDWR, S_IRUSR | S_IWUSR);
	if (ftruncate(shmId, sizeof(T_Data[dataLength])) == -1) {
		cout <<" ftruncate error: " << errno << endl;
	}
	cout << "done truncating." << endl;
	dataPtr = (T_Data*) mmap(NULL, sizeof(T_Data[dataLength]),
												PROT_READ | PROT_WRITE, MAP_SHARED, shmId, 0); 
	if (dataPtr == MAP_FAILED) {
		// 
		// Handle this better later on.
		//
		cout << "ERROR, MEMORY MAP FAILED." << endl;
		exit(1);
	}
	cout << "done mapping." << endl;
	return dataLength;
}


#endif

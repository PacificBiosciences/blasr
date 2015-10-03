#ifndef ALIGNMENT_MAPPING_SEMAPHORE_H_
#define ALIGNMENT_MAPPING_SEMAPHORE_H_
#include <vector>
#include <pthread.h>
#include <semaphore.h>

#ifndef __APPLE__
class MappingSemaphores {
  public:
	sem_t reader;
	sem_t writer;
	sem_t unaligned;
	sem_t hitCluster;
	MappingSemaphores& operator=(MappingSemaphores &rhs) {
		return *this;
	}

	void InitializeAll() {
		sem_init(&reader, 0, 1);
		sem_init(&writer, 0, 1);
		sem_init(&unaligned, 0, 1);
		sem_init(&hitCluster, 0, 1);
	}
};
#else
class MappingSemaphores {
  public:
	sem_t *reader;
	sem_t *writer;
	sem_t *unaligned;
	sem_t *hitCluster;
	MappingSemaphores& operator=(MappingSemaphores &rhs) {
		return *this;
	}

	void InitializeAll() {
		reader     = sem_open("/reader",     O_CREAT, 0644, 1);
		writer     = sem_open("/writer",     O_CREAT, 0644, 1);
		unaligned  = sem_open("/unaligned",  O_CREAT, 0644, 1);
		hitCluster = sem_open("/hitCluster", O_CREAT, 0644, 1);
	}
};
#endif

#endif

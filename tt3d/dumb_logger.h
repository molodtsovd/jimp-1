#ifndef DUMB_LOGGER
#define DUMB_LOGGER

#include <fstream>
#include <time.h>

namespace dumb{

class logger
{
	time_t t;
	char path[1024];
	FILE* f;
public:
	char buf[1024];
	void init(char* filename);
	logger& operator << (int);
};

}

void dumb::logger::init(char* filename)	
{ 
	int i = 0;
	while(filename[i])
	{
		path[i] = filename[i];
		++i;
	}
	path[i] = 0;
}

dumb::logger& dumb::logger::operator<<(int)
{
#pragma omp critical (logger)
	{
	f = fopen(path, "a");

	time(&t);
	tm* ht = localtime(&t); 
	
	fprintf(f,"\n%d/%d/%d %d:%02d:%02d	%s\n", ht->tm_mday, 1 + ht->tm_mon, 1900 + ht->tm_year, ht->tm_hour, ht->tm_min, ht->tm_sec, buf); 

	printf("\n%s\n", buf);
	
	fclose(f);
	}
	return *this;
}
#endif
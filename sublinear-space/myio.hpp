#ifndef _BINARY_GRAPH_READER_
#define _BINARY_GRAPH_READER_


#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <assert.h>

#include <unistd.h> // read
#include <sys/time.h> // gettimeofday

#include <fcntl.h> // open
#include <sys/stat.h>

#include <iostream>
#include <iomanip>


using namespace std;

template <class T> class BinaryGraphReader{
private:
  int fd;
  double iotime;
  T* buf;
  size_t buf_size;
  size_t cur_size;
  size_t p;
  T src, dst;
public:
  BinaryGraphReader(size_t bs);
  ~BinaryGraphReader();
  bool nextEdge();
  inline T getSrc(){return src;}
  inline T getDst(){return dst;}
  void setfd(int _fd);
  void reset();
  inline double gettime(){
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
  }
  inline double getiotime(){return iotime;}
};


template <class T> class BinaryWriter{
private:
  int fd;
  double iotime;
  T* buf;
  size_t buf_size;
  size_t p;
public:
  BinaryWriter(size_t bs);
  ~BinaryWriter();
  void setfd(int _fd);
  inline size_t bwwrite(const T* b, size_t count);
  inline size_t bwflush();
  inline double gettime(){
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
  }
  inline double getiotime(){return iotime;}
};

#include "myio.cpp"

#endif

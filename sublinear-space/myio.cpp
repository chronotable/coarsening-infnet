#include "myio.hpp"

///// BinaryGraphReader
template <class T> BinaryGraphReader<T>::BinaryGraphReader(size_t bs){
  assert((bs & 1) == 0);
  fd = -1;
  iotime = 0;
  buf_size = bs;
  try{
    buf = new T[buf_size];
  }catch(...){
    fprintf(stderr, "cannot new in %s\n", __func__);
    exit(1);
  }
  cur_size = 0;
  p = 0;
  src = dst = 0;
}

template <class T> BinaryGraphReader<T>::~BinaryGraphReader(){
  delete [] buf;
  cout << __func__ << ": iotime=" << setprecision(10) << iotime << endl;
}
template <class T> bool BinaryGraphReader<T>::nextEdge(){
  assert(fd >= 0);
  if(fd < 0)
    return false;
  if(p >= cur_size){
    double st = gettime();
    cur_size = read(fd, buf, buf_size * sizeof(T)) / sizeof(T);
    iotime += gettime() - st;
    if(cur_size <= 0) // -1 if error
      return false;
    p = 0;
  }
  src = buf[p];
  dst = buf[p+1];
  p += 2;
  return true;
}
template <class T> void BinaryGraphReader<T>::setfd(int _fd){
  fd = _fd;
  assert(lseek(fd, 0, SEEK_SET) == 0);
  cur_size = p = 0;
  src = dst = 0;
}
template <class T> void BinaryGraphReader<T>::reset(){
  assert(lseek(fd, 0, SEEK_SET) == 0);
  cur_size = p = 0;
  src = dst = 0;
}
/////




///// BinaryWriter
template <class T> BinaryWriter<T>::BinaryWriter(size_t bs){
  fd = -1;
  iotime = 0;
  buf_size = bs;
  try{
    buf = new T[buf_size];
  }catch(...){
    fprintf(stderr, "cannot new in %s\n", __func__);
    exit(1);
  }
  p = 0;
}
template <class T> BinaryWriter<T>::~BinaryWriter(){
  delete [] buf;
}
template <class T> void BinaryWriter<T>::setfd(int _fd){
  fd = _fd;
  assert(lseek(fd, 0, SEEK_SET) == 0);
  p = 0;
}
template <class T> inline size_t BinaryWriter<T>::bwwrite(const T* b, size_t count){
  assert(fd >= 0);
  for(size_t i = 0; i < count; i++){
    buf[p++] = b[i];
    if(p == buf_size){
      bwflush();
      p = 0;
    }
  }
  return count;
}
template <class T> inline size_t BinaryWriter<T>::bwflush(){
  if(p == 0)
    return 0;
  double st = gettime();
  const size_t unit = (1ll<<25);
  const size_t n = p / unit;
  const size_t m = p % unit;
  size_t ret = 0;
  {
    for(size_t i = 0; i < n; i++){
      ret += write(fd, buf + (i * unit), unit * sizeof(T)) / sizeof(T);
    }
    ret += write(fd, buf + (n * unit), m * sizeof(T)) / sizeof(T);
  }
  assert(ret == p);

  p = 0;
  iotime += gettime() - st;
  return ret;
}
/////







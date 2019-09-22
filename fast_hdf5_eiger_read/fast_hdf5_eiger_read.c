/* hello from 1990 - please compile with h5cc -o fast_hdf5_eiger_read.c  */

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "bitshuffle.h"

pthread_mutex_t hdf_mutex;

int n_jobs;
int job;

hid_t dataset, datasize;
hsize_t dims[3];

typedef struct chunk_t {
  char * chunk;
  size_t chunk_size;
  int index;
} chunk_t;

chunk_t next() {
  hsize_t offset[3];

  chunk_t retval;

  pthread_mutex_lock(&hdf_mutex);

  offset[0] = job;
  offset[1] = 0;
  offset[2] = 0;
  job += 1;

  if(job > n_jobs) {
    pthread_mutex_unlock(&hdf_mutex);
    retval.chunk_size = 0;
    retval.chunk = NULL;
    retval.index = 0;
    return retval;
  }

  hsize_t chunk_size;

  H5Dget_chunk_storage_size(dataset, offset, &chunk_size);
  retval.chunk_size = chunk_size;
  retval.index = offset[0];

  retval.chunk = (char *) malloc(chunk_size);

  uint32_t filter = 0;
  int status = H5DOread_chunk(dataset, H5P_DEFAULT, offset, &filter,
                              retval.chunk);
  pthread_mutex_unlock(&hdf_mutex);
  return retval;
}

void* worker(void * nonsense) {
  int size = dims[1] * dims[2];
  char* buffer = (char *) malloc(datasize * size);
  int32_t * longbuffer = (int32_t *) buffer;
  int16_t * shortbuffer = (int16_t *) buffer;

  while(1) {
    chunk_t chunk = next();
    if (chunk.chunk_size == 0) {
      free(buffer);
      return NULL;
    }
    bshuf_decompress_lz4((chunk.chunk)+12, (void *) buffer, size, datasize, 0);

    int64_t total = 0;

    for (size_t j = 0; j < size; j++) {
      if (datasize == 2) {
        if (shortbuffer[j] > 0) total += shortbuffer[j];
      } else {
        if (longbuffer[j] > 0) total += longbuffer[j];
      }
    }
    printf("Total counts for frame %d %lld\n", chunk.index, total);
    free(chunk.chunk);
  }
  return NULL;
}

int main(int argc,
         char ** argv)
{
  hid_t file, plist, space, datatype;
  H5D_layout_t layout;
  herr_t status;

  int * buffer = 0;
  char * chunk_buffer = 0;

  file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen(file, argv[2], H5P_DEFAULT);
  plist = H5Dget_create_plist(dataset);
  datatype = H5Dget_type(dataset);
  datasize = H5Tget_size(datatype);
  layout = H5Pget_layout(plist);

  printf("Element size %lld\n", datasize);

  if (layout != H5D_CHUNKED) {
    printf("You will not go to space, sorry\n");
    H5Pclose(plist);
    H5Dclose(dataset);
    H5Fclose(file);
    return 1;
  }

  space = H5Dget_space(dataset);

  printf("N dimensions %d\n", H5Sget_simple_extent_ndims(space));

  H5Sget_simple_extent_dims(space, dims, NULL);

  for (int j = 0; j < 3; j++) {
    printf("Dimension %d: %lld\n", j, dims[j]);
  }

  size_t size = dims[1] * dims[2];

  job = 0;
  n_jobs = dims[0];

  /* allocate and spin up threads */

  int n_threads = atoi(argv[3]);
  pthread_t * threads;

  pthread_mutex_init(&hdf_mutex, NULL);

  threads = (pthread_t *) malloc(sizeof(pthread_t) * n_threads);

  for (int j = 0; j < n_threads; j++) {
    pthread_create(&threads[j], NULL, worker, NULL);
  }

  for (int j = 0; j < n_threads; j++) {
    pthread_join(threads[j], NULL);
  }

  pthread_mutex_destroy(&hdf_mutex);

  H5Sclose(space);
  status = H5Pclose(plist);
  status = H5Dclose(dataset);
  status = H5Fclose(file);

  return 0;
}

/* hello from 1990 - please compile with h5cc -o fast_hdf5_eiger_read.c  */

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>
#include "bitshuffle.h"

/* global variables - necessary until I move the whole shebang into a class */

pthread_mutex_t hdf_mutex;

int n_jobs;
int job;

hid_t dataset, datasize;
hsize_t dims[3];

/* data structures to manage the work */

typedef struct chunk_t {
  char * chunk;
  size_t size;
  int index;
} chunk_t;

/* next() method will lock and read the next chunk from the file as needed */

chunk_t next() {
  uint32_t filter = 0;
  hsize_t offset[3], size;
  chunk_t chunk;

  pthread_mutex_lock(&hdf_mutex);

  offset[0] = job;
  offset[1] = 0;
  offset[2] = 0;
  job += 1;

  if(job > n_jobs) {
    pthread_mutex_unlock(&hdf_mutex);
    chunk.size = 0;
    chunk.chunk = NULL;
    chunk.index = 0;
    return chunk;
  }

  H5Dget_chunk_storage_size(dataset, offset, &size);
  chunk.index = offset[0];
  chunk.size = (int) size;
  chunk.chunk = (char *) malloc(size);

  H5DOread_chunk(dataset, H5P_DEFAULT, offset, &filter, chunk.chunk);
  pthread_mutex_unlock(&hdf_mutex);
  return chunk;
}

/* data structure to store the results */

pthread_mutex_t result_mutex;
uint64_t * histogram;
uint64_t histogram_bins;

void save_result(uint64_t * result) {
  pthread_mutex_lock(&result_mutex);

  for (int i = 0; i < histogram_bins; i++) {
    histogram[i] += result[i];
  }

  pthread_mutex_unlock(&result_mutex);
}

/* stupid interface to worker thread as it must be passed a void * and return
   the same */

void* worker(void* nonsense) {

  /* allocate frame storage - uses global information which is configured
     before threads spawned */

  int size = dims[1] * dims[2];
  char* buffer = (char *) malloc(datasize * size);

  /* allocate a buffer for the histogram */
  uint64_t * hbuffer = (uint64_t *) malloc(sizeof(uint64_t) * histogram_bins);

  /* alias the buffer as a short and long also, for ease of access */

  uint16_t * shortbuffer = (uint16_t *) buffer;
  unsigned char * charbuffer = (unsigned char *) buffer;

  /* while there is work to do, do work */

  while(1) {
    chunk_t chunk = next();
    if (chunk.size == 0) {
      free(buffer);
      free(hbuffer);
      return NULL;
    }

    /* reset histogram */
    for (int i = 0; i < histogram_bins; i++) {
      hbuffer[i] = 0;
    }

    /* decompress chunk - which starts 12 bytes in... */
    bshuf_decompress_lz4((chunk.chunk)+12, (void *) buffer, size, datasize, 0);

    /* accumulate histogram */

    for (size_t j = 0; j < size; j++) {
      if (datasize == 1) {
        hbuffer[charbuffer[j]] += 1;
      } else if (datasize == 2) {
        hbuffer[shortbuffer[j]] += 1;
      }
    }
    save_result(hbuffer);

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

  file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen(file, argv[2], H5P_DEFAULT);
  plist = H5Dget_create_plist(dataset);
  datatype = H5Dget_type(dataset);
  datasize = H5Tget_size(datatype);
  layout = H5Pget_layout(plist);

  printf("Element size %lld\n", datasize);

  if (datasize > 2) {
    printf("Not histogramming > 16 bit data\n");
    H5Pclose(plist);
    H5Dclose(dataset);
    H5Fclose(file);
    return 1;
  }

  if (layout != H5D_CHUNKED) {
    printf("You will not go to space, sorry\n");
    H5Pclose(plist);
    H5Dclose(dataset);
    H5Fclose(file);
    return 1;
  }

  if (datasize == 1) {
    histogram_bins = 1 + UCHAR_MAX;
  } else if (datasize == 2) {
    histogram_bins = 1 + USHRT_MAX;
  }

  histogram = (uint64_t *) malloc(sizeof(uint64_t) * histogram_bins);

  space = H5Dget_space(dataset);

  printf("N dimensions %d\n", H5Sget_simple_extent_ndims(space));

  H5Sget_simple_extent_dims(space, dims, NULL);

  for (int j = 0; j < 3; j++) {
    printf("Dimension %d: %lld\n", j, dims[j]);
  }

  /* set up global variables */

  job = 0;
  n_jobs = dims[0];

  /* useful static values
     UCHAR_MAX
     USHRT_MAX
     ULONG_MAX - I would not go bigger than here... */

  /* create a buffer here for the histogram */

  /* create a mutex for the histogram buffer */

  /* init with 0's */

  /* allocate and spin up threads */

  int n_threads = 1;

  if (argc > 3) {
    n_threads = atoi(argv[3]);
  }

  pthread_t * threads;

  pthread_mutex_init(&hdf_mutex, NULL);
  pthread_mutex_init(&result_mutex, NULL);

  threads = (pthread_t *) malloc(sizeof(pthread_t) * n_threads);

  for (int j = 0; j < n_threads; j++) {
    pthread_create(&threads[j], NULL, worker, NULL);
  }

  for (int j = 0; j < n_threads; j++) {
    pthread_join(threads[j], NULL);
  }

  pthread_mutex_destroy(&hdf_mutex);
  pthread_mutex_destroy(&result_mutex);

  /* print the results */
  for (int j = 0; j < histogram_bins; j++) {
    printf("Count %d total %lld\n", j, histogram[j]);
  }

  free(histogram);

  H5Sclose(space);
  status = H5Pclose(plist);
  status = H5Dclose(dataset);
  status = H5Fclose(file);

  return 0;
}

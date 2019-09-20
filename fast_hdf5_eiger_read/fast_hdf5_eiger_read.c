/* hello from 1990 - please compile with h5cc -o fast_hdf5_eiger_read.c  */

#include <time.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>
#include "bitshuffle.h"

int decompress(char * buffer, void * dest, size_t size, size_t element)
{
  // the data are 12 bytes in...
  return bshuf_decompress_lz4(buffer+12, (void *) dest, size, element, 0);
}
  
int main(int argc,
	 char ** argv)
{
  hid_t file, dataset, plist, space, memory, datatype, datasize;
  H5D_layout_t layout;
  herr_t status;

  hsize_t dims[3], offset[3], block[3];
  
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
  
  buffer = (int *) malloc (sizeof(int) * dims[1] * dims[2]);

  block[0] = 1;
  block[1] = dims[1];
  block[2] = dims[2];
  offset[0] = 0;
  offset[1] = 0;
  offset[2] = 0; 

  uint64_t total = 0;
  size_t size = dims[1] * dims[2];

  
  void * frame_buffer = malloc (datasize * size);

  double read_time = 0.0;
  double decompress_time = 0.0;

  clock_t start, end;
  
  for (int j = 0; j < dims[0]; j++) {
    hsize_t chunk_size;
    total += chunk_size;
    offset[0] = j;

    start = clock();
    H5Dget_chunk_storage_size(dataset, offset, &chunk_size);

    chunk_buffer = (char *) malloc (chunk_size);

    uint32_t filter = 0;
    status = H5DOread_chunk(dataset, H5P_DEFAULT, offset,
			    &filter, chunk_buffer);
    end = clock();
    read_time += (double) (end - start) / CLOCKS_PER_SEC;

    start = clock();
    decompress(chunk_buffer, frame_buffer, size, datasize);
    end = clock();
    decompress_time += (double) (end - start) / CLOCKS_PER_SEC;
    
    free(chunk_buffer);
    chunk_buffer = 0;
  }

  free(frame_buffer);
  
  printf("Read %lld images %lld bytes\n", dims[0], total);
  printf("Total read time: %.2fs\n", read_time);
  printf("Total decompress time: %.2fs\n", decompress_time);
  
  if (chunk_buffer) free(chunk_buffer);
  if (buffer) free(buffer);
  
  H5Sclose(space);
  status = H5Pclose(plist);
  status = H5Dclose(dataset);
  status = H5Fclose(file);

  return 0;
}

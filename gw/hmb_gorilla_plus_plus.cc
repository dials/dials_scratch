#include "dials/util/thread_pool.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>

struct Decompress {

  const char *buffer_;
  std::size_t size_;

  Decompress(const char *buffer, std::size_t size)
      : buffer_(buffer), size_(size) {}

  void operator()() const { std::cout << "Decompress" << std::endl; }
};

int main(int argc, char **argv) {
  std::size_t nthreads = 8;
  dials::util::ThreadPool pool(nthreads);

  hid_t file, dataset, plist, space, memory;
  H5D_layout_t layout;
  herr_t status;

  hsize_t dims[3], offset[3], block[3];

  int *buffer = 0;
  char *chunk_buffer = 0;

  file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen(file, argv[2], H5P_DEFAULT);
  plist = H5Dget_create_plist(dataset);
  layout = H5Pget_layout(plist);

  space = H5Dget_space(dataset);

  printf("N dimensions %d\n", H5Sget_simple_extent_ndims(space));

  H5Sget_simple_extent_dims(space, dims, NULL);

  for (int j = 0; j < 3; j++) {
    printf("Dimension %d: %d\n", j, dims[j]);
  }

  buffer = (int *)malloc(sizeof(int) * dims[1] * dims[2]);

  block[0] = 1;
  block[1] = dims[1];
  block[2] = dims[2];
  offset[0] = 0;
  offset[1] = 0;
  offset[2] = 0;

  status =
      H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, block, NULL);
  memory = H5Screate_simple(3, block, NULL);

  // status = H5Dread(dataset, H5T_NATIVE_INT, memory, space, H5P_DEFAULT,
  // buffer);

  for (int j = 0; j < dims[0]; j++) {
    hsize_t chunk_size;
    offset[0] = j;

    H5Dget_chunk_storage_size(dataset, offset, &chunk_size);

    printf("Chunk %d size: %d\n", j, chunk_size);

    chunk_buffer = (char *)malloc(chunk_size);

    uint32_t filter;
    H5DOread_chunk(dataset, H5P_DEFAULT, offset, &filter, chunk_buffer);

    pool.post(Decompress(chunk_buffer, chunk_size));

    free(chunk_buffer);
    chunk_buffer = 0;
  }

  pool.wait();

  switch (layout) {
  case H5D_COMPACT:
    printf("H5D_COMPACT\n");
    break;
  case H5D_CONTIGUOUS:
    printf("H5D_CONTIGUOUS\n");
    break;
  case H5D_CHUNKED:
    printf("H5D_CHUNKED\n");
    break;
  }

  if (chunk_buffer)
    free(chunk_buffer);
  if (buffer)
    free(buffer);

  H5Sclose(space);
  status = H5Pclose(plist);
  status = H5Dclose(dataset);
  status = H5Fclose(file);

  return 0;
}

/* hello from 1990 - please compile with h5cc hmb_gorilla.cc  */

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
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

  buffer = (int *)malloc(sizeof(int) * dims[1] * dims[2]);

  block[0] = 1;
  block[1] = dims[1];
  block[2] = dims[2];
  offset[0] = 0;
  offset[1] = 0;
  offset[2] = 0;

  uint64_t total = 0;

  for (int j = 0; j < dims[0]; j++) {
    hsize_t chunk_size;
    total += chunk_size;
    offset[0] = j;

    H5Dget_chunk_storage_size(dataset, offset, &chunk_size);

    chunk_buffer = (char *)malloc(chunk_size);

    uint32_t filter = 0;
    status =
        H5DOread_chunk(dataset, H5P_DEFAULT, offset, &filter, chunk_buffer);
    free(chunk_buffer);
    chunk_buffer = 0;
  }

  printf("Read %lld images %lld bytes\n", dims[0], total);

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

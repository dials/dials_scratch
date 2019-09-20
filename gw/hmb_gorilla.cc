#include <hdf5.h>
#include <stdio.h>

int main(int argc,
	 char ** argv)
{
  hid_t file, dataset, plist;
  H5D_layout_t layout;
  herr_t status;
  hsize_t dims[2];
  
  int * buffer;
  
  file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen(file, argv[2], H5P_DEFAULT);
  plist = H5Dget_create_plist(dataset);
  layout = H5Pget_layout(plist);
  
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

  status = H5Pclose(plist);
  status = H5Dclose(dataset);
  status = H5Fclose(file);

  return 0;
}

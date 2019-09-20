#include <hdf5.h>
#include <iostream>

int main(int argc,
	 char ** argv)
{
  hid_t file, dataset, layout;
  herr_t status;
  hsize_t dims[2];
  
  int * buffer;
  
  file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen(file, "/entry/data", H5P_DEFAULT);
  layout = H5Dget_create_plist(dataset);

  switch (layout) {
  case H5D_COMPACT:
    std::cout << "H5D_COMPACT" << std::endl;
    break;
  case H5D_CONTIGUOUS:
    std::cout << "H5D_CONTIGUOUS" << std::endl;
    break;
  case H5D_CHUNKED:
    std::cout << "H5D_CHUNKED" << std::endl;
  }

  status = H5Pclose(layout);
  status = H5Dclose(dataset);
  status = H5Fclose(file);

  return 0;
}

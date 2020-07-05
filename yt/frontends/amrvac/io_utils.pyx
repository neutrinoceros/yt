
# distutils: libraries = STD_LIBS
# distutils: include_dirs = LIB_DIR
cimport numpy as np
import numpy as np
import struct

def get_single_block_field_data(
    file_obj, # missing a type here (this is supposed to be a buffer)
    int block_offset,
    int field_per_block_size,
    int field_idx,
    np.ndarray[np.int64_t, ndim=1] output_shape,
    ):
    file_obj.seek(block_offset)

    # compute byte size of a single field
    fmt = "=" + field_per_block_size * 'd'

    #byte_size_field = struct.calcsize(fmt) # clean
    byte_size_field = 8 * field_per_block_size # fast

    # '1' means seek forward (as opposed to "from start of file")
    file_obj.seek(byte_size_field * field_idx, 1) 

    data = struct.unpack(fmt, file_obj.read(byte_size_field))
    data = np.reshape(data, output_shape, order='F')
    return data
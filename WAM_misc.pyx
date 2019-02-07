#from numpy cimport ndarray
import numpy as np
cimport numpy as np					# "cimport" is used to import special compile-time information about the numpy model(numpy.pxd)
							# .pxd files contain Cython declarations which are only meant for inclusion by Cython modules.
							# A pxd file is imported into a pyx module by using the cimport keyword.
cimport cython
from cython cimport parallel
import scipy.io as sio

DTYPE=np.float64
ctypedef np.float64_t DTYPE_t				# "ctypedef" assigns corresponding compile-time type to DTYPE_t.
							# For every type in the numpy module there's a corresponding compile-time type with a _t-suffix.
							# When to use DTYPE or DTYPE_t?

@cython.boundscheck(False)				# turn off bounds-checking for entire function
@cython.wraparound(False)				# turn off negative index wrapping for entire function
							# more explaination about these two options: http://cython.readthedocs.io/en/latest/src/tutorial/numpy.html

def std_mult_3d(object[DTYPE_t,ndim=3] buf1 not None,
                object[DTYPE_t,ndim=3] buf2 not None,
                object[DTYPE_t,ndim=3] output = None,
                int num_threads=2):
    cdef unsigned int x,y,z,inner,middle,outer

    if output is None:
        output = np.empty_like(buf1)
    outer=buf1.shape[0]
    middle=buf1.shape[1]
    inner=buf1.shape[2]

    with nogil:
        for x in parallel.prange(outer,num_threads=num_threads):
            for y in xrange(middle):
                for z in xrange(inner):
                    output[x,y,z]=buf1[x,y,z]*buf2[x,y,z]
    return output

#%% create empty array for track and time

def create_empty_array(count_time,divt,latitude,longitude,yearpart,years,datapathea):
# LG: variable $datapathea is not passed into this function. However, it can be used here?
    Sa_track_top  = np.zeros( (np.int(count_time*divt)+1,len(latitude),len(longitude)) )
    Sa_track_down = np.zeros( (np.int(count_time*divt)+1,len(latitude),len(longitude)) )
    Sa_time_top   = np.zeros( (np.int(count_time*divt)+1,len(latitude),len(longitude)) )
    Sa_time_down  = np.zeros( (np.int(count_time*divt)+1,len(latitude),len(longitude)) )
    Sa_dist_top   = np.zeros( (np.int(count_time*divt)+1,len(latitude),len(longitude)) )
    Sa_dist_down  = np.zeros( (np.int(count_time*divt)+1,len(latitude),len(longitude)) )

    if yearpart[0] == 359:
        sio.savemat( datapathea[0], {'Sa_track_top':Sa_track_top,'Sa_track_down':Sa_track_down}, do_compression=True )
        sio.savemat( datapathea[1], {'Sa_time_top':Sa_time_top  ,'Sa_time_down':Sa_time_down,'Sa_dist_top':Sa_dist_top,'Sa_dist_down':Sa_dist_down}  , do_compression=True )
    else:
        sio.savemat( datapathea[2], {'Sa_track_top':Sa_track_top,'Sa_track_down':Sa_track_down}, do_compression=True )
        sio.savemat( datapathea[3], {'Sa_time_top':Sa_time_top  ,'Sa_time_down':Sa_time_down,'Sa_dist_top':Sa_dist_top,'Sa_dist_down':Sa_dist_down}  , do_compression=True )

    return

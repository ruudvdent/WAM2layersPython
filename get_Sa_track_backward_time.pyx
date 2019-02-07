@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)						# False, Cython will adjust the remainder and quotient operators C types to match those of Python ints;
								# and raise a ZeroDivisionError when the right operand is 0.
								# http://cython.readthedocs.io/en/latest/src/reference/compilation.html?highlight=cdivision

cdef double fluxes_boundary(
    double flux_in1,
    double flux_in2,
    int pos) nogil:
    cdef double boundary

    boundary = 0.5*(flux_in1 + flux_in2)
    if boundary < 0:
        pos = pos-1
    else:
        pos = pos
    return boundary*pos

@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def get_Sa_track_backward_TIME(
    np.ndarray[np.float32_t,ndim=1] latitude not None,		# parameters and C-variables declared as an Extension type, may take the value None
    np.ndarray[np.float32_t,ndim=1] longitude not None,		# Cython does not check this for reasons of efficiency
    int count_time,						# not None can only be used in Python functions (declared with def not cdef)
    int divt,
    int timestep,
    int Kvf,
    np.ndarray[np.float32_t,ndim=2] Region,
    np.ndarray[np.float32_t,ndim=3] Fa_E_top not None,
    np.ndarray[np.float32_t,ndim=3] Fa_N_top not None,
    np.ndarray[np.float32_t,ndim=3] Fa_E_down not None,
    np.ndarray[np.float32_t,ndim=3] Fa_N_down not None,
    np.ndarray[np.float32_t,ndim=3] Fa_Vert not None,
    np.ndarray[np.float32_t,ndim=3] E not None,
    np.ndarray[np.float32_t,ndim=3] P,
    np.ndarray[np.float32_t,ndim=3] W_top not None,
    np.ndarray[np.float32_t,ndim=3] W_down not None,
    np.ndarray[np.float32_t,ndim=3] Sa_track_top_last not None,
    np.ndarray[np.float32_t,ndim=3] Sa_track_down_last not None,
    np.ndarray[np.float32_t,ndim=3] Sa_time_top_last not None,
    np.ndarray[np.float32_t,ndim=3] Sa_time_down_last not None,
    np.ndarray[np.float32_t,ndim=3] Sa_dist_top_last not None,
    np.ndarray[np.float32_t,ndim=3] Sa_dist_down_last not None,
    np.ndarray[np.float32_t,ndim=1] L_N_gridcell not None,
    np.ndarray[np.float32_t,ndim=1] L_S_gridcell not None,
    float L_EW_gridcell,
    int isglobal,
    int num_threads):

    cdef int x, y, z, outer, middle, inner, t, ti, nt, nx, ny

    # make P_region matrix
    nt = P.shape[0]
    ny = len(latitude)
    nx = len(longitude)
    cdef np.ndarray[np.float32_t,ndim=3] Region3D # to improve array lookups and assignments, type the contents of the narray objects
    cdef np.ndarray[np.float32_t,ndim=3] P_region # by using a "buffer" syntax which must be told the datatype and number of dimensions.

    Region3D = np.tile(np.reshape(Region,[1,ny,nx]),[nt,1,1])
    P_region = Region3D*P

    # Total moisture in the column
    cdef np.ndarray[np.float32_t,ndim=3] W
    W = W_top + W_down

    # separate the direction of the vertical flux and make it absolute
    cdef np.ndarray[np.float32_t,ndim=3] Fa_upward = np.zeros([nt,ny,nx],dtype=np.float32)
    Fa_upward[Fa_Vert <= 0 ]	= -1.0*Fa_Vert[Fa_Vert <= 0 ]
    cdef np.ndarray[np.float32_t,ndim=3] Fa_downward = np.zeros([nt,ny,nx],dtype=np.float32)
    Fa_downward[Fa_Vert >= 0 ]	= Fa_Vert[Fa_Vert >= 0 ]
    #Fa_upward			= np.abs(Fa_upward)

    # include the vertical dispersion
    if Kvf == 0:
        pass		# do nothing
    else:
        Fa_upward			= (1.+Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0]		= Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_downward			= (1.+Kvf) * Fa_downward
        Fa_downward[Fa_Vert <= 0]	= np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf

    # defining size of output
    cdef np.ndarray[np.float32_t,ndim=3] Sa_track_down = np.zeros(np.shape(W_down),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_track_top = np.zeros(np.shape(W_top),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_time_down = np.zeros(np.shape(W_down),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_time_top = np.zeros(np.shape(W_top),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_dist_down = np.zeros(np.shape(W_down),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_dist_top = np.zeros(np.shape(W_top),dtype=np.float32)

    # assign begin values of output == last (but first index) values of the previous time slot
    Sa_track_down[nt,:,:] = Sa_track_down_last
    Sa_track_top[nt,:,:] = Sa_track_top_last
    Sa_time_down[nt,:,:] = Sa_time_down_last
    Sa_time_top[nt,:,:] = Sa_time_top_last
    Sa_dist_down[nt,:,:] = Sa_dist_down_last
    Sa_dist_top[nt,:,:] = Sa_dist_top_last

    # defining sizes of tracked moisture
    cdef np.ndarray[np.float32_t,ndim=3] Sa_track_after_Fa_P_E_down = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_track_after_Fa_P_E_top = np.zeros([1,ny,nx],dtype=np.float32)

    cdef np.ndarray[np.float32_t,ndim=3] Sa_track_after_all_down = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_track_after_all_top = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_time_after_all_down = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_time_after_all_top = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_dist_after_all_down = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_dist_after_all_top = np.zeros([1,ny,nx],dtype=np.float32)

    cdef float my_Fa_E_top_WE,my_Fa_E_top_EW,my_Fa_W_top_WE,my_Fa_W_top_EW,my_Fa_N_top_SN,my_Fa_N_top_NS,my_Fa_S_top_SN,my_Fa_S_top_NS, \
    my_Fa_E_down_WE,my_Fa_E_down_EW,my_Fa_W_down_WE,my_Fa_W_down_EW,my_Fa_N_down_SN,my_Fa_N_down_NS,my_Fa_S_down_SN,my_Fa_S_down_NS

    cdef float Sa_track_after_Fa_down,Sa_track_after_Fa_top

    # define sizes of total moisture

    # define variables that find out what happens to the water
    cdef np.ndarray[np.float32_t,ndim=3] north_loss = np.zeros((np.int(count_time*divt),1,nx),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] south_loss = np.zeros((np.int(count_time*divt),1,nx),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] down_to_top = np.zeros(np.shape(P),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] top_to_down = np.zeros(np.shape(P),dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] water_lost = np.zeros(np.shape(P),dtype=np.float32)
#   cdef np.ndarray[np.float_t,ndim=3] water_lost_down = np.zeros(np.shape(P))
#   cdef np.ndarray[np.float_t,ndim=3] water_lost_top = np.zeros(np.shape(P))

    # defining sizes of timed moisture (1,107,240) [s]
    cdef np.ndarray[np.float32_t,ndim=3] Sa_time_after_Fa_P_E_down = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_time_after_Fa_P_E_top = np.zeros([1,ny,nx],dtype=np.float32)
    cdef float Sa_time_after_Fa_down,Sa_time_after_Fa_top

    # defining sizes of moisture distance (1,107,240) [m]
    cdef np.ndarray[np.float32_t,ndim=3] Sa_dist_after_Fa_P_E_down = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] Sa_dist_after_Fa_P_E_top = np.zeros([1,ny,nx],dtype=np.float32)
    cdef float Sa_dist_after_Fa_down,Sa_dist_after_Fa_top

    cdef np.ndarray[np.float32_t,ndim=3] di = np.zeros([1,ny,nx],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=3] xi
    cdef np.ndarray[np.float32_t,ndim=2] xi2d
    cdef float yi
    cdef int W_idx,E_idx

    ti = timestep/divt

    # delta_x is the mean of north and south boundaries [1,107,240]
    L_SN_gridcell = 0.5*(L_N_gridcell+L_S_gridcell)
    xi2d = np.tile(L_SN_gridcell,[1,nx])
    xi = np.reshape(xi2d,[1,ny,nx])
    # delta_y [1]
    yi = L_EW_gridcell

    # Sa calculation backward in time (96 ~ 1)
    for t in np.arange(np.int(count_time*divt),0,-1):
                # losses to the north and south
        	# delta_d is weighted mean of delta_x and delta_y [1,107,240]
        	# Fa_E and Fa_N can be either positive or negative
                #with nogil:
                    #for x in parallel.prange(outer,num_threads=num_threads):
        for y in range(ny):
            for x in range(nx):
                di[0,y,x] = (xi[0,y,x]*abs(Fa_E_down[t-1,y,x])+yi*abs(Fa_N_down[t-1,y,x]))/(abs(Fa_E_down[t-1,y,x])+abs(Fa_N_down[t-1,y,x]))
        	#di = 0.5*(xi+yi)

        for y in parallel.prange(1,ny-1,nogil=True,num_threads=num_threads,schedule='static'):
            for x in range(nx):
                if x == 0:
                    W_idx = nx-1
                    E_idx = x+1
                elif x == nx-1:
                    W_idx = x-1
                    E_idx = 0
                else:
                    W_idx = x-1
                    E_idx = x+1

                my_Fa_E_top_WE	= fluxes_boundary(Fa_E_top[t-1,y,x] ,Fa_E_top[t-1,y,E_idx] ,1)
                my_Fa_E_top_EW	= fluxes_boundary(Fa_E_top[t-1,y,x] ,Fa_E_top[t-1,y,E_idx] ,0)
                my_Fa_E_down_WE	= fluxes_boundary(Fa_E_down[t-1,y,x],Fa_E_down[t-1,y,E_idx],1)
                my_Fa_E_down_EW	= fluxes_boundary(Fa_E_down[t-1,y,x],Fa_E_down[t-1,y,E_idx],0)

                my_Fa_W_top_WE	= fluxes_boundary(Fa_E_top[t-1,y,x] ,Fa_E_top[t-1,y,W_idx] ,1)
                my_Fa_W_top_EW	= fluxes_boundary(Fa_E_top[t-1,y,x] ,Fa_E_top[t-1,y,W_idx] ,0)
                my_Fa_W_down_WE = fluxes_boundary(Fa_E_down[t-1,y,x],Fa_E_down[t-1,y,W_idx],1)
                my_Fa_W_down_EW = fluxes_boundary(Fa_E_down[t-1,y,x],Fa_E_down[t-1,y,W_idx],0)

                my_Fa_N_top_SN  = fluxes_boundary(Fa_N_top[t-1,y,x] ,Fa_N_top[t-1,y-1,x] ,1)
                my_Fa_N_top_NS  = fluxes_boundary(Fa_N_top[t-1,y,x] ,Fa_N_top[t-1,y-1,x] ,0)
                my_Fa_N_down_SN = fluxes_boundary(Fa_N_down[t-1,y,x],Fa_N_down[t-1,y-1,x],1)
                my_Fa_N_down_NS = fluxes_boundary(Fa_N_down[t-1,y,x],Fa_N_down[t-1,y-1,x],0)

                my_Fa_S_top_SN  = fluxes_boundary(Fa_N_top[t-1,y,x] ,Fa_N_top[t-1,y+1,x] ,1)
                my_Fa_S_top_NS  = fluxes_boundary(Fa_N_top[t-1,y,x] ,Fa_N_top[t-1,y+1,x] ,0)
                my_Fa_S_down_SN = fluxes_boundary(Fa_N_down[t-1,y,x],Fa_N_down[t-1,y+1,x],1)
                my_Fa_S_down_NS = fluxes_boundary(Fa_N_down[t-1,y,x],Fa_N_down[t-1,y+1,x],0)

                if y == 1:
                    north_loss[t-1,0,x]	= (my_Fa_N_top_NS*(Sa_track_top[t,y,x]/W_top[t,y,x])+my_Fa_N_down_NS*(Sa_track_down[t,y,x]/W_down[t,y,x]))
                if y == ny-2:
                    south_loss[t-1,0,x]	= (my_Fa_S_top_SN * (Sa_track_top[t,y,x] / W_top[t,y,x])+ my_Fa_S_down_SN * (Sa_track_down[t,y,x] / W_down[t,y,x]))

                # down: calculate with moisture fluxes
                Sa_track_after_Fa_down = (Sa_track_down[t,y,x]
                                          # Sa_E_down[0] is assign from W_down[t], so it should not be 0, however, 0 is in the first and last row
                                          # nan is caused by Sa_E_down[0], in other words, by W_down[t] which is 0 on the denominator
                                          # Need to solve this in Fluxes_and_Storage, before coming back and resuming diagnosis here (30/04/2017).
                                          + my_Fa_E_down_WE * (Sa_track_down[t,y,E_idx] / W_down[t,y,E_idx])
                                          - my_Fa_E_down_EW * (Sa_track_down[t,y,x]     / W_down[t,y,x])
                                          - my_Fa_W_down_WE * (Sa_track_down[t,y,x]     / W_down[t,y,x])
                                          + my_Fa_W_down_EW * (Sa_track_down[t,y,W_idx] / W_down[t,y,W_idx])
                                          + my_Fa_N_down_SN * (Sa_track_down[t,y-1,x]   / W_down[t,y-1,x])
                                          - my_Fa_N_down_NS * (Sa_track_down[t,y,x]     / W_down[t,y,x])
                                          - my_Fa_S_down_SN * (Sa_track_down[t,y,x]     / W_down[t,y,x])
                                          + my_Fa_S_down_NS * (Sa_track_down[t,y+1,x]   / W_down[t,y+1,x])
                                          - Fa_downward[t-1,y,x]  * (Sa_track_down[t,y,x]   / W_down[t,y,x])
                                          + Fa_upward[t-1,y,x]    * (Sa_track_top[t,y,x]    / W_top[t,y,x])
                                         )
                # top: calculate with moisture fluxes
                Sa_track_after_Fa_top = (Sa_track_top[t,y,x]
                                         + my_Fa_E_top_WE * (Sa_track_top[t,y,E_idx] / W_top[t,y,E_idx])
                                         - my_Fa_E_top_EW * (Sa_track_top[t,y,x]     / W_top[t,y,x])
                                         - my_Fa_W_top_WE * (Sa_track_top[t,y,x]     / W_top[t,y,x])
                                         + my_Fa_W_top_EW * (Sa_track_top[t,y,W_idx] / W_top[t,y,W_idx])
                                         + my_Fa_N_top_SN * (Sa_track_top[t,y-1,x]   / W_top[t,y-1,x])
                                         - my_Fa_N_top_NS * (Sa_track_top[t,y,x]     / W_top[t,y,x])
                                         - my_Fa_S_top_SN * (Sa_track_top[t,y,x]     / W_top[t,y,x])
                                         + my_Fa_S_top_NS * (Sa_track_top[t,y+1,x]   / W_top[t,y+1,x])
                                         + Fa_downward[t-1,y,x] * (Sa_track_down[t,y,x]  / W_down[t,y,x])
                                         - Fa_upward[t-1,y,x] * (Sa_track_top[t,y,x]     / W_top[t,y,x])
				        )
                # down: add precipitation and subtract evaporation
                Sa_track_after_Fa_P_E_down[0,y,x]   = (Sa_track_after_Fa_down
                                                    +P_region[t-1,y,x]*(W_down[t,y,x]/W[t,y,x])
                                                    -E[t-1,y,x]*(Sa_track_down[t,y,x]/W_down[t,y,x]))
                # top: add precipitation
                Sa_track_after_Fa_P_E_top[0,y,x]    = (Sa_track_after_Fa_top
                                                    +P_region[t-1,y,x]*(W_top[t,y,x]/W[t,y,x]))


                Sa_time_after_Fa_down = ((Sa_track_down[t,y,x]  * (ti +  Sa_time_down[t,y,x]) +
                                         my_Fa_E_down_WE * (ti +  Sa_time_down[t,y,E_idx]) *
                                         (Sa_track_down[t,y,E_idx] / W_down[t,y,E_idx]) - my_Fa_E_down_EW *
                                         (ti +  Sa_time_down[t,y,x]) * (Sa_track_down[t,y,x] / W_down[t,y,x]) - my_Fa_W_down_WE *
                                         (ti +  Sa_time_down[t,y,x]) * (Sa_track_down[t,y,x] / W_down[t,y,x]) + my_Fa_W_down_EW *
                                         (ti +  Sa_time_down[t,y,W_idx]) * (Sa_track_down[t,y,W_idx] / W_down[t,y,W_idx]) + my_Fa_N_down_SN *
                                         (ti +  Sa_time_down[t,y-1,x]) * (Sa_track_down[t,y-1,x] / W_down[t,y-1,x]) - my_Fa_N_down_NS *
                                         (ti +  Sa_time_down[t,y,x])   * (Sa_track_down[t,y,x] / W_down[t,y,x]) - my_Fa_S_down_SN *
                                         (ti +  Sa_time_down[t,y,x]) * (Sa_track_down[t,y,x]   / W_down[t,y,x]) + my_Fa_S_down_NS *
                                         (ti +  Sa_time_down[t,y+1,x]) * (Sa_track_down[t,y+1,x] / W_down[t,y+1,x]) - Fa_downward[t-1,y,x]  *
                                         (ti +  Sa_time_down[t,y,x])   * (Sa_track_down[t,y,x]  / W_down[t,y,x]) + Fa_upward[t-1,y,x] *
                                         (ti +  Sa_time_top[t,y,x])    * (Sa_track_top[t,y,x]   / W_top[t,y,x])) / Sa_track_after_Fa_down )
                # if Sa_track_after_Fa_down == 0, then Sa_time_after_Fa_down == inf.
                # Note here, 'inf' is number, but 'nan' literally means 'not a number'.
                # However, operation inf - inf yields a nan.
		# Following if statement does not work without GIL
                #if np.isnan(Sa_time_after_Fa_down):
                if Sa_time_after_Fa_down != Sa_time_after_Fa_down:
                    Sa_time_after_Fa_down = 0

                Sa_time_after_Fa_top = ((Sa_track_top[t,y,x] * (ti +  Sa_time_top[t,y,x]) +
                                        my_Fa_E_top_WE * (ti +  Sa_time_top[t,y,E_idx]) *
                                        (Sa_track_top[t,y,E_idx] / W_top[t,y,E_idx]) - my_Fa_E_top_EW *
                                        (ti +  Sa_time_top[t,y,x]) * (Sa_track_top[t,y,x] / W_top[t,y,x]) - my_Fa_W_top_WE *
                                        (ti +  Sa_time_top[t,y,x]) * (Sa_track_top[t,y,x] / W_top[t,y,x]) + my_Fa_W_top_EW *
                                        (ti +  Sa_time_top[t,y,W_idx]) * (Sa_track_top[t,y,W_idx] / W_top[t,y,W_idx]) + my_Fa_N_top_SN *
                                        (ti +  Sa_time_top[t,y-1,x]) * (Sa_track_top[t,y-1,x] / W_top[t,y-1,x]) - my_Fa_N_top_NS *
                                        (ti +  Sa_time_top[t,y,x]) * (Sa_track_top[t,y,x] / W_top[t,y,x]) - my_Fa_S_top_SN *
                                        (ti +  Sa_time_top[t,y,x]) * (Sa_track_top[t,y,x] / W_top[t,y,x]) + my_Fa_S_top_NS *
                                        (ti +  Sa_time_top[t,y+1,x]) * (Sa_track_top[t,y+1,x] / W_top[t,y+1,x]) + Fa_downward[t-1,y,x] *
                                        (ti +  Sa_time_down[t,y,x]) * (Sa_track_down[t,y,x] / W_down[t,y,x]) - Fa_upward[t-1,y,x] *
                                        (ti +  Sa_time_top[t,y,x]) * (Sa_track_top[t,y,x] / W_top[t,y,x])) / Sa_track_after_Fa_top)
                #if np.isnan(Sa_time_after_Fa_top):
                if Sa_time_after_Fa_top != Sa_time_after_Fa_top:
                    Sa_time_after_Fa_top = 0

                Sa_time_after_Fa_P_E_down[0,y,x] = ((Sa_track_after_Fa_down*Sa_time_after_Fa_down
                                                    +P_region[t-1,y,x]*ti/2.*(W_down[t,y,x]/W[t,y,x])
                                                    -E[t-1,y,x]*(ti+Sa_time_down[t,y,x])*(Sa_track_down[t,y,x]/W_down[t,y,x]))
                                                    /Sa_track_after_Fa_P_E_down[0,y,x])
                Sa_time_after_Fa_P_E_top[0,y,x] = ((Sa_track_after_Fa_top*Sa_time_after_Fa_top
                                                   +P_region[t-1,y,x]*ti/2*(W_top[t,y,x]/W[t,y,x]))
                                                   /Sa_track_after_Fa_P_E_top[0,y,x])
		    

                Sa_dist_after_Fa_down = ((Sa_track_down[t,y,x]*(               Sa_dist_down[t,y,x])
	                                 +my_Fa_E_down_WE * (xi[0,y,x] + Sa_dist_down[t,y,E_idx]) * (Sa_track_down[t,y,E_idx] / W_down[t,y,E_idx])
                                         -my_Fa_E_down_EW * (               Sa_dist_down[t,y,x])   * (Sa_track_down[t,y,x]   / W_down[t,y,x])
                                         -my_Fa_W_down_WE * (               Sa_dist_down[t,y,x])   * (Sa_track_down[t,y,x]   / W_down[t,y,x])
                                         +my_Fa_W_down_EW * (xi[0,y,x] + Sa_dist_down[t,y,W_idx]) * (Sa_track_down[t,y,W_idx] / W_down[t,y,W_idx])
                                         +my_Fa_N_down_SN * (yi           + Sa_dist_down[t,y-1,x]) * (Sa_track_down[t,y-1,x] / W_down[t,y-1,x])
                                         -my_Fa_N_down_NS * (               Sa_dist_down[t,y,x])   * (Sa_track_down[t,y,x]   / W_down[t,y,x])
                                         -my_Fa_S_down_SN * (               Sa_dist_down[t,y,x])   * (Sa_track_down[t,y,x]   / W_down[t,y,x])
                                         +my_Fa_S_down_NS * (yi           + Sa_dist_down[t,y+1,x]) * (Sa_track_down[t,y+1,x] / W_down[t,y+1,x])
                                         -Fa_downward[t-1,y,x]  * (               Sa_dist_down[t,y,x])   * (Sa_track_down[t,y,x]   / W_down[t,y,x])
                                         +Fa_upward[t-1,y,x]    * (               Sa_dist_top[t,y,x])    * (Sa_track_top[t,y,x]    / W_top[t,y,x]))
                                         /Sa_track_after_Fa_down)
                #if np.isnan(Sa_dist_after_Fa_down):
                if Sa_dist_after_Fa_down != Sa_dist_after_Fa_down:
                    Sa_dist_after_Fa_down = 0

                Sa_dist_after_Fa_top = ((Sa_track_top[t,y,x]  * (               Sa_dist_top[t,y,x])
                                        + my_Fa_E_top_WE * (xi[0,y,x] + Sa_dist_top[t,y,E_idx]) * (Sa_track_top[t,y,E_idx] / W_top[t,y,E_idx])
                                        - my_Fa_E_top_EW * (               Sa_dist_top[t,y,x])   * (Sa_track_top[t,y,x]   / W_top[t,y,x])
                                        - my_Fa_W_top_WE * (               Sa_dist_top[t,y,x])   * (Sa_track_top[t,y,x]   / W_top[t,y,x])
                                        + my_Fa_W_top_EW * (xi[0,y,x] + Sa_dist_top[t,y,W_idx]) * (Sa_track_top[t,y,W_idx] / W_top[t,y,W_idx])
                                        + my_Fa_N_top_SN * (yi           + Sa_dist_top[t,y-1,x]) * (Sa_track_top[t,y-1,x] / W_top[t,y-1,x])
                                        - my_Fa_N_top_NS * (               Sa_dist_top[t,y,x])   * (Sa_track_top[t,y,x]   / W_top[t,y,x])
                                        - my_Fa_S_top_SN * (               Sa_dist_top[t,y,x])   * (Sa_track_top[t,y,x]   / W_top[t,y,x])
                                        + my_Fa_S_top_NS * (yi           + Sa_dist_top[t,y+1,x]) * (Sa_track_top[t,y+1,x] / W_top[t,y+1,x])
                                        + Fa_downward[t-1,y,x] * (               Sa_dist_down[t,y,x])  * (Sa_track_down[t,y,x]  / W_down[t,y,x])
                                        - Fa_upward[t-1,y,x]   * (               Sa_dist_top[t,y,x])   * (Sa_track_top[t,y,x]   / W_top[t,y,x])
                                        ) / Sa_track_after_Fa_top)
                #if np.isnan(Sa_dist_after_Fa_top):
                if Sa_dist_after_Fa_top != Sa_dist_after_Fa_top:
                    Sa_dist_after_Fa_top = 0

                Sa_dist_after_Fa_P_E_down[0,y,x] = ((Sa_track_after_Fa_down * Sa_dist_after_Fa_down
        	                                    - E[t-1,y,x] * (Sa_dist_down[t,y,x]) * (Sa_track_down[t,y,x] / W_down[t,y,x])
                                                    ) / Sa_track_after_Fa_P_E_down[0,y,x])
                Sa_dist_after_Fa_P_E_top[0,y,x] = Sa_dist_after_Fa_top

        where_are_NaNs = np.isnan(Sa_time_after_Fa_P_E_down)
        Sa_time_after_Fa_P_E_down[where_are_NaNs] = 0
        where_are_NaNs = np.isnan(Sa_time_after_Fa_P_E_top)
        Sa_time_after_Fa_P_E_top[where_are_NaNs] = 0
        where_are_NaNs = np.isnan(Sa_dist_after_Fa_P_E_down)
        Sa_dist_after_Fa_P_E_down[where_are_NaNs] = 0

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t-1,:,:] = np.reshape(np.maximum(0,np.reshape(Sa_track_after_Fa_P_E_down,(np.size(Sa_track_after_Fa_P_E_down)))-np.reshape(W_down[t-1,:,:],(np.size(W_down[t-1,:,:])))),(ny,nx))
        top_to_down[t-1,:,:] = np.reshape(np.maximum(0,np.reshape(Sa_track_after_Fa_P_E_top,(np.size(Sa_track_after_Fa_P_E_top)))-np.reshape(W_top[t-1,:,:],(np.size(W_top[t-1,:,:])))),(ny,nx))


        for y in parallel.prange(ny,nogil=True,schedule='static',num_threads=num_threads):
            for x in range(nx):
                Sa_track_after_all_down[0,y,x] = Sa_track_after_Fa_P_E_down[0,y,x]-down_to_top[t-1,y,x]+top_to_down[t-1,y,x]
                Sa_track_after_all_top[0,y,x] = Sa_track_after_Fa_P_E_top[0,y,x]-top_to_down[t-1,y,x]+down_to_top[t-1,y,x]

                Sa_time_after_all_down[0,y,x] = ( (Sa_track_after_Fa_P_E_down[0,y,x] * Sa_time_after_Fa_P_E_down[0,y,x]
                                                 - down_to_top[t-1,y,x] * Sa_time_after_Fa_P_E_down[0,y,x]
                                                 + top_to_down[t-1,y,x] * Sa_time_after_Fa_P_E_top[0,y,x]
                                                 ) / Sa_track_after_all_down[0,y,x] )
                Sa_time_after_all_top[0,y,x] = ( (Sa_track_after_Fa_P_E_top[0,y,x] * Sa_time_after_Fa_P_E_top[0,y,x]
                                                - top_to_down[t-1,y,x] * Sa_time_after_Fa_P_E_top[0,y,x]
                                                + down_to_top[t-1,y,x] * Sa_time_after_Fa_P_E_down[0,y,x]
                                                ) / Sa_track_after_all_top[0,y,x] )

                Sa_dist_after_all_down[0,y,x] = ( (Sa_track_after_Fa_P_E_down[0,y,x] * Sa_dist_after_Fa_P_E_down[0,y,x]
                                                 - down_to_top[t-1,y,x] * Sa_dist_after_Fa_P_E_down[0,y,x]
                                                 + top_to_down[t-1,y,x] * Sa_dist_after_Fa_P_E_top[0,y,x]
                                                 ) / Sa_track_after_all_down[0,y,x] )
                Sa_dist_after_all_top[0,y,x] = ( (Sa_track_after_Fa_P_E_top[0,y,x] * Sa_dist_after_Fa_P_E_top[0,y,x]
                                                - top_to_down[t-1,y,x] * Sa_dist_after_Fa_P_E_top[0,y,x]
                                                + down_to_top[t-1,y,x] * Sa_dist_after_Fa_P_E_down[0,y,x]
                                                ) / Sa_track_after_all_top[0,y,x])

        where_are_NaNs = np.isnan(Sa_time_after_all_down)
        Sa_time_after_all_down[where_are_NaNs] = 0
        where_are_NaNs = np.isnan(Sa_time_after_all_top)
        Sa_time_after_all_top[where_are_NaNs] = 0
        where_are_NaNs = np.isnan(Sa_dist_after_all_down)
        Sa_dist_after_all_down[where_are_NaNs] = 0
        where_are_NaNs = np.isnan(Sa_dist_after_all_top)
        Sa_dist_after_all_top[where_are_NaNs] = 0

        # down and top: water lost to the system:
        #water_lost_down[t-1,:,:]	= np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down))) - np.reshape(W_down[t-1,:,:],
        #                                    (np.size(W_down[t-1,:,:])))), (ny,nx))
        #water_lost_top[t-1,:,:]	= np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top))) - np.reshape(W_top[t-1,:,:],
        #                                    (np.size(W_top[t-1,:,:])))), (ny,nx))
        water_lost[t-1,:,:]		= np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down))) - np.reshape(W_down[t-1,:,:],
                                            (np.size(W_down[t-1,:,:])))), (ny,nx)) + np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top))) - np.reshape(W_top[t-1,:,:],
                                            (np.size(W_top[t-1,:,:])))), (ny,nx)) #water_lost_down + water_lost_top

        # down: determine Sa_region of this next timestep 100% stable
        Sa_track_down[t-1,1:ny-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_down[t-1,1:ny-1,:], np.size(W_down[t-1,1:ny-1,:])), np.reshape(Sa_track_after_all_down[0,1:ny-1,:],
                                                np.size(Sa_track_after_all_down[0,1:ny-1,:])))), (len(latitude[1:ny-1]),nx))
        # top: determine Sa_region of this next timestep 100% stable
        Sa_track_top[t-1,1:ny-1,:]  = np.reshape(np.maximum(0,np.minimum(np.reshape(W_top[t-1,1:ny-1,:] , np.size(W_top[t-1,1:ny-1,:])) , np.reshape(Sa_track_after_all_top[0,1:ny-1,:] ,
                                                np.size(Sa_track_after_all_top[0,1:ny-1,:])))) , (len(latitude[1:ny-1]),nx))

        # down: determine Sa_region of this next timestep 100% stable
        Sa_time_down[t-1,1:ny-1,:] = Sa_time_after_all_down[0,1:ny-1,:]
        # top: determine Sa_region of this next timestep 100% stable
        Sa_time_top[t-1,1:ny-1,:]  = Sa_time_after_all_top[0,1:ny-1,:]

        # down: determine Sa_region of this next timestep 100% stable
        Sa_dist_down[t-1,1:ny-1,:] = Sa_dist_after_all_down[0,1:ny-1,:]
        # top: determine Sa_region of this next timestep 100% stable
        Sa_dist_top[t-1,1:ny-1,:]  = Sa_dist_after_all_top[0,1:ny-1,:]

    return Sa_dist_top,Sa_dist_down,Sa_time_top,Sa_time_down,Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost

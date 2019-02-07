@cython.boundscheck(False)
def get_Sa_track_backward(np.ndarray[np.float_t,ndim=1] latitude not None,
			  np.ndarray[np.float_t,ndim=1] longitude not None,
			  int count_time,
			  int divt,
			  int Kvf,
			  np.ndarray[np.float_t,ndim=2] Region not None,
			  np.ndarray[np.float_t,ndim=3] Fa_E_top not None,
			  np.ndarray[np.float_t,ndim=3] Fa_N_top not None,
			  np.ndarray[np.float_t,ndim=3] Fa_E_down not None,
			  np.ndarray[np.float_t,ndim=3] Fa_N_down not None,
			  np.ndarray[np.float_t,ndim=3] Fa_Vert not None,
			  np.ndarray[np.float_t,ndim=3] E not None,
			  np.ndarray[np.float_t,ndim=3] P not None,
			  np.ndarray[np.float_t,ndim=3] W_top not None,
			  np.ndarray[np.float_t,ndim=3] W_down not None,
			  np.ndarray[np.float_t,ndim=3] Sa_track_top_last not None,
			  np.ndarray[np.float_t,ndim=3] Sa_track_down_last not None,
			  int isglobal, int num_threads=2):

    cdef unsigned int x,y,z,inner,middle,outer
    cdef np.int Pmax=P.shape[0]
    cdef np.int nlat=len(latitude)
    cdef np.int nlon=len(longitude)
    cdef np.int Fa_Vert_s0 = Fa_Vert.shape[0]
    cdef np.int Fa_Vert_s1 = Fa_Vert.shape[1]
    cdef np.int Fa_Vert_s2 = Fa_Vert.shape[2]

    # make P_region matrix
    cdef np.ndarray[np.float_t,ndim=3] Region3D = np.tile(np.reshape(Region,[1,len(latitude),len(longitude)]),[len(P[:,0,0]),1,1])
    cdef np.ndarray[np.float_t,ndim=3] P_region = Region3D * P

    # Total moisture in the column
    cdef np.ndarray[np.float_t,ndim=3] W = W_top + W_down

    # separate the direction of the vertical flux and make it absolute
    cdef np.ndarray[np.float_t,ndim=3] Fa_upward = np.zeros(np.shape(Fa_Vert))
    Fa_upward[Fa_Vert <= 0 ] = Fa_Vert[Fa_Vert <= 0 ]
    cdef np.ndarray[np.float_t,ndim=3] Fa_downward = np.zeros(np.shape(Fa_Vert));
    Fa_downward[Fa_Vert >= 0 ] = Fa_Vert[Fa_Vert >= 0 ]
    Fa_upward = np.abs(Fa_upward)

    # include the vertical dispersion
    # LG: waiting for revisit
    if Kvf == 0:
        pass				# do nothing
    else:
        Fa_upward			= (1.+Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0]		= Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_downward			= (1.+Kvf) * Fa_downward
        Fa_downward[Fa_Vert <= 0]	= np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf

    # define the horizontal fluxes over the boundaries
    # fluxes over the eastern boundary
    cdef np.ndarray[np.float_t,ndim=3] Fa_E_top_boundary = np.zeros(np.shape(Fa_E_top))
    cdef np.ndarray[np.float_t,ndim=3] Fa_E_down_boundary = np.zeros(np.shape(Fa_E_down))
    outer = Fa_E_top_boundary.shape[0]
    middle = Fa_E_top_boundary.shape[1]
    inner = Fa_E_top_boundary.shape[2]
    with nogil:
        for x in parallel.prange(outer,num_threads=num_threads):
            for y in range(middle):
                for z in range(0,inner-1):
                    Fa_E_top_boundary[x,y,z] = 0.5 * (Fa_E_top[x,y,z] + Fa_E_top[x,y,z+1])
                    Fa_E_down_boundary[x,y,z] = 0.5 * (Fa_E_down[x,y,z] + Fa_E_down[x,y,z+1])
                if isglobal == 1:
                    Fa_E_top_boundary[x,y,inner] = 0.5 * (Fa_E_top[x,y,inner] + Fa_E_top[x,y,0])
                    Fa_E_down_boundary[x,y,inner] = 0.5 * (Fa_E_down[x,y,inner] + Fa_E_down[x,y,0])

    #Fa_E_top_boundary[:,:,:-1]		= 0.5 * (Fa_E_top[:,:,:-1] + Fa_E_top[:,:,1:])
    #Fa_E_down_boundary[:,:,:-1]		= 0.5 * (Fa_E_down[:,:,:-1] + Fa_E_down[:,:,1:])

    # find out where the positive and negative fluxes are
    cdef np.ndarray[np.float_t,ndim=3] Fa_E_top_pos = np.ones(np.shape(Fa_E_top))
    cdef np.ndarray[np.float_t,ndim=3] Fa_E_down_pos = np.ones(np.shape(Fa_E_down))
    Fa_E_top_pos[Fa_E_top_boundary < 0]		= 0
    Fa_E_down_pos[Fa_E_down_boundary < 0]	= 0
    cdef np.ndarray[np.float_t,ndim=3] Fa_E_top_neg = Fa_E_top_pos - 1
    cdef np.ndarray[np.float_t,ndim=3] Fa_E_down_neg = Fa_E_down_pos - 1

    # separate directions west-east (all positive numbers)
#    Fa_E_top_WE		= Fa_E_top_boundary  * Fa_E_top_pos
    cdef np.ndarray[np.float_t,ndim=3] Fa_E_top_WE, Fa_E_top_EW, Fa_E_down_WE, Fa_E_down_EW
    Fa_E_top_WE = std_mult_3d(Fa_E_top_boundary,Fa_E_top_pos)
#    Fa_E_top_EW		= Fa_E_top_boundary  * Fa_E_top_neg
    Fa_E_top_EW = std_mult_3d(Fa_E_top_boundary,Fa_E_top_neg)
#    Fa_E_down_WE	= Fa_E_down_boundary * Fa_E_down_pos
    Fa_E_down_WE = std_mult_3d(Fa_E_down_boundary,Fa_E_down_pos)
#    Fa_E_down_EW	= Fa_E_down_boundary * Fa_E_down_neg
    Fa_E_down_EW = std_mult_3d(Fa_E_down_boundary,Fa_E_down_neg)

    # fluxes over the western boundary
    cdef np.ndarray[np.float_t,ndim=3] Fa_W_top_WE = np.nan*np.zeros(np.shape(P))
    Fa_W_top_WE[:,:,1:]	= Fa_E_top_WE[:,:,:-1]
    Fa_W_top_WE[:,:,0]	= Fa_E_top_WE[:,:,-1]
    cdef np.ndarray[np.float_t,ndim=3] Fa_W_top_EW = np.nan*np.zeros(np.shape(P))
    Fa_W_top_EW[:,:,1:]	= Fa_E_top_EW[:,:,:-1]
    Fa_W_top_EW[:,:,0]	= Fa_E_top_EW[:,:,-1]
    cdef np.ndarray[np.float_t,ndim=3] Fa_W_down_WE = np.nan*np.zeros(np.shape(P))
    Fa_W_down_WE[:,:,1:]= Fa_E_down_WE[:,:,:-1]
    Fa_W_down_WE[:,:,0]	= Fa_E_down_WE[:,:,-1]
    cdef np.ndarray[np.float_t,ndim=3] Fa_W_down_EW = np.nan*np.zeros(np.shape(P))
    Fa_W_down_EW[:,:,1:]= Fa_E_down_EW[:,:,:-1]
    Fa_W_down_EW[:,:,0]	= Fa_E_down_EW[:,:,-1]

    # fluxes over the northern boundary
    cdef np.ndarray[np.float_t,ndim=3] Fa_N_top_boundary = np.nan*np.zeros(np.shape(Fa_N_top))
    cdef np.ndarray[np.float_t,ndim=3] Fa_N_down_boundary = np.nan*np.zeros(np.shape(Fa_N_down))
    outer = Fa_N_top_boundary.shape[0]
    middle = Fa_N_top_boundary.shape[1]
    inner = Fa_N_top_boundary.shape[2]
    with nogil:
        for x in parallel.prange(outer,num_threads=num_threads):
            for y in range(1,middle):
                for z in range(inner):
                    Fa_N_top_boundary[x,y,z]	= 0.5 * ( Fa_N_top[x,y-1,z] + Fa_N_top[x,y,z] )
                    Fa_N_down_boundary[x,y,z]	= 0.5 * ( Fa_N_down[x,y-1,z] + Fa_N_down[x,y,z] )

    # find out where the positive and negative fluxes are
    cdef np.ndarray[np.float_t,ndim=3] Fa_N_top_pos = np.ones(np.shape(Fa_N_top))
    cdef np.ndarray[np.float_t,ndim=3] Fa_N_down_pos = np.ones(np.shape(Fa_N_down))
    cdef np.ndarray[np.float_t,ndim=3] Fa_N_top_neg, Fa_N_down_neg
    Fa_N_top_pos[Fa_N_top_boundary < 0]		= 0
    Fa_N_down_pos[Fa_N_down_boundary < 0]	= 0
    Fa_N_top_neg				= Fa_N_top_pos - 1
    Fa_N_down_neg				= Fa_N_down_pos - 1

    # separate directions south-north (all positive numbers)
#    Fa_N_top_SN		= Fa_N_top_boundary  * Fa_N_top_pos
#    Fa_N_top_NS		= Fa_N_top_boundary  * Fa_N_top_neg
#    Fa_N_down_SN	= Fa_N_down_boundary * Fa_N_down_pos
#    Fa_N_down_NS	= Fa_N_down_boundary * Fa_N_down_neg
    cdef np.ndarray[np.float_t,ndim=3] Fa_N_top_SN,Fa_N_top_NS,Fa_N_down_SN,Fa_N_down_NS
    Fa_N_top_SN = std_mult_3d(Fa_N_top_boundary,Fa_N_top_pos)
    Fa_N_top_NS = std_mult_3d(Fa_N_top_boundary,Fa_N_top_neg)
    Fa_N_down_SN = std_mult_3d(Fa_N_down_boundary,Fa_N_down_pos)
    Fa_N_down_NS = std_mult_3d(Fa_N_down_boundary,Fa_N_down_neg)

    # fluxes over the southern boundary
    cdef np.ndarray[np.float_t,ndim=3] Fa_S_top_SN     		= np.nan*np.zeros(np.shape(P))
    Fa_S_top_SN[:,:-1,:]	= Fa_N_top_SN[:,1:,:]
    cdef np.ndarray[np.float_t,ndim=3] Fa_S_top_NS	      	= np.nan*np.zeros(np.shape(P))
    Fa_S_top_NS[:,:-1,:]	= Fa_N_top_NS[:,1:,:]
    cdef np.ndarray[np.float_t,ndim=3] Fa_S_down_SN		= np.nan*np.zeros(np.shape(P))
    Fa_S_down_SN[:,:-1,:]	= Fa_N_down_SN[:,1:,:]
    cdef np.ndarray[np.float_t,ndim=3] Fa_S_down_NS		= np.nan*np.zeros(np.shape(P))
    Fa_S_down_NS[:,:-1,:]	= Fa_N_down_NS[:,1:,:]

    # defining size of output
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_down	= np.zeros(np.shape(W_down))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_top	= np.zeros(np.shape(W_top))

    # assign begin values of output == last (but first index) values of the previous time slot
    Sa_track_down[-1,:,:]	= Sa_track_down_last
    Sa_track_top[-1,:,:]	= Sa_track_top_last

    # defining sizes of tracked moisture
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_after_Fa_down	= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_after_Fa_P_E_down	= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_E_down		= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_W_down		= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_N_down		= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_S_down		= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_after_Fa_top	= np.zeros(np.shape(Sa_track_top_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_after_Fa_P_E_top	= np.zeros(np.shape(Sa_track_top_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_E_top		= np.zeros(np.shape(Sa_track_top_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_W_top		= np.zeros(np.shape(Sa_track_top_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_N_top		= np.zeros(np.shape(Sa_track_top_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_track_S_top		= np.zeros(np.shape(Sa_track_top_last))

    # define sizes of total moisture
    cdef np.ndarray[np.float_t,ndim=3] Sa_E_down	= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_W_down	= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_N_down	= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_S_down	= np.zeros(np.shape(Sa_track_down_last))
    cdef np.ndarray[np.float_t,ndim=3] Sa_E_top	= np.zeros(np.shape(Sa_track_top_last) )
    cdef np.ndarray[np.float_t,ndim=3] Sa_W_top	= np.zeros(np.shape(Sa_track_top_last) )
    cdef np.ndarray[np.float_t,ndim=3] Sa_N_top	= np.zeros(np.shape(Sa_track_top_last) )
    cdef np.ndarray[np.float_t,ndim=3] Sa_S_top	= np.zeros(np.shape(Sa_track_top_last) )

    # define variables that find out what happens to the water
    cdef np.ndarray[np.float_t,ndim=3] north_loss		= np.zeros((np.int(count_time*divt),1,len(longitude)))
    cdef np.ndarray[np.float_t,ndim=3] south_loss		= np.zeros((np.int(count_time*divt),1,len(longitude)))
    cdef np.ndarray[np.float_t,ndim=3] down_to_top		= np.zeros(np.shape(P))
    cdef np.ndarray[np.float_t,ndim=3] top_to_down		= np.zeros(np.shape(P))
    cdef np.ndarray[np.float_t,ndim=3] water_lost		= np.zeros(np.shape(P))
    cdef np.ndarray[np.float_t,ndim=3] water_lost_down	        = np.zeros(np.shape(P))
    cdef np.ndarray[np.float_t,ndim=3] water_lost_top	        = np.zeros(np.shape(P))

    cdef np.ndarray[np.float_t,ndim=3] Sa_track_after_all_down, Sa_track_after_all_top

    # Sa calculation backward in time
    cdef int t
    for t in np.arange(np.int(count_time*divt),0,-1):

        # down: define values of total moisture
        Sa_E_down[0,:,:-1]	= W_down[t,:,1:]	# Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_E_down[0,:,-1]	= W_down[t,:,0]		# Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0,:,1:]	= W_down[t,:,:-1]	# Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_W_down[0,:,0]	= W_down[t,:,-1]	# Atmospheric storage of the cell to the west [m3]
        Sa_N_down[0,1:,:]	= W_down[t,0:-1,:]	# Atmospheric storage of the cell to the north [m3]
        Sa_S_down[0,:-1,:]	= W_down[t,1:,:]	# Atmospheric storage of the cell to the south [m3]

        # top: define values of total moisture
        Sa_E_top[0,:,:-1]	= W_top[t,:,1:]		# Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_E_top[0,:,-1]	= W_top[t,:,0]		# Atmospheric storage of the cell to the east [m3]
        Sa_W_top[0,:,1:]	= W_top[t,:,:-1]	# Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_W_top[0,:,0]		= W_top[t,:,-1]		# Atmospheric storage of the cell to the west [m3]
        Sa_N_top[0,1:,:]	= W_top[t,:-1,:]	# Atmospheric storage of the cell to the north [m3]
        Sa_S_top[0,:-1,:]	= W_top[t,1:,:]		# Atmospheric storage of the cell to the south [m3]

         # down: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_down[0,:,:-1]	= Sa_track_down[t,:,1:]	# Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_down[0,:,-1]	= Sa_track_down[t,:,0]	# Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0,:,1:]		= Sa_track_down[t,:,:-1]# Atmospheric storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_down[0,:,0]	= Sa_track_down[t,:,-1]	# Atmospheric storage of the cell to the west [m3]
        Sa_track_N_down[0,1:,:]		= Sa_track_down[t,:-1,:]# Atmospheric storage of the cell to the north [m3]
        Sa_track_S_down[0,:-1,:]	= Sa_track_down[t,1:,:]	# Atmospheric storage of the cell to the south [m3]

        # down: calculate with moisture fluxes
        outer = Sa_track_after_Fa_down.shape[1]
        inner = Sa_track_after_Fa_down.shape[2]
        with nogil:
            for x in parallel.prange(1,outer-1,num_threads=num_threads):
                for y in range(inner):
                    Sa_track_after_Fa_down[0,x,y] = (Sa_track_down[t,x,y]
                                                     + Fa_E_down_WE[t-1,x,y] * (Sa_track_E_down[0,x,y] / Sa_E_down[0,x,y])
                                                     - Fa_E_down_EW[t-1,x,y] * (Sa_track_down[t,x,y] / W_down[t,x,y])
                                                     - Fa_W_down_WE[t-1,x,y] * (Sa_track_down[t,x,y] / W_down[t,x,y])
                                                     + Fa_W_down_EW[t-1,x,y] * (Sa_track_W_down[0,x,y] / Sa_W_down[0,x,y])
                                                     + Fa_N_down_SN[t-1,x,y] * (Sa_track_N_down[0,x,y] / Sa_N_down[0,x,y])
                                                     - Fa_N_down_NS[t-1,x,y] * (Sa_track_down[t,x,y] / W_down[t,x,y])
                                                     - Fa_S_down_SN[t-1,x,y] * (Sa_track_down[t,x,y] / W_down[t,x,y])
                                                     + Fa_S_down_NS[t-1,x,y] * (Sa_track_S_down[0,x,y] / Sa_S_down[0,x,y])
                                                     - Fa_downward[t-1,x,y] * (Sa_track_down[t,x,y] / W_down[t,x,y])
                                                     + Fa_upward[t-1,x,y] * (Sa_track_top[t,x,y] / W_top[t,x,y]))

        # top: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_top[0,:,:-1]		= Sa_track_top[t,:,1:]	# Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_top[0,:,-1]	= Sa_track_top[t,:,0]	# Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_top[0,:,1:]		= Sa_track_top[t,:,:-1]	# Atmospheric tracked storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_top[0,:,0]	= Sa_track_top[t,:,-1]	# Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_top[0,1:,:]		= Sa_track_top[t,:-1,:]	# Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_top[0,:-1,:]		= Sa_track_top[t,1:,:]	# Atmospheric tracked storage of the cell to the south [m3]

        # top: calculate with moisture fluxes
        outer = Sa_track_after_Fa_top.shape[1]
        inner = Sa_track_after_Fa_top.shape[2]
        with nogil:
            for x in parallel.prange(1,outer-1,num_threads=num_threads):
                for y in range(inner):
                    Sa_track_after_Fa_top[0,x,y] = (Sa_track_top[t,x,y]
                                                       + Fa_E_top_WE[t-1,x,y] * (Sa_track_E_top[0,x,y] / Sa_E_top[0,x,y])
                                                       - Fa_E_top_EW[t-1,x,y] * (Sa_track_top[t,x,y] / W_top[t,x,y])
                                                       - Fa_W_top_WE[t-1,x,y] * (Sa_track_top[t,x,y] / W_top[t,x,y])
                                                       + Fa_W_top_EW[t-1,x,y] * (Sa_track_W_top[0,x,y] / Sa_W_top[0,x,y])
                                                       + Fa_N_top_SN[t-1,x,y] * (Sa_track_N_top[0,x,y] / Sa_N_top[0,x,y])
                                                       - Fa_N_top_NS[t-1,x,y] * (Sa_track_top[t,x,y]/ W_top[t,x,y])
                                                       - Fa_S_top_SN[t-1,x,y] * (Sa_track_top[t,x,y] / W_top[t,x,y])
                                                       + Fa_S_top_NS[t-1,x,y] * (Sa_track_S_top[0,x,y] / Sa_S_top[0,x,y])
                                                       + Fa_downward[t-1,x,y] * (Sa_track_down[t,x,y] / W_down[t,x,y])
                                                       - Fa_upward[t-1,x,y] * (Sa_track_top[t,x,y] / W_top[t,x,y]))

        # losses to the north and south
        inner = north_loss.shape[2]
        with nogil:
            for x in parallel.prange(inner,num_threads=num_threads):
                north_loss[t-1,0,x] = (Fa_N_top_NS[t-1,1,x]  * (Sa_track_top[t,1,x]  / W_top[t,1,x])  + Fa_N_down_NS[t-1,1,x]  * (Sa_track_down[t,1,x]  / W_down[t,1,x]) )
                south_loss[t-1,0,x] = (Fa_S_top_SN[t-1,-2,x] * (Sa_track_top[t,-2,x] / W_top[t,-2,x]) + Fa_S_down_SN[t-1,-2,x] * (Sa_track_down[t,-2,x] / W_down[t,-2,x]))

        # down: add precipitation and subtract evaporation
        outer = Sa_track_after_Fa_P_E_down.shape[1]
        inner = Sa_track_after_Fa_P_E_down.shape[2]
        with nogil:
            for x in parallel.prange(1,outer-1,num_threads=num_threads):
                for y in range(inner):
                    Sa_track_after_Fa_P_E_down[0,x,y] = (Sa_track_after_Fa_down[0,x,y]
                                                            + P_region[t-1,x,y] * (W_down[t,x,y] / W[t,x,y])
                                                            - E[t-1,x,y] * (Sa_track_down[t,x,y] / W_down[t,x,y]))
                    # top: add precipitation
                    Sa_track_after_Fa_P_E_top[0,x,y] = (Sa_track_after_Fa_top[0,x,y]
                                                           + P_region[t-1,x,y] * (W_top[t,x,y] / W[t,x,y]))

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t-1,:,:]	= np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_down, (np.size(Sa_track_after_Fa_P_E_down))) - np.reshape(W_down[t-1,:,:],
                                            (np.size(W_down[t-1,:,:])))), (len(latitude),len(longitude)))
        top_to_down[t-1,:,:]	= np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_top, (np.size(Sa_track_after_Fa_P_E_top))) - np.reshape(W_top[t-1,:,:],
                                            (np.size(W_top[t-1,:,:])))), (len(latitude),len(longitude)))
        Sa_track_after_all_down	= Sa_track_after_Fa_P_E_down - down_to_top[t-1,:,:] + top_to_down[t-1,:,:]
        Sa_track_after_all_top	= Sa_track_after_Fa_P_E_top - top_to_down[t-1,:,:] + down_to_top[t-1,:,:]

        # down and top: water lost to the system:
        water_lost_down[t-1,:,:]= np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down))) - np.reshape(W_down[t-1,:,:],
                                            (np.size(W_down[t-1,:,:])))), (len(latitude),len(longitude)))
        water_lost_top[t-1,:,:]	= np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top))) - np.reshape(W_top[t-1,:,:],
                                            (np.size(W_top[t-1,:,:])))), (len(latitude),len(longitude)))
        water_lost		= water_lost_down + water_lost_top

        # down: determine Sa_region of this next timestep 100% stable
        Sa_track_down[t-1,1:-1,:]= np.reshape(np.maximum(0,np.minimum(np.reshape(W_down[t-1,1:-1,:], np.size(W_down[t-1,1:-1,:])), np.reshape(Sa_track_after_all_down[0,1:-1,:],
                                                np.size(Sa_track_after_all_down[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
        # top: determine Sa_region of this next timestep 100% stable
        Sa_track_top[t-1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_top[t-1,1:-1,:], np.size(W_top[t-1,1:-1,:])), np.reshape(Sa_track_after_all_top[0,1:-1,:],
                                                np.size(Sa_track_after_all_top[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))

    return Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost

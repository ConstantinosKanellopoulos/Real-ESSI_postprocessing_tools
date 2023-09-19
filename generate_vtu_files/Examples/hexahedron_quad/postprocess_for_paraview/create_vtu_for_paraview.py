# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# create_vtu_for_paraview.py : python script that reads Real-ESSI output files (.feioutput), filters (if requested) the displacements, and creates .vtu files for ParaView
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# authors: Constantinos Kanellopoulos and Giulia Aguzzi(2023)
# Chair of Structural Dynamics and Earthquake Engineering, Chair of Structural Mechanics
# ETH, Zurich
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#! /usr/bin/env python3
import h5py
import scipy as sp
import time
import sys
import linecache
import numpy as np
from numpy  import *
import os 
import os.path
import scipy.io as sio
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import xarray as xr
import pyevtk as pe #needs pip3 (python3)
import matplotlib.pyplot as plt
from scipy import signal, interpolate


pwd = os.getcwd(); 


# *******************************
# ********** FUNCTIONS **********
# *******************************
def EXTRACT_GENERAL_DATA_SIMU(filename, pwd, total_number_of_elements): 
    ################################
    ###  GENERAL INFO SIMULATION ###
    ################################
    # MASTER PROCESS - PROCESSOR 0 
    os.chdir('../')
    master_filename   = filename+'.h5.feioutput'
    master_outputfile = h5py.File(master_filename,'r');
    master_groups     = master_outputfile['/']
    #[print(item) for item in master_groups.items()]  # Print Datasets and Groups within the Master Group 
    master_node_dataset = master_outputfile['/Model/Nodes/Partition'] # in the Model/Nodes/ of the 0 proc there's only partition no other info 
    master_elem_dataset = master_outputfile['/Model/Elements/Partition']
    print(' ******** INFO SIMUALATION ********')
    time_vec         = master_groups['time'][:]
    # print(time_vec)
    ntsteps          = len(time_vec)
    # print(ntsteps)
    nproc            = master_outputfile['Number_of_Processes_Used'][:] # number of cores-1 where the 1 represents the master 0 proc 
    nproc            = nproc[0]
    # print(nproc)
    print('Total Number of Timesteps:', ntsteps)
    print('Total Number of Processors:', nproc)
    nnodes           = master_node_dataset.size - 1 ; # Because the first element refers to processor 0
    # print(nnodes)
    nodes_procid     = master_node_dataset[:]         # Vector that shows, for each single node in which processor it's located. The size is equal to nnodes+1
    # print(nodes_procid)
    # print(len(nodes_procid))
    ncoords          = 3;                                         
    # nelems           = master_elem_dataset.size - 1; 
    nelems           = total_number_of_elements                                      
    # print(nelems)
    elems_procid     = master_elem_dataset[:nelems+1]
    # print(elems_procid)
    # print(len(elems_procid))
    print('Total Number of Nodes in the Whole Domain:', nnodes)
    print('Total Number of Elements in the Whole Domain:', nelems)
    # PROCESSOR 1 
    proc1_filename   = filename+'.h5.'+str(1).zfill(len(str(nproc)))+'.feioutput'   #h5.01.feioutput, why to complicate it
    proc1_outputfile = h5py.File(proc1_filename,'r');
    ndofs            = max(proc1_outputfile['/Model/Nodes/Number_of_DOFs'][:])      #### DELETE, when iterating in each processor should get also the dofs of each node
    print('Total Number of DOFs per node:', ndofs)
    print(' ********************************')
    #Close file 
    master_outputfile.close()
    proc1_outputfile.close()
    print(' ******* GENERAL INFO EXTRACTED *******')
    os.chdir(pwd)
    return time_vec, ntsteps, nproc, nnodes, nodes_procid, nelems, elems_procid, ncoords, ndofs


def EXTRACT_COORD_WAVEFIELD(filename,pwd,time_vec, ntsteps, nproc, nnodes, nodes_procid, ncoords, ndofs): 
    ##############################################
    ###  EXTRACT COORDINATES AND DISPLACEMENTS ###
    ##############################################
    #start = time.time()
    # Initialise Output Matrixes 
    os.chdir('../')
    coord = np.zeros((nnodes, ncoords)) # size is equal to nodes_procid-1 this is why I use "node_tags-1" later on 
    # print(coord)
    # print(coord.shape)
    ux = np.zeros((nnodes, ntsteps)); 
    uy = np.zeros((nnodes, ntsteps));
    uz = np.zeros((nnodes, ntsteps));
    # print(ux)
    # print(ux.shape)
    rx = np.zeros((nnodes, ntsteps)); 
    ry = np.zeros((nnodes, ntsteps));
    rz = np.zeros((nnodes, ntsteps));

    number_of_dofs = np.zeros((nnodes+1, 1));  #the first element will be always 0, so that number_of_dofs[1] refers to the 1st node

    # LOOP ON THE NUMBER OF PROCESSORS 
    for proc in range(1, nproc): 
        start = time.time()
        node_tags               = np.asarray(np.where(nodes_procid==proc))  # Find the index of the nodes that are localised in the processor "proc" e.g. proc = 1, node_tags returns the index of all the nodes in proc 1
        # asarray: Input data, in any form that can be converted to an array. This includes lists, lists of tuples, tuples, tuples of tuples, tuples of lists and ndarrays
        # print(node_tags)
        # print(len(node_tags[0]))
        # Extract Coordinates 
        proc_filename           = filename+'.h5.'+str(proc).zfill(len(str(nproc)))+'.feioutput'
        proc_outputfile         = h5py.File(proc_filename,'r');
        proc_nodes_idx_to_coord = proc_outputfile['/Model/Nodes/Index_to_Coordinates'][:]
        proc_nodes_coord        = proc_outputfile['/Model/Nodes/Coordinates'][:]
        proc_idx_extract        = proc_nodes_idx_to_coord[node_tags[0,:]]
        # print(node_tags[0,:])
        # print(proc_idx_extract)
        # Store Coordinates in a matrix 
        coord[node_tags[0,:]-1,0] = proc_nodes_coord[proc_idx_extract]
        coord[node_tags[0,:]-1,1] = proc_nodes_coord[proc_idx_extract+1]
        coord[node_tags[0,:]-1,2] = proc_nodes_coord[proc_idx_extract+2]
        # print(coord)
        # Extract Displacements 
        proc_nodes_idx_to_displ = proc_outputfile['/Model/Nodes/Index_to_Generalized_Displacements'][:]
        proc_nodes_displ        = proc_outputfile['/Model/Nodes/Generalized_Displacements'][:]
        proc_idx_extract        = proc_nodes_idx_to_displ[node_tags[0,:]]
        # print(proc_idx_extract)
        # Store Displacement Field in Matrices 
        ux[node_tags[0,:]-1, :] = proc_nodes_displ[proc_idx_extract,:]
        uy[node_tags[0,:]-1, :] = proc_nodes_displ[proc_idx_extract+1,:]
        uz[node_tags[0,:]-1, :] = proc_nodes_displ[proc_idx_extract+2,:]
        # print(ux)
        #rx[node_tags[0,:]-1, :] = proc_nodes_displ[proc_idx_extract+3,:]
        #ry[node_tags[0,:]-1, :] = proc_nodes_displ[proc_idx_extract+4,:]
        #rz[node_tags[0,:]-1, :] = proc_nodes_displ[proc_idx_extract+5,:]

        # extract the dofs per node
        proc_Number_of_DOFs = proc_outputfile['/Model/Nodes/Number_of_DOFs'][:]
        number_of_dofs[node_tags[0,:],0] = proc_Number_of_DOFs[node_tags[0,:]]
        number_of_dofs = number_of_dofs.astype(int)
        # print(number_of_dofs)
        # print(len(number_of_dofs))

        #Close file 
        proc_outputfile.close()
        print(proc) #to see at which iteration we are...
        end = time.time()
        print("Time consumed in working: ",end - start)


    if enable_relative_displacements == 1:
        # print(ux.shape)
        dispx0 = ux[:,0] #take the first column which has the initial displacements of each node (coming from the previous stage)
        dispy0 = uy[:,0]
        dispz0 = uz[:,0]
        # print(dispx0)
        # print(dispx0.shape) # (245,)
        dispx0 = np.reshape(dispx0, (nnodes,1)) #convert it to a 2d array (245,1) to match the dimensions for broadcasting
        dispy0 = np.reshape(dispy0, (nnodes,1))
        dispz0 = np.reshape(dispz0, (nnodes,1))
        # print(dispx0.shape)
        # print(ux)
        ux = ux - dispx0     # subtract the initial displacements
        uy = uy - dispy0
        uz = uz - dispz0
        # print(ux)


    # Store as xarray
    # Store as xarray - takes time... and needs RAM
    x  = xr.DataArray(coord[:,0], dims=['npts']); 
    y  = xr.DataArray(coord[:,1], dims=['npts']); 
    z  = xr.DataArray(coord[:,2], dims=['npts']); 
    ux = xr.DataArray(ux, dims=['npts','ntime']);  
    uy = xr.DataArray(uy, dims=['npts','ntime']);  
    uz = xr.DataArray(uz, dims=['npts','ntime']);  
    #rx = xr.DataArray(rx, dims=['npts','ntime']);  
    #ry = xr.DataArray(ry, dims=['npts','ntime']);  
    #rz = xr.DataArray(rz, dims=['npts','ntime']);  
    print(' ******* COORDINATES & FIELD DATA EXTRACTED *******')
    os.chdir(pwd)
    return  x, y, z, ux, uy, uz, number_of_dofs   #, rx, ry, rz


def EXTRACT_CONN_MATRIX(filename, pwd, time_vec, ntsteps, nproc, nelems, elems_procid):
    ####################################
    ###  EXTRACT CONNECTIVITY MATRIX ###
    ####################################
    os.chdir('../')
    nnodes_per_elem = 8; 
    ConnMat = np.zeros((nelems, nnodes_per_elem), dtype=int);
    # LOOP ON THE NUMBER OF PROCESSORS 
    for proc in range(1, nproc): 
        start = time.time()
        elems_tags              = np.asarray(np.where(elems_procid==proc))  # Find the index of the elements that are localised in the processor "proc" e.g. proc = 1, elems_tag returns the index of all the elements in proc 1
        # print(elems_tags)
        # Extract Coordinates 
        proc_filename           = filename+'.h5.'+str(proc).zfill(len(str(nproc)))+'.feioutput'
        proc_outputfile         = h5py.File(proc_filename,'r');
        proc_elems_idx_to_conn  = proc_outputfile['/Model/Elements/Index_to_Connectivity'][:]
        proc_elems_conn         = proc_outputfile['/Model/Elements/Connectivity'][:]


        if solid_elems_tags[0] != 0: #if true, it means that solid elements exist
            
            common_elements = intersect1d(elems_tags, solid_elems_tags)  #find common elements between elems_tags in i proc and the solid_elements
            # print(common_elements)
            # print(elems_tags[0,:])

            proc_idx_extract        = proc_elems_idx_to_conn[common_elements]
            # print(proc_idx_extract)
            # Store Coordinates in a matrix 
            ConnMat[common_elements-1,0]= proc_elems_conn[proc_idx_extract]
            ConnMat[common_elements-1,1]= proc_elems_conn[proc_idx_extract+1]
            ConnMat[common_elements-1,2]= proc_elems_conn[proc_idx_extract+2]
            ConnMat[common_elements-1,3]= proc_elems_conn[proc_idx_extract+3]
            ConnMat[common_elements-1,4]= proc_elems_conn[proc_idx_extract+4]
            ConnMat[common_elements-1,5]= proc_elems_conn[proc_idx_extract+5]
            ConnMat[common_elements-1,6]= proc_elems_conn[proc_idx_extract+6]
            ConnMat[common_elements-1,7]= proc_elems_conn[proc_idx_extract+7]


        if beam_elems_tags[0] != 0: #if true, it means that beam elements exist
            
            common_elements = intersect1d(elems_tags, beam_elems_tags)  #find common elements between elems_tags in i proc and the beam_elements
            # print(common_elements)
            # print(elems_tags[0,:])

            proc_idx_extract        = proc_elems_idx_to_conn[common_elements]
            # print(proc_idx_extract)
            # Store Coordinates in a matrix 
            ConnMat[common_elements-1,0]= proc_elems_conn[proc_idx_extract]
            ConnMat[common_elements-1,1]= proc_elems_conn[proc_idx_extract+1]


        if bonded_elems_tags[0] != 0: #if true, it means that bonded contact elements exist
            
            common_elements = intersect1d(elems_tags, bonded_elems_tags)  #find common elements between elems_tags in i proc and the bonded_elements
            # print(common_elements)
            # print(elems_tags[0,:])

            proc_idx_extract        = proc_elems_idx_to_conn[common_elements]
            # print(proc_idx_extract)
            # Store Coordinates in a matrix 
            ConnMat[common_elements-1,0]= proc_elems_conn[proc_idx_extract]
            ConnMat[common_elements-1,1]= proc_elems_conn[proc_idx_extract+1]

        if quad_elems_tags[0] != 0: #if true, it means that quad elements exist
            
            common_elements = intersect1d(elems_tags, quad_elems_tags)  #find common elements between elems_tags in i proc and the quad_elements
            # print(common_elements)
            # print(elems_tags[0,:])

            proc_idx_extract        = proc_elems_idx_to_conn[common_elements]
            # print(proc_idx_extract)
            # Store Coordinates in a matrix 
            ConnMat[common_elements-1,0]= proc_elems_conn[proc_idx_extract]
            ConnMat[common_elements-1,1]= proc_elems_conn[proc_idx_extract+1]
            ConnMat[common_elements-1,2]= proc_elems_conn[proc_idx_extract+2]
            ConnMat[common_elements-1,3]= proc_elems_conn[proc_idx_extract+3]

        proc_outputfile.close()
        print(proc)
        end = time.time()
        print("Time consumed in working: ",end - start)
    # print(ConnMat)
    # Store as Xarray 
    ConnMat = xr.DataArray(ConnMat, dims=['nelems', 'nnodes_per_elem']); 
    print(' ******* CONNECTIVITY MATRIX EXTRACTED *******')
    os.chdir(pwd)
    return ConnMat




# ************************************************************************************************************************
# ************************************************** BODY OF THE SCRIPT **************************************************
# ************************************************************************************************************************

np.set_printoptions(threshold=sys.maxsize)
###########################################################################################   INPUT DATA TO CHANGE - START 
###########################################################################################

total_number_of_elements = 154

#solid_elements
start_element = 1      # put 0 if no solid elements exist
end_element = 144      # put 0 if no solid elements exist
solid_elems_tags = np.arange(start_element, end_element+1, 1, dtype=int)
# solid_elems_tags = [*range(1,144+1)]  # * is argument-unpacking operator
# print(solid_elems_tags)

#beam_elements
start_element = 0      # put 0 if no solid elements exist
end_element = 0      # put 0 if no solid elements exist
beam_elems_tags = np.arange(start_element, end_element+1, 1, dtype=int)
# print(beam_elems_tags)

#bonded_contact_elements
start_element = 149      # put 0 if no solid elements exist
end_element = 154      # put 0 if no solid elements exist
bonded_elems_tags = np.arange(start_element, end_element+1, 1, dtype=int)
# print(bonded_elems_tags)

#quad_elements
start_element = 145      # put 0 if no quad elements exist
end_element = 148      # put 0 if no quad elements exist
quad_elems_tags = np.arange(start_element, end_element+1, 1, dtype=int)
# print(quad_elems_tags)

filename = "hexahedron_quad_DRM"   # give the name of the loading stage
enable_relative_displacements = 1   # 1: reset displacements to zero at the beginning of the requested loading stage, 0: keep displacments from previous loading stage

# filter related parameters
filter_OnOff = 1;               # 0: disable filter, 1: enable filter
Fmin = 5;                       # minumum pass frequency
Fmax = 10;                      # maximum pass frequency
interpolation_check_node = 246  # give a node of interest to see the result of filter and interpolation functions on it (3 figures will be plotted)
ts_interpolation = 0.0025       # time step of interpolation

###########################################################################################
###########################################################################################   INPUT DATA TO CHANGE - END


###  EXTRACT GENERAL INFO SIMU ###
time_vec, ntsteps, nproc, nnodes, nodes_procid, nelems, elems_procid, ncoords, ndofs = EXTRACT_GENERAL_DATA_SIMU(filename, pwd, total_number_of_elements); 
dt = np.mean(np.diff(time_vec))     # time_vec[i+1] - time_vec[i],  gives a mean time step to be used when filter is applied
# print(time_vec)
# print(dt)

###  EXTRACT WAVEFIELD ###
x_xr, y_xr, z_xr, ux_xr, uy_xr, uz_xr, number_of_dofs = EXTRACT_COORD_WAVEFIELD(filename, pwd, time_vec, ntsteps, nproc, nnodes, nodes_procid, ncoords, ndofs); 

###  EXTRACT CONNECTIVITY ###
ConnMat = EXTRACT_CONN_MATRIX(filename,pwd, time_vec, ntsteps, nproc, nelems, elems_procid)


# *** Variables necessary to rebuild the elements of the grid with Unstructured Grid ***  
conns    = ConnMat.values.ravel().astype('int64'); # Connectivity matrix of the whole mesh with 
conns    = conns - 1;                              # the indexing needs always to start from 0    --  WHY? that means node 1 becomes node 0
# print(ConnMat)
# print(conns)
# print(len(conns))

conns = conns[ conns != -1 ]  #remove all -1 from conns
# print(conns)
# print(len(conns))

offset   = np.zeros(nelems, dtype = np.int)        # number of nodes connected by every element -  2 because we have simple beam/line elements 
# print(offset)
celltype = np.zeros(nelems, dtype = np.int)

#solid_elements
if solid_elems_tags[0] != 0: #if true, it means that solid elements exist
    # offset[solid_elems_tags-1] = 8
    celltype[solid_elems_tags-1] = 12 # different index depending on the cell type - 12 for solid elemenents (VTK_HEXAHEDRON) (check the different types here: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) 
    # print(offset)
    # print(celltype)
#beam_elements
if beam_elems_tags[0] != 0: #if true, it means that beam elements exist
    # offset[beam_elems_tags-1] = 2
    celltype[beam_elems_tags-1] = 3
    # print(offset)
#bonded_contact_elements
if bonded_elems_tags[0] != 0: #if true, it means that bonded contact elements exist
    # offset[bonded_elems_tags-1] = 2
    celltype[bonded_elems_tags-1] = 3
    # print(offset)
    # print(len(offset))
if quad_elems_tags[0] != 0: #if true, it means that quad elements exist
    celltype[quad_elems_tags-1] = 9
    # print(offset)
    # print(len(offset))
# print(celltype)


# offset array needs to monotonically increase later in pe.hl.unstructuredGridToVTK, so depending on the element type we add the corresponding number of nodes (e.g. 8 for a 8-node brick element). This is needed as conns is an array
sum = 0
for i in range(0,nelems):
    if celltype[i] == 12:    # VTK_HEXAHEDRON
        offset[i] = sum + 8
    elif celltype[i] == 3:   # VTK_LINE
        offset[i] = sum + 2
    elif celltype[i] == 9:   # VTK_QUAD
        offset[i] = sum + 4
    sum = offset[i]

# print(offset)
# print(celltype)
# print(len(celltype))
# print(ConnMat.shape[1])


###  FILTER APPLICATION ###
node = interpolation_check_node
ts = ts_interpolation
plt.figure(1)
plt.plot(time_vec,ux_xr[node-1,:], '.-k', markersize=10, linewidth=3, label='original') #plot the last row of the matrix ux_xr

if filter_OnOff == 1: 
    # create interpolation function to help us create a uniform time history (evenly spaced time steps)
    f_ux = interpolate.interp1d(time_vec, ux_xr, axis=1)
    print("interp1d ux done")
    f_uy = interpolate.interp1d(time_vec, uy_xr, axis=1)
    print("interp1d uy done")
    f_uz = interpolate.interp1d(time_vec, uz_xr, axis=1)
    print("interp1d uz done")
    t_interpolated = np.arange(0, max(time_vec)+0.00001, ts) #create the interpolated time steps
    ntsteps = len(t_interpolated)
    print('New total number of time steps = ' + str(ntsteps) + ' with dt = ' + str(ts)) 
    ux_interpolated = f_ux(t_interpolated)   # produce the interpolated ux using the interpolation function returned by `interp1d`
    print("ux_interpolated created")
    uy_interpolated = f_uy(t_interpolated)
    print("uy_interpolated created")
    uz_interpolated = f_uz(t_interpolated)
    print("uz_interpolated created")
    # print(ux_interpolated)
    # print(type(ux_interpolated)) #it became numpy array
    # print(np.shape(ux_interpolated))
    plt.plot(t_interpolated,ux_interpolated[node-1,:], '.-g', markersize=5, linewidth=1, label='interpolated')
    plt.title('Displacement time history ux - node #' + str(node))
    plt.xlabel("Time [s]")
    plt.ylabel("Displacements [m]")
    # plt.show()

    ## filter the interpolated uniform time history (irrfilter - default)
    #filter by giulia
    BB,AA = signal.iirfilter(5, np.array([Fmin,Fmax])*2*ts) # 1 is the Nyquist frequency, so if f1=2hz and Nyquist=50Hz, Wn=(2*1)/50=0.04, see my notes p.3
    ux_filt=signal.filtfilt(BB, AA,ux_interpolated,axis=1)
    print("ux_filt done")
    uy_filt=signal.filtfilt(BB, AA,uy_interpolated,axis=1)
    print("uy_filt done")
    uz_filt=signal.filtfilt(BB, AA,uz_interpolated,axis=1)
    print("uz_filt done")
    plt.plot( t_interpolated, ux_filt[node-1,:], '-r', linewidth=2, label='filtered at '+str(Fmin)+'-'+str(Fmax)+' Hz' )
    plt.legend()
    plt.savefig("displacement_time_histories_ALL_"+str(Fmin)+"-"+str(Fmax)+"Hz.png")

    #plot the digital filter frequency response 
    w, h = signal.freqz(BB, AA)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    nyquist = (1/ts) / 2
    plt.semilogx((nyquist*w)/np.pi  , 20 * np.log10(abs(h)), '-m', linewidth=3) # w is normalized to the range [0,pi] and plots up to nyquist freq, we change it to hz
    ax.set_title('Digital filter frequency response - irrfilter (default)')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Amplitude [dB]')
    ax.axis((0.1, nyquist, -200, 10))
    ax.grid(which='both', axis='both')
    plt.savefig('filter_frequency_response_'+str(Fmin)+'-'+str(Fmax)+'Hz.png')

    # Perform Fast Fourier Transform (FFT) for the check node
    fft_result_original = np.fft.fft(ux_xr[node-1,:])
    fft_amp_original = np.abs(fft_result_original) / (len(time_vec) / 2)
    fft_freq_original = np.fft.fftfreq(len(time_vec), time_vec[1] - time_vec[0])

    fft_result_interpolated = np.fft.fft(ux_interpolated[node-1,:])
    fft_amp_interpolated = np.abs(fft_result_interpolated) / (len(t_interpolated) / 2)
    fft_freq_interpolated = np.fft.fftfreq(len(t_interpolated), t_interpolated[1] - t_interpolated[0])

    fft_result_filtered = np.fft.fft(ux_filt[node-1,:])
    fft_amp_filtered = np.abs(fft_result_filtered) / (len(t_interpolated) / 2)
    fft_freq_filtered = np.fft.fftfreq(len(t_interpolated), t_interpolated[1] - t_interpolated[0])

    # Plot the FFT
    plt.figure()
    plt.plot(fft_freq_original, fft_amp_original, '-k', linewidth=5, label='original')
    plt.plot(fft_freq_interpolated, fft_amp_interpolated, '-g', linewidth=2, label='interpolated')
    plt.plot(fft_freq_filtered, fft_amp_filtered, '-r', linewidth=3, label='filtered at '+str(Fmin)+'-'+str(Fmax)+' Hz' )
    plt.legend()
    plt.title('FFT ux - node #' + str(node))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.xlim(0, 30)
    # plt.grid(True)
    plt.savefig('FFT_ALL_'+str(Fmin)+'-'+str(Fmax)+'Hz.png')

    print('#### Filter Applied #####')



### SAVE VTU OUTPUT FILES ###
# folder_out = 'vtu_files_forParaview_filter/'
if filter_OnOff == 1:
    folder_out = 'vtu_files_for_ParaView_filter_'+str(Fmin)+'-'+str(Fmax)+'Hz_' + filename + '/'
else:
    folder_out = 'vtu_files_for_ParaView_' + filename + '/'
os.system(str('mkdir -pv '+folder_out))
# mem_size  = ntsteps * nnodes * 8 / (1024*1024*1024) * 3
# print('allocato vector:',mem_size,' Gb')
x = x_xr.values.ravel(); y = y_xr.values.ravel(); z = z_xr.values.ravel();

# Select output times 
tstart = 0; tend   = ntsteps; tstep  = 1; 

dframe = 1
for t in range(tstart,tend,tstep): # one way to speed up this loop would be to call values only for a limited number of time
  vtk_file = folder_out+'step' + str(t)
  #When using Xarray
  if t%dframe==0:
     print('azzerato');ct=0
     if filter_OnOff == 1: 
        dump1 = ux_filt[:,t:t+dframe]; 
        dump2 = uy_filt[:,t:t+dframe]; 
        dump3 = uz_filt[:,t:t+dframe];  
     else:
        dump1 = ux_xr[:,t:t+dframe].values; 
        dump2 = uy_xr[:,t:t+dframe].values; 
        dump3 = uz_xr[:,t:t+dframe].values;     
     filename = pe.hl.unstructuredGridToVTK(vtk_file, x, y, z, connectivity=conns, offsets=offset, cell_types=celltype, cellData=None,pointData= {"displacement [m]:" : (dump1[:,ct].ravel(), dump2[:,ct].ravel(), dump3[:,ct].ravel())}); 
     # print(dump1)
     # print(dump1[:,ct].ravel())
     ct+=1
  print('## saving VTK ##',filename)



























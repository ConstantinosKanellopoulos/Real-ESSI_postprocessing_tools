# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# general_postprocess.py : python script that reads Real-ESSI output files (.feioutput) and saves in the multiple_results.txt file the requested quantity of interest (e.g., acceleration/
#                          displacement time history of a node in x direction)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Useful reference: Lecture notes of Boris Jeremic (http://sokocalo.engr.ucdavis.edu/~jeremic/LectureNotes/Jeremic_et_al_CompMech_LectureNotes.pdf)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# authors: Nikolaos Psycharis and Constantinos Kanellopoulos (2022)
# Chair of Structural Dynamics and Earthquake Engineering
# ETH, Zurich
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# import python modules
import h5py
import scipy as sp
import time
import sys
import linecache
import numpy as np

#############
### INPUT ###
#############

nodes = [158]                                              # give the node(s) tags,    # 0 if empty
nodes_Physical_Groups = [0]                                # give the name of the physical group in quotes (e.g., 'top_of_the_building')    # 0 if empty
node_dofs = 0                                              # 0: x-direction,   1: y-direction,   2: z-direction   (For more details refer to Lecture notes of Boris Jeremic, at the section Node-specific output format)

accelerations = 1                                          # 0: no, 1: yes
displacements = 1                                          # 0: no, 1: yes

elements = [0]		                                       # give the element(s) tags,    # 0 if empty
elements_Physical_Groups = [0]                             # give the name of the physical group in quotes (e.g., 'soil_layer_1')    # 0 if empty
element_dofs = 0                                           # refer to Lecture notes of Boris Jeremic (go to section Element-specific output format)

gauss_elements = [0]                                       # give the element(s) tags,    # 0 if empty
gauss_elements_Physical_Groups = [0]                       # give the name of the physical group in quotes (e.g., 'soil_layer_1')    # 0 if empty
gauss_element_dofs = 0                                     # refer to Lecture notes of Boris Jeremic (go to section Element-gauss output format)

energy_nodes = [0]                                         # give the node(s) tags,    # 0 if empty
energy_nodes_Physical_Groups = [0]                         # give the name of the physical group in quotes (e.g., 'DRM_input_node')    # 0 if empty
energy_node_dofs = 0                                       # refer to Lecture notes of Boris Jeremic (go to section Energy output format)

energy_elements = [0]		                               # give the element(s) tags,    # 0 if empty
energy_elements_Physical_Groups = [0]                      # give the name of the physical group in quotes (e.g., 'soil_layer_1')    # 0 if empty
energy_element_dofs = 0		                               # refer to Lecture notes of Boris Jeremic (go to section Energy output format)

essi_master_filename = 'hexahedron_DRM.h5.feioutput'

############
### CODE ###
############

essi_master_h5file = h5py.File(essi_master_filename,'r')
partition_variable = essi_master_h5file['/Model/Nodes/'] #oti iparxei mesa sta nodes

if not isinstance(nodes,list): #checks if nodes is not a list (a single node)
    nodes = [nodes] #make it a list to be able to loop over its length
if not isinstance(nodes_Physical_Groups,list):
    nodes_Physical_Groups = [nodes_Physical_Groups]
if not isinstance(node_dofs,list):
    node_dofs = [node_dofs]
if not isinstance(elements,list):
    elements = [elements]
if not isinstance(elements_Physical_Groups,list):
    elements_Physical_Groups = [elements_Physical_Groups]
if not isinstance(element_dofs,list):
    element_dofs = [element_dofs]
if not isinstance(gauss_elements,list):
    gauss_elements = [gauss_elements]
if not isinstance(gauss_elements_Physical_Groups,list):
    gauss_elements_Physical_Groups= [gauss_elements_Physical_Groups]
if not isinstance(gauss_element_dofs,list):
    gauss_element_dofs = [gauss_element_dofs]
if not isinstance(energy_nodes,list):
    energy_nodes = [energy_nodes]
if not isinstance(energy_nodes_Physical_Groups,list):
    energy_nodes_Physical_Groups = [energy_nodes_Physical_Groups]
if not isinstance(energy_node_dofs,list):
    energy_node_dofs = [energy_node_dofs]
if not isinstance(energy_elements,list):
    energy_elements = [energy_elements]
if not isinstance(energy_elements_Physical_Groups,list):
    energy_elements_Physical_Groups = [energy_elements_Physical_Groups]
if not isinstance(energy_element_dofs,list):
    energy_element_dofs = [energy_element_dofs]


for i in range(0,len(nodes_Physical_Groups)): #tsekarei an exeis physical groups na kanei append sta nodes, elemenets...pou exw orisei manually stin arxi
	if nodes_Physical_Groups[i] != 0:
		nodes += essi_master_h5file['/Model/Physical_Groups/Physical_Node_Groups/'+nodes_Physical_Groups[i]]
for i in range(0,len(elements_Physical_Groups)):
	if elements_Physical_Groups[i] != 0:
		elements += essi_master_h5file['/Model/Physical_Groups/Physical_Element_Groups/'+elements_Physical_Groups[i]]
for i in range(0,len(gauss_elements_Physical_Groups)):
	if gauss_elements_Physical_Groups[i] != 0:
		gauss_elements += essi_master_h5file['/Model/Physical_Groups/Physical_Element_Groups/'+gauss_elements_Physical_Groups[i]]
for i in range(0,len(energy_nodes_Physical_Groups)):
	if energy_nodes_Physical_Groups[i] != 0:
		energy_nodes += essi_master_h5file['/Model/Physical_Groups/Physical_Node_Groups/'+energy_nodes_Physical_Groups[i]]
for i in range(0,len(energy_elements_Physical_Groups)):
	if energy_elements_Physical_Groups[i] != 0:
		energy_elements += essi_master_h5file['/Model/Physical_Groups/Physical_Element_Groups/'+energy_elements_Physical_Groups[i]]


nodes_sum = []
elements_sum = []
gauss_elements_sum = []
energy_elements_sum = []
energy_nodes_sum = []

if 'Partition' in partition_variable: #checks if you are in parallel

    partition_switch = 1  #to check for parallel
    partition_nodes = essi_master_h5file['/Model/Nodes/Partition']
    partition_elements = essi_master_h5file['/Model/Elements/Partition']
    total_nodes_partitions = max(partition_nodes)
    total_elements_partitions = max(partition_elements)
    total_partitions = max(total_nodes_partitions,total_elements_partitions) #number of processors
    open_partition_check = [] #variable to check if a partition should open
    essi_filename = []
    #initialize matrix
    nodes_p = [ [ 0 for i in range(len(nodes)) ] for j in range(total_partitions) ]
    elements_p = [ [ 0 for i in range(len(elements)) ] for j in range(total_partitions) ]
    gauss_elements_p = [ [ 0 for i in range(len(gauss_elements)) ] for j in range(total_partitions) ]
    energy_nodes_p = [ [ 0 for i in range(len(energy_nodes)) ] for j in range(total_partitions) ]
    energy_elements_p = [ [ 0 for i in range(len(energy_elements)) ] for j in range(total_partitions) ]

    for i in range(0,total_partitions):
    	if total_partitions > 98:
            if i > 98:
                essi_filename.append(essi_master_filename.replace("feioutput", (str(i+1) + '.feioutput')))
            elif i > 8:
                essi_filename.append(essi_master_filename.replace("feioutput", ("0"+str(i+1) + '.feioutput')))
            else:
                essi_filename.append(essi_master_filename.replace("feioutput", ("00"+str(i+1) + '.feioutput')))

        elif total_partitions > 8:
            if i > 8:
                essi_filename.append(essi_master_filename.replace("feioutput", (str(i+1) + '.feioutput')))
            else:
                essi_filename.append(essi_master_filename.replace("feioutput", ("0"+str(i+1) + '.feioutput')))

        else:   
            essi_filename.append(essi_master_filename.replace("feioutput", (str(i+1) + '.feioutput')))

        nodes_temp = 0  #to nodes_sum den xreiazetai mallon
        for j in range(0,len(nodes)):  #pame orizontia 
            if partition_nodes[nodes[j]] == i+1:
                nodes_p[i][j] = nodes[j]
                nodes_temp += nodes[j] # mia random timi apla megaluteri tou 0 gia na deikseis oti se auto to partition exeis nodes (pio katw sto open_partition_check) , tha mporouse na einai ison 1
        nodes_sum.append(nodes_temp) #mallon axristo

        elements_temp = 0
        for j in range(0,len(elements)):
            if partition_elements[elements[j]] == i+1:
                elements_p[i][j] = elements[j]
                elements_temp += elements[j]
        elements_sum.append(elements_temp)

        gauss_elements_temp = 0
        for j in range(0,len(gauss_elements)):
            if partition_elements[gauss_elements[j]] == i+1:
                gauss_elements_p[i][j] = gauss_elements[j]
                gauss_elements_temp += gauss_elements[j]
        gauss_elements_sum.append(gauss_elements_temp)

        energy_nodes_temp = 0
        for j in range(0,len(energy_nodes)):
            if partition_nodes[energy_nodes[j]] == i+1:
                energy_nodes_p[i][j] = energy_nodes[j]
                energy_nodes_temp += energy_nodes[j]
        energy_nodes_sum.append(energy_nodes_temp)

        energy_elements_temp = 0
        for j in range(0,len(energy_elements)):
            if partition_elements[energy_elements[j]] == i+1:
                energy_elements_p[i][j] = energy_elements[j]
                energy_elements_temp += energy_elements[j]
        energy_elements_sum.append(energy_elements_temp)

        open_partition_check.append(nodes_temp+elements_temp+gauss_elements_temp+energy_nodes_temp+energy_elements_temp)
else:

    partition_switch = 0

Generalized_Results = []
Names_Matrix = ['time']

if partition_switch == 1:
    for partitioni in range(0,total_partitions):  #pernaei apo ola ta partitions alla den ta anoigei ola
        if open_partition_check[partitioni] > 0:
            essi_h5file = h5py.File(essi_filename[partitioni],'r') #opens partition and open all indexes
            essi_Index_to_Generalized_Displacements = essi_h5file['/Model/Nodes/Index_to_Generalized_Displacements'][:]
            essi_Index_to_Element_Outputs = essi_h5file['/Model/Elements/Index_to_Element_Outputs'][:]
            essi_Index_to_Gauss_Outputs = essi_h5file['/Model/Elements/Index_to_Gauss_Outputs'][:]
            essi_Index_to_Nodal_Input_Energy = essi_h5file['/Energy/Index_to_Nodal_Input_Energy'][:]
            essi_Index_to_Energy_Element = essi_h5file['/Energy/Index_to_Energy_Element'][:]

            for nodei in range(0,len(nodes)):
                if nodes_p[partitioni][nodei] != 0:
                    if essi_Index_to_Generalized_Displacements[nodes_p[partitioni][nodei]] == -1: #if index is -1 print error
                        Names_Matrix.append(str(nodes_p[partitioni][nodei])+'_no disps/accels') 
                        Generalized_Results.append([]) #vazei ena keno stoixeio sti lista
                    else:
                        if accelerations == 1:
                            for dofi in range(0,len(node_dofs)):
                                essi_generalized_temp = essi_h5file['/Model/Nodes/Generalized_Accelerations'][essi_Index_to_Generalized_Displacements[nodes_p[partitioni][nodei]]+node_dofs[dofi],:]
                                Names_Matrix.append(str(nodes_p[partitioni][nodei])+'_acc'+str(node_dofs[dofi])+':')
                                Generalized_Results.append(essi_generalized_temp)
                        if displacements == 1:
                            for dofi in range(0,len(node_dofs)):
                                essi_generalized_temp = essi_h5file['/Model/Nodes/Generalized_Displacements'][essi_Index_to_Generalized_Displacements[nodes_p[partitioni][nodei]]+node_dofs[dofi],:]
                                Names_Matrix.append(str(nodes_p[partitioni][nodei])+'_disp'+str(node_dofs[dofi])+':')
                                Generalized_Results.append(essi_generalized_temp)

            for elementi in range(0,len(elements)):
                if elements_p[partitioni][elementi] != 0:
                    if essi_Index_to_Element_Outputs[elements_p[partitioni][elementi]] == -1:
                        Names_Matrix.append(str(elements_p[partitioni][elementi])+'_no elem. output')
                        Generalized_Results.append([])
                    else:
                        for dofi in range(0,len(element_dofs)):
                            essi_generalized_temp = essi_h5file['/Model/Elements/Element_Outputs'][essi_Index_to_Element_Outputs[elements_p[partitioni][elementi]]+element_dofs[dofi],:]
                            Names_Matrix.append(str(elements_p[partitioni][elementi])+'_elem.out'+str(element_dofs[dofi])+':')
                            Generalized_Results.append(essi_generalized_temp)

            for gausselementi in range(0,len(gauss_elements)):
                if gauss_elements_p[partitioni][gausselementi] != 0:
                    if essi_Index_to_Gauss_Outputs[gauss_elements_p[partitioni][gausselementi]] == -1:
                        Names_Matrix.append(str(gauss_elements_p[partitioni][gausselementi])+'_no gauss output')
                        Generalized_Results.append([])
                    else:
                        for dofi in range(0,len(gauss_element_dofs)):
                            essi_generalized_temp = essi_h5file['/Model/Elements/Gauss_Outputs'][essi_Index_to_Gauss_Outputs[gauss_elements_p[partitioni][gausselementi]]+gauss_element_dofs[dofi],:]
                            Names_Matrix.append(str(gauss_elements_p[partitioni][gausselementi])+'_gaus.out'+str(gauss_element_dofs[dofi])+':')
                            Generalized_Results.append(essi_generalized_temp)

            for energynodei in range(0,len(energy_nodes)):
                if energy_nodes_p[partitioni][energynodei] != 0:
                    if essi_Index_to_Nodal_Input_Energy[energy_nodes_p[partitioni][energynodei]] == -1:
                        Names_Matrix.append(str(energy_nodes_p[partitioni][energynodei])+'_no node energy')
                        Generalized_Results.append([])
                    else:
                        for dofi in range(0,len(energy_node_dofs)):
                            essi_generalized_temp = essi_h5file['/Energy/Nodal_Input_Energy'][essi_Index_to_Nodal_Input_Energy[energy_nodes_p[partitioni][energynodei]]+energy_node_dofs[dofi],:]
                            print(essi_generalized_temp)
                            Names_Matrix.append(str(energy_nodes_p[partitioni][energynodei])+'_enrgy.nod'+str(energy_node_dofs[dofi])+':')
                            Generalized_Results.append(essi_generalized_temp)

            for energyelementi in range(0,len(energy_elements)):
                if energy_elements_p[partitioni][energyelementi] != 0:
                    if essi_Index_to_Energy_Element[energy_elements_p[partitioni][energyelementi]] == -1:
                        Names_Matrix.append(str(energy_elements_p[partitioni][energyelementi])+'_no elem. energy')
                        Generalized_Results.append([])
                    else:
                        for dofi in range(0,len(energy_element_dofs)):
                            essi_generalized_temp = essi_h5file['/Energy/Energy_Element'][essi_Index_to_Energy_Element[energy_elements_p[partitioni][energyelementi]]+energy_element_dofs[dofi],:]
                            Names_Matrix.append(str(energy_elements_p[partitioni][energyelementi])+'_enrgy.elem'+str(energy_element_dofs[dofi])+':')
                            Generalized_Results.append(essi_generalized_temp)

else: #not parallel
    essi_h5file = essi_master_h5file
    essi_Index_to_Generalized_Displacements = essi_h5file['/Model/Nodes/Index_to_Generalized_Displacements'][:]
    essi_Index_to_Element_Outputs = essi_h5file['/Model/Elements/Index_to_Element_Outputs'][:]
    essi_Index_to_Gauss_Outputs = essi_h5file['/Model/Elements/Index_to_Gauss_Outputs'][:]
    essi_Index_to_Nodal_Input_Energy = essi_h5file['/Energy/Index_to_Nodal_Input_Energy'][:]
    essi_Index_to_Energy_Element = essi_h5file['/Energy/Index_to_Energy_Element'][:]

    for nodei in range(0,len(nodes)):
        if nodes[nodei] != 0:
            if essi_Index_to_Generalized_Displacements[nodes[nodei]] == -1:
                Names_Matrix.append(str(nodes[nodei])+'_no disps/accels')
                Generalized_Results.append([])
            else:
                if accelerations == 1:
                    for dofi in range(0,len(node_dofs)):
                        essi_generalized_temp = essi_h5file['/Model/Nodes/Generalized_Accelerations'][essi_Index_to_Generalized_Displacements[nodes[nodei]]+node_dofs[dofi],:]
                        Names_Matrix.append(str(nodes[nodei])+'_acc'+str(node_dofs[dofi])+':')
                        Generalized_Results.append(essi_generalized_temp)
                if displacements == 1:
                    for dofi in range(0,len(node_dofs)):
                        essi_generalized_temp = essi_h5file['/Model/Nodes/Generalized_Displacements'][essi_Index_to_Generalized_Displacements[nodes[nodei]]+node_dofs[dofi],:]
                        Names_Matrix.append(str(nodes[nodei])+'_disp'+str(node_dofs[dofi])+':')
                        Generalized_Results.append(essi_generalized_temp)

    for elementi in range(0,len(elements)):
        if elements[elementi] != 0:
            if essi_Index_to_Element_Outputs[elements[elementi]] == -1:
                Names_Matrix.append(str(elements[elementi])+'_no elem. output')
                Generalized_Results.append([])
            else:
                for dofi in range(0,len(element_dofs)):
                    essi_generalized_temp = essi_h5file['/Model/Elements/Element_Outputs'][essi_Index_to_Element_Outputs[elements[elementi]]+element_dofs[dofi],:]
                    Names_Matrix.append(str(elements[elementi])+'_elem.out'+str(element_dofs[dofi])+':')
                    Generalized_Results.append(essi_generalized_temp)

    for gausselementi in range(0,len(gauss_elements)):
        if gauss_elements[gausselementi] != 0:
            if essi_Index_to_Gauss_Outputs[gauss_elements[gausselementi]] == -1:
                Names_Matrix.append(str(gauss_elements[gausselementi])+'_no gauss output')
                Generalized_Results.append([])
            else:
                for dofi in range(0,len(gauss_element_dofs)):
                    essi_generalized_temp = essi_h5file['/Model/Elements/Gauss_Outputs'][essi_Index_to_Gauss_Outputs[gauss_elements[gausselementi]]+gauss_element_dofs[dofi],:]
                    Names_Matrix.append(str(gauss_elements[gausselementi])+'_gaus.out'+str(gauss_element_dofs[dofi])+':')
                    Generalized_Results.append(essi_generalized_temp)

    for energynodei in range(0,len(energy_nodes)):
        if energy_nodes[energynodei] != 0:
            if essi_Index_to_Nodal_Input_Energy[energy_nodes[energynodei]] == -1:
                Names_Matrix.append(str(energy_nodes[energynodei])+'_no node energy')
                Generalized_Results.append([])
            else:
                for dofi in range(0,len(energy_node_dofs)):
                    essi_generalized_temp = essi_h5file['/Energy/Nodal_Input_Energy'][essi_Index_to_Nodal_Input_Energy[energy_nodes[energynodei]]+energy_node_dofs[dofi],:]
                    Names_Matrix.append(str(energy_nodes[energynodei])+'_enrgy.nod'+str(energy_node_dofs[dofi])+':')
                    Generalized_Results.append(essi_generalized_temp)

    for energyelementi in range(0,len(energy_elements)):
        if energy_elements[energyelementi] != 0:
            if essi_Index_to_Energy_Element[energy_elements[energyelementi]] == -1:
                Names_Matrix.append(str(energy_elements[energyelementi])+'_no elem. energy')
                Generalized_Results.append([])
            else:
                for dofi in range(0,len(energy_element_dofs)):
                    essi_generalized_temp = essi_h5file['/Energy/Energy_Element'][essi_Index_to_Energy_Element[energy_elements[energyelementi]]+energy_element_dofs[dofi],:]
                    Names_Matrix.append(str(energy_elements[energyelementi])+'_enrgy.elem'+str(energy_element_dofs[dofi])+':')
                    Generalized_Results.append(essi_generalized_temp)

essi_Time = essi_h5file['time'][:]

#############
### PRINT ###
#############

filename = open('multiple_results.txt','w')
last = len(Generalized_Results[0])-1

filename.write(Names_Matrix[0]+'\t')
for i in range(0,len(essi_Time)):  #plots first line, time
    if i == last:
        filename.write(str(essi_Time[i]))
    else:
        filename.write(str(essi_Time[i])+'\t')
filename.write('\n')

namei = 1

for row in Generalized_Results:
    filename.write(Names_Matrix[namei]+'\t')
    for i in range(0,len(row)):
        if i == last:
            filename.write(str(row[i]))
        else:
            filename.write(str(row[i])+'\t')
    filename.write('\n')
    namei += 1

filename.close()



###############################################################


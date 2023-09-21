# Real-ESSI postprocessors



### Description

- [generate_vtu.py](./generate_vtu_files/generate_vtu.py) is a python script that extracts all nodal displacements from [Real-ESSI](http://real-essi.info/) output files (.feioutput), applies (if requested) a frequency filter, and finally generates .vtu files to visualize displacements in [ParaView](https://www.paraview.org/). For more information refer to [generate_vtu.py](./generate_vtu_files/generate_vtu.py).

This script was initially created for the needs of the paper: 
<p align="center">
<strong>" Dynamic Structure-Soil-Structure Interaction for Nuclear Power Plants "</strong>
</p>

An example figure from the paper after postpressing the output files with [generate_vtu.py](./generate_vtu_files/generate_vtu.py) follows :

<p align="center">
  <img src="https://github.com/ConstantinosKanellopoulos/images_for_my_repo/blob/master/NPP_reactor_tied_and_nonlinear_interfaces_filtered.png" width="500" height="300" >
</p>

<!-- <p align="center">
<strong>Figure: Snapshots of filtered displacement contours (10-13 Hz) for nonlinear (left) and tied (right) soil-NPP reactor building interfaces.</strong>
</p> -->

- [print_node_element_output.py](print_node-or-element_outputs/print_node_element_output.py) is a python script that reads [Real-ESSI](http://real-essi.info/) output files (.feioutput) and generates a multiple_results.txt file with the requested node/element output(s) (e.g. acceleration time history of the requested node in x direction). For more information refer to [print_node_element_output.py](print_node-or-element_outputs/print_node_element_output.py).

This script was initially created for the needs of the paper: 

<p align="center">
  <strong><a href="https://doi.org/10.1016/j.soildyn.2022.107366">Seismic resonant metamaterials for the protection of an elastic-plastic SDOF system against vertically propagating seismic shear waves (SH) in nonlinear soil</a></strong>
</p>


### Required software

- for [generate_vtu.py](./generate_vtu_files/generate_vtu.py) :
  - [python](https://www.python.org/)
  - [ParaView](https://www.paraview.org/)

- for [print_node_element_output.py](print_node-or-element_outputs/print_node_element_output.py)
  - [python](https://www.python.org/)


### How to run

- To run [generate_vtu.py](./generate_vtu_files/generate_vtu.py), simply open a terminal and type:

```bash
python generate_vtu.py
```
Note: The .feioutput files should be located one directory above (see in [Examples](./generate_vtu_files/Examples))

- To run [print_node_element_output.py](print_node-or-element_outputs/print_node_element_output.py), simply open a terminal and type:

```bash
python print_node_element_output.py
```
Note: The .feioutput files should be located in the same directory (see [Example](./print_node-or-element_outputs/Example))


### How to cite

- [generate_vtu.py](./generate_vtu_files/generate_vtu.py) :

- [print_node_element_output.py](print_node-or-element_outputs/print_node_element_output.py) :  [![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.soildyn.2022.107366-purple)](https://doi.org/10.1016/j.soildyn.2022.107366)

--

### Licence

[GNU General Public License v3.0](./COPYING)



### Acknowledgement

This research was part of the project '[INSPIRE](https://itn-inspire.eu/)' funded by the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement no. 813424. Additional support from [ETH Zurich](https://ethz.ch/en.html) is gratefully acknowledged.

<!-- <img align="center" src="https://github.com/ConstantinosKanellopoulos/images_for_my_repo/blob/master/logos.png"> -->

<p align="center">
  <img src="https://github.com/ConstantinosKanellopoulos/images_for_my_repo/blob/master/logos.png">
</p>



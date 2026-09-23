#!/usr/bin/env python
# coding: utf-8

# # Preparation of a FEMTIC inversion input files: data and mesh 
# 
# application à Annecy

# In[1]:


import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
femticpy_path = os.path.abspath(os.path.join(current_dir, '..', '..'))
sys.path.append(femticpy_path+'/src')
import femticPy


# # 1. Prepare inversion data file

# In[2]:


# the file will be written in the input_data directory
inversion = femticPy.DataGen(survey = 'test', outdir = './inversion')


# In[3]:


# Loading the data and the data coordinates
inversion.read_MTdata('./input_data/edi_files/edi_files')
inversion.read_MTdata_coordinates('./input_data', 'coords_annecy')


# In[4]:


## Data to be inverted for
inversion.invert_Z   = True
inversion.invert_VTF = True
inversion.invert_PT  = False


# In[5]:


# load the 3 files needed to create the mesh
inversion.topography = './input_data/topography.dat'
inversion.bathymetry = './input_data/bathymetry.dat'
inversion.coast_line = './input_data/coast_line.dat'


# In[6]:


# we center the data to a anchor point (center of the data set) which will also be the center of the future mesh
inversion.center_data()
inversion.anchor


# In[7]:


# we define the analysis domain, the extent of the total mesh including padding
inversion.analysis_domain = [[-30, 30],
                                [-30, 30],
                                [-60, 100]] 
inversion.plot_data_loc(plot_ids = False, zoom_core=False)
inversion.plot_coast_line()
# inversion.plot_topo_bathy()


# In[8]:


inversion.mt_coords


# In[9]:


#inversion.data_Z


# In[10]:


# Set up error floors
# by default, Cross-term scaling per pair
# 0: ZXX vs ZXY, 1-2: ZYX vs ZXY, 3: ZYY vs ZYX
inversion.error_floor_Z = [0.15, 0.05, 0.05, 0.15]  # [Zxx, Zxy, Zyx, Zyy]
inversion.error_floor_Tz = 0.03

# apply defined error values (equal to error floors) to diagonal components of Z and the tipper
inversion.data_Z = inversion.apply_defined_error(inversion.data_Z, data_type='Z',
                                comps=['ZXX.VAR', 'ZYY.VAR'],
                                err_val=[0.15, 0.15],
                                error_type=2)

inversion.data_VTF = inversion.apply_defined_error(inversion.data_VTF, data_type='VTF',
                                comps=['TXVAR.EXP','TYVAR.EXP'],
                                err_val=[0.03,0.03])

# In[11]:

# define frequency interval:
highfreq=1000
lowfreq=0.035

subsampling = 1
inversion.freqs[::subsampling]

# Add gaussian noise:
inversion.data_Z   = inversion.add_gaussian_noise(inversion.data_Z, data_type='Z')
inversion.data_VTF = inversion.add_gaussian_noise(inversion.data_VTF, data_type='VTF')


# In[12]:


inversion.write_observe(write = True, freq_bandwidth = [lowfreq,highfreq], subsampling = subsampling)


# # # 2. Prepare inversion control file

# # In[13]:


# # Define inversion parameters:

# #TRADE_OFF_PARAM
# TRADE_OFF_PARAM = [3]

# # Max number of iterations
# ITERATION = 30

# # NUM_THREADS
# NUM_THREADS = 4

# # CONVERGE
# CONVERGE = 0.2

# # output directory
# outdir = 'inversion'


# # In[14]:


# inversion.write_inversion_control(TRADE_OFF_PARAM = TRADE_OFF_PARAM,
#                                     ITERATION = ITERATION,
#                                     NUM_THREADS = NUM_THREADS,
#                                     CONVERGE = CONVERGE)


# # # 3. Prepare mesh files

# # In[15]:


# # output directory
# outdir = 'meshGen/files'


# # In[16]:


# # create meshGen object
# mesh01 = femticPy.MeshGen(inversion.survey,
#                          inversion.analysis_domain, 
#                          inversion.center, 
#                          inversion.mt_coords,
#                          inversion.nRx_Z,
#                          outdir)


# # In[17]:


# mesh01.sea = False
# mesh01.sea


# # In[18]:


# # Step 1: 

# # Define the mesh discretization, sphere centered on the mesh center 
#     # core area: 650km sphere radius (isotropic) / max cell size : 50km (anisotropic in z)
#     # padding: 2000km sphere radius (isotropic), max cell size: 500km (isotropic)

# mesh01.center = [0.0, 0.0, 0.0]
# mesh01.rotation = 0

# # mesh01.ellipsoids_control = [2,
# #                      [650.0, 50.0, 0.35, 0.0, 0.0],
# #                      [2000.0, 500.0, 0.0, 0.0, 0.0]]

# mesh01.ellipsoids_control = [2,
#                      [13.5, 1, 0.41, 0, 0],
#                      [40, 10, 0.0, 0.0, 0.0]]

# mesh01.coast_line = './input_data/coast_line.dat'

# # mesh01.ellipsoids_observing_sites =  [6, 
# #                                       1.0, 0.5,
# #                                       5.0, 1.0,
# #                                      50.0, 15.0,
# #                                       100.0, 25.0,
# #                                       150.0,50.0,
# #                                      200.0, 100.0]

# mesh01.ellipsoids_observing_sites =  [6, 
#                                       0.06, 0.03,
#                                       0.3, 0.06,
#                                       3., 1.1,
#                                       6., 1.5,
#                                       8., 2.5,
#                                       15., 7.]


# # In[19]:


# # Step 6: Refine
# # this step refines the tetras around the MT sites

# # mesh01.ellipsoids_obs_sites = [4,
# #                     [1.0, 0.5, 0.0] ,          
# #                     [10.0, 5.0, 0.0] ,           
# #                     [100.0, 25.0, 0.0],         
# #                     [200.0, 100.0, 0.0],]


# mesh01.ellipsoids_obs_sites = [4,
#                     [0.06, 0.03, 0.0] ,          
#                     [0.6, 0.3, 0.0] ,           
#                     [6., 1.5, 0.0],         
#                     [15., 7., 0.0],]


# # mesh01.ellipsoids_mtr = [2,
# #                     [150, 50.0, 0.0, 0.1, 0.1],
# #                     [2000.0, 500.0, 0.0, 0.0, 0.0]]

# mesh01.ellipsoids_mtr = [2,
#                     [15., 2., 0., 0., 0.],
#                     [30., 8., 0., 0., 0.]]


# # In[20]:


# # Step 7: Write mesh to FEMTIC format
# # the attr file will determine how tetras (elements) will be combined into cells (model parameters)

# # mesh01.ellipsoids_resistivity_attr = [2,
# #                     [400.0,       50.0,  0.0,  0.0],
# #                     [2000.0,  500.0,  0.0,  0.0]] 


# # mesh01.ellipsoids_resistivity_attr_sites = [3,                        
# #                     [10.0, 2.0],
# #                     [50.0, 10.0],
# #                     [150.0, 25.0]]

# mesh01.ellipsoids_resistivity_attr = [2,
#                      [15., 2., 0., 0.],
#                      [30., 8., 0., 0.]]

# mesh01.ellipsoids_resistivity_attr_sites = [3,        
#                     [0.6, 0.1] ,           
#                     [3., 0.6],         
#                     [10., 1.6]]


# # In[21]:


# # Define the starting model resistivities (in ohm.m)
# mesh01.resistivity_starting_model = 100


# # In[22]:


# mesh01.write_inputs()


# # In[23]:


# # number of regions (air, land, sea) 
# mesh01.region_attributes


# # # 4. Run `create_mesh.sh`

# # the `create_mesh.sh` script call the `meshGen.sh` pipeline, that sequentially runs all the steps required to create the mesh. Depending on the size of the mesh this can take long (~ 10 min). 
# # 
# # "Step 4" in `meshGen.sh` needs to be modified to your problem:

# # ```
# # #-----------------------"
# # #-----   STEP 4    -----"
# # #-----------------------"
# # 
# # echo "STEP 4... "
# # /scratch3/sei029/femtic/makeTetraMesh -stp 4
# # sed -i "$(( $(wc -l < output.poly) ))s/.*/3/" output.poly
# # echo "1   0.0 0.0 -40.0  10 1e9" >> output.poly
# # echo "2  -500.0 -500.0  1.0  20 1e9" >> output.poly
# # echo "3   0.0 0.0  40.0  30 1e9" >> output.poly
# # echo "   done!"
# # sleep 1
# # ```

# # In this step, 3 lines are added to the `output.poly` file. They correspond to:
# # 
# # `"region_attribute   x   y   z   region_number   1e9"`
# # 
# # x, y, z (in km) should indicate at which location a sample of each region can be found. z is positive downwards.
# # 
# # The region numbers have to match with the ones in `resistivity_attr.dat`
# # 
# # ```
# # 3
# # 10 1.0e+9 -1 1
# # 20 0.250  -1 1
# # 30 100.0   9 0
# # ```
# # 
# # 10 is for air, 20 for sea and 30 for earth (with associated resistivities)

# # "Step 5" in `meshGen.sh` calls a python script that fixes a bug that has been appearing sometimes. This is a temporary fix while understanding why a new region number (31 or 32) sometimes appears during step 5. 
# # 
# # ```
# # echo "fix output.1.ele if necessary... "
# # module load python
# # python fix_output1.py
# # sleep 1
# # echo "   done!"
# # ```

# # To check the size of the mesh, check the `resistivity_block_iter0.dat` file. The first number is the number of elements, the second is te number of cells.
# # 
# # ```
# # >head resistivity_block_iter0.dat 
# #    1453333    206736
# #          0         0
# # ....
# # ```

# # # 5. Visualise the mesh using paraview (before running the inversion)

# # Setting up the ellipsoids is cumbersome and takes some trials to obtain a mesh that looks ok, in terms of cells discretization and mesh size. 
# # 
# # If the shell script has run properly, the final files to be used for the inversion should be:
# # - resistivity_block_iter0.dat
# # - mesh.dat
# # 
# # 

# # The output files that can be copied over and are useful for visualization are:
# # 
# # - output.6.femtic.vtk
# #    - to check the model parameter discretization (the cells), use the blockSerial attribute of the .vtk file
# # - triangles_with_height.vtk
# # - coastLine_fine.vtk
# # - coastLine_rough.vtk
# # 
# # The last 3 files are in km, and a transform filter has to be applied in paraview to be visualized with `output.6.femtic.vtk` (which is m)

# # # 6. Running the inversion using `run_inversion.sh`

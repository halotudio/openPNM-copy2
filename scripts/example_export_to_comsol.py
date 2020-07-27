import sys
sys.path.insert(1, '/home/mehrez/Desktop/repositories/OpenPNM/')
import openpnm as op
import numpy as np

# work space and project
ws = op.Workspace()
proj = ws.new_project()

# network
np.random.seed(7)
net = op.network.Cubic(shape=[3, 3, 1], spacing=1e-4, project=proj)

# geometry
geo = op.geometry.StickAndBall(network=net,
                               pores=net.Ps,
                               throats=net.Ts)

# output network into a 2D Comsol file (network has to be a single layer!)
proj.export_data(filename='net_comsol_2d', filetype='Comsol', dimension='2D')

# output network into a 3D Comsol file
proj.export_data(filename='net_comsol_3d', filetype='Comsol', dimension='3D')

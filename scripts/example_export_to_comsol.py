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

# phase
phase = op.phases.Water(network=net)

# output results
proj.export_data(filename='comsol_test', filetype='Comsol', dimension='2D')

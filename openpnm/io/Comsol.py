import numpy as np
from openpnm.io import GenericIO
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class Comsol(GenericIO):
    r"""
    Writes files containing pores and throats of the considered network in a
    COMSOL object format.

    Notes
    -----
    The exported files contain COMSOL geometry objects, not meshes.
    This class exports in 2D only.

    """
    @classmethod
    def save(cls, network, phases=[], filename='', dimension=''):
        r"""
        Saves the network and geometry data from the given objects into the
        specified file. This exports in 2D only where throats and pores have
        rectangular and circular shapes, respectively.

        Parameters
        ----------
        network : OpenPNM Network Object
            The network containing the desired data

        phases : list of OpenPNM Phase Objects (optional, default is none)

        Notes
        -----
        This method only saves the network and geometry data, not any of the
        pore-scale models or other attributes.  To save an actual OpenPNM
        Project use the ``Workspace`` object.

        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        network = network[0]
        cn = network['throat.conns']
        p = network['pore.coords']
        p1 = p[cn[:, 0]]
        p2 = p[cn[:, 1]]

        # check network's dimenesion (1D, 2D or 3D)

        dim = np.all(p == p[0, :], axis=0)
        dim = np.where(~dim)[0]
        actual_dim = str(len(dim))+'D'

        if dimension not in ['2D', '3D', '2d', '3d']:
            logger.warning('Supported "dimension" is either "2D" or "3D"')
            dimension = '3D'

        if actual_dim in ['3D', '3d'] and dimension in ['2D', '2d']:
            logger.error('Cannot export a 3D network into a ' +
                         dimension + ' Comsol file')

        if actual_dim == '0D':
            logger.exception('Cannot export a "0D" network into a Comsol file')

        if dimension in ['2D', '2d']:
            # compute the rotation angle of throats
            try:
                dif_x = p2[:, dim[0]]-p1[:, dim[0]]
                dif_y = p2[:, dim[1]]-p1[:, dim[1]]
            # handle 1D cases
            except IndexError:
                if dim[0] == 0:
                    dif_x = p2[:, dim[0]]-p1[:, dim[0]]
                    dif_y = p2[:, dim[0]+1]-p1[:, dim[0]+1]
                else:
                    dif_x = p2[:, dim[0]-1]-p1[:, dim[0]-1]
                    dif_y = p2[:, dim[0]]-p1[:, dim[0]]
            # Avoid division by 0
            m = np.array([dif_x_i != 0 for dif_x_i in dif_x])
            r = np.zeros((len(dif_x)))
            r[~m] = np.inf
            r[m] = dif_y[m]/dif_x[m]
            angles = np.arctan(r)

            r_w = network['throat.diameter']
            c_c = network['pore.coords']
            c_r = network['pore.diameter']/2.0

            f = open(filename+'.mphtxt', 'w')
            header(file=f, Nt=network.Nt, Np=network.Np)
            rectangles(file=f, pores1=p1, pores2=p2, alphas=angles, widths=r_w)
            circles(file=f, centers=c_c, radii=c_r)
            f.close()

        elif dimension in ['3D', '3d']:
            s_c = network['pore.coords']
            s_r = network['pore.diameter']/2.0

            f = open(filename+'.mphtxt', 'w')
            header(file=f, Nt=0, Np=network.Np)
            spheres(file=f, centers=s_c, radii=s_r)
            f.close()


def header(file, Nt, Np):
    f = file

    f.write('# Geometry exported by OpenPNM'+2*'\n')

    f.write('# Major & minor version'+'\n')
    f.write('0 1'+2*'\n')

    f.write(str(Np+Nt)+' '+'# number of tags'+'\n')
    f.write('# Tags'+'\n')

    for r in range(1, Nt+1):
        tag = 't'+str(int(r))
        tag = str(len(tag))+' '+tag
        f.write(tag+'\n')

    for c in range(1, Np+1):
        tag = 'p'+str(int(c))
        tag = str(len(tag))+' '+tag
        f.write(tag+'\n')

    f.write('\n'+str(Np+Nt)+' '+'# number of types'+'\n')
    f.write('# Types'+'\n')

    for i in range(Np+Nt):
        f.write('3 obj'+'\n')

    f.write('\n')

    return


def rectangles(file, pores1, pores2, alphas, widths):
    f = file

    p1x = pores1[:, 0] + (widths/2)*np.sin(alphas)
    p1y = pores1[:, 1] - (widths/2)*np.cos(alphas)
    p2x = pores2[:, 0] + (widths/2)*np.sin(alphas)
    p2y = pores2[:, 1] - (widths/2)*np.cos(alphas)
    p3x = pores2[:, 0] - (widths/2)*np.sin(alphas)
    p3y = pores2[:, 1] + (widths/2)*np.cos(alphas)
    p4x = pores1[:, 0] - (widths/2)*np.sin(alphas)
    p4y = pores1[:, 1] + (widths/2)*np.cos(alphas)

    for r in range(len(pores1)):
        f.write('# --------- rectangle nbr '+str(r+1)+' ---------'+2*'\n')

        f.write('0 0 1'+'\n')
        f.write('5 Geom2 # class'+'\n')
        f.write('2 # version'+'\n')
        f.write('2 # type'+'\n')
        f.write('1 # voidsLabeled'+'\n')
        f.write('1e-010 # gtol'+'\n')
        f.write('0.0001 # resTol'+2*'\n')

        f.write('4 # number of vertices'+'\n')
        f.write('# Vertices'+'\n')
        f.write('# X Y dom tol'+'\n')
        f.write(str(p1x[r])+' '+str(p1y[r])+' '+'-1 NAN'+'\n')
        f.write(str(p2x[r])+' '+str(p2y[r])+' '+'-1 NAN'+'\n')
        f.write(str(p3x[r])+' '+str(p3y[r])+' '+'-1 NAN'+'\n')
        f.write(str(p4x[r])+' '+str(p4y[r])+' '+'-1 NAN'+2*'\n')

        f.write('4 # number of edges'+'\n')
        f.write('# Edges'+'\n')
        f.write('# vtx1 vtx2 s1 s2 up down mfd tol'+'\n')
        f.write('2 1 0 1 0 1 1 NAN'+'\n')
        f.write('3 2 0 1 0 1 2 NAN'+'\n')
        f.write('4 3 0 1 0 1 3 NAN'+'\n')
        f.write('1 4 0 1 0 1 4 NAN'+2*'\n')

        f.write('4 # number of manifolds'+'\n')
        f.write('# Manifolds'+2*'\n')

        f.write('# Manifold #0'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('1 0 # degrees'+'\n')
        f.write('2 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        f.write(str(p2x[r])+' '+str(p2y[r])+' '+'1'+'\n')
        f.write(str(p1x[r])+' '+str(p1y[r])+' '+'1'+2*'\n')

        f.write('# Manifold #1'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('1 0 # degrees'+'\n')
        f.write('2 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        f.write(str(p3x[r])+' '+str(p3y[r])+' '+'1'+'\n')
        f.write(str(p2x[r])+' '+str(p2y[r])+' '+'1'+2*'\n')

        f.write('# Manifold #2'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('1 0 # degrees'+'\n')
        f.write('2 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        f.write(str(p4x[r])+' '+str(p4y[r])+' '+'1'+'\n')
        f.write(str(p3x[r])+' '+str(p3y[r])+' '+'1'+2*'\n')

        f.write('# Manifold #3'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('1 0 # degrees'+'\n')
        f.write('2 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        f.write(str(p1x[r])+' '+str(p1y[r])+' '+'1'+'\n')
        f.write(str(p4x[r])+' '+str(p4y[r])+' '+'1'+2*'\n')

        f.write('# Attributes'+'\n')
        f.write('0 # nof attributes'+2*'\n')

    return


def circles(file, centers, radii):
    f = file

    p1x = centers[:, 0]-radii
    p1y = centers[:, 1]
    p2x = centers[:, 0]
    p2y = centers[:, 1]-radii
    p3x = centers[:, 0]+radii
    p3y = centers[:, 1]
    p4x = centers[:, 0]
    p4y = centers[:, 1]+radii

    for c in range(len(centers)):
        f.write('# --------- circle nbr '+str(c+1)+' ---------'+2*'\n')

        f.write('0 0 1'+'\n')
        f.write('5 Geom2 # class'+'\n')
        f.write('2 # version'+'\n')
        f.write('2 # type'+'\n')
        f.write('1 # voidsLabeled'+'\n')
        f.write('1e-010 # gtol'+'\n')
        f.write('0.0001 # resTol'+2*'\n')

        f.write('4 # number of vertices'+'\n')
        f.write('# Vertices'+'\n')
        f.write('# X Y dom tol'+'\n')
        # extreme left point of the circle
        f.write(str(p1x[c])+' '+str(p1y[c])+' '+'-1 NAN'+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x[c])+' '+str(p2y[c])+' '+'-1 NAN'+'\n')
        # extreme right point of the circle
        f.write(str(p3x[c])+' '+str(p3y[c])+' '+'-1 NAN'+'\n')
        # extreme top point of the circle
        f.write(str(p4x[c])+' '+str(p4y[c])+' '+'-1 NAN'+2*'\n')

        f.write('4 # number of edges'+'\n')
        f.write('# Edges'+'\n')
        f.write('# vtx1 vtx2 s1 s2 up down mfd tol'+'\n')
        f.write('2 1 0 1 0 1 1 NAN'+'\n')
        f.write('3 2 0 1 0 1 2 NAN'+'\n')
        f.write('4 3 0 1 0 1 3 NAN'+'\n')
        f.write('1 4 0 1 0 1 4 NAN'+2*'\n')

        f.write('4 # number of manifolds'+'\n')
        f.write('# Manifolds'+2*'\n')

        # bottom left quart of the circle
        f.write('# Manifold #0'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('2 0 # degrees'+'\n')
        f.write('3 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x[c])+' '+str(p2y[c])+' '+'1'+'\n')
        f.write(str(p1x[c])+' '+str(p2y[c])+' '+'0.70710678118654746'+'\n')
        # extreme left point of the circle
        f.write(str(p1x[c])+' '+str(p1y[c])+' '+'1'+2*'\n')

        # bottom right quart of the circle
        f.write('# Manifold #1'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('2 0 # degrees'+'\n')
        f.write('3 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        # extreme right point of the circle
        f.write(str(p3x[c])+' '+str(p3y[c])+' '+'1'+'\n')
        f.write(str(p3x[c])+' '+str(p2y[c])+' '+'0.70710678118654746'+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x[c])+' '+str(p2y[c])+' '+'1'+2*'\n')

        # top right quart of the circle
        f.write('# Manifold #2'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('2 0 # degrees'+'\n')
        f.write('3 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        # extreme top point of the circle
        f.write(str(p4x[c])+' '+str(p4y[c])+' '+'1'+'\n')
        f.write(str(p3x[c])+' '+str(p4y[c])+' '+'0.70710678118654746'+'\n')
        # extreme right point of the circle
        f.write(str(p3x[c])+' '+str(p3y[c])+' '+'1'+2*'\n')

        # top left quart of the circle
        f.write('# Manifold #3'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('2 0 # degrees'+'\n')
        f.write('3 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        # extreme left point of the circle
        f.write(str(p1x[c])+' '+str(p1y[c])+' '+'1'+'\n')
        f.write(str(p1x[c])+' '+str(p4y[c])+' '+'0.70710678118654746'+'\n')
        # extreme top point of the circle
        f.write(str(p4x[c])+' '+str(p4y[c])+' '+'1'+2*'\n')

        f.write('# Attributes'+'\n')
        f.write('0 # nof attributes'+2*'\n')

    return


def spheres(file, centers, radii):

    f = file

    p1x = centers[:, 0]-radii
    p1y = centers[:, 1]
    p1z = centers[:, 2]
    p2x = centers[:, 0]
    p2y = centers[:, 1]-radii
    p2z = centers[:, 2]
    p3x = centers[:, 0]
    p3y = centers[:, 1]
    p3z = centers[:, 2]-radii
    p4x = centers[:, 0]
    p4y = centers[:, 1]
    p4z = centers[:, 2]+radii
    p5x = centers[:, 0]
    p5y = centers[:, 1]+radii
    p5z = centers[:, 2]
    p6x = centers[:, 0]+radii
    p6y = centers[:, 1]
    p6z = centers[:, 2]

    for c in range(len(centers)):
        f.write('# --------- Object' + ' ' + str(c+1) + ' ' +
                '----------' + ' ' + '\n')

        f.write('0 0 1' + ' ' + '\n')
        f.write('5 Geom3 # class' + ' ' + '\n')
        f.write('2 # version' + ' ' + '\n')
        f.write('3 # type' + ' ' + '\n')
        f.write('1 # voidsLabeled' + ' ' + '\n')
        f.write('1e-010 # gtol' + ' ' + '\n')
        f.write('0.0001 # resTol' + ' ' + '\n')
        f.write('6 # number of vertices' + ' ' + '\n')
        f.write('# Vertices' + ' ' + '\n')
        f.write('# X Y Z dom tol' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '-1 1e-010' + ' ' + '\n')
        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '-1 1e-010' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '-1 1e-010' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '-1 1e-010' + ' ' + '\n')
        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) +
                ' ' + '-1 1e-010' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) +
                ' ' + '-1 1e-010' + ' ' + '\n')

        f.write('24 # number of parameter vertices' + ' ' + '\n')
        f.write('# Parameter vertices' + ' ' + '\n')
        f.write('# vtx s t fac mfd tol' + ' ' + '\n')
        f.write('1 0 0 -1 1 NAN' + ' ' + '\n')
        f.write('1 0 0 -1 2 NAN' + ' ' + '\n')
        f.write('1 0 0 -1 3 NAN' + ' ' + '\n')
        f.write('1 0 0 -1 4 NAN' + ' ' + '\n')
        f.write('2 1 0 -1 1 NAN' + ' ' + '\n')
        f.write('2 1 0 -1 2 NAN' + ' ' + '\n')
        f.write('2 0 0 -1 5 NAN' + ' ' + '\n')
        f.write('2 0 0 -1 6 NAN' + ' ' + '\n')
        f.write('3 0 1 -1 1 NAN' + ' ' + '\n')
        f.write('3 0 1 -1 3 NAN' + ' ' + '\n')
        f.write('3 0 1 -1 5 NAN' + ' ' + '\n')
        f.write('3 0 1 -1 7 NAN' + ' ' + '\n')
        f.write('4 0 1 -1 2 NAN' + ' ' + '\n')
        f.write('4 0 1 -1 4 NAN' + ' ' + '\n')
        f.write('4 0 1 -1 6 NAN' + ' ' + '\n')
        f.write('4 0 1 -1 8 NAN' + ' ' + '\n')
        f.write('5 1 0 -1 3 NAN' + ' ' + '\n')
        f.write('5 1 0 -1 4 NAN' + ' ' + '\n')
        f.write('5 0 0 -1 7 NAN' + ' ' + '\n')
        f.write('5 0 0 -1 8 NAN' + ' ' + '\n')
        f.write('6 1 0 -1 5 NAN' + ' ' + '\n')
        f.write('6 1 0 -1 6 NAN' + ' ' + '\n')
        f.write('6 1 0 -1 7 NAN' + ' ' + '\n')
        f.write('6 1 0 -1 8 NAN' + ' ' + '\n')

        f.write('12 # number of edges' + ' ' + '\n')
        f.write('# Edges' + ' ' + '\n')
        f.write('# vtx1 vtx2 s1 s2 dom mfd tol' + ' ' + '\n')
        f.write('1 2 0 1 -1 9 NAN' + ' ' + '\n')
        f.write('1 3 0 1 -1 10 NAN' + ' ' + '\n')
        f.write('1 4 0 1 -1 11 NAN' + ' ' + '\n')
        f.write('1 5 0 1 -1 12 NAN' + ' ' + '\n')
        f.write('2 3 0 1 -1 13 NAN' + ' ' + '\n')
        f.write('2 4 0 1 -1 14 NAN' + ' ' + '\n')
        f.write('2 6 0 1 -1 15 NAN' + ' ' + '\n')
        f.write('3 5 0 1 -1 16 NAN' + ' ' + '\n')
        f.write('3 6 0 1 -1 17 NAN' + ' ' + '\n')
        f.write('4 5 0 1 -1 18 NAN' + ' ' + '\n')
        f.write('4 6 0 1 -1 19 NAN' + ' ' + '\n')
        f.write('5 6 0 1 -1 20 NAN' + ' ' + '\n')
        f.write('24 # number of parameter edges' + ' ' + '\n')
        f.write('# Parameter edges' + ' ' + '\n')
        f.write('# edg v1 v2 s1 s2 up down mfdfac mfd tol' + ' ' + '\n')
        f.write('1 1 5 0 1 1 0 1 1 NAN' + ' ' + '\n')
        f.write('1 2 6 0 1 2 0 2 2 NAN' + ' ' + '\n')
        f.write('2 1 9 0 1 0 1 3 1 NAN' + ' ' + '\n')
        f.write('2 3 10 0 1 0 3 4 3 NAN' + ' ' + '\n')
        f.write('3 2 13 0 1 0 2 5 2 NAN' + ' ' + '\n')
        f.write('3 4 14 0 1 0 4 6 4 NAN' + ' ' + '\n')
        f.write('4 3 17 0 1 3 0 7 3 NAN' + ' ' + '\n')
        f.write('4 4 18 0 1 4 0 8 4 NAN' + ' ' + '\n')
        f.write('5 5 9 0 1 1 0 9 1 NAN' + ' ' + '\n')
        f.write('5 7 11 0 1 0 5 10 5 NAN' + ' ' + '\n')
        f.write('6 6 13 0 1 2 0 11 2 NAN' + ' ' + '\n')
        f.write('6 8 15 0 1 0 6 12 6 NAN' + ' ' + '\n')
        f.write('7 7 21 0 1 5 0 13 5 NAN' + ' ' + '\n')
        f.write('7 8 22 0 1 6 0 14 6 NAN' + ' ' + '\n')
        f.write('8 10 17 0 1 0 3 15 3 NAN' + ' ' + '\n')
        f.write('8 12 19 0 1 7 0 16 7 NAN' + ' ' + '\n')
        f.write('9 11 21 0 1 0 5 17 5 NAN' + ' ' + '\n')
        f.write('9 12 23 0 1 0 7 18 7 NAN' + ' ' + '\n')
        f.write('10 14 18 0 1 0 4 19 4 NAN' + ' ' + '\n')
        f.write('10 16 20 0 1 8 0 20 8 NAN' + ' ' + '\n')
        f.write('11 15 22 0 1 0 6 21 6 NAN' + ' ' + '\n')
        f.write('11 16 24 0 1 0 8 22 8 NAN' + ' ' + '\n')
        f.write('12 19 23 0 1 7 0 23 7 NAN' + ' ' + '\n')
        f.write('12 20 24 0 1 8 0 24 8 NAN' + ' ' + '\n')

        f.write('8 # number of faces' + ' ' + '\n')
        f.write('# Faces' + ' ' + '\n')
        f.write('# up down mfd tol' + ' ' + '\n')
        f.write('0 1 -1 NAN' + ' ' + '\n')
        f.write('0 1 2 NAN' + ' ' + '\n')
        f.write('0 1 3 NAN' + ' ' + '\n')
        f.write('0 1 -4 NAN' + ' ' + '\n')
        f.write('0 1 -5 NAN' + ' ' + '\n')
        f.write('0 1 6 NAN' + ' ' + '\n')
        f.write('0 1 7 NAN' + ' ' + '\n')
        f.write('0 1 -8 NAN' + ' ' + '\n')

        f.write('20 # number of 3D manifolds' + ' ' + '\n')
        f.write('# Manifolds' + ' ' + '\n')

        f.write('# Manifold #0' + ' ' + '\n')

        f.write('10 BezierSurf # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 2 # degrees' + ' ' + '\n')
        f.write('9 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p1x[c]) + ' ' + str(p1y[c]-radii[c]) + ' ' +
                str(p1z[c]) + ' ' + '0.70710678118654746'+'\n')
        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p3x[c]-radii[c]) + ' ' + str(p3y[c]) + ' ' +
                str(p3z[c]) + ' ' + '0.70710678118654746'+'\n')
        f.write(str(p3x[c]-radii[c]) + ' ' + str(p3y[c]-radii[c]) +
                ' ' + str(p3z[c]) + ' ' + '0.5' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]-radii[c]) + ' ' +
                str(p3z[c]) + ' ' + '0.70710678118654746'+'\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #1' + ' ' + '\n')

        f.write('10 BezierSurf # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 2 # degrees' + ' ' + '\n')
        f.write('9 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p1x[c]) + ' ' + str(p1y[c]-radii[c]) + ' ' +
                str(p1z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p4x[c]-radii[c]) + ' ' + str(p4y[c]) + ' ' +
                str(p4z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]-radii[c]) + ' ' + str(p4y[c]-radii[c]) +
                ' ' + str(p4z[c]) + ' ' + '0.5' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]-radii[c]) + ' ' +
                str(p4z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #2' + ' ' + '\n')

        f.write('10 BezierSurf # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 2 # degrees' + ' ' + '\n')
        f.write('9 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p1x[c]) + ' ' + str(p1y[c]+radii[c]) + ' ' +
                str(p1z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p3x[c]-radii[c]) + ' ' + str(p3y[c]) + ' ' +
                str(p3z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]-radii[c]) + ' ' + str(p3y[c]+radii[c]) +
                ' ' + str(p3z[c]) + ' ' + '0.5' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]+radii[c]) + ' ' +
                str(p3z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #3' + ' ' + '\n')

        f.write('10 BezierSurf # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 2 # degrees' + ' ' + '\n')
        f.write('9 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p1x[c]) + ' ' + str(p1y[c]+radii[c]) + ' ' +
                str(p1z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p4x[c]-radii[c]) + ' ' + str(p4y[c]) + ' ' +
                str(p4z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]-radii[c]) + ' ' + str(p4y[c]+radii[c]) +
                ' ' + str(p4z[c]) + ' ' + '0.5' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]+radii[c]) + ' ' +
                str(p4z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #4' + ' ' + '\n')

        f.write('10 BezierSurf # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 2 # degrees' + ' ' + '\n')
        f.write('9 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p2x[c]+radii[c]) + ' ' + str(p2y[c]) + ' ' +
                str(p2z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]-radii[c]) + ' ' +
                str(p3z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]+radii[c]) + ' ' + str(p3y[c]-radii[c]) +
                ' ' + str(p3z[c]) + ' ' + '0.5' + ' ' + '\n')
        f.write(str(p3x[c]+radii[c]) + ' ' + str(p3y[c]) + ' ' +
                str(p3z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #5' + ' ' + '\n')

        f.write('10 BezierSurf # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 2 # degrees' + ' ' + '\n')
        f.write('9 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p2x[c]+radii[c]) + ' ' + str(p2y[c]) + ' ' +
                str(p2z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]-radii[c]) + ' ' +
                str(p4z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]+radii[c]) + ' ' + str(p4y[c]-radii[c]) +
                ' ' + str(p4z[c]) + ' ' + '0.5' + ' ' + '\n')
        f.write(str(p4x[c]+radii[c]) + ' ' + str(p4y[c]) + ' ' +
                str(p4z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #6' + ' ' + '\n')

        f.write('10 BezierSurf # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 2 # degrees' + ' ' + '\n')
        f.write('9 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p5x[c]+radii[c]) + ' ' + str(p5y[c]) + ' ' +
                str(p5z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]+radii[c]) + ' ' +
                str(p3z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]+radii[c]) + ' ' + str(p3y[c]+radii[c]) +
                ' ' + str(p3z[c]) + ' ' + '0.5' + ' ' + '\n')
        f.write(str(p3x[c]+radii[c]) + ' ' + str(p3y[c]) + ' ' +
                str(p3z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #7' + ' ' + '\n')

        f.write('10 BezierSurf # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 2 # degrees' + ' ' + '\n')
        f.write('9 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p5x[c]+radii[c]) + ' ' + str(p5y[c]) + ' ' +
                str(p5z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]+radii[c]) + ' ' +
                str(p4z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]+radii[c]) + ' ' + str(p4y[c]+radii[c]) +
                ' ' + str(p4z[c]) + ' ' + '0.5' + ' ' + '\n')
        f.write(str(p4x[c]+radii[c]) + ' ' + str(p4y[c]) + ' ' +
                str(p4z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #8' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p1x[c]) + ' ' + str(p1y[c]-radii[c]) + ' ' +
                str(p1z[c]) + ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '1' + ' ' + '\n')

        f.write('# Manifold #9' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]-radii[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #10' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]+radii[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #11' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p1x[c]) + ' ' + str(p1y[c]+radii[c]) + ' ' + str(p1z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #12' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]-radii[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #13' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '1' + ' ' + '\n')
        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]+radii[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #14' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) + ' ' +
                '1' + ' ' + '\n')
        f.write(str(p2x[c]+radii[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #15' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) + ' ' +
                '1' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]+radii[c]) + ' ' + str(p3z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #16' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) + ' ' +
                '1' + ' ' + '\n')
        f.write(str(p3x[c]+radii[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #17' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) + ' ' +
                '1' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]+radii[c]) + ' ' + str(p4z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #18' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) + ' ' +
                '1' + ' ' + '\n')
        f.write(str(p4x[c]+radii[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('# Manifold #19' + ' ' + '\n')

        f.write('11 BezierCurve # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('3 # sdim' + ' ' + '\n')
        f.write('0 3 1 # transformation' + ' ' + '\n')
        f.write('2 0 # degrees' + ' ' + '\n')
        f.write('3 # number of control points' + ' ' + '\n')
        f.write('# control point coords and weights' + ' ' + '\n')

        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) + ' ' +
                '1' + ' ' + '\n')
        f.write(str(p5x[c]+radii[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) +
                ' ' + '0.70710678118654746' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) + ' ' +
                '1' + ' ' + '\n')

        f.write('24 # number of parameter curves' + ' ' + '\n')
        f.write('# Parameter curves' + ' ' + '\n')

        f.write('# Parameter curve #0' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')

        f.write('2 50 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 0 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #1' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0.020408163265306121 ' + ' ' +
                '0 0.040816326530612242 0 0.061224489795918366 ' + ' ' +
                '0 0.081632653061224483 0 0.1020408163265306 ' + ' ' +
                '0 0.12244897959183673 0 0.14285714285714285 ' + ' ' +
                '0 0.16326530612244897 0 0.18367346938775508 ' + ' ' +
                '0 0.2040816326530612 0 0.22448979591836732 ' + ' ' +
                '0 0.24489795918367346 0 0.26530612244897955 ' + ' ' +
                '0 0.2857142857142857 0 0.30612244897959179 ' + ' ' +
                '0 0.32653061224489793 0 0.34693877551020408 ' + ' ' +
                '0 0.36734693877551017 0 0.38775510204081631 ' + ' ' +
                '0 0.4081632653061224 0 0.42857142857142855 ' + ' ' +
                '0 0.44897959183673464 0 0.46938775510204078 ' + ' ' +
                '0 0.48979591836734693 0 0.51020408163265307 ' + ' ' +
                '0 0.53061224489795911 0 0.55102040816326525 ' + ' ' +
                '0 0.5714285714285714 0 0.59183673469387754 ' + ' ' +
                '0 0.61224489795918358 0 0.63265306122448972 ' + ' ' +
                '0 0.65306122448979587 0 0.67346938775510201 ' + ' ' +
                '0 0.69387755102040816 0 0.71428571428571419 ' + ' ' +
                '0 0.73469387755102034 0 0.75510204081632648 ' + ' ' +
                '0 0.77551020408163263 0 0.79591836734693866 ' + ' ' +
                '0 0.81632653061224481 0 0.83673469387755095 ' + ' ' +
                '0 0.8571428571428571 0 0.87755102040816324 ' + ' ' +
                '0 0.89795918367346927 0 0.91836734693877542 ' + ' ' +
                '0 0.93877551020408156 0 0.95918367346938771 ' + ' ' +
                '0 0.97959183673469385 0 ' + ' ' +
                '0.99999999999999989 0 # chain' + ' ' + '\n')

        f.write('# Parameter curve #2' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0 0.020408163265306121 ' + ' ' +
                '0 0.040816326530612242 0 0.061224489795918366 ' + ' ' +
                '0 0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #3' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 ' + ' ' +
                '0 0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 ' + ' ' +
                '0 0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #4' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #5' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0 0.020408163265306121 ' + ' ' +
                '0 0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #6' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 0 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #7' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 0 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #8' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 1 0 1 0.020408163265306121 1 ' + ' ' +
                '0.040816326530612242 1 0.061224489795918366 1 ' + ' ' +
                '0.081632653061224483 1 0.1020408163265306 1 ' + ' ' +
                '0.12244897959183673 1 0.14285714285714285 1 ' + ' ' +
                '0.16326530612244897 1 0.18367346938775508 1 ' + ' ' +
                '0.2040816326530612 1 0.22448979591836732 1 ' + ' ' +
                '0.24489795918367346 1 0.26530612244897955 1 ' + ' ' +
                '0.2857142857142857 1 0.30612244897959179 1 ' + ' ' +
                '0.32653061224489793 1 0.34693877551020408 1 ' + ' ' +
                '0.36734693877551017 1 0.38775510204081631 1 ' + ' ' +
                '0.4081632653061224 1 0.42857142857142855 1 ' + ' ' +
                '0.44897959183673464 1 0.46938775510204078 1 ' + ' ' +
                '0.48979591836734693 1 0.51020408163265307 1 ' + ' ' +
                '0.53061224489795911 1 0.55102040816326525 1 ' + ' ' +
                '0.5714285714285714 1 0.59183673469387754 1 ' + ' ' +
                '0.61224489795918358 1 0.63265306122448972 1 ' + ' ' +
                '0.65306122448979587 1 0.67346938775510201 1 ' + ' ' +
                '0.69387755102040816 1 0.71428571428571419 1 ' + ' ' +
                '0.73469387755102034 1 0.75510204081632648 1 ' + ' ' +
                '0.77551020408163263 1 0.79591836734693866 1 ' + ' ' +
                '0.81632653061224481 1 0.83673469387755095 1 ' + ' ' +
                '0.8571428571428571 1 0.87755102040816324 1 ' + ' ' +
                '0.89795918367346927 1 0.91836734693877542 1 ' + ' ' +
                '0.93877551020408156 1 0.95918367346938771 1 ' + ' ' +
                '0.97959183673469385 1 0.99999999999999989 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #9' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #10' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 1 0 1 0.020408163265306121 1 ' + ' ' +
                '0.040816326530612242 1 0.061224489795918366 1 ' + ' ' +
                '0.081632653061224483 1 0.1020408163265306 1 ' + ' ' +
                '0.12244897959183673 1 0.14285714285714285 1 ' + ' ' +
                '0.16326530612244897 1 0.18367346938775508 1 ' + ' ' +
                '0.2040816326530612 1 0.22448979591836732 1 ' + ' ' +
                '0.24489795918367346 1 0.26530612244897955 1 ' + ' ' +
                '0.2857142857142857 1 0.30612244897959179 1 ' + ' ' +
                '0.32653061224489793 1 0.34693877551020408 1 ' + ' ' +
                '0.36734693877551017 1 0.38775510204081631 1 ' + ' ' +
                '0.4081632653061224 1 0.42857142857142855 1 ' + ' ' +
                '0.44897959183673464 1 0.46938775510204078 1 ' + ' ' +
                '0.48979591836734693 1 0.51020408163265307 1 ' + ' ' +
                '0.53061224489795911 1 0.55102040816326525 1 ' + ' ' +
                '0.5714285714285714 1 0.59183673469387754 1 ' + ' ' +
                '0.61224489795918358 1 0.63265306122448972 1 ' + ' ' +
                '0.65306122448979587 1 0.67346938775510201 1 ' + ' ' +
                '0.69387755102040816 1 0.71428571428571419 1 ' + ' ' +
                '0.73469387755102034 1 0.75510204081632648 1 ' + ' ' +
                '0.77551020408163263 1 0.79591836734693866 1 ' + ' ' +
                '0.81632653061224481 1 0.83673469387755095 1 ' + ' ' +
                '0.8571428571428571 1 0.87755102040816324 1 ' + ' ' +
                '0.89795918367346927 1 0.91836734693877542 1 ' + ' ' +
                '0.93877551020408156 1 0.95918367346938771 1 ' + ' ' +
                '0.97959183673469385 1 0.99999999999999989 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #11' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #12' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 0 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #13' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 0 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #14' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 1 1 1 0.97959183673469385 1 ' + ' ' +
                '0.95918367346938771 1 0.93877551020408168 1 ' + ' ' +
                '0.91836734693877553 1 0.89795918367346939 1 ' + ' ' +
                '0.87755102040816324 1 0.85714285714285721 1 ' + ' ' +
                '0.83673469387755106 1 0.81632653061224492 1 ' + ' ' +
                '0.79591836734693877 1 0.77551020408163263 1 ' + ' ' +
                '0.75510204081632648 1 0.73469387755102045 1 ' + ' ' +
                '0.7142857142857143 1 0.69387755102040827 1 ' + ' ' +
                '0.67346938775510212 1 0.65306122448979598 1 ' + ' ' +
                '0.63265306122448983 1 0.61224489795918369 1 ' + ' ' +
                '0.59183673469387754 1 0.5714285714285714 1 ' + ' ' +
                '0.55102040816326536 1 0.53061224489795922 1 ' + ' ' +
                '0.51020408163265307 1 0.48979591836734693 1 ' + ' ' +
                '0.46938775510204089 1 0.44897959183673475 1 ' + ' ' +
                '0.4285714285714286 1 0.40816326530612246 1 ' + ' ' +
                '0.38775510204081642 1 0.36734693877551028 1 ' + ' ' +
                '0.34693877551020413 1 0.32653061224489799 1 ' + ' ' +
                '0.30612244897959184 1 0.28571428571428581 1 ' + ' ' +
                '0.26530612244897966 1 0.24489795918367352 1 ' + ' ' +
                '0.22448979591836737 1 0.20408163265306134 1 ' + ' ' +
                '0.18367346938775519 1 0.16326530612244905 1 ' + ' ' +
                '0.1428571428571429 1 0.12244897959183676 1 ' + ' ' +
                '0.10204081632653073 1 0.08163265306122458 1 ' + ' ' +
                '0.061224489795918435 1 0.04081632653061229 1 ' + ' ' +
                '0.020408163265306145 1 1.1102230246251565e-016 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #15' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 1 0 0.97959183673469385 0 ' + ' ' +
                '0.95918367346938771 0 0.93877551020408168 0 ' + ' ' +
                '0.91836734693877553 0 0.89795918367346939 0 ' + ' ' +
                '0.87755102040816324 0 0.85714285714285721 0 ' + ' ' +
                '0.83673469387755106 0 0.81632653061224492 0 ' + ' ' +
                '0.79591836734693877 0 0.77551020408163263 0 ' + ' ' +
                '0.75510204081632648 0 0.73469387755102045 0 ' + ' ' +
                '0.7142857142857143 0 0.69387755102040827 0 ' + ' ' +
                '0.67346938775510212 0 0.65306122448979598 0 ' + ' ' +
                '0.63265306122448983 0 0.61224489795918369 0 ' + ' ' +
                '0.59183673469387754 0 0.5714285714285714 0 ' + ' ' +
                '0.55102040816326536 0 0.53061224489795922 0 ' + ' ' +
                '0.51020408163265307 0 0.48979591836734693 0 ' + ' ' +
                '0.46938775510204089 0 0.44897959183673475 0 ' + ' ' +
                '0.4285714285714286 0 0.40816326530612246 0 ' + ' ' +
                '0.38775510204081642 0 0.36734693877551028 0 ' + ' ' +
                '0.34693877551020413 0 0.32653061224489799 0 ' + ' ' +
                '0.30612244897959184 0 0.28571428571428581 0 ' + ' ' +
                '0.26530612244897966 0 0.24489795918367352 0 ' + ' ' +
                '0.22448979591836737 0 0.20408163265306134 0 ' + ' ' +
                '0.18367346938775519 0 0.16326530612244905 0 ' + ' ' +
                '0.1428571428571429 0 0.12244897959183676 0 ' + ' ' +
                '0.10204081632653073 0 0.08163265306122458 0 ' + ' ' +
                '0.061224489795918435 0 0.04081632653061229 0 ' + ' ' +
                '0.020408163265306145 0 1.1102230246251565e-016 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #16' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 1 1 1 0.97959183673469385 1 ' + ' ' +
                '0.95918367346938771 1 0.93877551020408168 1 ' + ' ' +
                '0.91836734693877553 1 0.89795918367346939 1 ' + ' ' +
                '0.87755102040816324 1 0.85714285714285721 1 ' + ' ' +
                '0.83673469387755106 1 0.81632653061224492 1 ' + ' ' +
                '0.79591836734693877 1 0.77551020408163263 1 ' + ' ' +
                '0.75510204081632648 1 0.73469387755102045 1 ' + ' ' +
                '0.7142857142857143 1 0.69387755102040827 1 ' + ' ' +
                '0.67346938775510212 1 0.65306122448979598 1 ' + ' ' +
                '0.63265306122448983 1 0.61224489795918369 1 ' + ' ' +
                '0.59183673469387754 1 0.5714285714285714 1 ' + ' ' +
                '0.55102040816326536 1 0.53061224489795922 1 ' + ' ' +
                '0.51020408163265307 1 0.48979591836734693 1 ' + ' ' +
                '0.46938775510204089 1 0.44897959183673475 1 ' + ' ' +
                '0.4285714285714286 1 0.40816326530612246 1 ' + ' ' +
                '0.38775510204081642 1 0.36734693877551028 1 ' + ' ' +
                '0.34693877551020413 1 0.32653061224489799 1 ' + ' ' +
                '0.30612244897959184 1 0.28571428571428581 1 ' + ' ' +
                '0.26530612244897966 1 0.24489795918367352 1 ' + ' ' +
                '0.22448979591836737 1 0.20408163265306134 1 ' + ' ' +
                '0.18367346938775519 1 0.16326530612244905 1 ' + ' ' +
                '0.1428571428571429 1 0.12244897959183676 1 ' + ' ' +
                '0.10204081632653073 1 0.08163265306122458 1 ' + ' ' +
                '0.061224489795918435 1 0.04081632653061229 1 ' + ' ' +
                '0.020408163265306145 1 1.1102230246251565e-016 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #17' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 1 1 1 0.97959183673469385 1 ' + ' ' +
                '0.95918367346938771 1 0.93877551020408168 1 ' + ' ' +
                '0.91836734693877553 1 0.89795918367346939 1 ' + ' ' +
                '0.87755102040816324 1 0.85714285714285721 1 ' + ' ' +
                '0.83673469387755106 1 0.81632653061224492 1 ' + ' ' +
                '0.79591836734693877 1 0.77551020408163263 1 ' + ' ' +
                '0.75510204081632648 1 0.73469387755102045 1 ' + ' ' +
                '0.7142857142857143 1 0.69387755102040827 1 ' + ' ' +
                '0.67346938775510212 1 0.65306122448979598 1 ' + ' ' +
                '0.63265306122448983 1 0.61224489795918369 1 ' + ' ' +
                '0.59183673469387754 1 0.5714285714285714 1 ' + ' ' +
                '0.55102040816326536 1 0.53061224489795922 1 ' + ' ' +
                '0.51020408163265307 1 0.48979591836734693 1 ' + ' ' +
                '0.46938775510204089 1 0.44897959183673475 1 ' + ' ' +
                '0.4285714285714286 1 0.40816326530612246 1 ' + ' ' +
                '0.38775510204081642 1 0.36734693877551028 1 ' + ' ' +
                '0.34693877551020413 1 0.32653061224489799 1 ' + ' ' +
                '0.30612244897959184 1 0.28571428571428581 1 ' + ' ' +
                '0.26530612244897966 1 0.24489795918367352 1 ' + ' ' +
                '0.22448979591836737 1 0.20408163265306134 1 ' + ' ' +
                '0.18367346938775519 1 0.16326530612244905 1 ' + ' ' +
                '0.1428571428571429 1 0.12244897959183676 1 ' + ' ' +
                '0.10204081632653073 1 0.08163265306122458 1 ' + ' ' +
                '0.061224489795918435 1 0.04081632653061229 1 ' + ' ' +
                '0.020408163265306145 1 1.1102230246251565e-016 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #18' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 1 1 1 0.97959183673469385 1 ' + ' ' +
                '0.95918367346938771 1 0.93877551020408168 1 ' + ' ' +
                '0.91836734693877553 1 0.89795918367346939 1 ' + ' ' +
                '0.87755102040816324 1 0.85714285714285721 1 ' + ' ' +
                '0.83673469387755106 1 0.81632653061224492 1 ' + ' ' +
                '0.79591836734693877 1 0.77551020408163263 1 ' + ' ' +
                '0.75510204081632648 1 0.73469387755102045 1 ' + ' ' +
                '0.7142857142857143 1 0.69387755102040827 1 ' + ' ' +
                '0.67346938775510212 1 0.65306122448979598 1 ' + ' ' +
                '0.63265306122448983 1 0.61224489795918369 1 ' + ' ' +
                '0.59183673469387754 1 0.5714285714285714 1 ' + ' ' +
                '0.55102040816326536 1 0.53061224489795922 1 ' + ' ' +
                '0.51020408163265307 1 0.48979591836734693 1 ' + ' ' +
                '0.46938775510204089 1 0.44897959183673475 1 ' + ' ' +
                '0.4285714285714286 1 0.40816326530612246 1 ' + ' ' +
                '0.38775510204081642 1 0.36734693877551028 1 ' + ' ' +
                '0.34693877551020413 1 0.32653061224489799 1 ' + ' ' +
                '0.30612244897959184 1 0.28571428571428581 1 ' + ' ' +
                '0.26530612244897966 1 0.24489795918367352 1 ' + ' ' +
                '0.22448979591836737 1 0.20408163265306134 1 ' + ' ' +
                '0.18367346938775519 1 0.16326530612244905 1 ' + ' ' +
                '0.1428571428571429 1 0.12244897959183676 1 ' + ' ' +
                '0.10204081632653073 1 0.08163265306122458 1 ' + ' ' +
                '0.061224489795918435 1 0.04081632653061229 1 ' + ' ' +
                '0.020408163265306145 1 1.1102230246251565e-016 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #19' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 1 0 0.97959183673469385 0 ' + ' ' +
                '0.95918367346938771 0 0.93877551020408168 0 ' + ' ' +
                '0.91836734693877553 0 0.89795918367346939 0 ' + ' ' +
                '0.87755102040816324 0 0.85714285714285721 0 ' + ' ' +
                '0.83673469387755106 0 0.81632653061224492 0 ' + ' ' +
                '0.79591836734693877 0 0.77551020408163263 0 ' + ' ' +
                '0.75510204081632648 0 0.73469387755102045 0 ' + ' ' +
                '0.7142857142857143 0 0.69387755102040827 0 ' + ' ' +
                '0.67346938775510212 0 0.65306122448979598 0 ' + ' ' +
                '0.63265306122448983 0 0.61224489795918369 0 ' + ' ' +
                '0.59183673469387754 0 0.5714285714285714 0 ' + ' ' +
                '0.55102040816326536 0 0.53061224489795922 0 ' + ' ' +
                '0.51020408163265307 0 0.48979591836734693 0 ' + ' ' +
                '0.46938775510204089 0 0.44897959183673475 0 ' + ' ' +
                '0.4285714285714286 0 0.40816326530612246 0 ' + ' ' +
                '0.38775510204081642 0 0.36734693877551028 0 ' + ' ' +
                '0.34693877551020413 0 0.32653061224489799 0 ' + ' ' +
                '0.30612244897959184 0 0.28571428571428581 0 ' + ' ' +
                '0.26530612244897966 0 0.24489795918367352 0 ' + ' ' +
                '0.22448979591836737 0 0.20408163265306134 0 ' + ' ' +
                '0.18367346938775519 0 0.16326530612244905 0 ' + ' ' +
                '0.1428571428571429 0 0.12244897959183676 0 ' + ' ' +
                '0.10204081632653073 0 0.08163265306122458 0 ' + ' ' +
                '0.061224489795918435 0 0.04081632653061229 0 ' + ' ' +
                '0.020408163265306145 0 1.1102230246251565e-016 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #20' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 1 1 1 0.97959183673469385 1 ' + ' ' +
                '0.95918367346938771 1 0.93877551020408168 1 ' + ' ' +
                '0.91836734693877553 1 0.89795918367346939 1 ' + ' ' +
                '0.87755102040816324 1 0.85714285714285721 1 ' + ' ' +
                '0.83673469387755106 1 0.81632653061224492 1 ' + ' ' +
                '0.79591836734693877 1 0.77551020408163263 1 ' + ' ' +
                '0.75510204081632648 1 0.73469387755102045 1 ' + ' ' +
                '0.7142857142857143 1 0.69387755102040827 1 ' + ' ' +
                '0.67346938775510212 1 0.65306122448979598 1 ' + ' ' +
                '0.63265306122448983 1 0.61224489795918369 1 ' + ' ' +
                '0.59183673469387754 1 0.5714285714285714 1 ' + ' ' +
                '0.55102040816326536 1 0.53061224489795922 1 ' + ' ' +
                '0.51020408163265307 1 0.48979591836734693 1 ' + ' ' +
                '0.46938775510204089 1 0.44897959183673475 1 ' + ' ' +
                '0.4285714285714286 1 0.40816326530612246 1 ' + ' ' +
                '0.38775510204081642 1 0.36734693877551028 1 ' + ' ' +
                '0.34693877551020413 1 0.32653061224489799 1 ' + ' ' +
                '0.30612244897959184 1 0.28571428571428581 1 ' + ' ' +
                '0.26530612244897966 1 0.24489795918367352 1 ' + ' ' +
                '0.22448979591836737 1 0.20408163265306134 1 ' + ' ' +
                '0.18367346938775519 1 0.16326530612244905 1 ' + ' ' +
                '0.1428571428571429 1 0.12244897959183676 1 ' + ' ' +
                '0.10204081632653073 1 0.08163265306122458 1 ' + ' ' +
                '0.061224489795918435 1 0.04081632653061229 1 ' + ' ' +
                '0.020408163265306145 1 1.1102230246251565e-016 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #21' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 1 1 1 0.97959183673469385 1 ' + ' ' +
                '0.95918367346938771 1 0.93877551020408168 1 ' + ' ' +
                '0.91836734693877553 1 0.89795918367346939 1 ' + ' ' +
                '0.87755102040816324 1 0.85714285714285721 1 ' + ' ' +
                '0.83673469387755106 1 0.81632653061224492 1 ' + ' ' +
                '0.79591836734693877 1 0.77551020408163263 1 ' + ' ' +
                '0.75510204081632648 1 0.73469387755102045 1 ' + ' ' +
                '0.7142857142857143 1 0.69387755102040827 1 ' + ' ' +
                '0.67346938775510212 1 0.65306122448979598 1 ' + ' ' +
                '0.63265306122448983 1 0.61224489795918369 1 ' + ' ' +
                '0.59183673469387754 1 0.5714285714285714 1 ' + ' ' +
                '0.55102040816326536 1 0.53061224489795922 1 ' + ' ' +
                '0.51020408163265307 1 0.48979591836734693 1 ' + ' ' +
                '0.46938775510204089 1 0.44897959183673475 1 ' + ' ' +
                '0.4285714285714286 1 0.40816326530612246 1 ' + ' ' +
                '0.38775510204081642 1 0.36734693877551028 1 ' + ' ' +
                '0.34693877551020413 1 0.32653061224489799 1 ' + ' ' +
                '0.30612244897959184 1 0.28571428571428581 1 ' + ' ' +
                '0.26530612244897966 1 0.24489795918367352 1 ' + ' ' +
                '0.22448979591836737 1 0.20408163265306134 1 ' + ' ' +
                '0.18367346938775519 1 0.16326530612244905 1 ' + ' ' +
                '0.1428571428571429 1 0.12244897959183676 1 ' + ' ' +
                '0.10204081632653073 1 0.08163265306122458 1 ' + ' ' +
                '0.061224489795918435 1 0.04081632653061229 1 ' + ' ' +
                '0.020408163265306145 1 1.1102230246251565e-016 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #22' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 0 # chain' +
                ' ' + '\n')

        f.write('# Parameter curve #23' + ' ' + '\n')

        f.write('8 PolChain # class' + ' ' + '\n')
        f.write('0 0 # version' + ' ' + '\n')
        f.write('2 # sdim' + ' ' + '\n')
        f.write('0 2 1 # transformation' + ' ' + '\n')
        f.write('2 50 0 0 0.020408163265306121 0 ' + ' ' +
                '0.040816326530612242 0 0.061224489795918366 0 ' + ' ' +
                '0.081632653061224483 0 0.1020408163265306 0 ' + ' ' +
                '0.12244897959183673 0 0.14285714285714285 0 ' + ' ' +
                '0.16326530612244897 0 0.18367346938775508 0 ' + ' ' +
                '0.2040816326530612 0 0.22448979591836732 0 ' + ' ' +
                '0.24489795918367346 0 0.26530612244897955 0 ' + ' ' +
                '0.2857142857142857 0 0.30612244897959179 0 ' + ' ' +
                '0.32653061224489793 0 0.34693877551020408 0 ' + ' ' +
                '0.36734693877551017 0 0.38775510204081631 0 ' + ' ' +
                '0.4081632653061224 0 0.42857142857142855 0 ' + ' ' +
                '0.44897959183673464 0 0.46938775510204078 0 ' + ' ' +
                '0.48979591836734693 0 0.51020408163265307 0 ' + ' ' +
                '0.53061224489795911 0 0.55102040816326525 0 ' + ' ' +
                '0.5714285714285714 0 0.59183673469387754 0 ' + ' ' +
                '0.61224489795918358 0 0.63265306122448972 0 ' + ' ' +
                '0.65306122448979587 0 0.67346938775510201 0 ' + ' ' +
                '0.69387755102040816 0 0.71428571428571419 0 ' + ' ' +
                '0.73469387755102034 0 0.75510204081632648 0 ' + ' ' +
                '0.77551020408163263 0 0.79591836734693866 0 ' + ' ' +
                '0.81632653061224481 0 0.83673469387755095 0 ' + ' ' +
                '0.8571428571428571 0 0.87755102040816324 0 ' + ' ' +
                '0.89795918367346927 0 0.91836734693877542 0 ' + ' ' +
                '0.93877551020408156 0 0.95918367346938771 0 ' + ' ' +
                '0.97959183673469385 0 0.99999999999999989 0 # chain' +
                ' ' + '\n')

        f.write('# Attributes' + ' ' + '\n')
        f.write('0 # nof attributes' + ' ' + '\n')

    return


def cylinders(file, centers, radii):

    f = file

    p1x = centers[:, 0]-radii
    p1y = centers[:, 1]
    p1z = centers[:, 2]
    p2x = centers[:, 0]
    p2y = centers[:, 1]-radii
    p2z = centers[:, 2]
    p3x = centers[:, 0]
    p3y = centers[:, 1]
    p3z = centers[:, 2]-radii
    p4x = centers[:, 0]
    p4y = centers[:, 1]
    p4z = centers[:, 2]+radii
    p5x = centers[:, 0]
    p5y = centers[:, 1]+radii
    p5z = centers[:, 2]
    p6x = centers[:, 0]+radii
    p6y = centers[:, 1]
    p6z = centers[:, 2]

    p7x = centers[:, 0]+radii
    p7y = centers[:, 1]
    p7z = centers[:, 2]

    p8x = centers[:, 0]+radii
    p8y = centers[:, 1]
    p8z = centers[:, 2]

#0 -20 0 -1 NAN 
#0 0 -20 -1 NAN 
#0 0 20 -1 NAN 
#0 20 0 -1 NAN 
#
#80 -20 0 -1 NAN 
#80 0 -20 -1 NAN 
#80 0 20 -1 NAN 
#80 20 0 -1 NAN 

    for c in range(len(centers)):
        f.write('# --------- Object' + ' ' + str(c+1) + ' ' +
                '----------' + ' ' + '\n')

        f.write('0 0 1' + ' ' + '\n')
        f.write('5 Geom3 # class' + ' ' + '\n')
        f.write('2 # version' + ' ' + '\n')
        f.write('3 # type' + ' ' + '\n')
        f.write('1 # voidsLabeled' + ' ' + '\n')
        f.write('1e-010 # gtol' + ' ' + '\n')
        f.write('0.0001 # resTol' + ' ' + '\n')
        f.write('8 # number of vertices' + ' ' + '\n')
        f.write('# Vertices' + ' ' + '\n')
        f.write('# X Y Z dom tol' + ' ' + '\n')

        f.write(str(p1x[c]) + ' ' + str(p1y[c]) + ' ' + str(p1z[c]) +
                ' ' + '-1 NAN' + ' ' + '\n')
        f.write(str(p2x[c]) + ' ' + str(p2y[c]) + ' ' + str(p2z[c]) +
                ' ' + '-1 NAN' + ' ' + '\n')
        f.write(str(p3x[c]) + ' ' + str(p3y[c]) + ' ' + str(p3z[c]) +
                ' ' + '-1 NAN' + ' ' + '\n')
        f.write(str(p4x[c]) + ' ' + str(p4y[c]) + ' ' + str(p4z[c]) +
                ' ' + '-1 NAN' + ' ' + '\n')
        f.write(str(p5x[c]) + ' ' + str(p5y[c]) + ' ' + str(p5z[c]) +
                ' ' + '-1 NAN' + ' ' + '\n')
        f.write(str(p6x[c]) + ' ' + str(p6y[c]) + ' ' + str(p6z[c]) +
                ' ' + '-1 NAN' + ' ' + '\n')
        f.write(str(p7x[c]) + ' ' + str(p7y[c]) + ' ' + str(p7z[c]) +
                ' ' + '-1 NAN' + ' ' + '\n')
        f.write(str(p8x[c]) + ' ' + str(p8y[c]) + ' ' + str(p8z[c]) +
                ' ' + '-1 NAN' + ' ' + '\n')

    return

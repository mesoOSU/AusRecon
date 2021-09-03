# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 13:14:14 2021

@author: arger
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 12:49:17 2021

Heads up, this is really just a mock-up right now of one way to do this. 
Still up for debate.

@author: arger
"""

import numpy as np
import h5py
from scipy.spatial.transform import Rotation as R
import pyvista as pv
import math
#import meshpy
#import tetgen

# ===================== #
# HELPER FUNCTIONS (collapse and ignore)
# ===================== #

def h5disp(f, quiet=False):
    """ disply hdf5 data

    **h5disp** displays the groups, datasets, and attributes in an h5 file
    similar to the same command in Matlab

    Parameters
    ----------
    f :: h5py file class

    quiet :: True/False (False default), determines whether written to console

    Returns
    -------
    out :: text output string (optional)

    Example
    -------
    >>> import h5py
    >>> f = h5py.File('my_file.h5')
    >>> out = h5disp(f)
    """

    out = '\n\nHDF5 file: {0}\n----\n'.format(f.filename)

    mylist = []
    def func(name, obj):
        if isinstance(obj, h5py.Dataset):
            mylist.append(name)
    f.visititems(func) # get dataset and group names

    oldparent = ''
    tabs = ''
    rootGroupUsed = 0
    for item in mylist:
        newparent = f[item].parent.name
        if newparent == '/':
            tabs = ''
        else:
            tabs = '    ' * (len(f[mylist[0]].name.split('/')) - 1)

        if newparent != oldparent:
            if newparent == '/':
                if rootGroupUsed == 0:
                    out = out + '{0}Group \'{1}\'\n'.format(tabs, newparent)
                    rootGroupUsed = 1
                oldparent = newparent
            else:
                out = out + '{0}Group \'{1}\'\n'.format(tabs, newparent)
                oldparent = newparent


        out = out + "{0}    Dataset \'{1}\': shape={2}\n".format(\
                              tabs, f[item].name.split('/')[-1], f[item].shape)

        x = h5py.AttributeManager(f[item]).items()
        if len(x) != 0:
            out = out + "{0}        Attributes:\n".format(tabs)
            for a in x:
                out = out+ "{0}            \'{1}\': {2}\n".format(\
                                                              tabs, a[0], a[1])
    if quiet==False:
        print(out)

    return out

def Defaults(HDF5_location):
    DD = h5py.File('Default_Data.h5','r')
    try:
        dat =     np.asanyarray(DD[HDF5_location])
        DD.close()
        return(dat)
    except:
        print('Error: Requested Default data does not exist')
        DD.close()
        return()

# Make a static map of the contents od DD, so I can open and close it, and still
# know the addresses of the data I care about
DD = h5py.File('Default_Data.h5','r')
h5_map = h5disp(DD,quiet = True)
DD.close()

# ===================== #
#        OPTIONS
# ===================== #

fr_refine =   3#  % refinement level on FR
sp_refine =  10#  % refinement level on sphere
per_fiber = 100#  % points per fiber
pf_hkls = np.array([[1,1,1],[1,0,0],[1,1 ,0]])

wsopts = {
    'MakePoleFigures':  pf_hkls,
    'PointsPerFiber': per_fiber,
    'MakeFRL2IP':True,
    'MakeFRH1IP': True,
    'MakeSphL2IP':True,
    'MakeSphH1IP':True};

#cubic_fundamental_region_crd


# First off, recreate the Mesh structures from OdfPF, plus the Extended stuff
# Joel was writing.

def get_symmetry_elements(sym):
    """Some day, this will be able to return any symmetry group either from
    notation or from a set of generating functions. for now though, it just
    returns one of the three symmetry groups already defined in OdfPf
    """
    if np. any(np.array(['m3m','cub','m3_m','Cubic','cubic']) == sym):
        return(R(Defaults("Symmetries/Cubic")))
    elif np. any(np.array(['Hex','hex','hexagonal','Hexagonal','622']) == sym):
        return(R(Defaults("Symmetries/Hex")))
    elif np. any(np.array(['Ortho','Orthorhombic','ortho','orthorhombic','222']) == sym):
        return(R(Defaults("Symmetries/Hex")))
    else:
        return('DNE')
        

class FR_Mesh(object):
    def __init__(self,crd = 'default',con='default',eqv ='default',sym = 'm3m'):
# ==== Default settings ==== #
        if all([crd == 'default',con== 'default',eqv== 'default']):
            print("FR_Mesh was given no data, assuming Cubic Base Mesh")
            con = Defaults("Base_Meshes/Cubic/con")
            crd = Defaults("Base_Meshes/Cubic/crd")
            eqv = Defaults("Base_Meshes/Cubic/eqv")
            sym = 'm3m' # Lazy fix
        else:
            # if not ALL default, assert NONE are default (no half measures)
            assert any([crd == 'default',con== 'default',eqv== 'default']) == False, """
Either all or none of con,crd, and eqv must be defined when initializing FR_Mesh"""[1:]
# ==== Make sure crd,con, and eqv make sense  ==== #
# TODO: should add a check here to make sure they really can build an actual
# mesh and aren't nonsensical values
# TODO: since refine mesh will change these, I should either make it create a
# new FR_Mesh every time, OR I should move all these flips/asserts into internal
# functions that are called during initialization. leave for now.
        # arrange crd,con, and eqv into n x3 numpy arrays
        self.crd= np.asanyarray(crd).astype(float)
        if self.crd.shape[0] == 3 and not np.all(np.array(self.crd.shape)==3):
            self.crd = self.crd.T # should only trigger if Transpose needed
        self.con= np.asanyarray(con).astype(int)
        if self.con.shape[0] == 3 and not np.all(np.array(self.con.shape)==3):
            self.con = self.con.T # should only trigger if Transpose needed
        self.eqv= np.asanyarray(eqv).astype(int)
        if self.eqv.shape[0] == 3 and not np.all(np.array(self.eqv.shape)==3):
            self.eqv = self.eqv.T # should only trigger if Transpose needed
# ==== Define Symmetry elements using internal function ==== #
        self.set_symmetry(sym)

# ==== Functions that modify crd,con,eqv, and sym. === #

    def set_symmetry(self,sym):
        # redefine symmetry elements 
#TODO: ability to check that symmetry elements make sense (ie, actual complete
# group) not implemted yet. for now, just pass in name of sym group (eg 'm3m')
#TODO: might just move this outisde of the class and just make this a function call
        if type(sym) ==str:
            self.sym = get_symmetry_elements(sym)
        elif type(sym) == R:
            self.sym = sym
        else:
            self.sym =R(sym)
        assert type(get_symmetry_elements(sym)) == R ,"""
Not a valid option for sym. provide either a str-type name, a Rotation object,
or an n-by-4 list of numbers interpretable as quaternions."""
    def refine(self,refinement_level, relative = True):
#TODO: Mesh refinement is involved and might get used elsewhere as well, so
# this will be a call to an outside function.
        # or relative to current mesh
        self.con,self.crd,self.eqv = Refine_Mesh(self.con,self.crd,self.eqv)
        

# ==== Past here, initial object is made. Now build functions to fill in ==== #
# ====== other values that require calculation, do refinement, etc ====== #
        

        
    # SO, from here I am assuming not everyone always wants to calculate everything,
    # with that in mind, I'm adding functions to calc each individual thing in 
    # Joel's list of "things we care about", adding a CALC_ALL function,
    # and throwing in some if_then_else loops later to calc the things if they
    # don't exist, and also to recalc them if they don't line up (for instance,
    # if someone did a remesh but didn't recalculate values post mesh). Should 
    # Also probably make a warning trigger if data is updated AFTER the most recent
    # calculations of something. 
    def vol(self):
        if hasattr(self,'_volume') == False:
            # add assert that calc time follows last relevant update.
            self._volume = sum(self.l2ip)
        return(self._volume)

    def l2ip(self):
        if hasattr(self,'_l2ip') == False:
            # add assert that calc time follows last relevant update.
            l = len(self.crd) - len(self.con)
            self._l2ip = np.zeros([l,l])
            # add actual math here later
        return(self._l2ip)

    def _Calc_NpQpGradMatrix(self):
        [self._grad,self._MetricGij]=['do','math here']
        
    def grad(self,manual_input = 'None'):
        if manual_input != 'None':
            # Do some asserts here to make sure the input given isn't garbage
            self._grad = manual_input
        if hasattr(self,'_grad') == False:
            # add assert that calc time follows last relevant update.
            # Still confused on what this is.
            self._Calc_NpQpGradMatrix()
        return(self._l2ip)

    def metricGij(self,manual_input = 'None'):
        if manual_input != 'None':
            # Do some asserts here to make sure the input given isn't garbage
            self._grad = manual_input
        if hasattr(self,'_metricGij') == False:
            # add assert that calc time follows last relevant update.
            # Still confused on what this is.
            self._Calc_NpQpGradMatrix()
        return(self._metricGij)


    def h2ip(self):
        if hasattr(self,'_h2ip') == False:
            # add assert that calc time follows last relevant update.
            l = len(self.crd) - len(self.con)
            self._h2ip = np.zeros([l,l])
            # add actual math here later
        return(self._l2ip)
    
    def calc_nind(self):
        self.nind = len(self.crd) - len(self.con)
        
    def calc_ntot(self):
        self.ntot = len(self.crd)
        
class Cubic_Mesh(FR_Mesh):
    def __init__(self,crd = 'default',con='default'):
        FR_Mesh.__init__(self,crd,con)
        self.symm = "add cubic symm here"
        # Maybe add something here to fore symmetry to remain unchange? 
        # Perhaps to overwrite self.symm? might make problems later with inheritance



def Refine_Mesh(crd,con,eqv,level):
    assert level == 1, "this is a placeholder function. only works with level =1"
    return(crd,con,eqv)


# con = Defaults("Base_Meshes/Cubic/con")
# crd = Defaults("Base_Meshes/Cubic/crd")
# eqv = Defaults("Base_Meshes/Cubic/eqv")

# from pyvista import examples

# mesh = examples.load_hexbeam()
# bcpos = [(6.20, 3.00, 7.50),
#          (0.16, 0.13, 2.65),
#          (-0.28, 0.94, -0.21)]

# pl = pv.Plotter()
# pl.add_mesh(mesh, show_edges=True, color='white')
# pl.add_mesh(pv.PolyData(mesh.points), color='red',
#        point_size=10, render_points_as_spheres=True)
# pl.camera_position = bcpos
# pl.show()


# surf = pv.PolyData(crd,con)







# print("Testing mesh in Open3D...")
# mesh = o3dtut.get_knot_mesh()
# print(mesh)
# print('Vertices:')
# print(np.asarray(mesh.vertices))
# print('Triangles:')
# print(np.asarray(mesh.triangles))









def euler_to_quaternion(r):
    (yaw, pitch, roll) = (r[0], r[1], r[2])
    qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
    qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    return [qx, qy, qz, qw]

def quaternion_to_euler(q):
    (x, y, z, w) = (q[0], q[1], q[2], q[3])
    t0 = +2.0 * (w * x + y * z)
    t1 = +1.0 - 2.0 * (x * x + y * y)
    roll = math.atan2(t0, t1)
    t2 = +2.0 * (w * y - z * x)
    t2 = +1.0 if t2 > +1.0 else t2
    t2 = -1.0 if t2 < -1.0 else t2
    pitch = math.asin(t2)
    t3 = +2.0 * (w * z + x * y)
    t4 = +1.0 - 2.0 * (y * y + z * z)
    yaw = math.atan2(t3, t4)
    return [yaw, pitch, roll]


# import pyvista as pv
# import tetgen
# import numpy as np
# pv.set_plot_theme('document')

# sphere = pv.Sphere()
# tet = tetgen.TetGen(sphere)
# tet.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
# grid = tet.grid
# grid.plot(show_edges=True)



# # get cell centroids
# cells = grid.cells.reshape(-1, 5)[:, 1:]
# cell_center = grid.points[cells].mean(1)

# # extract cells below the 0 xy plane
# mask = cell_center[:, 2] < 0
# cell_ind = mask.nonzero()[0]
# subgrid = grid.extract_cells(cell_ind)

# # advanced plotting
# plotter = pv.Plotter()
# plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
# plotter.add_mesh(sphere, 'r', 'wireframe')
# plotter.add_legend([[' Input Mesh ', 'r'],
#                     [' Tesselated Mesh ', 'black']])
# plotter.show()

## THis is wrong, Spheres are their own thing# 
# class Sphere_Mesh(FR_Mesh):
#     def __init__(self,crd = 'default',con='default'):
#         FR_Mesh.__init__(self,crd,con)
#         self.symm = "add sphere symm here"
#         # Maybe add something here to fore symmetry to remain unchange? 
#         # Perhaps to overwrite self.symm? might make problems later with inheritance
#     def l2ip(self):
#         # have to overwrite with spherical-specific equation
#         if hasattr(self,'_l2ip') == False:
#             # add assert that calc time follows last relevant update.
#             l = len(self.crd) - len(self.con)
#             self._l2ip = np.zeros([l,l])
#             # add actual math here later
#         return(self._l2ip)
### THese classes also need print  functions and asserts.

#class symm()






# ========================================================================== #
# Goofing around and remembering how to do classes.     
    
class Person(object):
    def __init__(self, age,name):
        self.age = age
        self.name  = name
    def __str__(self):
        return"<"+self.name+","+str(self.age)+">"
    def set_age(self,age):
        self.age=age
    def get_age(self):
        return self.age
    def age_diff(self,other):
        diff =self.age - other.age
        return abs(diff)
class Married_Person(Person):
    def __init__(self,age,name,spouse, numper_of_children=0):
        Person.__init__(self,age,name)
        self.spouse_name = spouse.name
        self.numper_of_children = numper_of_children
        
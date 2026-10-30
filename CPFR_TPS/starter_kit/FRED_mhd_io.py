######################################################################
#
#         FRED interface for mhd input/output using SimpleITK
#
#         NB: important differences wrt to ITK conventions:
#         1) length units in FRED are cm
#         2) memory layout is column-major (F-ordering)
#         3) FRED uses "bounding box corner" instead of "centroid of first voxel"
#
######################################################################
import os
import numpy as np
import SimpleITK as sitk

def readMap(fname):
    img = sitk.ReadImage(fname)
    Map  = sitk.GetArrayFromImage(img)
    nn = np.array(img.GetSize(),dtype=np.int32)
    hs = np.array(img.GetSpacing(),dtype=np.float32)
    offset = np.array(img.GetOrigin(),dtype=np.float32)
    direction = np.array(img.GetDirection(),dtype=np.float32)
    R=direction.reshape(3,3) # rotation matrix
    x0 = offset - R.dot(0.5*hs)
    Map = Map.reshape(-1,order='A').reshape(nn,order='F') # change ordering, do not change memory
    hs /= 10 # mm -> cm
    x0 /= 10 # mm -> cm
    return nn,hs,x0,Map,direction
########################################################################
def writeMap(fname,hs,x0,Map,direction=None,splitFile=None):
    if direction is None: direction=(1,0,0,0,1,0,0,0,1)
    if splitFile is None: splitFile=False
    hs=np.array(hs)
    x0=np.array(x0)
    nn = Map.shape
    if Map.flags.f_contiguous:
        img = sitk.GetImageFromArray(Map.reshape(-1,order='A').reshape(nn[2],nn[1],nn[0],order='C'))
    else:
        img = sitk.GetImageFromArray(np.asfortranarray(Map).reshape(-1,order='A').reshape(nn[2],nn[1],nn[0],order='C'))
    img.SetSpacing(10*hs.astype(np.double))
    _dir=np.array(direction).astype(np.double)
    img.SetDirection(_dir)
    R=_dir.reshape(3,3)
    offset = x0 + R.dot(0.5*hs)
    img.SetOrigin(10*offset.astype(np.double))
    if fname.endswith('mhd') and not splitFile: # all in one mhd file
        sitk.WriteImage(img,'tmp.mha')
        os.rename('tmp.mha',fname)
    else:
        sitk.WriteImage(img,fname)
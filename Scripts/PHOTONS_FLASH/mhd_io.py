#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# read and write 3d maps in mhd format 
# Sep 2019 A.Schiavi and N. Krah
########################################################################
########################################################################
import numpy as np
import sys,os,re
########################################################################
MET_TYPES = ['MET_DOUBLE','MET_FLOAT','MET_INT','MET_LONG','MET_SHORT','MET_CHAR','MET_UCHAR','MET_UINT','MET_ULONG','MET_USHORT']
NUMPY_DTYPES = [np.float64,np.float32,np.int32,np.int64,np.int16,np.int8,np.uint8,np.uint32,np.uint64,np.uint16]
def mhd_read(fname):
    if fname[-4:] not in ['.mhd','.mha']:
        # print('fname',fname)
        if os.path.exists(fname+'.mhd'):
            fname = fname+'.mhd'
        elif os.path.exists(fname+'.mha'):
            fname = fname+'.mha'
        else:
            print('Error: file path extension is not mhd or mha')
            exit(-1)
    try:
        header = mhd_getHeader(fname)
        voxels=getMinimalMHD(fname,header)

        if getParam('BinaryData',header,1)!='True':
            print('BinaryData must be True');exit(1)

        if getParam('BinaryDataByteOrderMSB',header,1)!='False':
            print('BinaryDataByteOrderMSB must be False');exit(1)

        if getParam('CompressedData',header,1)!='False':
            print('CompressedData must be False');exit(1)

        ElementSpacing = getParam('ElementSpacing',header,3)
        if ElementSpacing==None:
            ElementSpacing=np.ones(3)
        else:
            try:
                ElementSpacing = np.array(list(map(float,ElementSpacing)))
            except:
                print('ElementSpacing is not a floating point list: ',ElementSpacing)

        TransformMatrix = getParam('TransformMatrix',header,9)
        if TransformMatrix==None:
            TransformMatrix=np.array([1,0,0,0,1,0,0,0,1])
        else:
            try:
                TransformMatrix = np.array(list(map(float,TransformMatrix)))
            except:
                print('TransformMatrix is not a floating point list: ',TransformMatrix)
        voxels['TransformMatrix']=TransformMatrix

        Offset = getParam('Offset',header,3)
        if Offset==None:
            Offset=ElementSpacing*0.5
        else:
            try:
                Offset = np.array(list(map(float,Offset)))
            except:
                print('Offset is not a floating point list: ',Offset)
        
        # from Offset to x0
        x0 = Offset
        Delta = ElementSpacing*0.5
        x0 -= Delta[0]*voxels['TransformMatrix'][0:3]
        x0 -= Delta[1]*voxels['TransformMatrix'][3:6]
        x0 -= Delta[2]*voxels['TransformMatrix'][6:9]

        voxels['hs']=ElementSpacing/10. # mm => cm
        voxels['x0']=x0/10. # mm => cm


        # CenterOfRotation = getParam('CenterOfRotation',header,3)
        # if CenterOfRotation==None:
        #     CenterOfRotation=np.zeros(3)
        # else:
        #     try:
        #         CenterOfRotation = np.array(list(map(float,CenterOfRotation)))
        #     except:
        #         print('CenterOfRotation is not a floating point list: ',CenterOfRotation)
        # voxels['CenterOfRotation']=CenterOfRotation

        return voxels
    except:
        print('IO error: cannot read mhd map from file:', fname)
    return None

########################################################################
def getParam(pname,lines,nvals):
    for l in lines:
        if l.strip()=='':
            continue
        tokens = l.split()
        # print(tokens)
        if len(tokens)<3 or tokens[1]!='=':
            print('syntax error in mhd file at line: ',l)
            exit(-1)
        if tokens[0]==pname:
            if nvals==1:
                return tokens[2]
            else:
                if len(tokens)!=2+nvals:
                    print('not enough values for parameter',pname)
                    exit(-1)
                return tokens[2:]
    return None
########################################################################
def mhd_getHeader(fname):
    header=[]
    try:
        fin = open(fname,'rb')
        while True:
            line = fin.readline().decode('utf-8').strip()
            # print(line)
            header.append(line)
            m = re.match(r'ElementDataFile\s*=\s*(.*)\s*$',line)
            if m:
                break
    except:
        print('IO error: cannot open input file or read the header', fname)
        return None    
    return header
########################################################################
def mhd_getData(fname,ElementDataFile):
    buffer=None
    if ElementDataFile=='LOCAL':
        try:
            fin = open(fname,'rb')
            while True:
                line = fin.readline().decode('utf-8').strip()
                # print(line)
                m = re.match(r'ElementDataFile\s*=\s*(.*)\s*$',line)
                if m:
                    if m.groups()[0]=='LOCAL':
                        buffer = fin.read()
                    break
        except:
            print('IO error: cannot open input file or read data from it', fname)
            return None
        return buffer
    else: # external file
        try:
            rawFilePath = os.path.join(os.path.dirname(fname),ElementDataFile)
            fin  = open(rawFilePath,'rb')
            buffer = fin.read()
        except:
            print('IO error: cannot open or read input file',ElementDataFile)
            return None
        return buffer
########################################################################
def getMinimalMHD(fname,lines):
    if getParam('ObjectType',lines,1)!='Image':
        print('Error: ObjectType != Image'); exit(-1)

    NDims = getParam('NDims',lines,1)
    if NDims==None:
        print('Error: NDims not found'); exit(-1)
    try:
        NDims = int(NDims)
    except:
        print('NDims is not an integer: ',NDims)
    
    DimSize = getParam('DimSize',lines,NDims)
    if DimSize==None:
        print('Error: DimSize not found'); exit(-1)
    try:
        DimSize = np.array(list(map(int,DimSize)))
    except:
        print('DimSize is not an integer list: ',DimSize)

    ElementType = getParam('ElementType',lines,1)
    if ElementType==None:
        print('Error: ElementType not found'); exit(-1)
    if ElementType not in MET_TYPES:
        print('ElementType ',ElementType,'is not one of',MET_TYPES)

    ElementDataFile = getParam('ElementDataFile',lines[-1:],1)
    if ElementDataFile==None:
        print('Error: ElementDataFile must be last tag of MetaImageHeader'); exit(-1)
    # print(ElementDataFile)

    voxels={}
    voxels['ElementDataFile']=ElementDataFile

    buffer = mhd_getData(fname,ElementDataFile)
    
    if buffer!=None:
        try:
            voxels['Map'] = np.frombuffer(buffer,dtype=NUMPY_DTYPES[MET_TYPES.index(ElementType)])
        except:
            print('IO error: cannot convert buffer to map')
    else:
        print('IO error: cannot open or read input file',ElementDataFile)
        exit(1)

    if voxels['Map'].size != np.prod(DimSize):
        print('Error: size of data in raw file does not match: read ',voxels['Map'].size,' elements instead of',np.prod(DimSize))
        exit(1)
    voxels['Map'] = np.reshape(voxels['Map'],DimSize,order='F')
    # voxels['Map'] = np.transpose(voxels['Map'])
    
    voxels['nn']=np.array(DimSize)
    return voxels


########################################################################
def mhd_write(fname,voxels):

    if len(fname)<5 or fname[-4:]!='.mhd': # add mhd extension if not present already
        fname=fname+'.mhd'

    Offset = voxels['x0'].copy()
    Delta = voxels['hs']/2.0
    Offset += Delta[0]*voxels['TransformMatrix'][0:3]
    Offset += Delta[1]*voxels['TransformMatrix'][3:6]
    Offset += Delta[2]*voxels['TransformMatrix'][6:9]
    # print('x0',voxels['x0'])
    # print('Offset',Offset)
    

    try:
        fout = open(fname,'w')
        print('ObjectType = Image',file=fout)
        print('NDims =',voxels['Map'].ndim,file=fout)
        print('DimSize =',' '.join(map(str,voxels['Map'].shape)),file=fout)
        print('ElementType =', MET_TYPES[NUMPY_DTYPES.index(voxels['Map'].dtype)],file=fout)
        print('ElementSpacing =',' '.join(map(str,voxels['hs']*10.)),file=fout)
        print('Offset =',' '.join(map(str,Offset*10.)),file=fout)
        print('BinaryData = True',file=fout)
        print('BinaryDataByteOrderMSB = False',file=fout)
        print('CompressedData = False',file=fout)
        
        # print('DimSize =',' '.join(map(str,voxels['Map'].shape)))

        for k in voxels.keys():
            if k in ['ElementDataFile','Map','nn','hs','x0']:
                continue
            print(k,'=',' '.join(map(str,voxels[k])),file=fout)
        # print('ElementDataFile =',voxels['ElementDataFile'],file=fout)
        print('ElementDataFile = LOCAL',file=fout)
        fout.close()
    except:
        print('IO error: cannot save mhd header to file:', fname)

    try:
        fout = open(fname,'ab')
        fout.write(voxels['Map'].tobytes(order='F'))
        fout.close()
    except:
        print('IO error: cannot save voxel map to file:', voxels['ElementDataFile'])

########################################################################
def packVoxels(nn,hs,x0,Map):
    voxels={}

    voxels['nn'] = nn.copy()
    voxels['hs'] = hs.copy()
    voxels['x0'] = x0.copy()
    voxels['Map']= Map.copy()
    voxels['TransformMatrix'] = np.array([1,0,0,0,1,0,0,0,1])

    return voxels
########################################################################
def unpackVoxels(voxels):
    return [voxels['nn'],voxels['hs'],voxels['x0'],voxels['Map']]
########################################################################
if __name__ == "__main__":
    if len(sys.argv)<2:
        print('usage: filename')
        exit(1)
    print('\n\nreading.....\n\n')
    voxels = mhd_read(sys.argv[1])
    print('dims=',voxels['nn'])
    print('spacing=',voxels['hs'])
    print('offset=',voxels['x0'])
    x0 = voxels['x0']
    hs = voxels['hs']
    nn = voxels['nn']
    print('x-spatial extent= [',x0[0],',',x0[0]+nn[0]*hs[0],'] => ',nn[0]*hs[0])
    print('y-spatial extent= [',x0[1],',',x0[1]+nn[1]*hs[1],'] => ',nn[1]*hs[1])
    print('z-spatial extent= [',x0[2],',',x0[2]+nn[2]*hs[2],'] => ',nn[2]*hs[2])
    print('datatype=',voxels['Map'].dtype )
    print('range=',np.min(voxels['Map']),np.max(voxels['Map']))
    print('sum=',np.sum(voxels['Map']))

    
    # Map32=Map.astype(np.float32)
    print('\n\nwriting.....\n\n')
    mhd_write('Map32.mhd',voxels)


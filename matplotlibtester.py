#!/usr/bin/env python3
import os,sys
import matplotlib
import numpy as np
########################################################################
# Nov 2024 aschiavi updated version
########################################################################
try:
    matplotlib.use(os.environ['FRED_MATPLOTLIB_BACKEND'])
except:
    pass

# Get the version of matplotlib as a tuple of integers
version = tuple(map(int, matplotlib.__version__.split(".")))
print('matplotlib version:',version)

# Check if the version is greater than 3.9.0
if version > (3, 9, 0):
    availBackends = matplotlib.backends.backend_registry.list_builtin(matplotlib.backends.BackendFilter.INTERACTIVE)
else:
    availBackends = matplotlib.rcsetup.interactive_bk

if len(sys.argv)>1:
    try:
        matplotlib.use(sys.argv[1])
    except:
        print('error: requested backend',sys.argv[1],'not found!!!!')
        print('Available interactive backends:',availBackends)
        exit(1)

print('Available interactive backends:',availBackends)
print('Current backend:', matplotlib.get_backend())
########################################################################
def on_key(event):
    exit(0)

import pylab as plt

fig, (ax1,ax2) = plt.subplots(2)

ax1.plot([0,1,2,3,4,5,6,7,8,9,10],[0,1,4,9,16,25,36,49,64,81,100],'r-o')
ax1.set_title(f'Matplotlib tester using {matplotlib.get_backend()} backend')

img=np.random.random((30,50))
ax2.imshow(img,cmap='jet')

plt.gcf().canvas.mpl_connect('key_press_event', on_key)
plt.show()

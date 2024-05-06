import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from num2tex import num2tex
import sys
plt.rcParams['axes.formatter.limits'] = -3, 3
# First create some toy data:

data=np.genfromtxt("temporary.raw",delimiter="")
print(data[0], data.shape)

xslices = 400
yslices = 400
costhslices = 1
cthbins=100
extent=(0,20,0,20)

time = data[0]


fig = plt.figure(figsize=(17, 7))
#plt.text("testxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",3,2)
#fig.suptitle("      t="+str(data[0]), fontsize=14)
#plt.tight_layout(pad=0.05, w_pad=0.001, h_pad=2.0)
datafin=np.zeros(xslices*yslices,dtype=float)
for i in range(xslices):
    for j in range(yslices):
        for k in range(0,1):
            #dataini[i*yslices+j]=data[0,1+8*i*yslices*costhslices+j*costhslices*8+8*k]
            datafin[i*yslices+j]=data[1+8*i*yslices*costhslices+j*costhslices*8+8*k]
datafin=datafin.reshape(xslices,yslices)
plt.subplot(2, 4, 1)
plt.imshow(datafin/cthbins,origin="lower",extent=extent)
#plt.text("test",3,2)
plt.colorbar()
plt.title("$\\rho_{ee}$",fontsize=14)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_visible(False)
frame1.axes.yaxis.set_visible(False)
#plt.clim(0,1.6)

datafin=np.zeros(xslices*yslices,dtype=float)
for i in range(xslices):
    for j in range(yslices):
        for k in range(0,1):
            #dataini[i*yslices+j]=data[0,1+8*i*yslices*costhslices+j*costhslices*8+8*k]
            datafin[i*yslices+j]=data[2+8*i*yslices*costhslices+j*costhslices*8+8*k]
datafin=datafin.reshape(xslices,yslices)
plt.subplot(2, 4, 2)
plt.imshow(datafin/cthbins,origin="lower",extent=extent)
plt.colorbar()
plt.title("$\\rho_{xx}$",fontsize=14)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_visible(False)
frame1.axes.yaxis.set_visible(False)

datafin=np.zeros(xslices*yslices,dtype=float)
for i in range(xslices):
    for j in range(yslices):
        for k in range(0,1):
            #dataini[i*yslices+j]=data[0,1+8*i*yslices*costhslices+j*costhslices*8+8*k]
            datafin[i*yslices+j]=data[5+8*i*yslices*costhslices+j*costhslices*8+8*k]
datafin=datafin.reshape(xslices,yslices)
plt.subplot(2, 4, 3)
plt.imshow(datafin/cthbins,origin="lower",extent=extent)
plt.colorbar()
plt.title("$\\bar{\\rho}_{ee}$",fontsize=14)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_visible(False)
frame1.axes.yaxis.set_visible(False)

datafin=np.zeros(xslices*yslices,dtype=float)
for i in range(xslices):
    for j in range(yslices):
        for k in range(0,1):
            #dataini[i*yslices+j]=data[0,1+8*i*yslices*costhslices+j*costhslices*8+8*k]
            datafin[i*yslices+j]=data[6+8*i*yslices*costhslices+j*costhslices*8+8*k]
datafin=datafin.reshape(xslices,yslices)
plt.subplot(2, 4, 4)
plt.imshow(datafin/cthbins,origin="lower",extent=extent)
plt.colorbar()
plt.title("$\\bar{\\rho}_{xx}$",fontsize=14)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_visible(False)
frame1.axes.yaxis.set_visible(False)

datafin=np.zeros(xslices*yslices,dtype=float)
for i in range(xslices):
    for j in range(yslices):
        for k in range(0,1):
            #dataini[i*yslices+j]=data[0,1+8*i*yslices*costhslices+j*costhslices*8+8*k]
            datafin[i*yslices+j]=data[3+8*i*yslices*costhslices+j*costhslices*8+8*k]
datafin=datafin.reshape(xslices,yslices)
plt.subplot(2, 4, 5)
plt.imshow(datafin/cthbins,origin="lower",extent=extent)
plt.colorbar()
plt.title("Re($\\rho_{ex}$)",fontsize=14)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_visible(False)
frame1.axes.yaxis.set_visible(False)


datafin=np.zeros(xslices*yslices,dtype=float)
for i in range(xslices):
    for j in range(yslices):
        for k in range(0,1):
            #dataini[i*yslices+j]=data[0,1+8*i*yslices*costhslices+j*costhslices*8+8*k]
            datafin[i*yslices+j]=data[4+8*i*yslices*costhslices+j*costhslices*8+8*k]
datafin=datafin.reshape(xslices,yslices)
plt.subplot(2, 4, 6)
plt.imshow(datafin/cthbins,origin="lower")
plt.colorbar()
plt.title("Im($\\rho_{ex}$)",fontsize=14)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_visible(False)
frame1.axes.yaxis.set_visible(False)

datafin=np.zeros(xslices*yslices,dtype=float)
for i in range(xslices):
    for j in range(yslices):
        for k in range(0,1):
            #dataini[i*yslices+j]=data[0,1+8*i*yslices*costhslices+j*costhslices*8+8*k]
            datafin[i*yslices+j]=data[7+8*i*yslices*costhslices+j*costhslices*8+8*k]
datafin=datafin.reshape(xslices,yslices)
plt.subplot(2, 4, 7)
plt.imshow(datafin/cthbins,origin="lower",extent=extent)
plt.colorbar()
plt.title("Re($\\bar{\\rho}_{ex}$)",fontsize=14)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_visible(False)
frame1.axes.yaxis.set_visible(False)

datafin=np.zeros(xslices*yslices,dtype=float)
for i in range(xslices):
    for j in range(yslices):
        for k in range(0,1):
            #dataini[i*yslices+j]=data[0,1+8*i*yslices*costhslices+j*costhslices*8+8*k]
            datafin[i*yslices+j]=data[8+8*i*yslices*costhslices+j*costhslices*8+8*k]
datafin=datafin.reshape(xslices,yslices)
plt.subplot(2, 4, 8)
plt.imshow(datafin/cthbins,origin="lower",extent=extent)
plt.colorbar()
plt.title("Im($\\bar{\\rho}_{ex}$)",fontsize=14)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_visible(False)
frame1.axes.yaxis.set_visible(False)

plt.subplots_adjust(left = 0.01,right = 0.96,bottom = 0.1,top = 0.9,wspace=-0.15)
time=str(time).split('e')
if(len(time)==1):
    fig.text(0.45,0.96,"$t=\ "+time[0]+"$ Seconds",fontsize=16)
else:
    fig.text(0.45,0.96,"$t =\ "+time[0]+"\\times 10^{"+str(time[1])+"}$ Seconds",fontsize=16)
#fig.text(0.5,0.5,"test",fontsize=16)

#plt.tight_layout()
#fname="image.pdf"
fname = 'image'+str(sys.argv[1])+'.pdf'
plt.savefig(fname)

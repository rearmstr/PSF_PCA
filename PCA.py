from matplotlib.patches import Circle, Wedge, Polygon,Rectangle
from matplotlib.collections import PatchCollection
import matplotlib
import pylab
import copy
import numpy as np


# Central coordinates of each CCD in mm
N7=["N7",16.908,191.670]
N6=["N6",16.908,127.780]
N5=["N5",16.908,63.890]
N4=["N4",16.908,0.]
N3=["N3",16.908,-63.890]
N2=["N2",16.908,-127.780]
N1=["N1",16.908,-191.670]

N13=["N13",50.724,159.725]
N12=["N12",50.724,95.835]
N11=["N11",50.724,31.945]
N10=["N10",50.724,-31.945]
N9=["N9",50.724,-95.835]
N8=["N8",50.724,-159.725]

N19=["N19",84.540,159.725]
N18=["N18",84.540,95.835]
N17=["N17",84.540,31.945]
N16=["N16",84.540,-31.945]
N15=["N15",84.540,-95.835]
N14=["N14",84.540,-159.725]
   
N24=["N24",118.356,127.780]
N23=["N23",118.356,63.890]
N22=["N22",118.356,0.]
N21=["N21",118.356,-63.890]
N20=["N20",118.356,-127.780]

N28=["N28",152.172,95.835]
N27=["N27",152.172,31.945]
N26=["N26",152.172,-31.945]
N25=["N25",152.172,-95.835]

N31=["N31",185.988,63.890]
N30=["N30",185.988,0.]
N29=["N29",185.988,-63.890]

S7=["S7",-16.908,191.670]
S6=["S6",-16.908,127.780]
S5=["S5",-16.908,63.890]
S4=["S4",-16.908,0.]
S3=["S3",-16.908,-63.890]
S2=["S2",-16.908,-127.780]
S1=["S1",-16.908,-191.670]

S13=["S13",-50.724,159.725]
S12=["S12",-50.724,95.835]
S11=["S11",-50.724,31.945]
S10=["S10",-50.724,-31.945]
S9=["S9",-50.724,-95.835]
S8=["S8",-50.724,-159.725]

S19=["S19",-84.540,159.725]
S18=["S18",-84.540,95.835]
S17=["S17",-84.540,31.945]
S16=["S16",-84.540,-31.945]
S15=["S15",-84.540,-95.835]
S14=["S14",-84.540,-159.725]

S24=["S24",-118.356,127.780]
S23=["S23",-118.356,63.890]
S22=["S22",-118.356,0.]
S21=["S21",-118.356,-63.890]
S20=["S20",-118.356,-127.780]

S28=["S28",-152.172,95.835]
S27=["S27",-152.172,31.945]
S26=["S26",-152.172,-31.945]
S25=["S25",-152.172,-95.835]

S31=["S31",-185.988,63.890]
S30=["S30",-185.988,0.]
S29=["S29",-185.988,-63.890]



# order of chips given in numeric order
ccdid = [S29,S30,S31,S25,S26,S27,S28,S20,S21,S22,S23,S24,S14,S15,S16,S17,S18,S19,S8,S9,S10,S11,S12,S13,S1,S2,S3,S4,S5,S6,S7,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31]

# defines the size of a chip.  One pixel=15 microns
xsize=2048*15e-6*1000
ysize=4096*15e-6*1000


# probably not labeled correclty, but it is the centers
# of the ccds
corners={}
for i,ext in enumerate(ccdid):
    xy=[]
    xy.append( (ext[1])-xsize/2)
    xy.append( (ext[2])-ysize/2)
    corners[i+1]=xy

patches=[]
for i,corner in corners.items():
    patches.append(Rectangle(corner,xsize,ysize))


def toFocal(ccd,x,y):
    xc=corners[ccd][0]
    yc=corners[ccd][1]
    return x*15e-6*1000+xc,y*15e-6*1000+yc


import pyfits



class PCAResult:
    def __init__(self,filename):
        self.filename=filename
        self.exps=[]
        self.bounds=[]
        self.skip61=True
        self.patches=[]
        for i,corner in corners.items():
            self.patches.append(Rectangle(corner,xsize,ysize))
        self.outline = PatchCollection(self.patches,facecolor='None')
        self.readFits()
        self.grid = PatchCollection(self.bounds, cmap=matplotlib.cm.jet,
                                    alpha=1,edgecolors=None,linewidths=0)

    def readFits(self):
        pyfile=pyfits.open(self.filename)
        self.exps=pyfile['exps'].data['exposure']
        self.ccd=pyfile['exps'].header['ccd']
        self.nvar=pyfile['exps'].header['nvar']
        self.nx=pyfile['exps'].header['nx']
        self.ny=pyfile['exps'].header['ny']
        self.type=pyfile['exps'].header['type']
        self.grid=[]
        self.xcell=[]
        self.ycell=[]
        for lx,ly,ux,uy in zip(pyfile['grid'].data['lower_x'],
                               pyfile['grid'].data['lower_y'],
                               pyfile['grid'].data['upper_x'],
                               pyfile['grid'].data['upper_y']):
            self.grid.append((lx,ly,ux,uy))

        if self.skip61:self.ccd-=1

        # build the bounding box of ccds
        for i in range(0,self.ccd*len(self.grid)):
            iccd=i/len(self.grid)+1
            icell=i%len(self.grid)
            xl=self.grid[icell][0]
            yl=self.grid[icell][1]
            xu=self.grid[icell][2]
            yu=self.grid[icell][3]
            

            if self.skip61:
                if iccd==61:iccd=62
            #convert to focal plane coordinates
            # for center of grid
            xc=(xu-xl)/2.+xl
            yc=(yu-yl)/2.+yl
            bx,by=toFocal(iccd,xc,yc)
            self.xcell.append(bx)
            self.ycell.append(by)

            # for grid boudaries
            bcx,bcy=toFocal(iccd,xl,yl)
            bcx2,bcy2=toFocal(iccd,xu,yu)

            corner=[]
            corner.append(bcx)
            corner.append(bcy)
            self.bounds.append(Rectangle(corner,bcx2-bcx,bcy2-bcy))

        self.vec=pyfile['vec'].data
        self.data=pyfile['data_mr'].data
        self.coeff=pyfile['coeff'].data
        self.singular=pyfile['singular'].data

    def getPC(self,var,pc):
        p=copy.deepcopy(self.grid)
        if var>=self.nvar:
            print "too high"
            return False
        #p.set_array(self.vec[var::self.nvar][pc])
        p.set_array(self.vec[var*self.nx*self.ny*self.ccd:
                           (var+1)*self.nx*self.ny*self.ccd][pc])
        
        return p

    def getPC2(self,ax,var1,var2,pc,scale=-np.sqrt(2)):
        #p=copy.deepcopy(self.grid)
        if var1>=self.nvar or var2>=self.nvar:
            print "too high"
            return False
        
        p1a=scale*self.vec[var1::self.nvar,pc]
        p2a=scale*self.vec[var2::self.nvar,pc]
        
        v1=[]
        v2=[]
        for p1,p2 in zip(p1a,p2a):
            theta=np.math.atan2(p2,p1)/2
            e=np.sqrt(p1**2+p2**2)
            ce1=e*np.cos(theta)
            ce2=e*np.sin(theta)
            v1.append(ce1)
            v2.append(ce2)
        
        return  ax.quiver(self.xcell,self.ycell,v1,v2,
                          angles='uv',scale=1./500,
                          units='xy',pivot='middle',width=1,headwidth=0.,
                          headlength=0., headaxislength=0.,color='blue')
            
            
fig=pylab.figure()
ax=fig.add_subplot(111)
r=PCAResult('/data2/home/rarmst/work/psf_testing/sv/pca/sjj.fits')
#pc0=r.getPC(0,0)
pe=r.getPC2(ax,0,1,15)



#
#ax.add_collection(pc0)
#pylab.xlabel(' Focal Plane X (mm)')
#pylab.ylabel(' Focal Plane Y (mm)')
#cb=pylab.colorbar(pc0)

ax.add_collection(r.outline)
pylab.axis([-250,250,-250,250])
fig.show()


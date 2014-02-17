from matplotlib.patches import Circle, Wedge, Polygon,Rectangle
from matplotlib.collections import PatchCollection
import matplotlib
import pylab
import copy
import numpy as np
from scipy import interpolate
import argparse
from scipy.interpolate import Rbf

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

# helper functions for fitting results

def definePXY(order,x,xmin,xmax):
    tmp=np.zeros(order+1)
    newx=(2.*x-xmin-xmax)/(xmax-xmin)
    tmp[0]=1
    if order>0: tmp[1] = newx
    for i in range(2,order+1):
        tmp[i] = ((2.*i-1.)*newx*tmp[i-1] - (i-1.)*tmp[i-2])/i
    return tmp

def setPRow(order,x,y, xmin,xmax,ymin,ymax):
    tmp=np.zeros((order+1)*(order+2)/2)
    px=definePXY(order,x,xmin,xmax)
    py=definePXY(order,y,ymin,ymax)
    pq=0
    for n in range(0,order+1):
        p=n
        q=n-p
        while True:
            if not q<=n:break
            tmp[pq]=px[p]*py[q]
            pq+=1
            p-=1
            q+=1
    return tmp


class PCAResult:
    def __init__(self,filename,rpc=-1):
        self.filename=filename
        self.exps=[]
        self.bounds=[]
        self.skip61=True
        self.patches=[]
        for i,corner in corners.items():
            self.patches.append(Rectangle(corner,xsize,ysize))
        self.outline = PatchCollection(self.patches,facecolor='None')
        self.ascii=False
        self.rpc=rpc
        self.readFits()
        #self.readAscii()
        self.gridr = PatchCollection(self.bounds, cmap=matplotlib.cm.jet,
                                    alpha=1,edgecolors=None,linewidths=0)


    def readFits(self):
        self.pyfile=pyfits.open(self.filename)

        try:
            self.exps=self.pyfile['exps'].data['exposure']
            self.nexp=len(self.exps)
            self.ccd=self.pyfile['exps'].header['ccd']
            self.nvar=self.pyfile['exps'].header['nvar']
            self.nx=self.pyfile['exps'].header['nx']
            self.ny=self.pyfile['exps'].header['ny']
            self.xmax=self.pyfile['exps'].header['xmax']
            self.ymax=self.pyfile['exps'].header['ymax']
            self.type=self.pyfile['exps'].header['type']
            self.order=self.pyfile['exps'].header['order']
            self.skip61=1#self.pyfile['exps'].header['skip61']

        except:
            print "Could not read all header information"
            return

        try:
            self.pyfile['grid']
        except:
            print "Could not read grid hdu"
            return


        self.grid=[]
        self.xcell=[]
        self.ycell=[]
        for lx,ly,ux,uy in zip(self.pyfile['grid'].data['lower_x'],
                               self.pyfile['grid'].data['lower_y'],
                               self.pyfile['grid'].data['upper_x'],
                               self.pyfile['grid'].data['upper_y']):
            self.grid.append((lx,ly,ux,uy))
        self.ncells=len(self.grid)
        if self.skip61 and self.ccd>61:self.ccd-=1
        self.totvar=self.nx*self.ny*self.ccd*self.nvar

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

        try:
            self.data_mr=self.pyfile['data_mr'].data.transpose()
            self.mean=self.pyfile['mean'].data.transpose()
        except:
            self.data_mr=0
            self.mean=0

        try:
            self.coeff=copy.deepcopy(self.pyfile['coeff'].data.transpose())
            self.singular=self.pyfile['singular'].data
        except:
            print "Could not find coeff,singular in file"
            return


        try:
            self.missing=self.pyfile['missing'].data.transpose()
        except:
            1
                        
        

        if not self.type=='fit':
            try:
                self.vec=self.pyfile['vec'].data.transpose()
                self.data=self.pyfile['data'].data.transpose()
                self.npc=self.vec.shape[0]
                if self.rpc<0 or self.rpc>self.npc:
                    self.dataR=self.pyfile['dataR'].data.transpose()
                else:
                    
                    self.S=np.zeros((self.rpc,self.rpc))
                    for i in range(0,self.rpc):
                        self.S[i,i]=self.singular[0][i]
                   
                    self.dataR=np.dot(self.coeff[0:,0:self.rpc],self.S)
                    self.dataR=np.dot(self.dataR,self.vec[0:self.rpc,0:])
                    for i in range(0,len(self.mean)):
                        self.dataR[:,i]+=self.mean[i]
                   

            except:
                print "Could not find vec or data in file"
                return
        else:

            # This currently takes a long time.  Need to figure out
            # how to speed things up
            try:
                self.cvec=self.pyfile['vec'].data.transpose()
                self.cdata=self.pyfile['data'].data.transpose()

            except:
                print "Could not find vec or data in file"
                return

            self.npc=self.cvec.shape[0]
            if self.rpc<0 : self.rpc=self.npc

            self.S=np.zeros((self.rpc,self.rpc))
            for i in range(0,self.rpc):
                self.S[i,i]=self.singular[0][i]

                
            self.dataR=np.dot(self.coeff[0:,0:self.rpc],self.S)
            self.dataR=np.dot(self.dataR,self.cvec[0:self.rpc,0:])

            
            self.fitorder=(self.order+1)*(self.order+2)/2
            self.fitvar=self.fitorder*self.nvar
            (self.npc,nvar)=self.cvec.shape

            self.data=np.zeros((self.nexp,self.totvar))
            for exp in range(0,self.nexp):
                for i in range(0,self.ccd*len(self.grid)):
                    iccd=i/self.ncells+1
                    icell=i%self.ncells
                    xl=self.grid[icell][0]
                    yl=self.grid[icell][1]
                    xu=self.grid[icell][2]
                    yu=self.grid[icell][3]

                    
                    a=setPRow(self.order,(xu-xl)/2.+xl,(yu-yl)/2.+yl,
                              xl,xu,yl,yu)
                    # grab the fit information
                    start=i*self.fitorder*self.nvar
                    end=start+self.fitorder*self.nvar
                    
                    xmat=self.cdata[exp,start:end]
                    xmat=xmat.reshape(self.fitorder,self.nvar)
                    b=a.dot(xmat)

                    # write out the values at the center of the grid
                    start=i*self.nvar
                    end=start+self.nvar
                    self.data[exp,start:end]=b

                    
            self.vec=np.zeros((self.npc,self.totvar))
            for pc in range(0,self.npc):
                for i in range(0,self.ccd*len(self.grid)):
                    iccd=i/self.ncells+1
                    icell=i%self.ncells
                    xl=self.grid[icell][0]
                    yl=self.grid[icell][1]
                    xu=self.grid[icell][2]
                    yu=self.grid[icell][3]
                    
                    a=setPRow(self.order,(xu-xl)/2.+xl,(yu-yl)/2.+yl,
                              xl,xu,yl,yu)
                    # grab the fit information
                    start=i*self.fitorder*self.nvar
                    end=start+self.fitorder*self.nvar
                    xmat=self.cvec[pc,start:end]
                    
                    xmat=xmat.reshape(self.fitorder,self.nvar)
                    
                    b=a.dot(xmat)

                    # write out the values at the center of the grid
                    start=i*self.nvar
                    end=start+self.nvar
                    self.vec[pc,start:end]=b

    def readAscii(self):
        self.nx=2
        self.ny=4
        self.nvar=2
        self.ccd=10
        self.nvar=self.nx*self.ny*self.nvar*self.ccd
        self.grid=[]
        self.xcell=[]
        self.ycell=[]
        grid_file=open('/data2/home/rarmst/work/psf_testing/sv/pca/gridXY')
        for i,line in enumerate(grid_file):

            coords=line.split()
            x=(float(coords[2])-float(coords[0]))/2.+float(coords[0])
            y=(float(coords[3])-float(coords[1]))/2.+float(coords[1])
            self.bounds.append([float(coords[0]),float(coords[2]),
                                float(coords[1]),float(coords[3])])
            iccd=int(i/(self.nx*self.ny))+1
            cell=i%(self.nx*self.ny)
            cellid=i
            
            if self.skip61:
                if iccd==61:iccd=62
            bx,by=toFocal(iccd,x,y)
            self.xcell.append(bx)
            self.ycell.append(by)

        self.vec=np.genfromtxt(self.filename+'_vec')
        self.coeff=np.genfromtxt(self.filename+'_coeff')
        self.ascii=True


    def getShapelet(self,exp,ccd,x,y,var=-1,interp=None):

        
        if not self.type=='fit' and not interp:
            binx=int(x/(self.xmax/self.nx))
            biny=int(y/(self.ymax/self.ny))
            bin=binx*self.ny+biny
            icell=(ccd-1)*self.nx*self.ny*self.nvar+bin*self.nvar

            return self.dataR[exp,icell:icell+self.nvar]
            #return [a+b for a,b in zip(self.dataR[exp,icell:icell+self.nvar],self.mean[icell:icell+self.nvar])]
        elif interp=='interp' and var>=0:
            None
            #data=self.dataR[exp,var::self.nvar]
            #rbf=Rbf(self.xcell,self.ycell,data,epsilon=2)
            
           
            
            #return rbf(x,y)
        else:
            binx=int(x/(self.xmax/self.nx))
            biny=int(y/(self.ymax/self.ny))
            # this is the x,y cell 
            bin=binx*self.ny+biny

            xl=self.grid[bin][0]
            yl=self.grid[bin][1]
            xu=self.grid[bin][2]
            yu=self.grid[bin][3]
            
            a=setPRow(self.order,x,y,xl,xu,yl,yu)

            # grab the fit information
            
            start=(ccd-1)*self.nx*self.ny*self.fitvar+bin*self.fitvar
            end=start+self.fitvar
            xmat=self.dataR[exp,start:end]
                    
            xmat=xmat.reshape(self.fitorder,self.nvar)
            
            b=a.dot(xmat)
            return b

    #expect two lists of two variables
    def getQuiver(self,ax,x,y,data1,data2,scale=-np.sqrt(2),qscale=500):

        v1=[]
        v2=[]
        for p1,p2 in data1:
            theta=np.math.atan2(p2,p1)/2
            e=np.sqrt(p1**2+p2**2)
#            print e,p1,p2
            ce1=e*np.cos(theta)
            ce2=e*np.sin(theta)
            v1.append(ce1)
            v2.append(ce2)
            
        q1=ax.quiver(x,y,v1,v2,
                     angles='uv',scale=1./qscale,
                     units='dots',pivot='middle',width=1,headwidth=0.,
                     headlength=0., headaxislength=0.,color='blue')

        v3=[]
        v4=[]
        for p1,p2 in data2:
            
            theta=np.math.atan2(p2,p1)/2
            e=np.sqrt(p1**2+p2**2)
            #print e,p1,p2
            ce1=e*np.cos(theta)
            ce2=e*np.sin(theta)
            v3.append(ce1)
            v4.append(ce2)

        q2=ax.quiver(x,y,v3,v4,
                     angles='uv',scale=1./qscale,
                     units='dots',pivot='middle',width=1,headwidth=0.,
                     headlength=0., headaxislength=0.,color='red')

        return (q1,q2)
            

    def getCumVar(self):
        var=[]
        sum=0
        for i in self.singular[0]:
            sum+=i**2
            var.append(sum)

        var/=sum
        return var
    
    def getPC(self,var,pc):
        p=copy.deepcopy(self.gridr)
        if var>=self.nvar:
            print "too high"
            return False
        if not self.ascii:
            p.set_array(self.vec[pc,var::self.nvar])
        else:
            p.set_array(self.vec[pc,var*self.nx*self.ny*self.ccd:
                                 (var+1)*self.nx*self.ny*self.ccd])
        
        return p

    def getExp(self,var,exp):
        p=copy.deepcopy(self.gridr)
        if var>=self.nvar:
            print "too high"
            return False
        
        p.set_array(self.data[exp,var::self.nvar])
        
        return p

    def getExp(self,var,exp):
        p=copy.deepcopy(self.gridr)
        if var>=self.nvar:
            print "too high"
            return False
        
        p.set_array(self.data[exp,var::self.nvar])
        
        return p


    def getExpR(self,var,exp):
        p=copy.deepcopy(self.gridr)
        if var>=self.nvar:
            print "too high"
            return False
        
        p.set_array(self.dataR[exp,var::self.nvar])
        
        return p

    def getExp2(self,var1,var2,exp,scale=-np.sqrt(2),qscale=500):

        if var1>=self.nvar or var2>=self.nvar:
            print "too high"
            return False
        
        p1a=scale*self.data[exp,var1::self.nvar]
        p2a=scale*self.data[exp,var2::self.nvar]
        missing=self.missing[exp]
        v1=[]
        v2=[]
        x1=[]
        y1=[]
        v2=[]
        v1m=[]
        v2m=[]
        x1m=[]
        y1m=[]
        for p1,p2,miss,x,y in zip(p1a,p2a,missing,self.xcell,self.ycell):

            theta=np.math.atan2(p2,p1)/2
            e=np.sqrt(p1**2+p2**2)
            #print e,p1,p2
            ce1=e*np.cos(theta)
            ce2=e*np.sin(theta)
            if not miss:
                v1.append(ce1)
                v2.append(ce2)
                x1.append(x)
                y1.append(y)
            else:
                v1m.append(ce1)
                v2m.append(ce2)
                x1m.append(x)
                y1m.append(y)
                
        
        q1=ax.quiver(x1,y1,v1,v2,
                     angles='uv',scale=1./qscale,
                     units='xy',pivot='middle',width=1,headwidth=0.,
                     headlength=0., headaxislength=0.,color='blue')
        q2=ax.quiver(x1m,y1m,v1m,v2m,
                     angles='uv',scale=1./qscale,
                     units='xy',pivot='middle',width=1,headwidth=0.,
                     headlength=0., headaxislength=0.,color='red')
        return (q1,q1)
        

    def getExp2R(self,var1,var2,exp,scale=-np.sqrt(2),qscale=500):

        if var1>=self.nvar or var2>=self.nvar:
            print "too high"
            return False
        
        p1a=scale*self.dataR[exp,var1::self.nvar]
        p2a=scale*self.dataR[exp,var2::self.nvar]
        missing=self.missing[exp]
        v1=[]
        v2=[]
        x1=[]
        y1=[]
        v2=[]
        v1m=[]
        v2m=[]
        x1m=[]
        y1m=[]
        for p1,p2,miss,x,y in zip(p1a,p2a,missing,self.xcell,self.ycell):

            theta=np.math.atan2(p2,p1)/2
            e=np.sqrt(p1**2+p2**2)
            #print e,p1,p2
            ce1=e*np.cos(theta)
            ce2=e*np.sin(theta)
            if not miss:
                v1.append(ce1)
                v2.append(ce2)
                x1.append(x)
                y1.append(y)
            else:
                v1m.append(ce1)
                v2m.append(ce2)
                x1m.append(x)
                y1m.append(y)
                
        
        q1=ax.quiver(x1,y1,v1,v2,
                     angles='uv',scale=1./qscale,
                     units='xy',pivot='middle',width=1,headwidth=0.,
                     headlength=0., headaxislength=0.,color='blue')
        q2=ax.quiver(x1m,y1m,v1m,v2m,
                     angles='uv',scale=1./qscale,
                     units='xy',pivot='middle',width=1,headwidth=0.,
                     headlength=0., headaxislength=0.,color='red')
        return (q1,q1)
    

    def getExp2Res(self,var1,var2,exp,scale=-np.sqrt(2),qscale=5000):

        if var1>=self.nvar or var2>=self.nvar:
            print "too high"
            return False
        
        p1ar=scale*self.dataR[exp,var1::self.nvar]
        p2ar=scale*self.dataR[exp,var2::self.nvar]
        p1a=scale*self.data[exp,var1::self.nvar]
        p2a=scale*self.data[exp,var2::self.nvar]
        v1=[]
        v2=[]
        for p1,p2,p1r,p2r in zip(p1a,p2a,p1ar,p2ar):

            r1=p1r-p1
            r2=p2r-p2
            
            theta=np.math.atan2(r2,r1)/2
            e=np.sqrt(r1**2+r2**2)
            ce1=e*np.cos(theta)
            ce2=e*np.sin(theta)
            v1.append(ce1)
            v2.append(ce1)
        
        return  ax.quiver(self.xcell,self.ycell,v1,v2,
                          angles='uv',scale=1./qscale,
                          units='xy',pivot='middle',width=1,headwidth=0.,
                          headlength=0., headaxislength=0.,color='blue')


    def get2Res(self,var1,var2):

        v1=[]
        v2=[]
        for exp in range(0,self.nexp):
            p1ar=self.dataR[exp,var1::self.nvar]
            p2ar=self.dataR[exp,var2::self.nvar]
            p1a=self.data[exp,var1::self.nvar]
            p2a=self.data[exp,var2::self.nvar]
            p1r=np.mean(p1ar-p1a)
            p2r=np.mean(p2ar-p2a)
            v1.extend(p1ar-p1a)
            v2.extend(p2ar-p2a)
                      
        return (v1,v2)

    def getMean2Res(self,var1,var2,scale=-np.sqrt(2),qscale=500):

        nvar=self.nx*self.ny*self.ccd
        
        pp1=np.zeros(nvar)
        pp2=np.zeros(nvar)

        for exp in range(0,self.nexp):
            p1ar=self.dataR[exp,var1::self.nvar]
            p2ar=self.dataR[exp,var2::self.nvar]
            p1a=self.data[exp,var1::self.nvar]
            p2a=self.data[exp,var2::self.nvar]
            p1r=np.mean(p1ar-p1a)
            p2r=np.mean(p2ar-p2a)
            pp1+=(p1ar-p1a)
            pp2+=(p1ar-p1a)
        pp1/=self.nexp
        pp2/=self.nexp

        v1=[]
        v2=[]

        for p1,p2 in zip(pp1,pp2):
            theta=np.math.atan2(p2,p1)/2
            e=np.sqrt(p1**2+p2**2)
            ce1=e*np.cos(theta)
            ce2=e*np.sin(theta)
            v1.append(ce1)
            v2.append(ce2)

        
        return  ax.quiver(self.xcell,self.ycell,v1,v2,
                          angles='uv',scale=1./qscale,
                          units='xy',pivot='middle',width=1,headwidth=0.,
                          headlength=0., headaxislength=0.,color='blue')
                      

    def getPC2(self,ax,var1,var2,pc,scale=-np.sqrt(2),qscale=500):

        if var1>=self.nvar or var2>=self.nvar:
            print "too high",self.nvar,var1,var2
            return False
        
        if self.ascii:
            p1a=scale*self.vec[pc,var1*self.nx*self.ny*self.ccd:
                               (var1+1)*self.nx*self.ny*self.ccd]
            p2a=scale*self.vec[pc,var2*self.nx*self.ny*self.ccd:
                               (var2+1)*self.nx*self.ny*self.ccd]
        else:
            p1a=scale*self.vec[pc,var1::self.nvar]
            p2a=scale*self.vec[pc,var2::self.nvar]
        
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
                          angles='uv',scale=1./qscale,
                          units='xy',pivot='middle',width=1,headwidth=0.,
                          headlength=0., headaxislength=0.,color='blue')


    def setNewBasis(self,other,tol=1e-6):
        if other.nvar!=self.nvar :
            return False

        #convert the inverse singular values to a matrix
        smat=np.zeros((other.npc,other.npc))
        nbad=0
        for i in range(other.npc):
            if other.singular[0,i]<tol:
                smat[i,i]=0.
                nbad+=1
            else:
                smat[i,i]=1./other.singular[0,i]
        tmp=np.dot(other.vec.transpose(),smat)

        self.coeff=np.dot(self.data_mr,tmp)

    # ignore the fact that if run without mean remove we don't have
    # the mean vector
    def getMean(self,var,rm_pc):
        p=copy.deepcopy(self.gridr)
        if var>=self.nvar:
            print "too high"
            return False
        
        p.set_array(self.mean[var::self.nvar].flatten())
        
        return p

    def getMean2(self,ax,var1,var2,rm=[],qscale=500,scale=-np.sqrt(2)):

        if var1>=self.nvar or var2>=self.nvar:
            print "too high"
            return False
        mean=copy.deepcopy(self.mean).flatten()
        if len(rm)>0:
            #smat=np.zeros((self.npc,self.npc))
            #for i in range(self.npc):
            #    if other.singular[0,i]<tol:
            #        smat[i,i]=0.
            #    else:
            #        smat[i,i]=1./other.singular[0,i]
            #tmp=np.dot(other.vec.transpose(),smat)
            #print mean
            mean=np.dot(mean,self.vec.transpose())
            for i in rm:
                mean[i]=0
            #print mean
            mean=np.dot(mean,self.vec)
           
        
        p1a=scale*mean[var1::self.nvar].flatten()
        p2a=scale*mean[var2::self.nvar].flatten()
        
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
                          angles='uv',scale=1./qscale,
                          units='xy',pivot='middle',width=1,headwidth=0.,
                          headlength=0., headaxislength=0.,color='blue')
        
if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description='Run single file')
    parser.add_argument('--file',
                        help='number of jobs to run concurrently')
    parser.add_argument('--print_pc',default=0,type=int,
                        help='print PCs')
    parser.add_argument('--print_pc2',default=0,type=int,
                        help='print PCs')
    parser.add_argument('--print_exp',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_expr',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_exp2',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_exp2c',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_exp2r',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_exp2res',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_exp2b',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_meanexp2res',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_res',default=0,type=int,
                        help='print exposures')
    parser.add_argument('--print_cumvar',default=0,type=int,
                        help='print cumulative variance')
    parser.add_argument('--npc',default=20,type=int,
                        help='number of PCs')
    parser.add_argument('--rpc',default=-1,type=int,
                        help='number of PCs used reconstruction')
    

    parser.add_argument('--var1',default=0,type=int,
                        help='1st variable to plot')
    parser.add_argument('--var2',default=1,type=int,
                        help='2nd variable to plot')

    parser.add_argument('--scale',default=500,type=float,
                        help='scale of tick mark')
    parser.add_argument('--suffix',default='',
                        help='suffix at end of plot names')
 
    args = parser.parse_args()
     

    fig=pylab.figure()
    ax=fig.add_subplot(111)
    r=PCAResult(args.file,args.rpc)
    #r2=PCAResult('/data2/home/rarmst/work/psf_testing/sv/pca/test_new3c.fits')
#pc=r.getMean(1)
#ax.add_collection(pc)
#cb=pylab.colorbar(pc)
#pylab.axis([-250,250,-250,250])
#fig.show()

#a=r.getCumVar()
#print a

#pe=r.getMean2(ax,0,1,[0,1])
#    pylab.axis([-250,250,-250,250])
#    ax.add_collection(r.outline)

     
#r2.setNewBasis(r)

    if args.print_res:
        #n, bins, patches = ax.hist(x, 50, normed=1, histtype='stepfilled')
        ax.clear()
        (r1,r2)=r.get2Res(args.var1,args.var2)
        ax.hist(r1,100)
        ax.hist(r2,100)
        #ax.plot(r1)
        #ax.plot(r2)
        print np.mean(r1),np.std(r1)
        print np.mean(r2),np.std(r2)
        fig.show()
    
    if args.print_exp2:
        for iexp in range(0,r.nexp):
            ax.clear()
            (q1,q2)=r.getExp2(args.var1,args.var2,iexp,-np.sqrt(2),args.scale)
            pylab.axis([-250,250,-250,250])
            ax.add_collection(r.outline)
            fig.savefig('%s_var%d_%d%s.png'%(r.exps[iexp],args.var1,args.var2,args.suffix))

    if args.print_exp2r:
        for iexp in range(0,r.nexp):
            ax.clear()
            exp1=r.getExp2R(args.var1,args.var2,iexp,-np.sqrt(2),args.scale)
            pylab.axis([-250,250,-250,250])
            ax.add_collection(r.outline)
            fig.savefig('r_%s_var%d_%d%s.png'%(r.exps[iexp],args.var1,args.var2,args.suffix))

    if args.print_meanexp2res:
        ax.clear()
        exp1=r.getMean2Res(args.var1,args.var2)
        pylab.axis([-250,250,-250,250])
        ax.add_collection(r.outline)
        fig.show()

    if args.print_exp2res:
        for iexp in range(0,r.nexp):
            ax.clear()
            exp1=r.getExp2Res(args.var1,args.var2,iexp)
            pylab.axis([-250,250,-250,250])
            ax.add_collection(r.outline)
            fig.savefig('res_%s_var%d_%d%s.png'%(r.exps[iexp],args.var1,args.var2,args.suffix))

   
    if args.print_exp2c:
        for iexp in range(0,r.nexp):
            ax.clear()
            exp1=r.getExp2(args.var1,args.var2,iexp)
            exp1r=r.getExp2R(args.var1,args.var2,iexp)
            exp1r[0].set_color('red')
            pylab.axis([-250,250,-250,250])
            ax.add_collection(r.outline)
            fig.savefig('c_%s_var%d_%d%s.png'%(r.exps[iexp],args.var1,args.var2,args.suffix))

    if args.print_exp:
        for iexp in range(0,r.nexp):
            print r.exps[iexp]
            ax.clear()
            fig.clear()
            ax=fig.add_subplot(111)
                        
            exp1=r.getExp(args.var1,iexp)
            ax.add_collection(exp1)
            pylab.axis([-250,250,-250,250])
            cb=pylab.colorbar(exp1)
            ax.add_collection(r.outline)
            pylab.xlabel(' Focal Plane X (mm)')
            pylab.ylabel(' Focal Plane Y (mm)')
            fig.savefig('%s_var%d%s.png'%(r.exps[iexp],args.var1,args.suffix))

    if args.print_expr:
        for iexp in range(0,r.nexp):
            print r.exps[iexp]
            ax.clear()
            fig.clear()
            ax=fig.add_subplot(111)
                        
            exp1=r.getExpR(args.var1,iexp)
            ax.add_collection(exp1)
            pylab.axis([-250,250,-250,250])
            cb=pylab.colorbar(exp1)
            ax.add_collection(r.outline)
            pylab.xlabel(' Focal Plane X (mm)')
            pylab.ylabel(' Focal Plane Y (mm)')
            fig.savefig('r_%s_var%d%s.png'%(r.exps[iexp],args.var1,args.suffix))

    if args.print_pc2:
        for ipc in range(0,args.npc):
            ax.clear()
            pe=r.getPC2(ax,args.var1,args.var2,ipc)
            pylab.axis([-250,250,-250,250])
            ax.add_collection(r.outline)
            pylab.xlabel(' Focal Plane X (mm)')
            pylab.ylabel(' Focal Plane Y (mm)')
            fig.savefig('pc%d_var%d_var%d%s.png'%(ipc,args.var1,args.var2,args.suffix))

    if args.print_pc:
        for ipc in range(0,args.npc):
            fig.clear()
            ax=fig.add_subplot(111)
            ax.clear()
            pc=r.getPC(args.var1,ipc)
            ax.add_collection(pc)
            
            cb=pylab.colorbar(pc)
            pylab.axis([-250,250,-250,250])
            pylab.xlabel(' Focal Plane X (mm)')
            pylab.ylabel(' Focal Plane Y (mm)')
            fig.savefig('pc%d_var%d%s.png'%(ipc,args.var1,args.suffix))

    if args.print_cumvar:

        fig.clear()
        ax=fig.add_subplot(111)
        var=r.getCumVar()[0:r.rpc]
        ax.plot(var,'o-')
        pylab.ylabel('Cumulative Variance')
        pylab.xlabel('PC')
        fig.savefig('var%s.png'%(args.suffix))

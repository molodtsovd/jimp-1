import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
#from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate

# get number of lines without loading the file into RAM
def getNumLines(fname):
    import subprocess
    return int(subprocess.check_output('wc -l {}'.format(fname), shell=True).split()[0])

def plotSD(fname):

  plt.figure()

  #ID,RGridE,RGridN,RElev,SGridE,SGridN,SElev,WaveID,TimeObs,TimeErr,TimeMod

  rx, sx, waveID, Tobs, Terr, Tcal = np.loadtxt(fname,delimiter=',', skiprows=1, 
                                          usecols=(1,4,7,8,9,10), unpack=True)

  Nt = len(rx)

  sx0 = sx[0];
  waveID0 = waveID[0];
  i0 = 0;
  for i in range(Nt):
    if i == Nt-1 or sx[i] != sx0 or (i!=0 and abs(rx[i]-rx[i-1]) > 0.15) or waveID[i] != waveID0 :
#      h1 = plt.errorbar(rx[i0:i], Tobs[i0:i], yerr=Terr[i0:i], fmt=' ', color="red", label="$t_{obs}$")
#      h1, = plt.plot(rx[i0:i], Tobs[i0:i], '-', color="red", label="$t_{obs}$")
      h1 = plt.fill_between(rx[i0:i], Tobs[i0:i]-Terr[i0:i], Tobs[i0:i]+Terr[i0:i], color="red", label="$t_{obs}$")
      h2, = plt.plot(rx[i0:i], Tcal[i0:i], '-', color="blue", label="$t_{cal}$")
	  
      i0 = i;
      sx0 = sx[i];
      waveID0 = waveID[i];

  plt.legend(handles=[h1,h2], loc='upper right outside', bbox_to_anchor=(1,1))
  plt.gca().set_xlabel("y [km]")
  plt.gca().set_ylabel("t [s]")
  tit = "$\delta t$=%g s" % (np.linalg.norm(Tobs - Tcal) / np.sqrt(Nt))
  plt.title(tit)
  plt.gca().set_xlim([min(rx), max(rx)])
  plt.show()

def plotV(v, h, x=0, y=0, z=0, v0=1., v1=9., dv=1., prop='$v_p$ [km/s]'):


  cm='jet'
    

  fig = plt.figure()
  yz = fig.add_subplot(221)
  zx = fig.add_subplot(222)#, sharey=yz)
  xy = fig.add_subplot(223)#, sharex=yz)

  xy.set_xlabel("y [km]")
  xy.set_ylabel("x [km]")
  ext = (0, h * np.shape(v)[1], 0, h * np.shape(v)[2])
  xy.imshow(v[z,:,:].T, extent=ext, vmin=v0, vmax=v1, cmap=cm)
  xy.set_aspect('auto')


  zx.set_xlabel("x [km]")
  zx.set_ylabel("z [km]")
  ext = (0, h * np.shape(v)[2], -h * np.shape(v)[0], 0)
  zx.imshow(v[:,y,:], extent=ext, vmin=v0, vmax=v1, cmap=cm)
  zx.set_aspect('auto')
  

  yz.set_xlabel("y [km]")
  yz.set_ylabel("z [km]")
  ext = (0, h * np.shape(v)[1], -h * np.shape(v)[0], 0)
  yzh = yz.imshow(v[:,:,x], extent=ext, vmin=v0, vmax=v1, cmap=cm)
  yz.set_aspect('auto')
  
  
  bns=np.linspace(v0, v1)
  tix=np.arange(v0, v1+dv, dv)
  fig.colorbar(yzh, ax=[xy,zx,yz], extend='both', boundaries=bns, ticks=tix, 
                    orientation='vertical', aspect=50, label=prop)

#  ax = fig.add_subplot(224, projection='3d')
#
#  nx = np.shape(v)[2]
#  ny = np.shape(v)[1]
#  nz = np.shape(v)[0]
#  
#  xx, yy = np.meshgrid(np.linspace(0, h*nx, nx), np.linspace(0, h*ny, ny))
#  ax.contourf(xx, yy, v[z,:,:], zdir='z', offset = -h*z, cmap=cm)
#  
#  xx, zz = np.meshgrid(np.linspace(0, h*nx, nx), np.linspace(-h*nz, 0, nz))
#  ax.contourf(xx, zz, v[:,y,:], zdir='y', offset = h*y, cmap=cm)
#  
#  yy, zz = np.meshgrid(np.linspace(0, h*ny, ny), np.linspace(-h*nz, 0, nz))
#  ax.contourf(yy, zz, v[:,:,x], zdir='z', offset = h*x, cmap=cm)
  
  return yz

def plotV2d(v, h, v0=1., v1=9., dv=1., prop='$v_p$ [km/s]'):


  cm='jet'
    

  fig = plt.figure()
  yz = fig.add_subplot(111)

  yz.set_xlabel("y [km]")
  yz.set_ylabel("z [km]")
  ext = (0, h * np.shape(v)[1], -h * np.shape(v)[0], 0)
  yzh = yz.imshow(v, extent=ext, vmin=v0, vmax=v1, cmap=cm)
  yz.set_aspect('auto')
  
  
  bns=np.linspace(v0, v1)
  tix=np.arange(v0, v1+dv, dv)
  fig.colorbar(yzh, extend='both', boundaries=bns, ticks=tix, 
                    orientation='vertical', aspect=50, label=prop)
  
#  sx,sz=np.genfromtxt('test/t.sd', usecols=(4,6), unpack=True, skip_header=1,delimiter=',')
#  yz.plot(sx,sz,'^',color='red')
      
  return yz

def plotH(H, Hfd, ax, style='-', width=1.5, color='black'):
    y = Hfd[0,:]
    z = np.zeros(np.shape(y))
    for i in range(len(H)):
        z -= Hfd[i+1,:]
        spl = interpolate.PchipInterpolator(y, z, extrapolate=True)
#        if i > 1:
#        ax.plot(y, spl(y), linestyle=style, linewidth=width, color=color)
        ax.plot([p.y for p in H[i]], spl([p.y for p in H[i]]), 
                '.', linewidth=width, color=color)
    ax.plot([p.y for p in H[4][9:27]], spl([p.y for p in H[4][9:27]]), 
            '.', linewidth=width, color='#00ff00')
  
    
    
    
def plotTemp(t, h, tit=''):

  cm='jet'
    
  fig = plt.figure()
  yz = fig.add_subplot(111)

  yz.set_xlabel("y [km]")
  yz.set_ylabel("z [km]")
  ext = (0, h * np.shape(t)[1], -h * np.shape(t)[0], 0)
  yzh = yz.imshow(t[:,:,0], extent=ext, cmap=cm)
  yz.set_aspect('auto')

  fig.colorbar(yzh, orientation='vertical', aspect=50, label='$T$ K')
  
  plt.title(tit)
   
  return yz    


def plotGD_2D(fname, fname2=None):

  plt.figure()
  
  # X,Y,Z,GOBS,ERR,GCAL

  y, g_obs, err, g_cal = np.loadtxt(fname,delimiter=',', skiprows=1, 
                                          usecols=(1,3,4,5), unpack=True)

  g_obs -= np.average(g_obs)
  g_cal -= np.average(g_cal)

  N = len(y)

  h1 = plt.errorbar(y, g_obs,fmt=' ', yerr=err, capsize=2, color="blue", label="$\Delta g_{obs}$")
#  h1 = plt.fill_between(y, g_obs-err, g_obs+err, color="red", label="$\Delta g_{obs}$")
  h2, = plt.plot(y, g_cal, color="red", label="$\Delta g_{cal}$")

  plt.legend(handles=[h1,h2], loc='upper right')
  plt.gca().set_xlabel("y [km]")
  plt.gca().set_ylabel("$\Delta g$ [mGal]")
  tit = "$\delta \Delta g$=%g mGal" % (np.linalg.norm(g_obs - g_cal) / np.sqrt(N))
  plt.title(tit)
  plt.gca().set_xlim([min(y), max(y)])
  plt.show()

  if fname2 != None:
    # compare with fname2
      g_obs2 = np.loadtxt(fname2,delimiter=',', skiprows=1, 
                                              usecols=3, unpack=True)
      h3, = plt.plot(y, g_obs2, color="green", label="$\Delta g_{obs} 2D$")
      h4, = plt.plot(y, g_obs2-g_obs, color="black", label="error")
      plt.legend(handles=[h1,h2,h3,h4], loc='upper right')


def plotElev(fname, bath, ax):
    
  y_bath, z_bath = np.loadtxt(bath,delimiter=',', skiprows=1, unpack=True)
  spl = interpolate.interp1d(y_bath, z_bath, fill_value="extrapolate")
  
  y, z_obs, z_cal = np.loadtxt(fname,delimiter=',', skiprows=1, unpack=True)
  ax.plot(y, spl(y) + z_cal - z_obs, color='red')


def plotCov(C, mode=None, tit=None, vmax=None):
    M = np.size(C,0)
    fig, ax = plt.subplots()
    if vmax == None:
        vmax=np.max(abs(C))
    im = plt.imshow(C, vmin=-vmax, vmax=vmax, cmap='seismic')
    fig.colorbar(im)

    m = [35,35,35,6,7,6,18]
    if mode == 'all_ox':
        m = [35,35,35,8,8,7,18]
    tix=[sum(m[0:i])+0.5*m[i] for i in range(len(m))]
    labels=['$H_0$','$H_1$','$H_2$','$C_0$','$C_1$','$C_2$', 'Srp']
    tixlen=0
    
    if mode=='h':
        m = [35,35,35]
        tix=[sum(m[0:i])+0.5*m[i] for i in range(len(m))]
        labels=['$H_0$','$H_1$','$H_2$']
    if mode=='c':
        m = [6,7,6,18]
        tix=[sum(m[0:i])+0.5*m[i] for i in range(len(m))]
        print(tix)
        labels = ['$C_0$','$C_1$','$C_2$', 'Srp']
    if mode=='ox':
        # oxides
        m = [8,8,7]
        tix = range(M)
        labels = ['H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
                  'H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
                  'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O']
        tixlen=3.5
    if mode=='ph':
        # phases
        m = [6,7,6]
        tix = range(M)
        labels = ['Opx','Fsp','Fsp','Cpx','Bio','q',
                  'Opx','Fsp','Fsp','Fsp','Cpx','Bio','q',
                  'Opx','Fsp','Fsp','Fsp','Sp','Cpx']  
        tixlen=3.5
    
    plt.yticks(tix, labels)
    plt.xticks(tix, labels, rotation=90)
    for i in range(len(m)-1):
        a = [-0.5,M-0.5], [sum(m[0:i+1])-.5,sum(m[0:i+1])-.5]
        color = 'gray'
        width = 1.5
        plt.plot(a[0], a[1], linewidth=width, color=color)
        plt.plot(a[1], a[0], linewidth=width, color=color)
    ax.tick_params(length=tixlen)
    if tit != None:
        plt.title(tit)



def plotHistZ(finPrefix, Nlayers, R0, R1, dx, y):
    Nx = (R1[0]-R0[0])/dx
    x = np.linspace(R0[0], R1[0], Nx)
    print(x)
        
    fig, axes = plt.subplots(Nlayers,1)
    for i in range(Nlayers):
        hist = np.genfromtxt('%s%d'%(finPrefix,i))  
        zBin = np.linspace(R0[1], R1[1], np.size(hist,0))
        spl = interpolate.interp2d([y], zBin, hist)
        z = zBin
        if i < Nlayers-1:
            axes[i].set_xticklabels([])
        else:
            axes[i].set_xlabel("y [km]")
        axes[i].set_ylabel("z [km]")
        ext = (R0[0], R1[0], R1[1], R0[1])
        im = axes[i].imshow(np.flipud(spl(x,z)), extent=ext)
        axes[i].set_aspect('auto')
        
    fig.colorbar(im, ax=axes, orientation='vertical')
    return axes

def plotRMS(fin, period, chains, targetRMS=1.):
    T = np.zeros(len(chains))
    for k in chains:
        T[k] = np.genfromtxt(fin+str(k), usecols=3, max_rows=1)
    Tmax = np.log(max(T))
    Tmin = np.log(min(T))
    cmap = plt.get_cmap('jet_r')
    fig = plt.figure()
    im = plt.imshow(T[:, np.newaxis], cmap=cmap, norm=mpl.colors.LogNorm())
    plt.clf()
    axes = fig.subplots(3,1)
    fig.colorbar(im, ax=axes, label='Temperature')
    handles = []
    for k in chains:
        rms = np.genfromtxt(fin+str(k))
        itn = np.arange(0, np.size(rms,0),period, dtype=int)
        for i in range(0,3):
            h, = axes[i].loglog(itn, rms[0::period,i],'-',
                     color=cmap((np.log(rms[0,3])-Tmin)/(Tmax-Tmin)), label=str(k)) #label=('T=%g'%T[k]))
            if i == 0:
                handles.append(h)

    for i in range(0,3):
        xlim = axes[i].get_xlim()
        axes[i].loglog([xlim[0],xlim[1]], [targetRMS,targetRMS],'--',color='black',linewidth=1.)
        axes[i].set_xlim(xlim)
        axes[i].set_ylabel('RMS misfit')
        axes[i].grid(True,which='both')
        if i < 2:
            axes[i].set_xticklabels([])

    axes[0].set_yticks([1,10])
    axes[0].set_title('Traveltimes')
    axes[1].set_title('Gravity')
    axes[2].set_title('Elevation')
    axes[2].set_xlabel("Sample")
#    fig.legend(handles=handles,loc='center right',fontsize='x-small',ncol=3)
    
    
def plotSwap(fin):
    t = np.genfromtxt(fin)
    plt.figure()
    for i in range(np.size(t,0)):
        if t[i,1] == t[i,2] :#or (t[i,1] != 1. and t[i,2] != 1.):
            t[i,1] = t[i,2] = np.NaN
    plt.loglog(t[:,0],t[:,1],'.',color='blue')
    plt.loglog(t[:,0],t[:,2],'.',color='red')

#def plotModels(fin, fin_true=None, mode='ph'):
#    period = 8
#    m = np.cumsum([0,35,35,35,6,7,6,18], dtype=int)
#    Clabels = ['Opx','Fsp','Fsp','Cpx','Bio','q',
#              'Opx','Fsp','Fsp','Fsp','Cpx','Bio','q',
#              'Opx','Fsp','Fsp','Fsp','Sp','Cpx']  
#    if mode=='ox':
#        m = np.cumsum([0,35,35,35,8,8,7,18], dtype=int)
#        Clabels = ['H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
#                  'H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
#                  'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O']
#
#    Nsamp = int(getNumLines(fin) / period)
#    
#    cmap = plt.get_cmap('jet')
#    Nmark = len(Line2D.filled_markers)
#    fig = plt.figure()
#    axes = fig.subplots(3,1)
#    handles = []
#    for i in range(Nsamp):       
#        H = np.genfromtxt(fin, skip_header=i*period, max_rows=3)
#        C = np.empty(0)
#        for k in range(3):
#            a = np.genfromtxt(fin, skip_header=i*period+3+k, max_rows=1)
#            C = np.concatenate([C,a])        
#        Srp = np.genfromtxt(fin, skip_header=i*period+6, max_rows=1)
#
#        for k in range(3):
#            xi = range(m[k],m[k+1])
##            xc = np.linspace(m[k],m[k+1]-1,300)
##            spl = interpolate.PchipInterpolator(xi, H[k,:])           
##            axes[0].plot(xc, spl(xc), '-',color=cmap(i/Nsamp))
#            axes[0].plot(xi, H[k,:], marker=Line2D.filled_markers[i%Nmark],
#                fillstyle='none',color=cmap(i/Nsamp))
#        axes[1].plot(C, marker=Line2D.filled_markers[i%Nmark],
#            fillstyle='none',color=cmap(i/Nsamp))
#        h, = axes[2].plot(Srp, marker=Line2D.filled_markers[i%Nmark],
#            fillstyle='none',color=cmap(i/Nsamp), label=str(i))
#        handles.append(h)
#
#    if fin_true != None:        
#        H = np.genfromtxt(fin_true, max_rows=3)
#        C = np.empty(0)
#        for k in range(3):
#            a = np.genfromtxt(fin_true, skip_header=3+k, max_rows=1)
#            C = np.concatenate([C,a])        
#        Srp = np.genfromtxt(fin_true, skip_header=6, max_rows=1)
#
#        for k in range(3):
#            xi = range(m[k],m[k+1])
##            xc = np.linspace(m[k],m[k+1]-1,300)
##            spl = interpolate.PchipInterpolator(xi, H[k,:])           
##            axes[0].plot(xc, spl(xc), '-',color=cmap(i/Nsamp))
#            axes[0].plot(xi, H[k,:], '--',color='black')
#        axes[1].plot(C,'--',color='black')
#        h, = axes[2].plot(Srp, '--',color='black', label='True')
#        handles.append(h)
#        
#    tix=[m[1]*0.5, m[2]-(m[1]-m[0])*0.5, m[3]-(m[2]-m[1])*0.5]
#    labels=['Layer 0','Layer 1','Layer 2']
#    axes[0].set_xticks(tix)
#    axes[0].set_xticklabels(labels)
#    axes[0].tick_params(length=0)
#    ylim = axes[0].get_ylim()
#    for i in range(1,3):
#        axes[0].plot([m[i]-0.5,m[i]-0.5], [ylim[0], ylim[1]],color='black',linewidth=0.5)
#    axes[0].set_ylim(ylim)
#    axes[0].set_ylabel('Layer thickness [km]')
#
#    axes[1].set_xticks(range(m[6]-m[3]))
#    axes[1].set_xticklabels(Clabels)
#    ylim = axes[1].get_ylim()
#    for i in range(1,3):
#        axes[1].plot([m[3+i]-m[3]-0.5,m[3+i]-m[3]-0.5], [ylim[0], ylim[1]],color='black',linewidth=0.5)
#    axes[1].set_ylim(ylim)
#    axes[1].set_ylabel('vol%')
#    
#    axes[2].set_xticklabels([])
#    axes[2].tick_params(length=0)
#    axes[2].set_xlabel('Mantle serpentinization')
#    axes[2].set_ylabel('wt%')
#
#    fig.legend(handles=handles,loc='center right',fontsize='x-small')
#    
#def histAll(fin, Nbin, fin_true=None, mode='ph', 
#            h01=[0.,30.], c01=[0.,100.], s01=[0.,100.]):
#    period = 8
#    m = np.cumsum([0,35,35,35,6,7,6,18], dtype=int)
#    Clabels = ['Opx','Fsp','Fsp','Cpx','Bio','q',
#              'Opx','Fsp','Fsp','Fsp','Cpx','Bio','q',
#              'Opx','Fsp','Fsp','Fsp','Sp','Cpx']  
#    if mode=='ox':
#        m = np.cumsum([0,35,35,35,8,8,7,18], dtype=int)
#        Clabels = ['H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
#                  'H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
#                  'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O']
#    
#    dh = (h01[1] - h01[0]) / Nbin
#    dc = (c01[1] - c01[0]) / Nbin
#    ds = (s01[1] - s01[0]) / Nbin
#    v0 = np.concatenate([h01[0]*np.ones(3), c01[0]*np.ones(3), s01[0]*np.ones(1)])  
#    dv = np.concatenate([dh*np.ones(3), dc*np.ones(3), ds*np.ones(1)])  
#
#    hist = np.zeros((Nbin, m[-1]))      
#    
#    Nsamp = getNumLines(fin) / period
#    print('Number of samples:', Nsamp)
#    Nsamp = int(Nsamp)
#    
#    f = open(fin, 'r')
#    buf = f.readlines()
#    f.close()
#
#    for i in range(Nsamp):
#        for j in range(len(m)-1):
#            v = np.fromstring(buf[i*period+j][:], sep=' ')
#            for k in range(len(v)):
#                n = int((v[k] - v0[j]) / dv[j])
#                if n < np.size(hist,0):
#                    hist[n, m[j]+k] += 1
#    hist /= Nsamp
#    hist = np.flipud(hist)
#            
#    fig = plt.figure()
#    axes = fig.subplots(3,1)
#    vmax = np.max(hist[:,m[6]:m[-1]])
#    axes[0].imshow(hist[:,m[0]:m[3]], vmax=vmax,
#            extent=[0,m[3],h01[0],h01[1]], aspect='auto')
#    axes[1].imshow(hist[:,m[3]:m[6]], vmax=vmax,
#            extent=[0,m[6]-m[3],c01[0],c01[1]], aspect='auto')
#    im=axes[2].imshow(hist[:,m[6]:m[-1]], vmax=vmax,
#            extent=[0,m[-1]-m[6],s01[0],s01[1]], aspect='auto')
#
#    if fin_true != None:        
#        H = np.genfromtxt(fin_true, max_rows=3)
#        C = np.empty(0)
#        for k in range(3):
#            a = np.genfromtxt(fin_true, skip_header=3+k, max_rows=1)
#            C = np.concatenate([C,a])        
#        Srp = np.genfromtxt(fin_true, skip_header=6, max_rows=1)
#
#        for k in range(3):
#            xi = np.arange(m[k],m[k+1])
#            axes[0].plot(xi+0.5, H[k,:], '--',color='white')
#        axes[1].plot(np.arange(0.5,m[6]-m[3]), C,'--',color='white')
#        axes[2].plot(np.arange(0.5,m[-1]-m[6]), Srp, '--',color='white')
#        
#    tix=[m[1]*0.5, m[2]-(m[1]-m[0])*0.5, m[3]-(m[2]-m[1])*0.5]
#    labels=['Layer 0','Layer 1','Layer 2']
#    axes[0].set_xticks(tix)
#    axes[0].set_xticklabels(labels)
#    axes[0].tick_params(length=0)
#    ylim = axes[0].get_ylim()
#    for i in range(1,3):
#        axes[0].plot([m[i],m[i]], [ylim[0], ylim[1]],color='white',linewidth=0.5)
#    axes[0].set_ylim(ylim)
#    axes[0].set_ylabel('Layer thickness [km]')
#
#    axes[1].set_xticks(np.arange(0.5,m[6]-m[3]))
#    axes[1].set_xticklabels(Clabels)
#    ylim = axes[1].get_ylim()
#    for i in range(1,3):
#        axes[1].plot([m[3+i]-m[3],m[3+i]-m[3]], [ylim[0], ylim[1]],color='white',linewidth=0.5)
#    axes[1].set_ylim(ylim)
#    axes[1].set_ylabel('vol%')
#    
#    axes[2].set_xticklabels([])
#    axes[2].tick_params(length=0)
#    axes[2].set_xlabel('Mantle serpentinization')
#    axes[2].set_ylabel('wt%')
#    fig.colorbar(im, ax=axes)

def histAll(fin, Nbin, fin_true=None, mode='ph', 
            h01=[0.,30.], c01=[0.,100.], s01=[0.,100.]):
    period = 8
    m = np.cumsum([0,35,35,35,6,7,6,18], dtype=int)
    Clabels = ['Opx','Fsp','Fsp','Cpx','Bio','q',
              'Opx','Fsp','Fsp','Fsp','Cpx','Bio','q',
              'Opx','Fsp','Fsp','Fsp','Sp','Cpx']  
    if mode=='ox':
        m = np.cumsum([0,35,35,35,8,8,7,18], dtype=int)
        Clabels = ['H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
                  'H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
                  'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O']
    
    dh = (h01[1] - h01[0]) / Nbin
    dc = (c01[1] - c01[0]) / Nbin
    ds = (s01[1] - s01[0]) / Nbin
    v0 = np.concatenate([h01[0]*np.ones(3), c01[0]*np.ones(3), s01[0]*np.ones(1)])  
    dv = np.concatenate([dh*np.ones(3), dc*np.ones(3), ds*np.ones(1)])  

    hist = np.zeros((Nbin, m[-1]))      
    
    Nsamp = getNumLines(fin) / period
    print('Number of samples:', Nsamp)
    Nsamp = int(Nsamp)
    
    f = open(fin, 'r')
    buf = f.readlines()
    f.close()

    for i in range(Nsamp):
        for j in range(len(m)-1):
            v = np.fromstring(buf[i*period+j][:], sep=' ')
            for k in range(len(v)):
                n = int((v[k] - v0[j]) / dv[j])
                if n < np.size(hist,0):
                    hist[n, m[j]+k] += 1
    hist /= Nsamp
    hist = np.flipud(hist)
            
    fig = plt.figure()
    axes = fig.subplots(3,3)
    
    gs = axes[2, 0].get_subplotspec().get_gridspec()
    # remove the underlying axes
    for ax in axes[2, :]:
        ax.remove() 
    ax2 = fig.add_subplot(gs[2, :])
    
    vmax = np.max(hist[:,m[6]:m[-1]])
    for i in range(3):
        axes[0,i].imshow(hist[:,m[i]:m[i+1]], vmax=vmax,
                extent=[0,m[i+1]-m[i],h01[0],h01[1]], aspect='auto')
        axes[1,i].imshow(hist[:,m[i+3]:m[i+4]], vmax=vmax,
                extent=[0,m[i+4]-m[i+3],c01[0],c01[1]], aspect='auto')
    im=ax2.imshow(hist[:,m[6]:m[-1]], vmax=vmax,
            extent=[0,m[-1]-m[6],s01[0],s01[1]], aspect='auto')

    if fin_true != None:        
        H = np.genfromtxt(fin_true, max_rows=3)
        C = []
        for k in range(3):
            a = np.genfromtxt(fin_true, skip_header=3+k, max_rows=1)
            C.append(a)        
        Srp = np.genfromtxt(fin_true, skip_header=6, max_rows=1)

        for k in range(3):
            axes[0,k].plot(np.arange(m[k+1]-m[k])+0.5, H[k,:], '--',color='red')
            axes[1,k].plot(np.arange(m[k+4]-m[k+3])+0.5, C[k],'--',color='red')
        ax2.plot(np.arange(0.5,m[-1]-m[6]), Srp, '--',color='red')
    
    for ax in axes[0,:]:
        ax.set_xticklabels([])
    axes[0,0].set_ylabel('Layer thickness [km]')
    for i in range(3):
        axes[0,i].set_title('Layer %d'%i)
        axes[1,i].set_title('Layer %d'%i)
        axes[1,i].set_xticks(np.arange(0.5,m[i+4]-m[i+3]))
        axes[1,i].set_xticklabels(Clabels[m[i+3]-m[3]:m[i+4]-m[3]], fontsize=13)
    for i in range(1,3):
        axes[0,i].set_yticklabels([])
        axes[1,i].set_yticklabels([])
    axes[1,0].set_ylabel('mol%')
    ax2.set_xticklabels([])
    ax2.tick_params(length=0)
    ax2.set_xlabel('Mantle serpentinization')
    ax2.set_ylabel('wt%')
    fig.subplots_adjust(wspace=0.03)
    plt.colorbar(im, ax=np.append(axes[0:2,0:3],ax2))

def plotModels(fin, fin_true=None, mode='ph'):
    period = 8
    m = np.cumsum([0,35,35,35,6,7,6,18], dtype=int)
    Clabels = ['Opx','Fsp','Fsp','Cpx','Bio','q',
              'Opx','Fsp','Fsp','Fsp','Cpx','Bio','q',
              'Opx','Fsp','Fsp','Fsp','Sp','Cpx']  
    if mode=='ox':
        m = np.cumsum([0,35,35,35,8,8,7,18], dtype=int)
        Clabels = ['H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
                  'H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
                  'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O']

    Nsamp = int(getNumLines(fin) / period)

    cmap = plt.get_cmap('jet')
    Nmark = len(Line2D.filled_markers)
    fig = plt.figure()
    axes = fig.subplots(3,3)    
    gs = axes[2, 0].get_subplotspec().get_gridspec()
    # remove the underlying axes
    for ax in axes[2, :]:
        ax.remove() 
    ax2 = fig.add_subplot(gs[2, :])
    handles = []

    rank = [26,7,3,24,5,21]
    rms =[3.92618476,1.13189222,1.26774333,1.31058307,1.3700102,1.3863569]
    Nsamp = len(rank)
    for q in range(Nsamp):
        i = rank[q]
#    for i in range(Nsamp):
        H = np.genfromtxt(fin, skip_header=i*period, max_rows=3)
        C = []
        for k in range(3):
            a = np.genfromtxt(fin, skip_header=i*period+3+k, max_rows=1)
            C.append(a)
        Srp = np.genfromtxt(fin, skip_header=i*period+6, max_rows=1)

        if q == 0:
            for k in range(3):
                axes[0,k].plot(np.arange(m[1])+0.5, H[k], '-', color=cmap(q/Nsamp))
                axes[1,k].plot(np.arange(m[k+4]-m[k+3])+0.5, C[k], '-', color=cmap(q/Nsamp))
            h, = ax2.plot(np.arange(0.5,m[-1]-m[6]), Srp, '-', 
                          color=cmap(q/Nsamp), label='Total mean\n RMS=%.2f'%rms[q])
            handles.append(h)
        else:
            for k in range(3):
                axes[0,k].plot(np.arange(m[1])+0.5, H[k], 
                    marker=Line2D.filled_markers[q%Nmark],
                    fillstyle='none',color=cmap(q/Nsamp))
                axes[1,k].plot(np.arange(m[k+4]-m[k+3])+0.5, C[k], 
                    marker=Line2D.filled_markers[q%Nmark],
                    fillstyle='none',color=cmap(q/Nsamp))
            h, = ax2.plot(np.arange(0.5,m[-1]-m[6]), Srp, 
                     marker=Line2D.filled_markers[q%Nmark],
                     fillstyle='none',color=cmap(q/Nsamp), label='GM comp. %d\n RMS=%.2f'%(i,rms[q]))
            handles.append(h)
        
    if fin_true != None:        
        H = np.genfromtxt(fin_true, max_rows=3)
        C = []
        for k in range(3):
            a = np.genfromtxt(fin_true, skip_header=3+k, max_rows=1)
            C.append(a)        
        Srp = np.genfromtxt(fin_true, skip_header=6, max_rows=1)

        for k in range(3):
            axes[0,k].plot(np.arange(m[k+1]-m[k])+0.5, H[k,:], '--',color='black')
            axes[1,k].plot(np.arange(m[k+4]-m[k+3])+0.5, C[k],'--',color='black')
        h,=ax2.plot(np.arange(0.5,m[-1]-m[6]), Srp, '--',color='black', label='True')
        handles.append(h)
    
    for ax in axes[0,:]:
        ax.set_xticklabels([])
    axes[0,0].set_ylabel('Layer thickness [km]')
    for i in range(3):
        axes[0,i].set_title('Layer %d'%i)
        axes[1,i].set_title('Layer %d'%i)
        axes[1,i].set_xticks(np.arange(0.5,m[i+4]-m[i+3]))
        axes[1,i].set_xticklabels(Clabels[m[i+3]-m[3]:m[i+4]-m[3]])
    for i in range(1,3):
        axes[0,i].set_yticklabels([])
        axes[1,i].set_yticklabels([])
    axes[1,0].set_ylabel('vol%')
    ax2.set_xticklabels([])
    ax2.tick_params(length=0)
    ax2.set_xlabel('Mantle serpentinization')
    ax2.set_ylabel('wt%')
    fig.subplots_adjust(wspace=0.03)

    fig.legend(handles=handles,loc='center right',fontsize='small')
#
#def plotModels(fin, fin_true=None, mode='ph'):
#    period = 8
#    m = np.cumsum([0,35,35,35,6,7,6,18], dtype=int)
#    Clabels = ['Opx','Fsp','Fsp','Cpx','Bio','q',
#              'Opx','Fsp','Fsp','Fsp','Cpx','Bio','q',
#              'Opx','Fsp','Fsp','Fsp','Sp','Cpx']  
#    if mode=='ox':
#        m = np.cumsum([0,35,35,35,8,8,7,18], dtype=int)
#        Clabels = ['H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
#                  'H$_2$O','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','SiO$_2$',
#                  'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O']
#
#    Nsamp = int(getNumLines(fin) / period)
#    
#    cmap = plt.get_cmap('jet')
#    Nmark = len(Line2D.filled_markers)
#    fig = plt.figure()
#    axes = fig.subplots(3,1)
#    handles = []
#    for i in range(Nsamp):       
#        H = np.genfromtxt(fin, skip_header=i*period, max_rows=3)
#        C = np.empty(0)
#        for k in range(3):
#            a = np.genfromtxt(fin, skip_header=i*period+3+k, max_rows=1)
#            C = np.concatenate([C,a])        
#        Srp = np.genfromtxt(fin, skip_header=i*period+6, max_rows=1)
#
#        for k in range(3):
#            xi = range(m[k],m[k+1])
##            xc = np.linspace(m[k],m[k+1]-1,300)
##            spl = interpolate.PchipInterpolator(xi, H[k,:])           
##            axes[0].plot(xc, spl(xc), '-',color=cmap(i/Nsamp))
#            axes[0].plot(xi, H[k,:], marker=Line2D.filled_markers[i%Nmark],
#                fillstyle='none',color=cmap(i/Nsamp))
#        axes[1].plot(C, marker=Line2D.filled_markers[i%Nmark],
#            fillstyle='none',color=cmap(i/Nsamp))
#        h, = axes[2].plot(Srp, marker=Line2D.filled_markers[i%Nmark],
#            fillstyle='none',color=cmap(i/Nsamp), label=str(i))
#        handles.append(h)
#
#    if fin_true != None:        
#        H = np.genfromtxt(fin_true, max_rows=3)
#        C = np.empty(0)
#        for k in range(3):
#            a = np.genfromtxt(fin_true, skip_header=3+k, max_rows=1)
#            C = np.concatenate([C,a])        
#        Srp = np.genfromtxt(fin_true, skip_header=6, max_rows=1)
#
#        for k in range(3):
#            xi = range(m[k],m[k+1])
##            xc = np.linspace(m[k],m[k+1]-1,300)
##            spl = interpolate.PchipInterpolator(xi, H[k,:])           
##            axes[0].plot(xc, spl(xc), '-',color=cmap(i/Nsamp))
#            axes[0].plot(xi, H[k,:], '--',color='black')
#        axes[1].plot(C,'--',color='black')
#        h, = axes[2].plot(Srp, '--',color='black', label='True')
#        handles.append(h)
#        
#    tix=[m[1]*0.5, m[2]-(m[1]-m[0])*0.5, m[3]-(m[2]-m[1])*0.5]
#    labels=['Layer 0','Layer 1','Layer 2']
#    axes[0].set_xticks(tix)
#    axes[0].set_xticklabels(labels)
#    axes[0].tick_params(length=0)
#    ylim = axes[0].get_ylim()
#    for i in range(1,3):
#        axes[0].plot([m[i]-0.5,m[i]-0.5], [ylim[0], ylim[1]],color='black',linewidth=0.5)
#    axes[0].set_ylim(ylim)
#    axes[0].set_ylabel('Layer thickness [km]')
#
#    axes[1].set_xticks(range(m[6]-m[3]))
#    axes[1].set_xticklabels(Clabels)
#    ylim = axes[1].get_ylim()
#    for i in range(1,3):
#        axes[1].plot([m[3+i]-m[3]-0.5,m[3+i]-m[3]-0.5], [ylim[0], ylim[1]],color='black',linewidth=0.5)
#    axes[1].set_ylim(ylim)
#    axes[1].set_ylabel('vol%')
#    
#    axes[2].set_xticklabels([])
#    axes[2].tick_params(length=0)
#    axes[2].set_xlabel('Mantle serpentinization')
#    axes[2].set_ylabel('wt%')
#
#    fig.legend(handles=handles,loc='center right',fontsize='x-small')
#
# 
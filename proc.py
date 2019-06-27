import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import sklearn.mixture as mix
from compro import material
import ownplt

# get number of lines without loading the file into RAM
def getNumLines(fname):
    import subprocess
    return int(subprocess.check_output('wc -l {}'.format(fname), shell=True).split()[0])


def saveModelVec(v, fname, mode='ph'):
    m = np.cumsum([0,35,35,35,6,7,6,18], dtype=int)
    if mode == 'ox':
        m = np.cumsum([0,35,35,35,8,8,7,18], dtype=int)
    f = open(fname, 'a')
    for i in range(len(m)-1):
        np.savetxt(f, [v[m[i]:m[i+1]]], fmt='%.2f')
    f.write('\n')
    f.close()

def filtRMS(fin, fout, fRMSin, fRMSout, blockLen, maxRMS, T=None):

    period = blockLen+1
    
    Nsamp = getNumLines(fin) / period
    print('Number of samples:', Nsamp)
    Nsamp = int(Nsamp)
    
    f = open(fin, 'r')
    buf = f.readlines()
    f.close()

    rms = np.genfromtxt(fRMSin)

    f = open(fout, 'w')
    frms = open(fRMSout, 'w')
    offset = 0
    accepted = 0
    for i in range(np.size(rms,0)):#range(Nsamp):
        if np.all(rms[i, :3] < maxRMS) and (T==None or rms[i,3]==T):
            for j in range(blockLen):
                p = np.fromstring(buf[offset+j][:], sep=' ')
                np.savetxt(f, [p], fmt='%.2f')
            f.write('\n')
            np.savetxt(frms, [rms[i,:]], fmt='%g')
            accepted += 1
        offset += period
        if i%1000 == 0:
            print('%d / %d' % (i, Nsamp))
    print('Samples with RMS < %g: %d' % (maxRMS, accepted))
    f.close()            
    frms.close()

def filtRMScat(prefin, fout, prefRMS, prefRMSout, N, blockLen, maxRMS, T=None):
    os.system('rm %s' % fout)
    os.system('rm %s' % prefRMSout)
    for i in range(N):
        filtRMS(prefin+str(i), fout+str(i), prefRMS+str(i), 
                prefRMSout+str(i), blockLen, maxRMS, T)
        os.system('cat %s >> %s' % (fout+str(i), fout))
        os.system('rm %s' % fout+str(i))
        os.system('cat %s >> %s' % (prefRMSout+str(i), prefRMSout))
        os.system('rm %s' % prefRMSout+str(i))


def ph2ox(fin, fout, refComp, Noxides, Nphases, offset, period, 
          blockLen, phaseBlockOffset, phaseBlockLen):

    mat = list(refComp)
    for i in range(len(refComp)):
        mat[i] = material()
        mat[i].load(refComp[i], Nphases[i], Noxides[i])

    Nsamp = (getNumLines(fin) - offset) / period
    print('Number of samples:', Nsamp)
    Nsamp = int(Nsamp)
    
    f = open(fin, 'r')
    buf = f.readlines()
    f.close()

    f = open(fout, 'w')
    for i in range(Nsamp):
        for j in range(blockLen):
            p = np.fromstring(buf[offset+j][:], sep=' ')
            if j >= phaseBlockOffset and j < phaseBlockOffset + phaseBlockLen:
                mat[j-phaseBlockOffset].phacomtxt[:,1] = p.T
                p = mat[j-phaseBlockOffset].getOxMolPerc()
            np.savetxt(f, [p], fmt='%.2f')
        offset += period
        f.write('\n')
        if i%1000 == 0:
            print('%d / %d' % (i, Nsamp))
    f.close()            


def meanAll(fin, fout, offset, period, blockLen):

    Nsamp = (getNumLines(fin) - offset) / period
    print('Number of samples:', Nsamp)
    Nsamp = int(Nsamp)
    
    f = open(fin, 'r')
    buf = f.readlines()
    f.close()
    
    p = list(range(blockLen))
    off = offset
    for i in range(blockLen):
        off = offset + i
        p[i] = np.zeros_like(np.fromstring(buf[off][:], sep=' '))
        for j in range(Nsamp):
            p[i] += np.fromstring(buf[off][:], sep=' ')
            off += period
            if j%1000 == 0:
                print('%d %d / %d' % (i, j, Nsamp))
        p[i] /= Nsamp
    
    f = open(fout, 'w')
    for i in range(blockLen):
        np.savetxt(f, [p[i]], fmt='%.2f')
    f.close()
    
def covAll(f_samp, offset, period, blockOffset, blockLen, mode=None, vmax=None):
    Nsamp = (getNumLines(f_samp) - offset) / period
    print('Number of samples:', Nsamp)
    Nsamp = int(Nsamp)
    
    # get length of model vector
    M = sum(len(np.genfromtxt(f_samp, skip_header=i, max_rows=1)) \
            for i in range(blockOffset,blockOffset+blockLen))

    f = open(f_samp, 'r')
    buf = f.readlines()
    f.close()

    # build matrix of samples
    p = np.zeros((Nsamp,M))
    off = offset + blockOffset
    for i in range(Nsamp):
        a = np.empty(0)
        for j in range(blockLen):
            a = np.append(a, np.fromstring(buf[off+j][:], sep=' '))
        p[i,:] = a
        off += period
        if i%1000 == 0:
            print('%d / %d' % (i, Nsamp))

    C = np.cov(p, rowvar=False)

    ownplt.plotCov(C, mode, vmax=vmax)
    
    std = np.sqrt(np.diag(C))
    print('std = ', std)
    return std

def meanComp(fin, fout, offset, period, blockOffset, eps):

    Nsamp = (getNumLines(fin) - offset) / period
    print('Number of samples:', Nsamp)
    Nsamp = int(Nsamp)
    
    f = open(fin, 'r')
    buf = f.readlines()
    f.close()
    
    offset += blockOffset
    p = np.zeros_like(np.fromstring(buf[offset][:], sep=' '))
    for i in range(Nsamp):
        c = np.fromstring(buf[offset][:], sep=' ')
        c[c < eps] = eps
        p += np.log(c)
        offset += period
        if i%1000 == 0:
            print('%d / %d' % (i, Nsamp))
    p = np.exp(p / Nsamp)
    p *= 100./sum(p)
    np.savetxt(fout, [p], fmt='%.2f')
    
def covComp(fin, offset, period, blockOffset):
    return

def histZ(fin, foutPrefix, z01, Nz, offset, period, blockOffset, blockLen):
# NOTE: number and horizontal position of nodes in all layers must be the same!
# zBin must be descending
    # bathymetry
    y0, z0 = np.genfromtxt("test/rapids4/in.txt", max_rows=2)
    z0 = np.flip(z0,0)
    # basement top
    y1, z1 = np.genfromtxt("test/rapids4/in.txt", skip_header=6, max_rows=2)
    z1 = np.flip(z1,0)
# NOTE: this fixes a bug in init
    spl = interpolate.PchipInterpolator(y0, z0)
    H0 = z1 + spl(y1)
#    
    zBin = np.linspace(z01[0], z01[1], Nz)
    
    Nx = len(H0)
    hist = list(range(blockLen))
    for i in range(blockLen):
        hist[i] = np.zeros((Nz, Nx))
    h = zBin[0] - zBin[1]
    
    Nsamp = (getNumLines(fin) - offset) / period
    print('Number of samples:', Nsamp)
    Nsamp = int(Nsamp)
    
    f = open(fin, 'r')
    buf = f.readlines()
    f.close()

    for i in range(Nsamp):
        z = -H0.copy()
        for j in range(blockOffset, blockLen):
            z -= np.fromstring(buf[offset+j][:], sep=' ')
#            if np.any(z < z01[1]):
#                print('error')
#                return
            for k in range(Nx):
                hist[j][int((z01[0] - z[k]) / h)][k] += 1
        offset += period
            
    fig, axes = plt.subplots(3,1)
    for i in range(blockLen):
        hist[i] /= Nsamp
        np.savetxt('%s%d'%(foutPrefix,i), hist[i])
        im = axes[i].imshow(hist[i])
        axes[i].set_aspect('auto')
    fig.colorbar(im, ax=axes, orientation='horizontal')
 

def multimodal(f_samp, Ncomp_range, f_rms=None, blockLen=7, period=8, mode='ph'):
    Nsamp = getNumLines(f_samp) / period
    print('Number of samples:', Nsamp)
    Nsamp = int(Nsamp)

    # get length of model vector
    M = sum(len(np.genfromtxt(f_samp, skip_header=i, max_rows=1)) for i in range(blockLen))
    
    f = open(f_samp, 'r')
    buf = f.readlines()
    f.close()

    # build matrix of samples
    samp = np.zeros([Nsamp, M])
    off = 0
    for i in range(Nsamp):
        a = np.empty(0)
        for j in range(blockLen):
            a = np.append(a, np.fromstring(buf[off+j][:], sep=' '))
        samp[i,:] = a
        off += period
        if i%1000 == 0:
            print('%d / %d' % (i, Nsamp))
#    # weights
#    rms_s, rms_g, rms_e, T = np.genfromtxt(f_rms, unpack=True)
#    phi = 0.5 * (rms_s**2 * 88381 + rms_g**2 * 201 + rms_e**2 * 462)
#    phi_s = 0.5 * (rms_s**2 + rms_g**2 + rms_e**2) / (1./88381 + 1./201 + 1./462)

    # estimate optimal number of components by BIC
    minBIC = np.infty
    gmm = None
    bic = np.zeros(max(Ncomp_range)+1)
    for i in Ncomp_range:
        trial_gmm = mix.GaussianMixture(n_components=i, #tol=1.e-8, reg_covar=1.e-11, 
                                  covariance_type='full', verbose=1)    
        trial_gmm.fit(samp)
        bic[i] = trial_gmm.bic(samp)
        print('Ncomp=%d/%d, BIC=%g' % (i, len(Ncomp_range), bic[i]))
        if bic[i] < minBIC:
            minBIC = bic[i]
            gmm = trial_gmm
#    plt.figure()
#    plt.plot(Ncomp_range, bic[Ncomp_range],'.-')
#    plt.gca().set_xlabel("Number of components")
#    plt.gca().set_ylabel("BIC")
#    plt.title('Optimal Ncomp=%d' % gmm.n_components)
    np.savetxt('bic', np.stack([Ncomp_range, bic[Ncomp_range]]).T)
    f_means = f_samp + '_means'
    os.system('rm %s' % f_means)        
    for i in range(gmm.n_components):
        # save means
        saveModelVec(gmm.means_[i], f_means, mode=mode)
        
#        # plot covariances
#        m = np.cumsum([0,35,35,35,6,7,6,18], dtype=int)
#        ownplt.plotCov(gmm.covariances_[i][m[0]:m[7],m[0]:m[7]],
#                       tit='GMM component '+str(i))
        
   

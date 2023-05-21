import numpy as np
from numpy.linalg import eigvalsh
from scipy.interpolate import CubicSpline

def diag_spline(dists, matcomp):
    spline = CubicSpline(dists, matcomp)
    return spline

ndists=56

E=np.zeros((ndists,34,34))
Etest=np.zeros((ndists))
dists=np.zeros((ndists))

for i in range(34):
    for j in range(34):
        if i+1 < 10:
            if j+1 < 10:
                Etest = np.loadtxt('s0'+str(i+1)+'-0'+str(j+1)+'.dat')
            else:
                Etest = np.loadtxt('s0'+str(i+1)+'-'+str(j+1)+'.dat')
        else:
            if j+1 < 10:
                Etest = np.loadtxt('s'+str(i+1)+'-0'+str(j+1)+'.dat')
            else:
                Etest = np.loadtxt('s'+str(i+1)+'-'+str(j+1)+'.dat')
        for k in range(ndists):
            E[k][i][j] = Etest[k,1]
            dists[k] = Etest[k,0]

npoints=500
Esp=np.zeros((npoints,2))
xs=np.linspace(2.0,10.0,npoints)


for i in range(34):
    for j in range(34):
        matcomp = diag_spline(dists, E[:,i,j])
        for k in range(npoints):
            Esp[k,0] = xs[k]
            Esp[k,1] = matcomp(xs[k])
                    #if i==1:
                    #    print(E[k][i][j])
        if i+1 < 10:
            if j+1 < 10:
                np.savetxt('splines-python-test/s0'+str(i+1)+'-0'+str(j+1)+'.dat',Esp,fmt='%6f')
            else:
                np.savetxt('splines-python-test/s0'+str(i+1)+'-'+str(j+1)+'.dat',Esp,fmt='%6f')
        else:
            if j+1 < 10:
                np.savetxt('splines-python-test/s'+str(i+1)+'-0'+str(j+1)+'.dat',Esp,fmt='%6f')
            else:
                np.savetxt('splines-python-test/s'+str(i+1)+'-'+str(j+1)+'.dat',Esp,fmt='%6f')
        
print(E[10].shape)
eigenvalues, eigenvectors = np.linalg.eig(E[10])

print(eigenvalues)
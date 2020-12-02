
import numpy as np
from numpy import linalg as LA


#Function that computes binned data.
def binnedData(minx,maxx,interval,dataForBins,dataToBin):
    data=[]
    bins=np.arange(minx, maxx, interval)
    mean=[]
    std=[]
    stderr=[]
    totsum=0
    for i in range(len(bins)): 
        #Getting data with this conditions doesn't include nan values
        data.append(dataToBin[(dataForBins >= minx + i*interval ) & (dataForBins < minx + (i+1)*interval )] )
        totsum=totsum+len(data[i])
        mean.append(np.nanmean(data[i]))
        std.append(np.nanstd(data[i]))
        stderr.append( std[i]/np.sum(~np.isnan(data[i]))**0.5 ) 
    stderr=np.asarray(stderr)
    return bins, data, mean, std, stderr, totsum    


def get_coserr(path,target,stepsDir): 
    #vector from path to target
    vector =  target-path
    norm=np.sqrt(vector[:,0]**2+vector[:,1]**2+vector[:,2]**2)
    #report files where the patches encounter target  before the last time point
    if np.sum(norm[0:-1]==0)>0:
        print(file)

    #make matrix with distance to target
    norm_matrix = np.tile(norm,(3,1))
    unit_vector = vector/np.transpose(norm_matrix)

    #cos(err) as dot product of unit vectors (len(tracks)-1)
    coserr=np.sum(stepsDir*unit_vector[0:-1],1)
    return norm[0:-1] , np.arccos(coserr) , coserr

def get_plane_vecs5(path):
    
    plane_vecs = []
    for k in range(len(path)):
        if k >= 2 or k <= len(path) - 3:
            xs=path[k-2:k+3,0]
            ys=path[k-2:k+3,1]
            zs=path[k-2:k+3,2]
            
        if k<2:
            xs=path[k:k+3,0]
            ys=path[k:k+3,1]
            zs=path[k:k+3,2]
            
        if k > len(path)-3:
            xs=path[k-2:k+1,0]
            ys=path[k-2:k+1,1]
            zs=path[k-2:k+1,2]    
            
        # do fit
        tmp_A = []
        tmp_b = []
        for i in range(len(xs)):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_b.append(zs[i])
        b = np.matrix(tmp_b).T
        A = np.matrix(tmp_A)
        fit = (A.T * A).I * A.T * b
        errors = b - A * fit
        residual = np.linalg.norm(errors)
        #print "%f x + %f y + %f = z" % (fit[0], fit[1], fit[2])
        plane_vecs.append([fit[0,0], fit[1,0], -1])        
    
    plane_vecs= np.asarray(plane_vecs)    
    #print plane_vecs    
    return plane_vecs


def get_plane_vecs4(path):
    
    plane_vecs = []
    for k in range(len(path)):
        if k>0 or k < len(path) - 2:
            xs=path[k-1:k+3,0]
            ys=path[k-1:k+3,1]
            zs=path[k-1:k+3,2]
            
        if k==0:
            xs=path[0:3,0]
            ys=path[0:3,1]
            zs=path[0:3,2]
            
        if k >= len(path)-2:
            xs=path[k-2:k+1,0]
            ys=path[k-2:k+1,1]
            zs=path[k-2:k+1,2]    
            
        # do fit
        tmp_A = []
        tmp_b = []
        for i in range(len(xs)):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_b.append(zs[i])
        b = np.matrix(tmp_b).T
        A = np.matrix(tmp_A)
        fit = (A.T * A).I * A.T * b
        errors = b - A * fit
        residual = np.linalg.norm(errors)
        #print "%f x + %f y + %f = z" % (fit[0], fit[1], fit[2])
        plane_vecs.append([fit[0,0], fit[1,0], -1])        
    
    plane_vecs= np.asarray(plane_vecs)    
    #print plane_vecs    
    return plane_vecs

def get_plane_vecs(path):
    
    plane_vecs = []
    for k in range(len(path)):
        if k>0 or k < len(path) - 1:
            xs=path[k-1:k+2,0]
            ys=path[k-1:k+2,1]
            zs=path[k-1:k+2,2]
            
        if k==0:
            xs=path[0:3,0]
            ys=path[0:3,1]
            zs=path[0:3,2]
            
        if k==len(path)-1:
            xs=path[-3:,0]
            ys=path[-3:,1]
            zs=path[-3:,2]    
            
        # do fit
        tmp_A = []
        tmp_b = []
        for i in range(len(xs)):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_b.append(zs[i])
        b = np.matrix(tmp_b).T
        A = np.matrix(tmp_A)
        fit = (A.T * A).I * A.T * b
        errors = b - A * fit
        residual = np.linalg.norm(errors)
        #print "%f x + %f y + %f = z" % (fit[0], fit[1], fit[2])
        plane_vecs.append([fit[0,0], fit[1,0], -1])        
    
    plane_vecs= np.asarray(plane_vecs)    
    #print plane_vecs    
    return plane_vecs
        

def get_coserr_plane(path,target,stepsDir): 
    #vector from path to target
    vectorTarg =  target-path
    normVecTarg=np.sqrt(vectorTarg[:,0]**2+vectorTarg[:,1]**2+vectorTarg[:,2]**2)
    
    #report files where the patches encounter target  before the last time point
    if np.sum(normVecTarg[0:-1]==0)>0:
        print(file)

    #make matrix with distance to target
    #norm_matrixTarg = np.tile(normVecTarg,(3,1))
    #unit_vectorTarg = vectorTarg/np.transpose(norm_matrixTarg)

    #find ortogonal vectors to tangent plane for each position on path
    vecsPlane = get_plane_vecs(path)
    
    #project vectors on plane
    
    #normalize ortogonals vectors to plane
    normVecP = np.sqrt(vecsPlane[:,0]**2+vecsPlane[:,1]**2+vecsPlane[:,2]**2)
    norm_matrixVecP = np.tile(normVecP,(3,1))
    vecPunit= vecsPlane/np.transpose(norm_matrixVecP)
    
    #get projection of step directions on plane
    cos_VecPUnit_x_stepsDir = np.sum(stepsDir*vecPunit[0:-1],1)
    cos_VecPUnit_x_stepsDir_matrix = np.tile(cos_VecPUnit_x_stepsDir,(3,1))
    stepsDirProjVecP = vecPunit[0:-1]*np.transpose(cos_VecPUnit_x_stepsDir_matrix)    
    stepsDirProjPlane= stepsDir - stepsDirProjVecP
    
    #get projection of target vectors on plane
    cos_VecPUnit_x_vecTarg = np.sum(vectorTarg*vecPunit,1)
    cos_VecPUnit_x_vecTarg_matrix = np.tile(cos_VecPUnit_x_vecTarg,(3,1))
    vecTargProjVecP = vecPunit*np.transpose(cos_VecPUnit_x_vecTarg_matrix) 
    vecTargProjPlane = vectorTarg - vecTargProjVecP
    
    #normalize projection of step directions on plane
    normstepsDirProjPlane = np.sqrt(stepsDirProjPlane[:,0]**2+stepsDirProjPlane[:,1]**2+stepsDirProjPlane[:,2]**2)
    norm_matrixstepsDirProjPlane = np.tile(normstepsDirProjPlane,(3,1))
    stepsDirProjPlaneunit= stepsDirProjPlane/np.transpose(norm_matrixstepsDirProjPlane)
    
    #normalize projection of target vector on plane
    normvecTargProjPlane = np.sqrt(vecTargProjPlane[:,0]**2+vecTargProjPlane[:,1]**2+vecTargProjPlane[:,2]**2)
    norm_matrixvecTargProjPlane = np.tile(normvecTargProjPlane,(3,1))
    vecTargProjPlaneunit= vecTargProjPlane/np.transpose(norm_matrixvecTargProjPlane)
    
    #cos(err) as dot product of unit vectors on plane (len(tracks)-1)
    coserr=np.sum(stepsDirProjPlaneunit*vecTargProjPlaneunit[0:-1],1)
    return normVecTarg[0:-1] , np.arccos(coserr) , coserr

def sum_squared_displacements(path, maxstep):
    #get subpaths such that succesive steps have step lengt less than maxstep
    #skip timepoints where there are not patch. 
    subpaths=[]
   
    #get first position that contains data
    for i in range(0,len(path)):
        if ~np.isnan(path[i,0]):
            break
    i_ini=i
    
    
    for i in range(i_ini+1,len(path)):
        #if current and previous point have data and the step lenght is long
        if ~np.isnan(path[i,0]) and ~np.isnan(path[i-1,0]) and LA.norm(path[i,:]-path[i-1,:]) > maxstep :
            subpaths.append([i_ini,i-1])
            i_ini=i
            
        #if current point is nan and the previous point has data close subpath 
        if np.isnan(path[i,0]) and ~np.isnan(path[i-1,0]):
            subpaths.append([i_ini,i-1])
            
        #if the current point has data but the previous point was nan, set i_ini      
        if ~np.isnan(path[i,0]) and np.isnan(path[i-1,0]):
            i_ini=i
    
    #add last subpath 
    if ~np.isnan(path[len(path)-1,0]):       
        subpaths.append([i_ini,len(path)-1])
        
    sum_norm_disp = np.zeros(len(path)-1)
    sum_squared_disp = np.zeros(len(path)-1)
    sum_squared_disp_squared = np.zeros(len(path)-1)    
    #number of each intervals computed
    n_squared_disp = np.zeros(len(path)-1)
    
    for subp in range(len(subpaths)):
        #get msds in subpath
        sublen = subpaths[subp][1] - subpaths[subp][0]
        for interval in range(1,sublen+1):
            #compute along the path
            #for pos in range( first, frist + sublen - interval + 1)
            for pos in range( subpaths[subp][0], subpaths[subp][0] + sublen - interval + 1):
                norm_disp = LA.norm( path[pos,:] - path[pos+interval,:])
                sum_norm_disp[interval-1] += norm_disp
                sum_squared_disp[interval-1] += norm_disp**2
                sum_squared_disp_squared[interval-1] += norm_disp**4
                n_squared_disp[interval-1] += 1
    return sum_norm_disp, sum_squared_disp, sum_squared_disp_squared, n_squared_disp
    

    
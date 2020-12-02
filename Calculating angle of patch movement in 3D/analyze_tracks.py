
import numpy as np
from os import listdir
from os.path import isfile, join
import re
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib 

from functions_tracks import get_coserr
from functions_tracks import get_coserr_plane
from functions_tracks import binnedData
from functions_tracks import sum_squared_displacements

plt.ion()

#1 is a, 2 is alpha

stepsizes1=np.empty(0)
stepsizes2=np.empty(0)

distances=np.empty(0)
coserrs1=np.empty(0)
coserrs2=np.empty(0)
errs1=np.empty(0)
errs2=np.empty(0)


coserrs1plane=np.empty(0)
coserrs2plane=np.empty(0)

errs1plane=np.empty(0)
errs2plane=np.empty(0)

costurns1=np.empty(0)
costurns2=np.empty(0)
turns1=np.empty(0)
turns2=np.empty(0)

siz=100
sum_norm_disps1=np.zeros(siz) 
sum_squared_disps1=np.zeros(siz) 
sum_squared_disps_squared1=np.zeros(siz)
n_squared_disps1=np.zeros(siz)

sum_norm_disps2=np.zeros(siz) 
sum_squared_disps2=np.zeros(siz) 
sum_squared_disps_squared2=np.zeros(siz)
n_squared_disps2=np.zeros(siz)


'''
path="./tracks_same_sex_3D_katie/"
all = [f for f in listdir(path) if re.match(r'.*\.csv',f)]
celltype="same"
co="red"
colora="red"
coloralpha="blue"
files=all
sampling=1.5; #min
'''
#wt same
'''
path="./paths/"
same = [f for f in listdir(path) if re.match(r'.*whi5\.csv',f)]
#all = [f for f in listdir(path) if re.match(r'.*\.csv',f)]
celltype="wt same type"
co="blue"
files = same
colora="black"
coloralpha="green"
'''

#wt opposite

path="./paths/"
opposite = [f for f in listdir(path) if re.match(r'.*opp\.csv',f)]
#all = [f for f in listdir(path) if re.match(r'.*\.csv',f)]
celltype="wt opposite type"
co="red"
files = opposite
colora="magenta"
coloralpha="cyan"


for i,file in enumerate(files):
#for i,file in enumerate (opposite[0:2]):

    #read tracks from 2 cells
    tracks =np.genfromtxt(path+file,delimiter=',')
    
    #step size (len(tracks)-1)
    stepsizefile = np.sqrt((tracks[0:-1,[0,3]]-tracks[1:,[0,3]])**2 +(tracks[0:-1,[1,4]]-tracks[1:,[1,4]])**2 + (tracks[0:-1,[2,5]]-tracks[1:,[2,5]])**2)
    #add step size to the whole pool 
    stepsizes1 = np.concatenate((stepsizes1,stepsizefile[:,0]))
    stepsizes2 = np.concatenate((stepsizes2,stepsizefile[:,1]))

    #vector steps (len(tracks)-1)
    vector_steps = tracks[1:,:] - tracks[0:-1,:]
    vector_steps1 = vector_steps[:,0:3]
    vector_steps2 = vector_steps[:,3:6]
    matrix_stepsize1 = np.tile(stepsizefile[:,0],(3,1))
    matrix_stepsize2 = np.tile(stepsizefile[:,1],(3,1))
    unit_vector_steps1 = vector_steps1/np.transpose(matrix_stepsize1)
    unit_vector_steps2 = vector_steps2/np.transpose(matrix_stepsize2)
   
    #TURNING ANGLES
    costurn1=np.sum(unit_vector_steps1[0:-1]*unit_vector_steps1[1:],1)   
    turn1=np.arccos(costurn1)
    costurn2=np.sum(unit_vector_steps2[0:-1]*unit_vector_steps2[1:],1)   
    turn2=np.arccos(costurn2)
    
    costurns1=np.concatenate((costurns1,costurn1))
    costurns2=np.concatenate((costurns2,costurn2))
    turns1=np.concatenate((turns1,turn1))
    turns2=np.concatenate((turns2,turn2))
    
        
    #COS ERR TO OTHER PATCH
    distance, err1, coserr1 = get_coserr(tracks[:,0:3],tracks[:,3:],unit_vector_steps1)
    distance, err2, coserr2 = get_coserr(tracks[:,3:],tracks[:,0:3],unit_vector_steps2)
    #add distances from file without the last value to the whole pool (len(tracks)-1)
    distances=np.concatenate((distances,distance))
    #Add coserrs to whole pool
    errs1=np.concatenate((errs1,err1[:]))
    errs2=np.concatenate((errs2,err2[:]))
    coserrs1=np.concatenate((coserrs1,coserr1[:]))
    coserrs2=np.concatenate((coserrs2,coserr2[:]))
    
    #COS ERR ON PLANE TO OTHER PATCH
    distance, err1plane,coserr1plane = get_coserr_plane(tracks[:,0:3],tracks[:,3:],unit_vector_steps1)
    distance, err2plane,coserr2plane = get_coserr_plane(tracks[:,3:],tracks[:,0:3],unit_vector_steps2)
    #add distances from file without the last value to the whole pool (len(tracks)-1)
    #distances=np.concatenate((distances,distance))
    #Add coserrs to whole pool
    coserrs1plane=np.concatenate((coserrs1plane,coserr1plane[:]))
    coserrs2plane=np.concatenate((coserrs2plane,coserr2plane[:]))
    errs1plane=np.concatenate((errs1plane,err1plane[:]))
    errs2plane=np.concatenate((errs2plane,err2plane[:]))
    
    #Sums of displacement, squared displacement
    
    maxstepsize=3.5
    sum_norm_disp1, sum_squared_disp1, sum_squared_disp_squared1, n_squared_disp1 = \
    sum_squared_displacements(tracks[:,0:3], maxstepsize)    
    sum_norm_disps1[0:len(sum_norm_disp1)] +=  sum_norm_disp1   
    sum_squared_disps1[0:len(sum_squared_disp1)] +=  sum_squared_disp1  
    sum_squared_disps_squared1[0:len(sum_squared_disp_squared1)] +=  sum_squared_disp_squared1 
    n_squared_disps1[0:len(n_squared_disp1)] +=  n_squared_disp1
    
    sum_norm_disp2, sum_squared_disp2, sum_squared_disp_squared2, n_squared_disp2 = \
    sum_squared_displacements(tracks[:,3:6], maxstepsize)    
    sum_norm_disps2[0:len(sum_norm_disp2)] +=  sum_norm_disp2   
    sum_squared_disps2[0:len(sum_squared_disp2)] +=  sum_squared_disp2  
    sum_squared_disps_squared2[0:len(sum_squared_disp_squared2)] +=  sum_squared_disp_squared2 
    n_squared_disps2[0:len(n_squared_disp2)] +=  n_squared_disp2
    

#==============    
    
matplotlib.rcParams.update({'font.size': 12})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 2



#MSD computation
msd1 = sum_squared_disps1/n_squared_disps1
stddevmsd1 = np.sqrt(sum_squared_disps_squared1/n_squared_disps1 - msd1**2) 
stderrmsd1 = stddevmsd1/n_squared_disps1**0.5
msd2 = sum_squared_disps2/n_squared_disps2
stddevmsd2 = np.sqrt(sum_squared_disps_squared2/n_squared_disps2 - msd2**2) 
stderrmsd2 = stddevmsd2/n_squared_disps2**0.5

#FIGURE MSD
#plt.close()
'''
plt.figure(1)
msd1=np.concatenate(([0],msd1))
msd2=np.concatenate(([0],msd2))
stderrmsd1=np.concatenate(([0],stderrmsd1))
stderrmsd2=np.concatenate(([0],stderrmsd2))

plt.errorbar(sampling*np.arange(len(msd1)), msd1, 2*stderrmsd1 , label="a "+celltype, color =colora, capsize=2)
plt.errorbar(sampling*np.arange(len(msd2)), msd2, 2*stderrmsd2 , label="alpha "+celltype, color =coloralpha, capsize=2)
#plt.xlim(0,50*10)
#plt.ylim(0,10000)
plt.title("MSD max step size "+str(maxstepsize)+"um")
plt.xlabel('Interval (min)')
plt.ylabel(r'$MSD \pm 2SE(\mu m^2)$')
plt.ylim(0,10)
plt.xlim(0,30)
plt.legend(loc=2)
#fig.set_size_inches(5, 4)
'''

#HISTOGRAM  TURNS
'''
plt.figure(2)
#bins=np.arange(-1, 1+0.05, 0.05)
bins=np.arange(0, np.pi + np.pi/100, np.pi/100)
plt.title(celltype)
turns0pi=np.concatenate((turns1[~np.isnan(turns1)],turns2[~np.isnan(turns2)]))
np.save("turns0pi",turns0pi)
#plt.hist(costurns2[~np.isnan(costurns2)] ,color='blue',label="a",alpha=0.6)
plt.hist(turns0pi,bins,color='blue',label="same")
plt.xlabel("turn(0-pi)")
plt.ylabel("Frequency")
plt.legend()
plt.show()
'''

#HISTOGRAM COS TURNS
'''
plt.figure(2)
bins=np.arange(-1, 1+0.05, 0.05)
plt.title(celltype)
costurns=np.concatenate((costurns1,costurns2))
plt.hist(costurns[~np.isnan(costurns)],bins,color='blue',label="same")
plt.xlabel("cos turn")
plt.ylabel("Frequency")
plt.legend()
plt.show()
'''

#MEAN COS ERR TO OTHER PATCH VS DISTANCE OTHER PATCH
"""
plt.figure(1)
minx=0
maxx=12
interval=1
bins_distances, data_coserrs, mean_coserrs, std_coserrs, stderr_coserrs, totsum_coserr  = \
binnedData(minx,maxx,interval,np.concatenate((distances, distances)),np.concatenate(( coserrs1, coserrs2) ))
plt.errorbar(bins_distances + interval/2.0 , 180/np.pi*np.asarray(mean_coserrs),180/np.pi*2*stderr_coserrs, color='green',label='opposite')
plt.xlabel(r'Separation between patches ($\mu$m)')
plt.ylabel('mean cos(err)')
plt.legend()
plt.ylim(-0.3, 0.7)
plt.xlim(0, 10)
plt.show()
"""

#MEAN  ERR TO OTHER PATCH PROJECTED ON TANGENT PLANE VS DISTANCE OTHER PATCH
plt.figure(1)
minx=0
maxx=12
interval=1
bins_distances, data, mean_data, stddev, stderr, totsum  = \
binnedData(minx,maxx,interval,np.concatenate((distances, distances)),np.concatenate(( errs1plane, errs2plane) ))
plt.errorbar(bins_distances + interval/2.0 ,180/np.pi*np.asarray(mean_data),180/np.pi*2*stderr, color=co,label=celltype, capsize=5)
plt.xlabel(r'Separation between patches ($\mu$m)')
plt.ylabel('mean err plane')
plt.legend()
plt.ylim(30, 120)
plt.xlim(1, 9)
plt.show()


#MEAN COS ERR TO OTHER PATCH VS STEP SIZE
'''
interval=1.8
bins_stepsizes, data_coserrs, mean_coserrs, std_coserrs, stderr_coserrs, totsum_coserr  = \
binnedData(0,5.4,interval,np.concatenate((stepsizes1, stepsizes2)),np.concatenate(( coserrs1, coserrs2) ))
plt.errorbar(bins_stepsizes + interval/2 , mean_coserrs, 2*stderr_coserrs, color='red',label='opposite',capsize=5)
plt.xlabel(r'Step size ($\mu$m)')
plt.ylabel('mean cos(err)')
plt.legend()
plt.show()
'''


#HISTOGRAM STEP SIZES
'''
plt.figure(3)
bins=np.arange(0.0, 6, 0.1)
nonanvalues=~np.isnan(np.concatenate((stepsizes1,stepsizes2)))
stepsizes = np.concatenate((stepsizes1,stepsizes2))[nonanvalues]
np.save("steps",stepsizes)
plt.hist(stepsizes,bins)
plt.axvline(x=np.nanmean(np.concatenate((stepsizes1,stepsizes2))),color='black')
plt.xlabel(r'Step size ($\mu$m)')
plt.xlim(0,6)
#plt.ylim(0,160)
plt.ylabel('Number')
plt.show()
'''

#plot step size vs distance to other patch
#plt.scatter(np.concatenate((distances,distances)), np.concatenate((stepsizes1,stepsizes2)) )
#plt.scatter(distances,stepsizes1,color='blue',label='a')
#plt.scatter(distances,stepsizes2,color='red',label='alpha')
#plt.ylabel(r'Step size ($\mu$m)')
#plt.xlabel(r'Distance to other patch ($\mu$m)')
#plt.legend()
#plt.title('Opposite mating type')
#plt.show()

#plot step size vs step length*cos(err) 
#plt.scatter(stepsizes1,coserrs1,color='blue',label='a')
#plt.scatter(stepsizes2,coserrs2,color='red',label='alpha')
#plt.xlabel(r'Step size ($\mu$m)')
#plt.ylabel('cos(err)')
#plt.legend()
#plt.title('Same mating type')
#plt.show()

#plot length*cos(err) vs distance
#plt.scatter(distances,coserrs1,color='blue',label='a')
#plt.scatter(distances,coserrs2,color='red',label='alpha')
#plt.xlabel(r'Distance to other patch ($\mu$m)')
#plt.ylabel('cos(err)')
#plt.legend()
#plt.title('Opposite mating type')
#plt.show()


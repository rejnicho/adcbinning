import numpy as np 
import matplotlib.pyplot as plt
import pickle as pkl
import scipy.interpolate
import time
from scipy.signal import savgol_filter
import copy

def checkbins(pcounts, pbins):
    issueindex = []
    for x in range(len(pbins)-1):
        if pbins[x] != pbins[x-1] +1:
            issueindex.append(x)
            if len(issueindex) != 0: #take the last index in the list and begin the lists there
                lastindex = issueindex[-1]
                bins = pbins[lastindex:]
                counts = pcounts[lastindex:]
            else: 
                counts = pcounts
                bins = pbins 
    return counts, bin

def checkstartingbin(start, bins): #verifies data exists for starting bin 
    if start < bins[0]:
        print("binning cannot start at specified value")
        print("binning can start at", bins[0]) 
        raise KeyboardInterrupt

def load_dataset(amp): #load and prepares the data
    superhist = "histdataset1.pkl"
    with open(superhist, "rb") as f:
        ds = pkl.load(f)
    countarr = ds[amp] #strictly the compiled im_arr of all the images 
    mbins = np.arange(26500, 26500+len(countarr)) #built into the pkl file
    mincounts = 150
    scounts = countarr > mincounts #bin indexes of sufficient counts
    pcounts = countarr[scounts]
    pbins = mbins[scounts]
    ##check that bins have no gaps; it will cause issues
    counts, bins = checkbins(pcounts, pbins)
    return counts, bins

def combined_filter_dataset(counts, bins, typ): #typ is padding, none, or positive 
    if typ == 'padding': #use the padded method
        smoothed = savgol_filter(counts, 33, 3, mode='constant', cval=0) #use the padding
        positive = smoothed >= 0 #require that the result of smoothed >= 0 
        positivesmoothed = smoothed[positive]
        positivebins = bins[positive]
        filtered, filtbins = checkbins(positivesmoothed, positivebins) #use same func as before 
        cumulativesum = np.cumsum(filtered)
        cs = scipy.interpolate.CubicSpline(filtbins+1, cumulativesum) #have right edge
        validcounts = counts[len(bins)-len(filtbins):]
    elif typ = 'None': #use no padding (will likely return an error)
        filtered = savgol_filter(counts, 33, 3) #apply filter
        cumulsum = np.cumsum(filtered) #make the cumulative sum
        cs = scipy.interpolate.CubicSpline(bins +1, cumulsum) #spline fitting for the cumulative values 
        validcounts = counts
        filtbins = bins
    elif typ == 'positive': #use positive outputs only from filter 
        smoothed = savgol_filter(counts, 33, 3) #apply filter
        #require result to be positive 
        positive = smoothed > 0
        positivesmoothed = smoothed[positive]
        positivebins = bins[positive]
        # remove any jumps in the dataset 
        filtered, filtbins = checkbins(positivesmoothed, positivebins) #use same function as before
        #now fit the spline to the cumulative sum of the filtered values greater than zero 
        cumulativesum = np.cumsum(filtered)
        cs = scipy.interpolate.CubicSpline(filtbins +1, cumulativesum)
        validcounts = counts[len(bins)-len(filtbins):] 
    
    return cs, validcounts, filtbins, filtered 

def makenarrow(bcvar, right, left):
    width = right - left
    width *= (1-bcvar) #reduce the width of the bin 
    right = width + left 
    bcvar *= 0.85 #make the change amount reduced by 85% 
    return right, bcvar
    
def makewide(bcvar, right, left): #function increases the bin width 
    width = right - left
    width *= (1+bcvar) #increase the width of the bin 
    right = width + left 
    bcvar *= 0.85 
    return right, bcvar

def makeedges(counts, bins, spline, un, start, end, redrid): #redrid whether to reduce convergence based on last bin 
    st = time.time() #time the algorithm  
    checkstartingbin(start, bins)
    edges = [start] #edges start at leftmost bin
    residualcumsum = spline(start) #the cumsum of the bins not included in the binning 
    firstcount = np.argwhere(bins == start) #start using values where specified     
    rightgoals = np.cumsum(counts[firstcount[0][0]:]) #the cumsum of the counts to be achieved; ie the right edge
    goals = rightgoals + residualcumsum #add in the ped of the spline discarded 
    lun = un #left convergence
    run = un #right convergence
    residuals =[]
    en = np.argwhere(bins == end)
    for a in range(en[0][0]):
    #for a in range(len(goals)): 
        left = edges[-1]
        if left > 102999.5: #prevent interpolation beyond what we have for bins
            break
        right = left + 1 #standard bin width = 1
        cv = spline(right) #current value, based on current right edge
        bcvar = 0.25 #bin change variable (start with 25% change)  
        x = 0
        b = 0 
        while len(edges) != a+2: 
            if goals[a]-lun < cv < goals[a] + run: #specify the requirement
                edges.append(right)
                #the residual is made this way 
                resid = goals[a] - cv
                residuals.append(resid)
                
                if redrid == True: ### Name some variable to activate this 
                    if resid > 0: #then the bin is too small 
                        lun = (un/2)
                        run = un
                    else: #the bin is too wide
                        lun = un
                        run = (un/2) 
                        
            elif cv < goals[a]: #bin is too small, widen right edge
                right, bcvar = makewide(bcvar, right, left)
                cv = spline(right)
                if x >1500: #reset bcvar 
                    if b <5: 
                        bcvar = 0.25 #gets stuck w too small, increase to help get closer
                        b +=1
                    else: 
                        print("index", a, "edge cannot be computed")
                        raise KeyboardInterrupt
                    
            elif cv > goals[a]:  #bin is too big, reduce right edge
                right, bcvar = makenarrow(bcvar, right, left)
                cv = spline(right)
                if x >1500: #reset bcvar
                    if b <5: #only let bcvar reset so many times
                        bcvar = 0.25 #gets stuck w too small, increase to help get closer
                        b+=1 
                    else: 
                        print("index", a, "edge cannot be computed")
                        raise KeyboardInterrupt          
            x+= 1 #count iterations 
    end = time.time()
    print("the computation time is", end-st)
    return edges, residuals

def plotinl(edges, start): #here binsnum is binnumbers 
    inl = getinl(edges, start)
    plt.plot(inlbins, inl, color='black', label='inl')
    plt.title('inl by bin number')
    plt.xlabel("bin number")
    plt.ylabel("inl")
    plt.legend()
    plt.grid()
    plt.show()
    return inl 

def plotdnl(edges, start):
    dnl = getdnl(edges, start)
    plt.scatter(np.arange(start, start+len(dnl)), dnl, s=5, label='dnl')
    plt.title('dnl by bin')
    plt.xlabel("bin")
    plt.ylabel("dnl")
    plt.grid()
    plt.show()
    return dnl

def getinl(edges, start): 
    inlbins = np.arange(start, start+len(edges)-1)
    idealmids = 0.5 + inlbins
    adcmids = [sum(i) for i in zip(edges[1:], edges[:-1])] #add together elements
    adcmids[:] = [x / 2 for x in adcmids] #find their mean, verified 3 entires, sufficient
    arraymids = np.array(adcmids)*(-1) 
    inl = [sum(i) for i in zip(idealmids, arraymids)]
    return inl

def getdnl(edges, start): 
    widths = edges[1:] - edges[:-1]
    dnl = np.array(widths) -1  #ideal width =1
    return dnl


#TODO efficiency:
# Can we avoid vomputing subrank for all pockets but focus only on subpockets when needed ?

import numpy as np
from siteFerret import global_module
from siteFerret.global_module import Error
from siteFerret.functions import Pdist_C,getIndex,getEntrance
import re



rmin=global_module.rmin
rmin_entrance = global_module.rmin_entrance

ML = {"IsolationForest":1,"Extended IsolationForest":2}

minsize = int(global_module.count_Threshold/10)
# minsize =10
# shiftThreshold=4
shiftThreshold = global_module.shift_Threshold


# Hard coding classes (20 Canonical aminoacids)
hClass = {}
hClass['ILE'] = 1 #Very hydrophobic (very low polarity)
hClass['VAL'] = 1 
hClass['LEU'] = 2 #Hydrophobic
hClass['ALA'] = 3 #Quite hydrophobic
hClass['CYS'] = 3
hClass['PHE'] = 3
hClass['MET'] = 3
hClass['GLY'] = 4 #~neutral
hClass['PRO'] = 5 #hydrophilic
hClass['SER'] = 5
hClass['THR'] = 5
hClass['TRP'] = 5
hClass['TYR'] = 5 
hClass['ASP'] = 6 #Very hidrophilic
hClass['GLU'] = 6
hClass['HIS'] = 6

#Alternative nomenclature for HIS
hClass['HID'] = 6
hClass['HIE'] = 6
hClass['HIP'] = 6
#########

hClass['LYS'] = 6 
hClass['ASN'] = 6
hClass['GLN'] = 6 
hClass['ARG'] = 6 

likeWater={}
likeWater[1] = -3
likeWater[2] = -2
likeWater[3] = -1
likeWater[4] = 0
likeWater[5] = 1
likeWater[6] = 2


resClass = {}
resClass['ALA'] = 1
resClass['ARG'] = 2 
resClass['ASN'] = 3
resClass['ASP'] = 4 
resClass['CYS'] = 5
resClass['GLN'] = 6
resClass['GLU'] = 7
resClass['GLY'] = 8 
resClass['HIS'] = 9
#Alternative nomenclature for HIS
resClass['HID'] = 9
resClass['HIE'] = 9
resClass['HIP'] = 9
#########
resClass['ILE'] = 10 
resClass['LEU'] = 11
resClass['LYS'] = 12 
resClass['MET'] = 13
resClass['PHE'] = 14
resClass['PRO'] = 15 
resClass['SER'] = 16
resClass['THR'] = 17
resClass['TRP'] = 18
resClass['TYR'] = 19
resClass['VAL'] = 20

def loadDict(fileName):
    import pickle
    dictList = []
    with open(fileName, 'rb') as dbfile:
        while True:
            try:
                db = pickle.load(dbfile)
            except EOFError:
                break
            dictList.append(db)
    return dictList

def concatenateRanks(a,b):
    ''' Concatenate and then like np unique but preserving order'''
    c = np.concatenate([a,b])
    ind = np.unique(c, return_index=True)[1]
    return np.array([c[index] for index in sorted(ind)])

def orderPreservedUnique(a):
    u,ind = np.unique(a,return_index = True)
    return u[np.argsort(ind)]

def weightSamples(weigth_list,original_samples,isP=False):
    #1. RENORMALIZE WEIGHTS
    VS_min = 0.2
    OS_min = 0.5
    if(isP):
        print("Large pockets weights")
        minscore = VS_min
    else:
        print("Small pockets weights")
        minscore =0.5*(VS_min + OS_min)
    
    weight_list = np.array(weigth_list)
    weight_list = (weight_list-minscore)/(1-minscore)

    print(weight_list[:100])
    #2. REWEIGHTING
    import random
    if (isinstance(original_samples, np.ndarray)):
        original_samples=original_samples.tolist()
    return_samples=original_samples.copy()
    for k,s in enumerate(original_samples):
        r = random.random()
        if ((r<=weigth_list[k])):
            # return_samples = np.vstack((return_samples, s))
            return_samples +=[s]
        else:
            pass
    
    return np.array(return_samples)
                


def getRes(resMap,atomList,getCoord=False,getAtomNumber=False):
        """
        Returns residue info of given atoms list (parsed like atom indexes)
        atomList are checked for redundancies (expected, since the spheres can be tangent to identical atoms)
        NEW: mask where the residues where within 4 AA of ligand coord, if given.
        """

        dummyChain=None
        resList = np.asarray(atomList)
        resList = np.unique(resList)
        # resMap= self.protein.resMap
        if(getAtomNumber):
            resInfo=  [resMap[r]['atomNumber'] for r in resList]
            return resInfo 
        #check if chain provided 
        if (not('resChain' in resMap[0])):
            dummyChain = 'A'
            if(getCoord):
                resInfo = [(resMap[r]['resName'],resMap[r]['resNum'],dummyChain,tuple(resMap[r]['coord']+[resMap[r]['radius']])) for r in resList]
            else:
                resInfo = [(resMap[r]['resName'],resMap[r]['resNum'],dummyChain) for r in resList]
        else:
            if(getCoord):
                resInfo = [(resMap[r]['resName'],resMap[r]['resNum'],resMap[r]['resChain'],tuple(resMap[r]['coord']+[resMap[r]['radius']])) for r in resList]
                # resInfo = [(resMap[r]['resAtom'],resMap[r]['resName'],resMap[r]['resNum'],np.append(np.asarray(resMap[r]['coord']),resMap[r]['radius'])) for r in resList]
            else:
                resInfo = [(resMap[r]['resName'],resMap[r]['resNum'],resMap[r]['resChain']) for r in resList]
        
        return resInfo

def hydroCounter(res):
    countClass = np.zeros(6)
    # n=len(f['res'])
    n=0
    content = set()
    for r in res:
        # key = (r[0],r[1])
        resname = r[0]
        resid = r[1]
        resChain = r[2]
        # print(resid,resname)
        if (resid,resname,resChain) in content:
            continue
        content.add((resid,resname,resChain))

        if resname in hClass:
            countClass[hClass[resname]-1]+=1
            # print("found", resname, pClass[resname])
        
    norm = len(content)
    # print('Normalization =',norm)
    if norm ==0:
        normCount = countClass*0 #might happen with ligand contact since the trheshold was 4 angstrom instead than 5 (Walter's idea..)
    else:
        normCount = np.round(countClass/norm *100,2)
    
    return normCount # already as a percentage..

def resCounter(res):
    countClass = np.zeros(20)
    # n=len(f['res'])
    n=0
    content = set()
    for r in res:
        # key = (r[0],r[1])
        resname = r[0]
        resid = r[1]
        resChain = r[2]
        # print(resid,resname)
        if (resid,resname,resChain) in content:
            continue
        content.add((resid,resname,resChain))

        if resname in resClass:
            countClass[resClass[resname]-1]+=1
            # print("found", resname, pClass[resname])
        
    norm = len(content)
    # print('Normalization =',norm)
    if norm ==0:
        normCount = countClass*0 #might happen with ligand contact since the trheshold was 4 angstrom instead than 5 (Walter's idea..)
    else:
        normCount = np.round(countClass/norm *100,2)
    
    return normCount # already as a percentage..

def densityHydro(res):
    D = 5  
    N = len(res)
    
    coord = np.array([r[3] for r in res])

    d,_f=Pdist_C(coord[:,0:3],coord[:,0:3])
    index = np.where(d<D)[0]

    counter = np.ones(N)

    counterHPL = np.array([1 if likeWater[hClass[r[0]]]>0 else 0 for r in res]) #Hydrophilic
    counterHP = np.array([1 if likeWater[hClass[r[0]]]<0 else 0 for r in res]) #Hydrophobic

    for k in index:
        # indexesMAP.append(getIndex(1,k,norm)) #CHECK
        resnameA = res[getIndex(1,k,N)[0]][0]
        resnameB = res[getIndex(1,k,N)[1]][0]
        if((likeWater[hClass[resnameA]]>0)and(likeWater[hClass[resnameB]]>0)):
            counterHPL[getIndex(1,k,N)[0]]+=1
            counterHPL[getIndex(1,k,N)[1]]+=1
        elif((likeWater[hClass[resnameA]]<0)and(likeWater[hClass[resnameB]]<0)):
            counterHP[getIndex(1,k,N)[0]]+=1
            counterHP[getIndex(1,k,N)[1]]+=1
       
        counter[getIndex(1,k,N)[0]]+=1
        counter[getIndex(1,k,N)[1]]+=1

    #Like this is the actual density in a sphere of D=5. Is more a surrogate of g(r) with no need to define area
    #If GLY are present is less dense but this is ok.
    # REMARK: is important in this sense to use same D as the one used to define matching scores.

    density_hydrophilic = sum(counterHPL)/sum(counter)*100 
    density_hydrophobic = sum(counterHP)/sum(counter) *100
    return density_hydrophilic,density_hydrophobic



def buildGeomFeatures(pocket,imSub=False):
    """
    Return the feature list of a single pocket.
    IMPORTANT: Respect the ordering of the training features
    """
    

    err = Error()


    # +++++++++++++ GEOMETRICAL/CLUSTERING FEATURES ++++++++++++++++++

    nn = pocket['node'].count #NORMALIZATION FACTOR
    
    volume,area,err = pocket['node'].NSVolume(accurate=True)
    if(err.value==1):
        #Trying with ARVO
        volume,area,err = pocket['node'].volume()
        print("Using ARVO to compute volume..")
        if(err.value==1):
            return err
    
    HW = np.pi**(1./3.)*(6*volume)**(2./3.)/area
    depthList = np.asarray([d for d,r,count in pocket['persistence']])
    aV = [a[1:] for a in pocket['aggregations']] #aggregation vector
    aVr = list(filter(lambda x:(max(x[2],x[3])/min(x[2],x[3]) <= 5) and min(x[2],x[3]) >=minsize, aV))#FILTERED AGGREGATIONS FOR RELEVANT EVENTS 
    bV = pocket['btlnk_info']

    isCavity = 0
    if(imSub):
        #Thus is a "subpocket" and I have to hand build this raw feature
        clusteredEntrances = [(pocket['depth'],pocket['rtop'],pocket['rtop'],1)]
        entrances = [(pocket['node'].coord,pocket['depth'],pocket['rtop'],pocket['rtop'],1)]
        n_entrances = 1
        if(clusteredEntrances[0][1]<rmin_entrance):
            # print("here")
            # input()
            # print(clusteredEntrances[0][1])
            isCavity=1
    else:
        entrances=getEntrance(pocket['large_mouths'],pwDist=True) # already contains rmin_entrance threshold..
        clusteredEntrances = [e[1:] for e in entrances]
        n_entrances = len(clusteredEntrances)
    
        if(n_entrances == 0):
            isCavity = 1 #abuse of notation, could also be a deep cleft..
            entrances=getEntrance(pocket['mouths'],pwDist=True)  
            clusteredEntrances = [e[1:] for e in entrances]
            n_entrances = len(clusteredEntrances)


    pocket['entrances'] = entrances #Save entrances in the dictionary --> useful for debugging

    # +++ Entrance scores +++
    
    av_radius_entrances = sum([r_eff for d,r_eff,r_av,cluster_size in clusteredEntrances])/n_entrances
    #I think it make sense to weight this with the entrance cluster score..
    av_entrance_depth = sum([d*cluster_size for d,r_eff,r_av,cluster_size in clusteredEntrances])/sum([cluster_size for d,r_eff,r_av,cluster_size in clusteredEntrances])
    
    # +++++++++++ geometrical scores ++++++++
    ram =  np.std(depthList)
    tr = abs(av_entrance_depth - np.average(depthList)) - ram #"Protrusion" =  ramification on superficial features
    av_entrance_depth = av_entrance_depth/shiftThreshold # 1 is the minimum

    # +++++ aggregation and clustering  Scores +++++
    score = (sum(depthList)+len(aV))/nn *100
    aggr_score = len(aV)/nn *100
    relevantAggr_score = len(aVr)/nn *100
    #COMMENT: so far I'm not using radii
    pers_score=sum(depthList)/nn *100 #how many aligned sferes over total 

    if(aVr):
        average_relevantAggrRadius = sum([(a[2]+a[3])*a[1] for a in aVr])/sum([a[2]+a[3] for a in aVr]) #weigth more important aggregations
    else:
        average_relevantAggrRadius = 0
    # ++++ bottleneck scores +++
    n_bottl = len(bV)
    if(bV):
        average_r_bottl = sum(r for r,n in bV)/n_bottl
    else:
        average_r_bottl=0
    
    n_bottl = len(bV)/nn * 100
  
    Xg=[nn,volume,HW,n_bottl,average_r_bottl,pers_score,aggr_score,relevantAggr_score,average_relevantAggrRadius,
    n_entrances,av_radius_entrances,av_entrance_depth,score,volume/nn, tr,ram, isCavity]

   
    return err,Xg

def buildChemFeaturesNoDensity(pocket,resMap):
    """
    Return the feature list of a single pocket.
    IMPORTANT: Respect the ordering of the training features
    """
    

    err = Error()

    # +++++++++++++ CHEMICAL FEATURES +++++++++++++++
    
    resData = getRes(resMap,pocket['node'].getAtomsIds())
    # res = [(r[2],r[3],r[4]) for r in resData]
    countHydroClasses = hydroCounter(resData)
    hydrophilic_score = sum(countHydroClasses[4:6])
    hydrophobic_score = sum(countHydroClasses[0:3])
    countRes = resCounter(resData)
  

    Xc=[countRes[0],countRes[1],countRes[2],countRes[3],countRes[4],countRes[5],countRes[6],countRes[7],countRes[8],
        countRes[9],countRes[10],countRes[11],countRes[12],countRes[13],countRes[14],countRes[15],countRes[16],countRes[17],countRes[18],
        countRes[19],hydrophilic_score, hydrophobic_score]
   

    return err,Xc


def buildChemFeaturesDensity(pocket,resMap):
    """
    Return the feature list of a single pocket.
    IMPORTANT: Respect the ordering of the training features
    """
    

    err = Error()

    # +++++++++++++ CHEMICAL FEATURES +++++++++++++++
    
    resData = getRes(resMap,pocket['node'].getAtomsIds(),getCoord=True)
    # res = [(r[2],r[3],r[4]) for r in resData]
    densityHydrophilic,densityHydrophobic = densityHydro(resData)
    countRes = resCounter(resData)
    

    Xc=[countRes[0],countRes[1],countRes[2],countRes[3],countRes[4],countRes[5],countRes[6],countRes[7],countRes[8],
        countRes[9],countRes[10],countRes[11],countRes[12],countRes[13],countRes[14],countRes[15],countRes[16],countRes[17],countRes[18],
        countRes[19],densityHydrophilic, densityHydrophobic]
   

    return err,Xc


def buildTestFeaturesUnique(pocket,resMap,imSub=False):
    """
    Return the feature list of a single pocket.
    IMPORTANT: Respect the ordering of the training features
    """

    err = Error()

    # +++++++++++++ CHEMICAL FEATURES +++++++++++++++
    
    resData = getRes(resMap,pocket['node'].getAtomsIds(),getCoord=True)
    # res = [(r[2],r[3],r[4]) for r in resData]
    countHydroClasses = hydroCounter(resData)
    hydrophilic_score = sum(countHydroClasses[4:6])
    hydrophobic_score = sum(countHydroClasses[0:3])
    countRes = resCounter(resData)
    # densityHydrophilic,densityHydrophobic = densityHydro(resData) 
    #could be redundant with hydroScores above (this is conceptually better though)


    # +++++++++++++ GEOMETRICAL/CLUSTERING FEATURES ++++++++++++++++++

    nn = pocket['node'].count #NORMALIZATION FACTOR
    
    volume,area,err = pocket['node'].NSVolume(accurate=True)
    if(err.value==1):
        #Trying with ARVO
        volume,area,err = pocket['node'].volume()
        print("Using ARVO to compute volume..")
        if(err.value==1):
            return err
    
    HW = np.pi**(1./3.)*(6*volume)**(2./3.)/area
    depthList = np.asarray([d for d,r,count in pocket['persistence']])
    aV = [a[1:] for a in pocket['aggregations']] #aggregation vector
    aVr = list(filter(lambda x:(max(x[2],x[3])/min(x[2],x[3]) <= 5) and min(x[2],x[3]) >=minsize, aV))#FILTERED AGGREGATIONS FOR RELEVANT EVENTS 
    bV = pocket['btlnk_info']

    isCavity = 0
    if(imSub):
        #Thus is a "subpocket" and I have to hand build this raw feature
        clusteredEntrances = [(pocket['depth'],pocket['rtop'],pocket['rtop'],1)]
        entrances = [(pocket['node'].coord,pocket['depth'],pocket['rtop'],pocket['rtop'],1)]
        n_entrances = 1
        if(clusteredEntrances[0][1]<rmin_entrance):
            # print(clusteredEntrances[0][1])
            isCavity=1
    else:
        entrances=getEntrance(pocket['large_mouths'],pwDist=True)
        clusteredEntrances = [e[1:] for e in entrances]
        n_entrances = len(clusteredEntrances)
    
        if(n_entrances == 0):
            isCavity = 1 #abuse of notation, could also be a deep cleft..
            entrances=getEntrance(pocket['mouths'],pwDist=True)  
            clusteredEntrances = [e[1:] for e in entrances]
            n_entrances = len(clusteredEntrances)


    pocket['entrances'] = entrances #Save entrances in the dictionary --> useful for debugging

    # +++ Entrance scores +++
    
    av_radius_entrances = sum([r_eff for d,r_eff,r_av,cluster_size in clusteredEntrances])/n_entrances
    av_entrance_depth = sum([d*cluster_size for d,r_eff,r_av,cluster_size in clusteredEntrances])/sum([cluster_size for d,r_eff,r_av,cluster_size in clusteredEntrances])
    
    # +++++++++++ geometrical scores ++++++++
    ram =  np.std(depthList)
    tr = abs(av_entrance_depth - np.average(depthList)) - ram #"Protrusion" =  ramification on superficial features
    av_entrance_depth = av_entrance_depth/shiftThreshold # 1 is the minimum

    # +++++ aggregation and clustering  Scores +++++
    score = (sum(depthList)+len(aV))/nn *100
    aggr_score = len(aV)/nn *100
    relevantAggr_score = len(aVr)/nn *100

    pers_score=sum(depthList)/nn *100 #how many aligned sferes over total 

    if(aVr):
        average_relevantAggrRadius = sum([(a[2]+a[3])*a[1] for a in aVr])/sum([a[2]+a[3] for a in aVr]) #weigth more important aggregations
    else:
        average_relevantAggrRadius = 0
    # ++++ bottleneck scores +++
    n_bottl = len(bV)
    if(bV):
        average_r_bottl = sum(r for r,n in bV)/n_bottl
    else:
        average_r_bottl=0
    
    n_bottl = len(bV)/nn * 100
  
    Xt=[nn,volume,HW,n_bottl,average_r_bottl,pers_score,aggr_score,relevantAggr_score,average_relevantAggrRadius,
    n_entrances,av_radius_entrances,av_entrance_depth,score,volume/nn, tr,ram, isCavity,countRes[0],countRes[1],countRes[2],countRes[3],countRes[4],countRes[5],countRes[6],countRes[7],countRes[8],
        countRes[9],countRes[10],countRes[11],countRes[12],countRes[13],countRes[14],countRes[15],countRes[16],countRes[17],countRes[18],
        countRes[19],hydrophilic_score, hydrophobic_score]#,densityHydrophilic,densityHydrophobic]
   

    return err,Xt

##################################################

class MLmodel(object):
    """
    To be flexible in the integration of other ML methods
    """
    def __init__(self, externalObj, mType = 1):
        self.mType=mType
        self.myself = externalObj
        self.ntrees = None
        self.sample = None
        if(mType == ML['IsolationForest']):
            self.ntrees = self.myself.get_params(False)["n_estimators"]
            self.sample = self.myself.get_params(False)["max_samples"]
            self.bootstrap =self.myself.get_params(False)["bootstrap"]
            self.n_features = self.myself.n_features_
            self.exlevel = 0 #by definition from EIF perspective..
        elif(mType == ML['Extended IsolationForest']):
            self.ntrees = self.myself.ntrees
            self.sample = self.myself.sample
            self.bootstrap =True #I think..
            self.exlevel = self.myself.exlevel
    def getScore(self,Xt):
        Xt = np.array(Xt)
        if(self.mType== ML['IsolationForest']):
            # STANDARD ISOLATION FOREST
            return -self.myself.score_samples(Xt)
        elif(self.mType==ML['Extended IsolationForest']):
            #EXTENDED ISOLATION FOREST
            return self.myself.compute_paths(Xt)


class Scoring(object):
    #Can be trained or loaded with an already trained model 
    def __init__(self,modelL=None,modelS=None,modelChem=None,modelChemS = None,modelP=None,modelPall=None):
        """
        Allows to load a pretrained model
        """
        self._modelL = modelL
        self._modelS = modelS
        self._modelChem= modelChem
        self._modelChemS= modelChemS

        self._modelP = modelP
        self._modelPall = modelPall

        self._Sscore=None
        self._Lscore = None    

        self._chemScore = None
        self._chemScoreS = None

        self._pScore = None
        self._pScoreS = None


    def load(self,filename_model,modelType,unique=False):
        err = Error()
        import pickle
        filename_model = global_module.trainingData+filename_model
        print("Loading ML trained model")
        print(filename_model)
        print("Checking for large pocket and small pocket interpreters")
        if(unique):
            try:
                inFile=open(filename_model+"_L.pkl",'rb')
                self._modelL= MLmodel(pickle.load(inFile),modelType)
                inFile.close()
                inFile=open(filename_model+"_S.pkl",'rb')
                self._modelS= MLmodel(pickle.load(inFile),modelType)
                inFile.close()
            except:
                err.value=2
                err.info = "Cannot load trained model(s) in "+global_module.trainingData
                return err
            if(modelType==ML["IsolationForest"]):
                # Isolation Forest
                print("Standard Isolation Forest")
                nTrees = self._modelL.ntrees
                sample = self._modelL.sample
                bootstrap = self._modelL.bootstrap
                n_features = self._modelL.n_features
                if(sample=="auto"):
                    sample =256
                print("Unique tree with geometry and chemistry")
                print("Number of trees = %d, n features = %d, sample size = %d, bootstrap = %s"%(nTrees,n_features,sample,bootstrap))
        
        else:
            
            try:
                print("Loading geometry")
                inFile=open(filename_model+"_geometryL.pkl",'rb')
                self._modelL= MLmodel(pickle.load(inFile),modelType)
                inFile.close()
                inFile=open(filename_model+"_geometryS.pkl",'rb')
                self._modelS= MLmodel(pickle.load(inFile),modelType)
                inFile.close()
                
                print("Loading chemistry")
                inFile=open(filename_model+"_chemistryL.pkl",'rb')
                self._modelChem= MLmodel(pickle.load(inFile),modelType)
                inFile.close()
                inFile=open(filename_model+"_chemistryS.pkl",'rb')
                self._modelChemALL= MLmodel(pickle.load(inFile),modelType)
                inFile.close()
            except:
                err.value=2
                err.info = "Cannot load trained model(s)"
                return err
                ### PRINT INFO ###

            if(modelType==ML["IsolationForest"]):
                # Isolation Forest
                print("Standard Isolation Forest")
                nTrees = self._modelL.ntrees
                sample = self._modelL.sample
                bootstrap = self._modelL.bootstrap
                n_features = self._modelL.n_features
                if(sample=="auto"):
                    sample =256
                print("Geometry")
                print("Number of trees = %d, n features = %d, sample size = %d, bootstrap = %s"%(nTrees,n_features,sample,bootstrap))

                nTrees = self._modelChem.ntrees
                sample = self._modelChem.sample
                bootstrap = self._modelChem.bootstrap
                n_features = self._modelChem.n_features
                if(sample=="auto"):
                    sample =256
                print("Chemistry")
                print("Number of trees = %d, n features = %d, sample size = %d, bootstrap = %s"%(nTrees,n_features,sample,bootstrap))

                ############



        if(modelType==ML["Extended IsolationForest"]):
            # Isolation Forestsample
            print("Extended Isolation Forest")
            nTrees = self._modelL.ntrees
            sample = self._modelL.sample
            exlevel = self._modelL.exlevel
            print("Number of trees = %d, sample size = %d, extension level= %d"%(nTrees,sample,exlevel))
        return err

        

    
    def _getScores(self,featListGeom=[],featListChem=[]):
        if(featListGeom):
            Lscore = self._modelL.getScore(featListGeom)
            Sscore = self._modelS.getScore(featListGeom)
            self._Lscore = Lscore
            self._Sscore=Sscore
        if(featListChem):
            chemScore = self._modelChem.getScore(featListChem)
            chemScoreS = self._modelChemALL.getScore(featListChem)
            self._chemScore = chemScore
            self._chemScoreS = chemScoreS
        return 


    def getScoresUnique(self,featList):
        LscoreU = self._modelL.getScore(featList)
        SscoreU = self._modelS.getScore(featList)
        return LscoreU,SscoreU

    def resetRank(self):
        self._Lscore = None
        self._Sscore = None
        self._chemScore = None
        self._chemScoreS = None



    def getRanking(self,featListGeom,featListChem):
        
        if self._Lscore is None:
                self._getScores(featListGeom,featListChem)
        else:
            pass
        Lscore = self._Lscore
        Sscore = self._Sscore
        
        chemScore = self._chemScore
        chemScoreS = self._chemScoreS
        avScoreL = 0.5 * (Lscore + chemScore)
        
        # Sscore = 0.5*(Lscore + Sscore) #Alternative inspired by previous analysis (globalRank) NO better S forest
        avScoreS = 0.5 * (Sscore + chemScoreS)
        
        mainRank  = np.argsort(avScoreL)
        subrank = np.argsort(avScoreS)

        return mainRank,subrank,np.sort(np.round(avScoreL,3)),np.sort(np.round(avScoreS,3))


    
    
    def getRankingOnlyGeom(self,featListGeom):

        if self._Lscore is None:
            self._getScores(featListGeom)
        else:
            pass
        Lscore = self._Lscore
        Sscore = self._Sscore
        
        # Sscore = 0.5*(Lscore + Sscore) 

        mainRank  = np.argsort(Lscore)
        subrank = np.argsort(Sscore)

        return mainRank,subrank,np.sort(np.round(Lscore,3)),np.sort(np.round(Sscore,3))
    
    
    def getRankingOnlyChem(self,featListChem):
        
        if self._chemScore is None:
            self._getScores(featListChem=featListChem)
        else:
            pass
        chemScore = self._chemScore
        chemScoreS = self._chemScoreS
        mainRank  = np.argsort(chemScore)
        subrank = np.argsort(chemScoreS)

        return mainRank,subrank,np.sort(np.round(chemScore,3)),np.sort(np.round(chemScoreS,3))
    
    

    def getRankingUnique(self,featList):
        
        LscoreUnique,SscoreUnique = self.getScoresUnique(featList)
        mainRank  = np.argsort(LscoreUnique)
        subrank = np.argsort(SscoreUnique)

        return mainRank,subrank,np.sort(np.round(LscoreUnique,3)),np.sort(np.round(SscoreUnique,3))

    
    

def getRank(pList,resMap,name,unique):
    """
    Ranking function
    """
    err = Error()

    score = Scoring()
    if(re.match(".*IF\d*",name)):
        m = ML["IsolationForest"]
    elif(re.match(".*EIF\d*",name)):
            m = ML["Extended IsolationForest"]
    else: 
        print("Model not provided or naming scheme not respected. EXIT")
        print(name)
        err.value=2
        err.info = "Model not provided or naming scheme not respected"
        
        exit()
    
    if(unique):
        err = score.load(name,modelType=m,unique=True)
    else:
        err = score.load(name,modelType=m)

    if(err.value==2):
        print("<<ERROR>> Cannot load trained model\n Aborting.")
        exit()

    #**************************************
   
    featListGeom,featListChem,map_direct,map_reverse = getFeatNoDensity(pList,resMap)
    return_featList = featListGeom
    if(unique):
        print('UNIQUE CHEM GEOM FOREST')
        rankedIndexes,rankedIndexesSub,numericalScore,numericalScoreSub = score.getRankingUnique([fg + fc for fg,fc in zip(featListGeom,featListChem)])
    else:
        rankedIndexes,rankedIndexesSub,numericalScore,numericalScoreSub = score.getRanking(featListGeom,featListChem)
        
    return (rankedIndexes,numericalScore),(rankedIndexesSub,numericalScoreSub),return_featList,map_direct,map_reverse

def getFeatNoDensity(pList,resMap):
    '''
    Output: List of list. First index features per pocket, second index spans the features
    '''
    featListGeom = [] 
    featListChem = []
    mapd={} 
    mapr={}
    s = 0
    for parentPindex,p in enumerate(pList):
        # print(parentPindex)
        err,fg= buildGeomFeatures(p)
        err,fc= buildChemFeaturesNoDensity(p,resMap)
        if(err.value==1):
            #Maybe in these scenarios assume 0 volume? (but is the same, it will cause the pocket to be out from ranking)
            print("Skipping pocket, since unable to estimate the volume")
            continue
        featListGeom.append(fg)
        featListChem.append(fc)
        mapd[s] = (parentPindex,None)
        mapr [(parentPindex,None)]=s
        s+=1
        for subIndex,sub in enumerate(p['subpockets']):
            err,fg = buildGeomFeatures(sub,imSub=True)
            err,fc= buildChemFeaturesNoDensity(sub,resMap)
            if(err.value==1):
                print("Skipping sub-pocket, since unable to estimate the volume")
                continue
            featListGeom.append(fg)
            featListChem.append(fc)
            mapd[s] = (parentPindex,subIndex)
            mapr [(parentPindex,subIndex)]=s
            s+=1
        
    return featListGeom,featListChem,mapd,mapr


def getFeatOnlyDensity(pList,resMap):
    '''
    Output: List of list. First index features per pocket, second index spans the features
    '''
    featListGeom = [] 
    featListChem = []
    mapd={} 
    mapr={}
    s = 0
    for parentPindex,p in enumerate(pList):
        # print(parentPindex)
        err,fg= buildGeomFeatures(p)
        err,fc= buildChemFeaturesDensity(p,resMap)
        if(err.value==1):
            #Maybe in these scenarios assume 0 volume? (but is the same, it will cause the pocket to be out from ranking)
            print("Skipping pocket, since unable to estimate the volume")
            continue
        featListGeom.append(fg)
        featListChem.append(fc)
        mapd[s] = (parentPindex,None)
        mapr [(parentPindex,None)]=s
        s+=1
        for subIndex,sub in enumerate(p['subpockets']):
            err,fg = buildGeomFeatures(sub,imSub=True)
            err,fc= buildChemFeaturesDensity(sub,resMap)
            if(err.value==1):
                print("Skipping sub-pocket, since unable to estimate the volume")
                continue
            featListGeom.append(fg)
            featListChem.append(fc)
            mapd[s] = (parentPindex,subIndex)
            mapr [(parentPindex,subIndex)]=s
            s+=1
        
    return featListGeom,featListChem,mapd,mapr


def getFeatUnique(pList,resMap):
    '''
    Output: List of list. First index features per pocket, second index spans the features
    '''
    featList=[]
    mapd={} 
    mapr={}
    s = 0
    for parentPindex,p in enumerate(pList):
        # print(parentPindex)
        err,f = buildTestFeaturesUnique(p,resMap)
        if(err.value==1):
            #Maybe in these scenarios assume 0 volume? (but is the same, it will cause the pocket to be out from ranking)
            print("Skipping pocket, since unable to estimate the volume")
            continue
        featList.append(f)
    
        mapd[s] = (parentPindex,None)
        mapr [(parentPindex,None)]=s
        s+=1
        for subIndex,sub in enumerate(p['subpockets']):
            err,f = buildTestFeaturesUnique(sub,resMap,imSub=True)
            if(err.value==1):
                print("Skipping sub-pocket, since unable to estimate the volume")
                continue
            featList.append(f)
            
            mapd[s] = (parentPindex,subIndex)
            mapr [(parentPindex,subIndex)]=s
            s+=1
        
    return featList,mapd,mapr


def rerankVol_withSub(originalPlist,rankingForest,scoresIF,map_direct,map_reverse,keep=10,cutoff=1):
    '''
    New ranking method which uses forest for top10 and volume among them. But among top10 subpockets are not counted as filling the ranking.
    (r index not progressing)
    '''
    r=0
    r_in=0
    selfContained = set()
    alreadyDoneSub=set()
    pListTop = []
    vols =[]
    scores=[]

    num_sub=0
    
    while ((r<10) and (r_in < rankingForest.size)and (scoresIF[rankingForest[r_in]]<=cutoff)):
        pi,si = map_direct[rankingForest[r_in]]
        
        if(si is not None):
    # ------------------- SUBPOCKET--------------------------- 
            if(pi in selfContained):
                #skip subpocket of a master pocket ahead in the ranking
                #going to next element (positions do not progress in the r
                r_in+=1
                continue # r index does not advance
            n_subs=0
            pocket = originalPlist[pi]['subpockets'][si] #IS A SUBPOCKET
            pListTop.append(pocket)
            volume,_A,_err = pocket['node'].volume()
            vols.append(volume)
            scores.append(scoresIF[map_reverse[(pi,si)]])

            if(len(originalPlist[pi]['subpockets'])==1):
                #Single subpocket--> master pocket in black list
                selfContained.add(pi)
            else:
                alreadyDoneSub.add((pi,si))
        else:
            # --------- PARENT POCKET (or with no subpockets) --------------------
            if(pi in selfContained): #SKIP PARENT POCKET OF A SINGLE SUBPOCKET ALREADY EXPRESSED IN THE RANKING
                r_in+=1
                continue
            selfContained.add(pi) #to filter out subpockets already expressed by the master pocket        
            subs = originalPlist[pi]['subpockets']
            n_subs = len(subs)
            if(n_subs==0):
                #pocket with no subpockets
                pocket = originalPlist[pi]
                pListTop.append(pocket)
                volume,_A,_err = pocket['node'].volume()
                vols.append(volume)
                scores.append(scoresIF[map_reverse[(pi,si)]])
            for i in range(n_subs):
                if((pi,i)in alreadyDoneSub):
                    # I think this scenario is very unlikely
                    pass
                else:
                    num_sub+=1
                    subpocket = subs[i]
                    pListTop.append(subpocket)
                    volume,_A,_err = subpocket['node'].volume()
                    vols.append(volume)
                    scores.append(scoresIF[map_reverse[(pi,i)]])
        r+=1
        r_in+=1
        
    print('Original number of pockets considered',len(vols))
    print('number of subpockets opened:',num_sub)
    # print(vols)
    ind = np.argsort(vols)[::-1]
    pListTop = np.array(pListTop)[ind]
    vols = np.array(vols)[ind]
    print('sorted volumes:',vols)

    return pListTop[:keep],scores


def rerankVol_noSub(originalPlist,rankingForest,scoresIF,map_direct,map_reverse,keep=10,cutoff=1):

    r=0
    r_in=0
    selfContained = set()
    pListTop = []
    vols =[]
    scores=[]
    while ((r<10) and (r_in < rankingForest.size) and (scoresIF[rankingForest[r_in]]<=cutoff)):
        pi,si = map_direct[rankingForest[r_in]]
        if(si is not None):
    # ------------------- SUBPOCKET--------------------------- 
            if(pi in selfContained):
                #skip subpocket of a master pocket ahead in the ranking
                #going to next element (positions do not progress in the r
                r_in+=1
                continue # r index does not advance
            pocket = originalPlist[pi]['subpockets'][si] #IS A SUBPOCKET
            pListTop.append(pocket)
            volume,_A,_err = pocket['node'].volume()
            vols.append(volume) 
            scores.append(scoresIF[map_reverse[(pi,si)]])
        else:
            # --------- PARENT POCKET (or with no subpockets) --------------------
            selfContained.add(pi) #to filter out subpockets already expressed by the master pocket
            pocket = originalPlist[pi]
            pListTop.append(pocket)
            volume,_A,_err = pocket['node'].volume()
            vols.append(volume)
            scores.append(scoresIF[map_reverse[(pi,si)]])
            

        r+=1
        r_in+=1

    print('Original number of pockets considered keeping only master pockets or single clusters',len(vols))
    ind = np.argsort(vols)[::-1]
    vols = np.array(vols)[ind]
    pListTop = np.array(pListTop)[ind]
    return pListTop[:keep],scores

################# LOOP.PY ###############
########### Author: Luca Gagliardi -- ConceptLAB - IIT Genova
#
# Usage: python3 loop.py.<optional:folder_name containing structures in pqr (path relative to working directory)>
# 1. By default, if the optional folder name is not given inline, the script will look for the "structures" folder
# 2. Default clustering parameters can be changed by providing the config.txt file in the working directory. A template is provided.
# Runs SiteFerret over all structures provided in the folder with default or given clustering parameters.
# Ouput: pn_atm for each structure where n is the ranking position (max n=10). pn_atm contains the protein atoms contacted (within 5 Angstrom) by the pocket.


import sys
import os
import subprocess
import numpy as np

from siteFerret import *
from siteFerret import global_module

################## INITIALIZATIONS ###############################
model="IF" 


errFile = open("errorLog.txt",'w+')
errFile.write("## ERROR LOG FILE ##\n")

runPath = global_module.exPath
if not os.path.exists(runPath):
    os.makedirs(runPath)
subprocess.call('cp '+global_module.pathTo_NS_ex+'* '+runPath+"/", shell=True)
argv = sys.argv[1:]
if argv:
    structureFolderName= argv[0]
else:
    structureFolderName = "structures"
err = initFolders(pdbName=None,structureFolder=structureFolderName)
err.handle(errFile)


logFile = open("logFile.txt", 'w')
logFile.write("******************************** LOGFILE ********************************** \n")

logFile.write(" LOOP MODE: pockets are stored for the given clustering parameters and numbering given by the Isolation Forest ranking.\n No feature vector is saved. Use getData.py to generate pocket features of succesful pockets. \n")

main_crono =Crono()
meso_crono = Crono()

t = main_crono.init()
logFile.write(t)
print(t)


#######################3

## READING AND SETTING INPUTS

config = ReadConfig()
try:
    confFile=open("config.txt",'r')
    err=config.get(confFile)
except FileNotFoundError:
    print("Config file not found, using defaut parameters")
    logFile.write("\nConfig file not found, using defaut parameters")
    err=Error()
    err.info = "Config file not found, using defaut parameters"
    err.value = 1
err.handle(errFile)
gamma = config.gamma
beta = config.beta
rp_max = config.rp_max

print("gamma = %.1f\tbeta = %.1f \tmaximum probe radius = %.1f"%(gamma,beta,rp_max))

_input = ReadInput()
structures = _input.getStructureIterator()

    
print(structures)
n_structures = len(structures)

logFile.write("\nNumber of structures : %d\n" 
%(len(structures)))
print("\n ++Number of structures : %d \n" 
%(len(structures)))
logFile.write("\nStructures analysed: " )
for s in structures:
        logFile.write("\n")
        for k, v in s.items():
            logFile.write(str(k) + ' -->'+ str(v)+'\t')
            


######## RUNNING CLUSTERING ALGOTITHM IN LOOP AND STORING RESULTS. NOTE: Feature vector is not stored
print("\n\n ** STARTING CORE ALGORITHM **\n")
logFile.write("\n\n ** STARTING CORE ALGORITHM ** \n")

from siteFerret.train_classifier import Scoring,getFeatNoDensity
from siteFerret.functions import saveResSimple


#   *********** LOAD TRAINED MODEL ****************

print("\n++++TEST MODE: assessing performance of predictor on structures +++++\n")#generic, could be training as new samples
logFile.write("\n++++TEST MODE: assessing performance of predictor on structures +++++\n")
# +++++++++++++++++++ Load trained classifier +++++++++++++++++++++

scoreIF = Scoring()



err = scoreIF.load(model,modelType=1,unique=False) #change here and below to test different classifiers..
err.handle(errFile)

#LOOP OVER STRUCTURES
skipped_structures = 0
clustering = NS_clustering()
for s in range(n_structures):
    print("\n\n ******** Current analysis: " + structures[s]['pqr'] + " ********** \n" )
    logFile.write("\n\n ******** Current analysis: " + structures[s]['pqr'] + " ********** \n")
    proteinFile = structures[s]['pqr']

    # INITIALIZATION : loading structure data...
    err = clustering.init(structure_name=proteinFile)

    if(err.value==2):
        print("A MAJOR WARNING WAS PRODUCED")
        print(err.info + "\nSkipping the structure. Cannot load " + proteinFile)
        logFile.write("\n A MAJOR WARNING WAS PRODUCED\n")
        err.info = err.info + "\nSkipping the structure. Cannot load " + proteinFile
        err.value = 1
        err.handle(errFile)
        skipped_structures+=1
        continue
    err.handle(errFile) 


    meso_crono.init()
    

    try:
        err,pList=clustering.build(gamma,beta,rpMAX = 3.0)
    except Exception:
        err.value = 2
        print("\n Critical problem 1 ")
        err.info = err.info + "An unhandled exception was produced\n "#Skipping gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max) +"of current protein="+proteinFile
        logFile.write("\n A CRITICAL WARNING WAS PRODUCED\n")
        err.handle(errFile)
    if(err.value==1):
        err.info= err.info + "\n Skipping current triplet due to failure in clusters construsction \n Skipping gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max) +"of current protein="+proteinFile
        # err.value = 1 #prevents exit, in this case we want to skip
        err.handle(errFile)
        #NOT EASY TO ACCOUNT FOR THIS WHEN AVERAGING MANY RUNS.. BETTER FULL SKIP OF THE STRUCTURE?
        continue 
    elif(err.value==2):
        print("\n Critical problem 2")
        err.info = err.info + "An unhandled exception was produced\n"#Skipping gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max) +"of current protein="+proteinFile
        logFile.write("\n A CRITICAL WARNING WAS PRODUCED\n")
        err.handle(errFile)
        break

    # ************ RANKING *********************************8
    
    featListGeom,featListChem,map_direct,map_reverse = getFeatNoDensity(pList,clustering.get_protein().resMap)
    
    # featList = [fg + fc for fg,fc in zip(featListGeom,featListChem)]
    scoreIF.resetRank()
    rankedIndexesIF,rankedIndexesSubIF,_numericalScoreIF,_numericalScoreSubIF = scoreIF.getRanking(featListGeom,featListChem)

    rankIF = (rankedIndexesIF,rankedIndexesSubIF)

    
    save_path = 'output/'+structures[s]['pqr']+'_pockets'

    rankMain = rankIF[0]
    rankSub =  rankIF[1]# For some ranking modes it coincides with the above..
    #ligand-structure loop
    r=0 # rank position
    r_in=0 # running over ranked indexes array
    selfContained = set()

    while ((r<10) and (r_in < rankMain.size)):
        pi,si = map_direct[rankMain[r_in]] #pocket and subpocket original index
        # print(r,r_in)
        if(si is not None):
            # ------------------- SUBPOCKET --------------------
            if(pi in selfContained):
                #going to next element (positions do not progress in the ranking) since sub contained in master pocket higher in ranking
                r_in+=1
                continue # r index does not advance

            pocket = pList[pi]['subpockets'][si]['node'] #IS A SUBPOCKET
                    
            if(len(pList[pi]['subpockets'])==1):
            #Single subpocket--> master pocket in black list
                selfContained.add(pi)

            if not os.path.exists(save_path):
                os.makedirs(save_path)
            saveResSimple(r+1,save_path,pocket,clustering.get_protein().resMap)

        else:
            # --------- PARENT POCKET (or with no subpockets) --------------------
            if(pi in selfContained): #SKIP PARENT POCKET OF A SINGLE SUBPOCKET ALREADY EXPRESSED IN THE RANKING
                r_in+=1
                continue
            selfContained.add(pi) #to filter out subpockets already expressed by the master pocket

            pocket = pList[pi]['node']
            
            subs = pList[pi]['subpockets']
            n_subs = len(subs)
            
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            saveResSimple(r+1,save_path,pocket,clustering.get_protein().resMap)
            
            internal_rank=[] # container for score of subpockets analysed 
            
            indS = []# just to avoid an extra "if" clause
            if(n_subs>1):
                # Save subpockets only when more than one

                for sub_i in range(n_subs):
                    internal_rank.append(np.where(rankSub==map_reverse[(pi,sub_i)])[0][0])
                    indS = np.argsort(internal_rank) #relative rank among subpockets
                for ns,sub_i in enumerate(indS):
                    subpocket = subs[sub_i]['node']
                    subpocket.load_protein(clustering.get_protein())
                    
                    if not os.path.exists(save_path+'/sub'):
                        os.makedirs(save_path+'/sub')
                    saveResSimple(r+1,save_path+'/sub',subpocket,clustering.get_protein().resMap,nsub=ns)
        r+=1
        r_in+=1
        
    advancement =(100*np.round((s+1)/n_structures,3),s+1-skipped_structures,n_structures)
    print("ADVANCEMENT: %d %% of structures"%advancement[0])
    t=meso_crono.get()
    print('\n',t)
    logFile.write(t+"\n")
    logFile.write("ADVANCEMENT: %d %% of structures"% advancement[0])
    logFile.flush()
    errFile.flush()
   


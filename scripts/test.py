
############################################ TEST.PY ##################################################
#
########### Author: Luca Gagliardi -- ConceptLAB - IIT Genova
#
# Usage: python3 test.py. 
#  1. A structure folder must be present containing all structures (pqr) and ligands (xyz). A ligand-structure map file must be present in the folder.
#  2. The input.prm file must be in the working directory. A sample is provided in the script folder.
# Pre-requisite: Build the content of the structure folder using the lfetch.py script: https://github.com/concept-lab/MOAD_ligandFinder 
#       Runs SiteFerret over all structures provided in the structures folder on all clustering parameters provided by input file
# Ouput: 1. Statistics of performance of the ranking with all parameters explored in terms of top1, top3, top10, etc...

import sys
import os
import subprocess
import numpy as np

from siteFerret import *
from siteFerret import global_module
from siteFerret.train_classifier import Scoring,getFeatNoDensity
from siteFerret.functions import save_rankingStats



# ********************** CHANGE HERE IF NEEDED*************************88

modelName = "IF" # <-- will look for <modelName>_geometryL.pkl and <modelName>_chemistryL.pkl, and <modelName>_geometryS.pkl and <modelName>_chemistryS.pkl 
                   # If the option unique is given (see scoreIF.load() ), only <modelName>_L.pkl and <modelName>_S.pkl will be used (no average between geometry and chemistry)
PATH = "./" #Insert None if using default path of the SiteFerret module (default trained model in siteFerret/trainedModels/)

# NOTE: Predisposed to use in alternative "unique" forest trained on geometrical and chemical features together 
# (distintion still kept on "large" (main ranker) and "small" (only for sub-pockets within a pocket) Isolation Forests)


excludePeptide = False
onlyPeptides = False

##################### EVALUATION METRIC'S THRESHOLDS ###
VS_threshold = 0.2
OV_threshold = 0.5


#**********************
#*************


np.seterr( invalid='ignore') #avoids warning when no match causes 0/0 division when gathering stats
amino = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HID','HIE','HIP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}




print("THRESHOLDS:\nOV=%.1f\tVS=%.1f"%(OV_threshold,VS_threshold))

######################

class ContinueI(Exception):
    pass


################## INITIALIZATIONS ###############################
##      INITIALIZING ERROR CLASS AND FILE
errFile = open("errorLog.txt",'w+')
errFile.write("## ERROR LOG FILE ##\n")

################## SET UP FOLDERS #
runPath = 'temp'
if not os.path.exists(runPath):
    os.makedirs(runPath)
subprocess.call('cp '+global_module.pathTo_NS_ex+'* '+runPath+"/", shell=True) 
err = initFolders()
err.handle(errFile)



################ READING AND SETTING INPUTS

_input = ReadInput()
try:
    inFile= open("input.prm",'r')
except FileNotFoundError as info:
    print("Input file not found. You can provide pqr filename inline for simple run.\n Aborting.")
    err = Error()
    err.info=(str(type(info))+str(info.args))
    err.value=2
    err.write(errFile)
    exit()
    
err,(gammas,betas,radii) = _input.get(inFile)
err.handle(errFile)


structures = _input.getStructureIterator()

######################## INFO ################

main_crono =Crono()
meso_crono = Crono()
t = main_crono.init()
print(t)

logFile = open("logFile.txt", 'w')
logFile.write("******************************** LOGFILE ********************************** \n")

print("\n++++ TEST MODE: storing successful pockets and gathering stats +++++\n")
logFile.write("\n+++++ TEST MODE: Gathering stats on ranking performance +++++\n")
logFile.write(t)
    
print("gammas array=",gammas)
print("betas array=",betas)
print("radii array =",radii)
print(structures)
n_structures = len(structures)

logFile.write("\nNumber of structures : %d \n Number of parameters to check x structure: %d = n_gamma=%d X n_beta=%d X n_radii=%d\n" 
%(len(structures),gammas.size * betas.size * radii.size,gammas.size,betas.size,radii.size))
print("\n ++Number of structures : %d \n Number of parameters to check x structure: %d = n_gamma=%d X n_beta=%d X n_radii=%d ++ \n" 
%(len(structures),gammas.size * betas.size * radii.size,gammas.size,betas.size,radii.size))
logFile.write("\ngammas: " )
logFile.writelines(str(gammas))
logFile.write("\nBetas: " )
logFile.writelines(str(betas))
logFile.write("\nRadii: " )
logFile.writelines(str(radii))
logFile.write("\nStructures analysed: " )
for s in structures:
        logFile.write("\n")
        for k, v in s.items():
            logFile.write(str(k) + ' -->'+ str(v)+'\t')
if(onlyPeptides):
    print("CAREFUL: only peptides are analyzed")
    input('continue?')
if(excludePeptide):
    print("CAREFUL: peptides excluded from the anaysis")
    input('continue?')

######################## MAIN ##########################################3

print("\n\n ** STARTING CORE ALGORITHM **\n")
logFile.write("\n\n ** STARTING CORE ALGORITHM ** \n")

#   *********** LOAD TRAINED MODEL ****************

print("\n++++TEST MODE: assessing performance of predictor on structures +++++\n")#generic, could be training as new samples
logFile.write("\n++++TEST MODE: assessing performance of predictor on structures +++++\n")

scoreIF = Scoring() 
# names = ["IF_geomChem","IF_onlyGeom","IF_onlyChem","VolumeScore"]
names = ["IF_geomChem"]

err = scoreIF.load(modelName,modelType=1,path=PATH,unique=False) #change here and below to test different classifiers..
err.handle(errFile)

#############################

############ SET UP COUNTERS ######################

n_models = 1 # Ranking modes..

hitTop1 = np.zeros((gammas.size,betas.size,radii.size,n_models))
hitTop3 = np.zeros((gammas.size,betas.size,radii.size,n_models))
hitTop10 = np.zeros((gammas.size,betas.size,radii.size,n_models))
# hitMatrix = np.zeros((gammas.size,betas.size,radii.size,n_models,10))

hitTopWithSub = np.zeros((gammas.size,betas.size,radii.size,n_models))
hitTopOnlySub = np.zeros((gammas.size,betas.size,radii.size,n_models))

OS_log = np.zeros((gammas.size,betas.size,radii.size,n_models))
VS_log = np.zeros((gammas.size,betas.size,radii.size,n_models))
volume_log = np.zeros((gammas.size,betas.size,radii.size,n_models))

norm = np.zeros((gammas.size,betas.size,radii.size,n_models))
singleHit = np.zeros((gammas.size,betas.size,radii.size,n_models)) #counter of how many times pocket with 0 subpockets hit

nAnalysed = 0
nPockets=np.empty((gammas.size,betas.size,radii.size))
sum_nPockets=np.zeros((gammas.size,betas.size,radii.size))
nSubs=np.zeros((gammas.size,betas.size,radii.size,n_models))

nohitMap =[]
for i in range(n_models):
    nohitMap.append({})

ff = [open("TEST_noTopHit"+modelName+".out",'w')]


#LOOP OVER STRUCTURES
skipped_structures = 0
clustering = NS_clustering()
analStructures = 0 
for s in range(n_structures):            
    proteinFile = structures[s]['pqr']
    ligands = structures[s]['ligands']
    isPeptide = []
    for ln in ligands:
        if (bool(set(ln.split('_')) & amino)):
            isPeptide.append(True)
        else:
            isPeptide.append(False)
    
    ligands = [(ln,isP) for ln,isP in zip(ligands,isPeptide)]
    ligands_coord=[]

#####       INITIALIZE CLUSTERING ALGORITHM 

    err = clustering.init(structure_name=proteinFile)
    if(err.value==2):
        print("\nSkipping the structure. Cannot load " + proteinFile)
        logFile.write("\n A MAJOR WARNING WAS PRODUCED\n")
        err.info = err.info + "\nSkipping the structure. Cannot load " + proteinFile
        err.value = 1
        err.handle(errFile)
        skipped_structures+=1
        continue
    err.handle(errFile) 

    #Prepare container with ligand coordinates
    for ligand_name in ligands:
        #NEW OPTION --> exclude Peptide
        if((ligand_name[1] == True) and excludePeptide):
            print("SKIPPING PEPTIDE "+ligand_name[0])
            logFile.write("\n SKIPPING PEPTIDE "+ligand_name[0])
            err.info = err.info + "\n" + "Skipping "+ligand_name[0]+" of structure "+proteinFile+": PEPTIDE"
            err.value = 1
            err.handle(errFile)
            continue
        elif(onlyPeptides):
            if(ligand_name[1] == True):
                print('Peptide, keeping..')
            else:
                continue


        try:
            coord = readLigands(global_module.pdbFolder_path+ligand_name[0],proteinCoord=clustering.get_protein().atoms)
            if (len(coord)==0):
                print("Skipping "+ligand_name[0]+" of structure "+proteinFile+ " : No ligands heavy atoms within 5A from protein")
                err.value = 1
                err.info = "Skipping "+ligand_name[0]+" of structure "+proteinFile + " : No ligands heavy atoms within 5A from protein!"
                err.handle(errFile)

                continue
            ligands_coord.append({"name":ligand_name[0],"coord":coord, "isPeptide":ligand_name[1]})
            
        except NameError:
            err.value = 1
            err.info = "Ligand " + ligand_name[0] +" could not be loaded. \n Skipping, it won't be considered among available comparesons. "
            err.handle(errFile)
            print("SKIPPING: Ligand  "+ligand_name[0]+ " could not be loaded")
            continue
    if ((not ligands_coord)):
        err.value = 1 
        err.info = "Skipping the structure " + proteinFile+". No valid ligands found."
        err.handle(errFile)
        if(not onlyPeptides):
            print("SKIPPING the structure: no valid ligands found")
            logFile.write("\n SKIPPING: no  valid ligands found..")
        skipped_structures+=1
        #skip current iteration: go to next structure
        continue

    
    print("\n\n ******** Current analysis: " + structures[s]['pqr'] + " ********** \n" )
    logFile.write("\n\n ******** Current analysis: " + structures[s]['pqr'] + " ********** \n")
    print("Ligands: ",ligands)
    
    # ++++++++++++++++++ LOOP OVER ALL PARAMETERS ++++++++++++++++++
    meso_crono.init()
    
    for i,gamma in enumerate(gammas):
        for j,beta in enumerate(betas):
            for k,rp_max in enumerate(radii):
                
                ######## RUNNING CLUSTERING ALGOTITHM
                
                try:
                    err,pList=clustering.build(gamma,beta,rpMAX = rp_max)
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
                # ** ANALYSE CANDIDATE POCKETS FOUND **
                
                numberPockets = len(pList)

                featListGeom,featListChem,map_direct,map_reverse = getFeatNoDensity(pList,clustering.get_protein().resMap)
                #getFeatUnique(pList,clustering.get_protein().resMap) in case unique forest is used
                
                scoreIF.resetRank()
                rankedIndexesIF,rankedIndexesSubIF,numericalScoreIF,numericalScoreSubIF = scoreIF.getRanking(featListGeom,featListChem)
                # getRankingUnique(featList) in case unique forest is used

                # rankedIndexesIF,rankedIndexesSubIF,_numericalScoreIF,_numericalScoreSubIF = scoreIF.getRankingOnlyChem(featListChem)
                # rankIF_onlyChem = (rankedIndexesIF,rankedIndexesSubIF)

                # rankedIndexesIF,rankedIndexesSubIF,_numericalScoreIF,_numericalScoreSubIF = scoreIF.getRankingOnlyGeom(featListGeom)
                # rankIF_onlyGeom = (rankedIndexesIF,rankedIndexesSubIF)

                # rankedIndexesIF,rankedIndexesSubIF,numericalScoreIF,numericalScoreSubIF = getRanking_volume(featListGeom)
                # rankIF_volume = (rankedIndexesIF,rankedIndexesSubIF) #actually identical ranking for the 2 categories..
                # score_volume = (numericalScoreIF,numericalScoreSubIF)

                rankIF_geomChem = (rankedIndexesIF,rankedIndexesSubIF)
                score_geomChem = (numericalScoreIF,numericalScoreSubIF)

                ranks = [rankIF_geomChem]
                scores =[score_geomChem]

                OS_kept = np.zeros(n_models)
                VS_kept = np.zeros(n_models)
                volume = np.zeros(n_models)
                hitNsubs=np.zeros(n_models)
                for m in range(n_models):
                    rankMain = ranks[m][0]
                    rankSub =  ranks[m][1]
                    scoreMain = scores[m][0]
                    scoreSub = scores[m][1]

                    r_in=0 # running over ranked indexes array
                    selfContained = set()
                    top10 = 0
                    hittenLigands = set() #-->keep tracks of found ligands, avoid to double count..
                    ligandMap = np.zeros(len(ligands_coord))
                    r = 0 # Number of empty places preceeding.. note--> now the ranking is not anymore completely independent from the number of ligands encountered before
                    numberEmptyPlaces = 0
                    while ((r<10) and (r_in < rankMain.size)):
                        pi,si = map_direct[rankMain[r_in]] #pocket and subpocket original index
                        gotHit = False
                        if(si is not None):
                            # ------------------- SUBPOCKET  (TREATED AS MASTER FROM USER PERSPECTIVE)--------------------
                            if(pi in selfContained):
                                #going to next element (positions do not progress in the ranking) since sub contained in master pocket higher in ranking
                                r_in+=1
                                continue # r index does not advance

                            if(len(pList[pi]['subpockets'])==1):
                            #Single subpocket--> master pocket in black list
                                selfContained.add(pi)

                            n_subs=0
                            pocket = pList[pi]['subpockets'][si]['node'] #IS A SUBPOCKET
                            pocket.load_protein(clustering.get_protein())
                            
                            for ln,l in enumerate(ligands_coord):
                                if ln in hittenLigands:
                                    continue
                                hit,OS,VS = pocket.matchScore(l['coord'],coverageTh=VS_threshold,matchScoreTh=OV_threshold)
                                if hit:
                                    hittenLigands.add(ln)
                                    gotHit = True
                                    ligandMap[ln] = 1
                                    # SUBPOCKET HIT TREATED AS SINGLE POCKET
                                    # hitMatrix[i,j,k,m,:] += hitCounter(r)

                                    hitNsubs[m] += n_subs
                                    OS_kept[m] += OS
                                    VS_kept[m] += VS
                                    vol,_A,_err = pocket.volume()
                                    volume[m]+=vol
                                    hitTop10[i,j,k,m] +=1
                                    top10  += 1
                                    score=scoreMain[r_in]
                                    
                                    singleHit[i,j,k,m] +=1 #counter of how many times pocket with 0 or 1 subpockets hit

                                    if(r<3):#within top3
                                        if(r==0):
                                            hitTop1[i,j,k,m] += 1
                                        hitTop3[i,j,k,m] += 1
                        else:
                            # --------- PARENT POCKET (or with no subpockets) --------------------

                            if(pi in selfContained): #SKIP PARENT POCKET OF A SINGLE SUBPOCKET ALREADY EXPRESSED IN THE RANKING
                                r_in+=1
                                continue
                            selfContained.add(pi) #to filter out subpockets already expressed by the master pocket
                            pocket = pList[pi]['node']
                            pocket.load_protein(clustering.get_protein())
                            
                            for ln,l in enumerate(ligands_coord):
                                if ln in hittenLigands:
                                    continue    
        
                                hit,OS,VS = pocket.matchScore(l['coord'],coverageTh=VS_threshold,matchScoreTh=OV_threshold)
                                subs = pList[pi]['subpockets']
                                n_subs = len(subs)
                                score=scoreMain[r_in]

                                if hit:
                                    hittenLigands.add(ln)
                                    gotHit = True
                                    ligandMap[ln] = 1

                                    # hitMatrix[i,j,k,m,:] += hitCounter(r)
                                    OS_kept[m] += OS
                                    VS_kept[m] += VS
                                    vol,_A,_err = pocket.volume()
                                    volume[m]+=vol

                                    hitTop10[i,j,k,m] +=1
                                    top10  += 1
                                    
                                    if(r<3):#within top
                                        if(r==0):
                                            hitTop1[i,j,k,m] += 1
                                        hitTop3[i,j,k,m] += 1
                                    if(n_subs == 0 ):
                                        singleHit[i,j,k,m] +=1 #counter of how many times pocket with 0 subpockets hit

                                # STUDY OF SUBPOCKET HITS
                                # A SINGLE SUBPOCKET IN PRACTICE WOULD BE PROMOTED TO MASTER

                                internal_rank=[] # container for score of subpockets analysed 
                                indS = []
                                for sub_i in range(n_subs):
                                    try:
                                        internal_rank.append(np.where(rankSub==map_reverse[(pi,sub_i)])[0][0]) #list with subranks of inquired pockets
                                    except KeyError:
                                        #The subpocket has been skipped for some reasons (Nan in feature vector or cannot compute V..)
                                        print("Tried to find skipped subpocket")
                                        internal_rank.append(rankSub.size) # end of ranking
                                        scoreSub = np.append(scoreSub,1)
                                indS = np.argsort(internal_rank) #relative rank among subpockets
                                # print(indS)
                                hitSub = False
                                for ns,sub_i in enumerate(indS):
                                    if(ns>2):
                                        # print('\n **Breaking: maximum number of sub to consider reached')
                                        break
                                    subpocket = subs[sub_i]['node']
                                    subpocket.load_protein(clustering.get_protein())
                                    
                                    hitSub,OS_sub,VS_sub = subpocket.matchScore(l['coord'],coverageTh=VS_threshold,matchScoreTh=OV_threshold)

                                    if (hitSub): 
                                        #ONLY TOP 3 SUBS MATTER FOR PROMOTION OF MASTER POCKET
                                        subscore = scoreSub[internal_rank[sub_i]]
                                        if(ns < 3):
                                            hitTopWithSub[i,j,k,m] += 1 # counter indicating if the one hitting also hit within top3 pocket (if more than 1)
                                        if(n_subs>1):
                                            if(not gotHit):
                                                #NOTE: only sub counted fro scenario with more than 1 subpockets. If 1 if correctly absorbed in master
                                                hitTopOnlySub[i,j,k,m] += 1
                                        # COUNTER PROMOTING SUB HIT (even if single subpocket)
                                        if(not gotHit):
                                            hittenLigands.add(ln)
                                            gotHit= True
                                            ligandMap[ln] = 1
                                            # hitMatrix[i,j,k,m,:] += hitCounter(r)
                                            hitTop10[i,j,k,m] +=1
                                            top10  += 1
                                            if(r<3):
                                                if(r==0):
                                                    hitTop1[i,j,k,m] += 1
                                                hitTop3[i,j,k,m] += 1
                                                
                                            OS_kept[m] += OS_sub
                                            VS_kept[m] += VS_sub
                                            #volume[m],_A,err = subpocket.volume() 
                                            vol,_A,_err = pocket.volume()
                                            volume[m]+=vol
                                        break # out from subs loop                                         
                            if gotHit:
                                #OF MASTER POCKET
                                hitNsubs[m] += n_subs
                        if gotHit: 
                            if(top10 ==len(ligands_coord)):
                                print('all ligands found breaking')
                                #prevents to count different pockets hitting same ligand
                                break       
                        else:
                            numberEmptyPlaces+=1   
                        r= numberEmptyPlaces 
                                
                        r_in+=1

                        
                    # OUT OF RANKED POCKETS LOOP

                    if((r==10) or  (r_in == rankMain.size) ):
                        if(top10<len(ligands_coord)):
                            # print('Missed ligands:')
                            for ln,isMissed in enumerate(ligandMap):
                                if(isMissed==0):
                                    # print(structures[s]['ligands'][ln])
                                    if((gamma,beta,rp_max)in nohitMap[m]):
                                        nohitMap[m][gamma,beta,rp_max].append(structures[s]['pqr']+'_'+structures[s]['ligands'][ln])
                                    else:
                                        nohitMap[m][gamma,beta,rp_max]=[structures[s]['pqr']+'_'+structures[s]['ligands'][ln]]
                        
                    #OUT FROM MODELS LOOP
                    
                OS_log[i,j,k,:]+=OS_kept
                VS_log[i,j,k,:]+=VS_kept
                volume_log [i,j,k,:]+=volume
                norm[i,j,k,:] +=len(ligands_coord)
                nSubs[i,j,k,:] +=hitNsubs

                #out of ligand loop    
                nPockets[i,j,k] = numberPockets
                
                
    #OUT OF ALL PARAMETER LOOPS

    analStructures += 1.
    nAnalysed += len(ligands_coord) #n structure ligand pairs 
    sum_nPockets+=nPockets    
    
    avHitTop1 = np.nan_to_num(hitTop1/norm)           
    avHitTop3 = np.nan_to_num(hitTop3/norm)
    avHitTop10 = np.nan_to_num(hitTop10/norm)
    normForsub = hitTop10-singleHit #count when subs are succesfull among all hits which are not "single" hits that is with 0 sub-pockets
    #norm - singleHit #corrNorm: when master with no sub or sigle subpocket hit

    avHitTopwithSub = np.nan_to_num(hitTopWithSub/normForsub)

    avSingleHit = singleHit/norm
    avHitTopOnlySub = np.nan_to_num(hitTopOnlySub/norm)

    average_OS_log=np.nan_to_num(OS_log/hitTop10)#average only over success
    average_VS_log=np.nan_to_num(VS_log/hitTop10)
    average_volume_log = np.nan_to_num(volume_log/hitTop10)
    average_nSubs = np.nan_to_num(nSubs/hitTop10) # average number of subpockets OF HITTING POCKETS



    # SAVING METHOD.. 
    for m in range(n_models):
        save_rankingStats(m,gammas,betas,radii,avHitTop1,avHitTop3,avHitTop10,
        avHitTopwithSub,avHitTopOnlySub,avSingleHit,  #respectively top3 sub,subhit when no masterHit,top1 sub,single pocket or sub hit, sub hit more than 3rd
        average_OS_log,average_VS_log,average_volume_log,sum_nPockets/analStructures,average_nSubs,
        nAnalysed,names,date=Crono().init())
        


    advancement =(100*np.round((s+1)/n_structures,3),s+1-skipped_structures,n_structures)
    print("ADVANCEMENT: %d %% of structures"%advancement[0])
    t=meso_crono.get()
    print('\n',t)
    logFile.write(t+"\n")
    logFile.write("ADVANCEMENT: %d %% of structures"% advancement[0])
    logFile.flush()
    errFile.flush()
#Summary of structures-Ligand not classified into top 10
for m in range(n_models):
    for gamma in gammas:
        for beta in betas:
            for rp_max in radii:
                if((gamma,beta,rp_max)in nohitMap[m]):
                    ff[m].write("\n%.2f\t%.2f\t%.2f" %(gamma,beta,rp_max))
                    ff[m].write(repr(nohitMap[m][gamma,beta,rp_max]).replace("[","\n").replace("]","").replace(", ","\n"))
    ff[m].close()

t= main_crono.get()
print("\n ---------------- FINISHED -------------- \n")
print(t)  

logFile.write("\n ----------- FINISHED ---------- \n")
logFile.write("\n\n"+t)

n_warnings = Error.n_errors

if(n_warnings>0):
    logFile.write("\n\n <INFO> "+str(Error.n_errors)+" Warnings were produced.")
err.handle(errFile)
errFile.close()
logFile.close()

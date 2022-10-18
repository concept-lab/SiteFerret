############################################ GETDATA.PY ##################################################
#
########### Author: Luca Gagliardi -- ConceptLAB - IIT Genova
#
# Usage: python3 getData.py. 
#  1. A structure folder must be present containing all structures (pqr) and ligands (xyz). A ligand-structure map file must be present in the folder.
#  2. The input.prm file must be in the working directory. A sample is provided in the script folder.
# Pre-requisite: Build the content of the structure folder using the lfetch.py script: https://github.com/concept-lab/MOAD_ligandFinder 
# Runs SiteFerret over all structures provided in the structures folder on all the clustering parameters provided by input file
# Ouput: 1. (main) feature vector --> pickle file containing features for all matching pockets found (according to the metrcis. Threshold for success can be modified)
#                                     The clustering parameters are indicated as well as the exact score obtained. Furthermore, we keep trace if the ligand is a peptide.        
#        2. Statistics over the performance for each parameter.

import os
import subprocess
import numpy as np

from siteFerret import *
from siteFerret import global_module
from siteFerret.functions import getEntrance,writeBinary,save_avHit


amino = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HID','HIE','HIP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
np.seterr( invalid='ignore') #avoids warning when no match causes 0/0 division when gathering stats


################## INITIALIZATIONS ###############################
##      INITIALIZING ERROR CLASS AND FILE
errFile = open("errorLog.txt",'w+')
errFile.write("## ERROR LOG FILE ##\n")

### SET UP FOLDERS 

runPath = global_module.exPath 
if not os.path.exists(runPath):
    os.makedirs(runPath)
    subprocess.call('cp '+global_module.pathTo_NS_ex+'* '+runPath+"/", shell=True)
err = initFolders()
err.handle(errFile)
####
####### OUTPUT LOG FILE
logFile = open("logFile.txt", 'w')
logFile.write("******************************** LOGFILE ********************************** \n")

## SETTING TIME COUNTERS 
main_crono =Crono()
meso_crono = Crono()
local_crono = Crono()
t = main_crono.init()
logFile.write(t)
print(t)


###
###########     READING INPUT FILE ################
_input = ReadInput()
try:
    inFile= open("input.prm",'r')
except FileNotFoundError as info:
    print("Input file not found. Aborting.")
    err = Error()
    err.info=(str(type(info))+str(info.args))
    err.value=2
    err.write(errFile)
    exit()
    
err,(gammas,betas,radii) = _input.get(inFile)
err.handle(errFile)
structures = _input.getStructureIterator() # Structure--Ligand(s) map
n_structures = len(structures)

######################## INFO ################

print("\n++++ANALYSIS MODE: storing successful pockets and gathering stats +++++\n")
logFile.write("\n+++++ ANALYSIS MODE: storing successful pockets and gathering stats +++++\n")

print("gamma array=",gammas)
print("beta array=",betas)
print("radius array =",radii)
logFile.write("\nNumber of structures : %d \n Number of parameters to check x structure: %d = n_gamma=%d X n_beta=%d X n_radii=%d\n" 
%(len(structures),gammas.size * betas.size * radii.size,gammas.size,betas.size,radii.size))
print("\n ++Number of structures : %d \n Number of parameters to check x structure: %d = n_gamma=%d X n_beta=%d X n_radii=%d ++ \n" 
%(len(structures),gammas.size * betas.size * radii.size,gammas.size,betas.size,radii.size))
logFile.write("\nGammas: " )
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

####################3

## ANALYSIS RESULTS CONTAINERS
# Ligand match measures
Hit=np.zeros((gammas.size,betas.size,radii.size))
norm = np.zeros((gammas.size,betas.size,radii.size))
OS_log = np.zeros((gammas.size,betas.size,radii.size))
VS_log = np.zeros((gammas.size,betas.size,radii.size))

nohitMap ={}

# Containers for stats on all pockets ans subpockets found (independently from match)

processed_structures = np.zeros((gammas.size,betas.size,radii.size)) #counts processed structure x parameter (robust to skipping)

nPockets=np.empty((gammas.size,betas.size,radii.size))
sum_nPockets=np.zeros((gammas.size,betas.size,radii.size))
sum_norm_nPockets=np.zeros((gammas.size,betas.size,radii.size))

nSubs = np.empty((gammas.size,betas.size,radii.size))
sum_nSubs=np.zeros((gammas.size,betas.size,radii.size))
sum_norm_nSubs = np.zeros((gammas.size,betas.size,radii.size))

nSpheres =np.empty((gammas.size,betas.size,radii.size))
sum_nSpheres = np.zeros((gammas.size,betas.size,radii.size))
sum_norm_nSpheres = np.zeros((gammas.size,betas.size,radii.size))

nSpheres_xP=np.empty((gammas.size,betas.size,radii.size))
sum_nSpheres_xP = np.zeros((gammas.size,betas.size,radii.size))
sum_norm_nSpheres_xP = np.zeros((gammas.size,betas.size,radii.size))

Volume = np.empty((gammas.size,betas.size,radii.size))
sum_Volume = np.zeros((gammas.size,betas.size,radii.size))
norm_Volume = np.zeros((gammas.size,betas.size,radii.size))
#######################
#####################

### ++ File containing list of failed hits per parameter ++
ff = open("ANALYSIS_failureList.txt",'w')
ff.write("#gamma\tbeta\trp_max")

# for every triplet, list of structure-ligands with no hit.

clustering = NS_clustering()
# if(not lightMode):
pstream = writeBinary("features_data.pkl")
nAnalysed = 0
skipped_structures = 0
nPeptides =0 
for s in range(n_structures):
    saveP_path = 'output/'+structures[s]['pqr']+'_Pres'
    #MAIN LOOP OVER STRUCTURES
    print("\n\n ******** Current analysis: " + structures[s]['pqr'] + " ********** \n" )
    logFile.write("\n\n ******** Current analysis: " + structures[s]['pqr'] + " ********** \n")
    proteinFile = structures[s]['pqr']
    ligands = structures[s]['ligands']
    
    # print("ligands:", ligands)
    isPeptide = []
    for ln in ligands:
        if (bool(set(ln.split('_')) & amino)):
            isPeptide.append(True)
        else:
            isPeptide.append(False)
    # print("Loading ligand(s) coordinates")

    ligands = [(ln,isP) for ln,isP in zip(ligands,isPeptide)]

    ligands_coord=[]
    localHitFile=[]

    # INITIALIZATION : loading structure data...
    err = clustering.init(structure_name=proteinFile)
    if(err.value==2):
        logFile.write("\n A MAJOR WARNING WAS PRODUCED\n")
        print("SKIPPING THE STRUCTURE: check errorLog")
        err.info = err.info + "\nSkipping the structure. Cannot load " + proteinFile
        err.value = 1
        err.handle(errFile)
        skipped_structures+=1
        continue
    err.handle(errFile) 

    for ligand_name in ligands:
        # print(ligand_name)    
        try:
            coord = readLigands(global_module.pdbFolder_path+ligand_name[0],proteinCoord=clustering.get_protein().atoms)
            if (len(coord)==0):
                print("Skipping "+ligand_name[0]+" of structure "+proteinFile+ " : No ligands heavy atoms within 5A from protein")
                err.value = 1
                err.info = "Skipping "+ligand_name[0]+" of structure "+proteinFile + " : No ligands heavy atoms within 5A from protein!"
                err.handle(errFile)

                continue
            ligands_coord.append({"name":ligand_name[0],"coord":coord, "isPeptide":ligand_name[1]})
            if(ligand_name[1]==True):
                nPeptides+=1
        except NameError:
            err.value = 1
            err.info = "Ligand " + ligand_name[0] +" could not be loaded. \n Skipping, it won't be considered among available comparesons. "
            err.handle(errFile)
            print("SKIPPING: Ligand  "+ligand_name[0]+ " could not be loaded")
            continue
    if not ligands_coord:
        err.value = 1 
        err.info = "Skipping the structure " + proteinFile+". No ligands found."
        err.handle(errFile)
        print("SKIPPING the structure: no ligands found")
        logFile.write("\n SKIPPING: no ligands found..")
        skipped_structures+=1
        #skip current iteration: go to next structure
        continue
    # print(proteinFile)


    meso_crono.init()
    for i,gamma in enumerate(gammas):
        for j,beta in enumerate(betas):
            for k,rp_max in enumerate(radii):
                # print("\nComputing: gamma= %.2f, beta= %.2f ,rp_max = %.2f " %(gamma,beta,rp_max))
                # logFile.write("\n\nComputing: gamma= %.1f, beta= %.1f ,rp_max = %.1f \n"%(gamma,beta,rp_max))
                ###            **CLUSTERING STEP**
                local_crono.init()
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
                # err.handle(errFile) 
                ## ------------------------------

                #NEW: sorting from larger to smaller ensures in the rare cases of parity that potentially "worst" VS is returned
                #(NOTE that OS wins anyway. This might only influence outcome if OS identical)
                pList=sorted(pList, key = lambda x: x['node'].count)[::-1]

                #I think that the scale does not matter for accTriang = false
                # clustering.pocketMode(grid_scale=2,grid_selfInt=3,maxProbes_selfInt=200,accTriang=False)#setup volume calculation for pocket clusters

                # EVALUATION OF LIGANDS
                # hit is maximum 1 per ligand-structure pair. Then I keep best score if more than 1 hit
                
                
                # try:
                    #PARAMETERS ITERATOR
                for ln,l in enumerate(ligands_coord):
                    #LIGAND ITERATOR
                    hit = 0
                    OS_max= 0
                    VS_kept = 0
                    
                    # --> If more than one hit per ligand, keep the best according to OS score <--
                    # try:
                    volume_t = 0
                    nspheres = 0
                    nsubs = 0
                    pCounter=0
                    for pn,p in enumerate(pList):
                        #POCKET ITERATOR
                        nspheres += p['node'].count
                        nsubs += len(p['subpockets']) #0 if lightMode on, since not built..
                    
                        volume,area,err = p['node'].NSVolume(accurate=False)

                        # print("P"+str(pn) +" Volume= "+str(volume)+"A= "+str(area))
                        if(err.value==1):
                            volume,area,err = p['node'].volume()
                            # print(volume,area)
                            if (err.value==1):
                                #I expect this never happens..
                                print("SKipping pocket"+str(pn)+": Cannot compute volume. Parameters: gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max))
                                err.info="\n"+"Structure= "+proteinFile+" Cannot compute volume of pocket"+str(pn)+"parameters: gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max)+ "\nSKIPPING pocket!\n"
                                err.handle(errFile)
                                continue

                        volume_t += volume
                        pCounter+=1 
                    
                        # if(onlySub and p['subpockets'] and (not lightMode)):
                        # else: 
                        p['node'].load_protein(clustering.get_protein())
                        # try:
                        gotHit,OS,VS = p['node'].matchScore(l['coord']) #success test performed internally. OS,VS>0 only if success 
                        if(gotHit):
                            #RECOMPUTE ACCURATE VOLUME AND AREA FOR SUCCESSFUL POCKETS
                            
                            volume,area,err = p['node'].NSVolume(accurate=True)
                            
                            #If it failed with NS and did not with ARVO, an accurate volume calculation exists, and was stored
                            #if triang=false suceeded for some misterious reason, I use Arvo here. I don't expect this happens..
                            if(err.value==1):
                                print("Error in volume computation of pocket "+str(pn)+" for parameters: gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max))
                                print("Trying with analytic approach..")
                                volume,area,err = p['node'].volume()
                                print(volume,area,err)
                    
                            # print("p"+str(pn)+" accurate area and vol:"+str(volume)+","+str(area))

                            resData = clustering.getRes(p['node'].getAtomsIds(),getCoord=True,ligandCoord=l['coord'])
                            entrances_sizeBased = getEntrance(p['mouths'],pwDist=True)
                            entrances_rBased = getEntrance(p['large_mouths'],pwDist=True)
                            quality_correlator={"name":structures[s]['pqr'],"parameters":(gamma,beta,rp_max),"score":(OS,VS),"size":p['node'].count,
                            "volume":volume,"area":area, "res":resData,"entrances_sizeBased":[e[1:] for e in entrances_sizeBased],
                            "entrances_rBased":[e[1:] for e in entrances_rBased],
                            "persistence":p['persistence'],"top":(p['mouths'][0][1],p['mouths'][0][0].r),"bottlenecks":p['btlnk_info'],
                            "aggregations":[a[1:] for a in p['aggregations']],"n_subs":len(p['subpockets']),"isub":False,"isPeptide":l["isPeptide"]}
                            pstream.storeDict(quality_correlator) #write persitence in python pickle format
                            # if((OS>OS_max)and(gotHit)):
                            if(OS>OS_max):
                                #SAVE ONLY BEST MATCH PER LIGAND FOR STATISTICS
                                hit =1
                                # print("Previous score for current ligand (0=no previous hit in the pocket list -->with smaller score)",OS_max)
                                OS_max=OS
                                VS_kept =VS
                        OS_sub_max = 0
                        for sn,sub in enumerate(p['subpockets']):
                            sub['node'].load_protein(clustering.get_protein())
                            gotHit,OS_sub,VS_sub = sub['node'].matchScore(l['coord'])
                            # if(OS>0):
                            if(gotHit):
                                #More efficient: compute volume only if needed
                                # print('succesfull subpocket')
                                s_volume,s_area,err = sub['node'].NSVolume(accurate=True)
                                # print("sub"+str(sn)+" accurate area and vol:"+str(s_volume)+","+str(s_area))
                                if(err.value==1):
                                    # if(ln==0):
                                    #     err.info= "Error in volume computation of subpocket"+str(sn)+"in p"+str(pn)+" for parameters: gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max)+":  "+err.info + "\nTrying with ARVO module computing volume of overlapping spheres"
                                    #     err.handle(errFile)
                                        # print("Error in volume computation of subpocket"+str(sn)+" in p"+str(pn)+" for parameters: gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max))
                                        # print("Trying with analytic approach..")
                                    s_volume,s_area,err = sub['node'].volume()
                                # print(volume,area)
                                    if (err.value==1):
                                        print("SKIPPING saving of structure= "+proteinFile+": Cannot compute volume of subpocket"+str(sn)+" in p"+str(pn)+" for parameters: gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max))
                                        err.info="\n SKIPPING storage :Cannot compute volume of subpocket"+str(sn)+"in p"+ str(pn) +" for parameters: gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max) + "of structure= "+proteinFile
                                        err.handle(errFile)
                                        continue
                                    # else:
                                    #     print("ok!")
                                resData = clustering.getRes(sub['node'].getAtomsIds(),getCoord=True,ligandCoord=l['coord'])
                                quality_correlator={"name":structures[s]['pqr'],"parameters":(gamma,beta,rp_max),"score":(OS_sub,VS_sub),"size":sub['node'].count,
                                "volume":s_volume,"area":s_area,"res":resData,"entrances_sizeBased":[(sub['depth'],sub['rtop'],sub['rtop'],1)],"entrances_rBased":[(sub['depth'],sub['rtop'],sub['rtop'],1)],
                                "persistence":sub['persistence'],"top":(sub['depth'],sub['rtop']),"bottlenecks":sub['btlnk_info'],
                                "aggregations":[a[1:] for a in sub['aggregations']],"n_subs":0,"isub":True,"isPeptide":l["isPeptide"]}
                                pstream.storeDict(quality_correlator) #write persitence in python pickle format
                                
                    volume_t=volume_t/pCounter              
                    # except ContinueI:
                    #     #Raise continueI to go here..
                        
                    #     print("\nPROBLEM: SKIPPING COMPARESON WITH CURRENT LIGAND")
                    #     logFile.write("\n A MAJOR WARNING WAS PRODUCED\n")
                    #     err.info = err.info + "\n Skipping current Ligand: "+ l['name']
                    #     err.value = 1 #prevents exit, in this case we want to skip
                    #     err.handle(errFile)
                    #     continue
                        
                    if(hit==0):
                        if((gamma,beta,rp_max)in nohitMap):
                            nohitMap[gamma,beta,rp_max].append(structures[s]['pqr']+'_'+structures[s]['ligands'][ln])
                        else:
                            nohitMap[gamma,beta,rp_max] = [structures[s]['pqr']+'_'+structures[s]['ligands'][ln]]  
                    

                    Hit[i,j,k] +=hit #max 1 hit X structure-ligand
                    OS_log[i,j,k]+=OS_max
                    VS_log[i,j,k]+=VS_kept
                    norm[i,j,k] +=1 #doing so I correctly count only analysed structures (ligand+protein = 1 structure)   

                    # print(norm)
                                                            
                # --------------- **Out from ligand loop** -------------------------

                # THIS IS NOT LINKED TO LIGANDS ANALYSED BUT IS A GENERAL PROPERTY OF POCKETS FOUND 
                # for each triplet of parameters x structure
                nPockets[i,j,k] = len(pList)
                nSubs[i,j,k] = nsubs
                nSpheres[i,j,k] = nspheres
                nSpheres_xP[i,j,k] = nspheres/len(pList) #n spheres x pocket
                Volume[i,j,k] =volume_t
                
                # except Exception:
                #     print("\nProblem, skipping current triplet due to failure in pocket analysis\n")
                #     logFile.write("\n A MAJOR WARNING WAS PRODUCED\n")
                #     err.info = err.info + "\n Skipping current triplet due to failure in pocket analysis \n Skipping gamma= %.2f, beta= %.2f ,rp_max = %.2f," %(gamma,beta,rp_max) +"of current protein="+proteinFile
                #     err.value = 1 #prevents exit, in this case we want to skip
                #     err.handle(errFile)
                #     #Skip current parameters triplet
                #     continue

                processed_structures[i,j,k]+=1



    #END OF 1st NESTED CYCLE: all parameters over given structure
    print("\n\n ** -- Done -- *** ")
    nAnalysed += len(ligands_coord) #n structure ligand pairs
    
    advancement =(100*np.round((s+1)/n_structures,3),s+1-skipped_structures,n_structures)
    print("ADVANCEMENT: %d %% of structures"%advancement[0])
    

    ########### current averages
    avHit = np.nan_to_num(Hit/norm)
    average_OS_log=np.nan_to_num(OS_log/Hit)#average only over success
    average_VS_log=np.nan_to_num(VS_log/Hit)
    ###
    sum_nPockets +=nPockets
    nPockets = nPockets/np.amin(nPockets)#number of pockets for each parameter divided by minimum number of pockets found in the structure for each param
    sum_norm_nPockets +=nPockets

    sum_nSubs +=nSubs #number of subpockets WITHOUT counting parent pocket (a simple pocket with no sub - pocket is not counted here)
    nSubs = np.nan_to_num(nSubs/np.amin(nSubs)) 
    sum_norm_nSubs +=nSubs

    sum_nSpheres += nSpheres #total probe spheres 
    nSpheres = nSpheres/np.amin(nSpheres) 
    sum_norm_nSpheres += nSpheres

    sum_nSpheres_xP += nSpheres_xP #probe spheres x pocket for all parameters over single structure
    nSpheres_xP =  np.nan_to_num(nSpheres_xP/np.amin(nSpheres_xP))
    sum_norm_nSpheres_xP += nSpheres_xP

    sum_Volume +=Volume
    Volume = np.nan_to_num(Volume/np.amin(Volume))
    norm_Volume +=Volume 
    ###


    # processed_structures = s+1-skipped_structures
    save_avHit(gammas,betas,radii,average_OS_log,average_VS_log,avHit,
    sum_nPockets/(processed_structures),sum_norm_nPockets/(processed_structures),sum_nSubs/(processed_structures),
    sum_norm_nSubs/(processed_structures),sum_Volume/processed_structures,norm_Volume/processed_structures,
    sum_nSpheres_xP/(processed_structures),sum_norm_nSpheres_xP/(processed_structures),sum_nSpheres/(processed_structures),sum_norm_nSpheres/(processed_structures),
    advancement,nAnalysed,date=Crono().init())

    # FLUSH BUFFER 
    t=meso_crono.get()
    print(t)
    logFile.write(t+"\n")
    logFile.write("ADVANCEMENT: %d %% of structures"% advancement[0])
    logFile.flush()
    errFile.flush()

    
############## END OF MAIN CYCLE ######################
# ######### write structures which failed for given parameters
for gamma in gammas:
        for beta in betas:
            for rp_max in radii:
                if((gamma,beta,rp_max)in nohitMap):
                    ff.write("\n%.2f\t%.2f\t%.2f" %(gamma,beta,rp_max))
                    ff.write(repr(nohitMap[gamma,beta,rp_max]).replace("[","\n").replace("]","").replace(", ","\n"))
ff.close()
##################

pstream.end()

err.handle(errFile)
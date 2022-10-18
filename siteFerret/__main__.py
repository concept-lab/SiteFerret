#########################           SITEFERRET  ###########################################
############ Clustering of SES virtual probes for pocket generation and ranking via Isolation Forest  ############
#Author:								Luca Gagliardi 
#Copiright:					    © 2022 Istituto Italiano di Tecnologia   
# Original Raw project: PICKPOCKET: https://github.com/lucagl/pickPocket 
############################################################### 

from siteFerret import *
from siteFerret import global_module
import sys
import os
import subprocess
import re

import getopt #to manage user options..


default_keep = 10
default_rankingThreshold = 0.56
argv = sys.argv[1:]

amino = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HID','HIE','HIP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}


rankingThreshold = default_rankingThreshold
keep = default_keep
try:
    opts, args = getopt.getopt(argv,"h",["range=","maxScore=","help"])
except getopt.GetoptError:
    print ('<<ERROR>> Uncorrect formatting of options or unavaible option')
    sys.exit(2)

for opt, arg in opts:
    if opt in ["-h","--help"]:
        print('Usage:\npython3 -m siteFerret <options> <structureName>\n The structure must be in pqr format.\nOptions:')
        print('--range=<max number of pockets in ranking>: Default = %d'%default_keep)
        print('--maxScore=<(0 ,1]>: Default =%.2f, a smaller maximum score reduces the number of returned pockets'%default_rankingThreshold)
        input('\n')
        sys.exit()
    elif opt == '--range':
        keep = int(arg)
        print("OPTION: max ranked pockets= %d"%keep)
    elif opt == '--maxScore':
        rankingThreshold = float(arg)
        print("OPTION: max score Isolation Forest= %.2f"%rankingThreshold)

if not args:
    exit("Please provide structure name")
else:
    pass
# INITIALIZATIONS #
### Create temp folder (if not existing)###
runPath = 'temp'
if not os.path.exists(runPath):
    os.makedirs(runPath)

subprocess.call('cp '+global_module.pathTo_NS_ex+'* '+runPath+"/", shell=True)

inputFile  = args[0] 
match = re.match('([\w]*)',inputFile) #accepts both <pqrname> and <pqrname.pqr>
inputFile = match.group(1)

proteinFile=inputFile

errFile = open("errorLog.txt",'w+')
errFile.write("## ERROR LOG FILE ##\n")

logFile = open("logFile.txt", 'w')
logFile.write("******************************** LOGFILE ********************************** \n")

print()
main_crono =Crono()

t = main_crono.init()
logFile.write(t)
print(t)


## READING AND SETTING INPUTS

isAnalysis=False
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
err=initFolders(inputFile)
err.handle(errFile)


logFile.write("Protein structure name= "+proteinFile+"\n")
logFile.write("Clustering parameters are:\n gamma= %.1f beta= %.1f \n Maximum probe radius= %.1f, minimum probe radius= %.1f, probe increment= %.2f" 
% (gamma,beta,rp_max,global_module.R_WATER,global_module.delta_rp))
print("Clustering parameters are:\n gamma= %.1f beta= %.1f \n Maximum probe radius= %.1f, minimum probe radius= %.1f, probe increment= %.2f\n" 
% (gamma,beta,rp_max,global_module.R_WATER,global_module.delta_rp))
print("Protein structure name= "+proteinFile+"\n")


####################### MAIN #################################################

clustering = NS_clustering()
err = clustering.init(proteinFile)
err.handle(errFile)
err,pList=clustering.build(gamma,beta,rpMAX = rp_max)
err.handle(errFile)

if(pList):
    err,pList_ranked = clustering.printInfo(saveSpheres=True,keep=keep,rankingThreshold = rankingThreshold)
    logFile.write("\n Detailed info in protein output file("+global_module.pdbFolder_path+") and status.txt file("+global_module.runFolder_path)
    err.handle(errFile)
else:
    print("NO pockets found with current parameters!")
###############

############# BUILD TRIANGULATION OF STRUCTURE
if(global_module.accTriang):
    clustering.VMD_accTriang()

logFile.write("\n\n Number of (large) pockets found: "+str(len(pList)) + " of which %d were kept in ranking"%keep)


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

##################

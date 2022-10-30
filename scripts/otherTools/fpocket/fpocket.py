import subprocess
import os
import numpy as np
import re
import glob
from  C_functs import Pdist_C,getIndex
from time import time, strftime, localtime
import sys
np.seterr( invalid='ignore')

excludePeptide = False
onlyPeptides = False

VS_threshold = 0.2
OS_threshold = 0.5

print("Volume score threshold: ", VS_threshold)
print("Overlap score threshold: ", OS_threshold)

runFolder = os.path.abspath('structures')+'/'
readFolder = os.path.abspath('structures')+'/output/'
#readFolder = runFolder
outFolder = 'results/'
if os.path.isdir(outFolder):
    pass
else:
    print('creating output directory')
    subprocess.run(['mkdir',outFolder])

# NOTES:

#  - Important: For maximum coherency we use PQR file to compare vicinity of protein coordinates to ligand coordinates: 
#      Doing so we exactly compare to the ligand coordinates used by our clustering algorithm

#  - Numbering of pocket_atm files follow the ranking of Fpocket



amino = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HID','HIE','HIP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}

def getCM(verts,weights):
    print(verts)
    print(weights)

    CM = np.average(verts,weights=weights,axis=0)

    return CM
#########

def CMhit(pvert,ligand_coords,distanceTh=4):
    verts = pvert[0]
    weights = pvert[1]
    CM =getCM(verts,weights)
    d,_flag = Pdist_C(ligand_coords,CM.reshape((-1,1)).T)
    # print(d)
    d = np.round(d)
    # print(d)

    success = np.any(d<=distanceTh)

    return success



def readPocket_verts(pname):
    comment =['#', 'CRYST[0-9]?','HEADER','HETATM']
    remark = ['REMARK']
    termination = ['TER', 'END', '\n']
    skip = comment+remark+termination
    skip = '(?:% s)' % '|'.join(skip)

    header = slice(0,4)
    number = slice(7,11)
    pol = slice(13,16) #polarity
    resname = slice(17,20) #dummy resname
    ranking = slice(25,26)
    coordIndx = slice(32,38)
    coordIndy = slice(40,46)
    coordIndz = slice(48,54)
    dummyCharge = slice(58,62)
    radius = slice(67,71)
    coord = []
    weights =[]
    try:
        inFile = open(pname,'r')
        # print(pname+'_atm.pdb')
    except Exception:
        raise NameError("Cannot load PDB file")
    for line in inFile:
        if(re.match(skip,line)): 
            continue
        if (line[header]=='ATOM'):
            coord.append([float(line[coordIndx]),float(line[coordIndy]),float(line[coordIndz])])
            weights.append(float(line[radius]))
        else:
            pass
    
    coord = np.array(coord)

    return coord,weights

def readPocket(pname):

    comment =['#', 'CRYST[0-9]?','HEADER','HETATM'] #IMPORTANT: the need to include HETATM shows weakness of fpocket, 
    # which, being related to pdb, does not filter out hetatm.. --> thus is performing voronoi on everything!!
    remark = ['REMARK']
    termination = ['TER', 'END', '\n']
    skip = comment+remark+termination
    skip = '(?:% s)' % '|'.join(skip)

    # resInd = 5
    # chargeInd = 
    nameInd = slice(17,20) #Resnmae
    atomInd = slice(12,16) #Atom name
    coordIndx = slice(30,38)
    coordIndy = slice(38,46)
    coordIndz = slice(46,54)
# ATOM   1660   CG LYS A 209      30.153  40.690  26.319  1.00 16.04           C 0
    resMap =[]
    try:
        inFile = open(pname,'r')
        # print(pname+'_atm.pdb')
    except Exception:
        raise NameError("Cannot load PDB file")
    for line in inFile:
        if(re.match(skip,line)): 
            continue
        # print(line) 
        if(line[0:4]=='ATOM'):
            content = {'resName':line[nameInd],'resAtom': line[atomInd],
                    'coord':[float(line[coordIndx]),float(line[coordIndy]),float(line[coordIndz])]}
            resMap.append(content)

    coords = np.array([p['coord'] for p in resMap])

    return coords

def readLigand(lname,proteinCoord=[]):
    '''
    Return ligand atoms closer than 5\AA to protein atoms
    '''
    inFile = open(lname+".xyz",'r')
    ligand_coord = np.loadtxt(inFile)
    if len(proteinCoord)>0:
        d,flag = Pdist_C(ligand_coord[:,0:3],proteinCoord[:,0:3])
        index = np.where(d<=5)[0]
        lindex,_pindex=getIndex(flag,index,proteinCoord.shape[0])
        # print(lindex)
        # print(np.unique(lindex))
        # print(ligand_coord[np.unique(lindex)])
        ligand_coord = ligand_coord[np.unique(lindex)]
    else:
        pass

    return ligand_coord


def get_volume(lines,pnumber):
    
    volume = 0
    s=0
    while s < len(lines):
        line = lines[s]
        if re.match("(Pocket %d)"%pnumber,line):
            # print(line)
            volume = float(lines[s + 7].split()[2])
            # print(lines[s + 7])
            # print(volume)
            break
        s+=1
    if (volume==0):
        raise Exception("Cannot find volume")
    return volume

def get_proteinAtoms(structure):
    '''
    Use PQR to filter out undesired atoms such as H or in any case keep atoms meaningful to the SES
    '''
    try:
        # print(structure+'.pdb')
        inFile = open(structure+'.pqr','r')
    except Exception:
        raise NameError("Cannot load PQR file")
    # try:
    #     # print(structure+'.pdb')
    #     _check = open(structure+'.pdb','r')
    # except Exception:
    #     raise NameError("Cannot load PDB file")
    comment =['#', 'CRYST[0-9]?']
    remark = ['REMARK']
    termination = ['TER', 'END', '\n']
    skip = comment+remark+termination
    skip = '(?:% s)' % '|'.join(skip)
    for line in inFile: 
        if(re.match(skip,line)): 
            pass 
        else:
            linegNOChain=re.match("(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)",line)
            linegChain = re.match("(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+([\w0-9]+)\s*(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)",line)
            break

    if(linegChain):
        # print("PQR contains CHAIN_ID")
        isChainID=1                                                        #resID
        matchPattern = "(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+([\w0-9]+)\s*(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)"
    elif(linegNOChain):
        # print("PQR does NOT contain CHAIN_ID")
        isChainID =0                                        # resID
        matchPattern = "(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)"
    else:
        raise NameError("Incorrect pqr file formatting")
    if(isChainID):
        resInd = 5
        chargeInd = 9
        rInd = 10
    else:
        resInd = 4
        chargeInd = 8
        rInd = 9
    nameInd = 3
    atomInd = 2
    coordInd = resInd +1
    
    inFile.seek(0)

    resMap=[]
    for line in inFile:
        if(re.match(skip,line)): 
            pass 
        else: 
            mline=re.match(matchPattern,line).groups()
            content = {'resName':mline[nameInd],'resNum':mline[resInd],'resAtom': mline[atomInd],
                    'charge':mline[chargeInd],'coord':list(map(float, mline[coordInd:coordInd+3])),'radius':float(mline[rInd])}
            resMap.append(content)

    
    coords = np.array([p['coord'] for p in resMap])
    return coords

def getStructureIterator(mapFile):
        
    comment = ['^#','\n']
    comment = '(?:% s)' % '|'.join(comment)
    # mapFile = 'ligandMap.txt'
    structures=[]
    try:
        inFile = open(mapFile,'r')
    except Exception:
        raise NameError("Cannot load mapFile")

    content = inFile.readlines()

    structures = []
    s=0
    while s <(len(content)):
        current_line=content[s]
        if(re.match(comment,current_line)):
            s+=1
            continue
        current_line = current_line.split()
        n_ligands = int(current_line[0])
        name = current_line[1]
        ligand_names=[]
        for i in range(1,n_ligands+1):
            following_line=content[s+i].split()[0]
            if(re.match(comment,following_line)):
                #skipping ligand
                continue
            ligand_names.append(following_line)
        structures.append({'pdb':name,'ligands': ligand_names})
        s+=n_ligands +1
            
    return structures


def matchScore(pocket_atm,ligand_coordinates,distanceTh=5, matchScoreTh=0.5, coverageTh=0.2):
     
    n_ligand = ligand_coordinates.shape[0]
    n_pocket = pocket_atm.shape[0]

    # print("Number of atoms ligand=%d\n Number of atoms pocket=%d" %(n_ligand,n_pocket))#CHECK

    d = np.empty(n_ligand*n_pocket)

    d,_flag = Pdist_C(ligand_coordinates,pocket_atm)
    
    rowmap = np.ones(n_pocket,bool)*False
    n_inLig=0
    for i in range(n_ligand):
        rw = d[i*n_pocket:(i+1)*n_pocket]<=distanceTh #= row: distance relations between atom ligand and all pocket atoms"
        rowmap = np.logical_or(rw,rowmap)
        n_inLig+=np.any(rw) #at least one per row = one hit
    n_inP= np.sum(rowmap)#can be seen as a mask where a true entry is a hit, a pocket atom within the distance threshold


    matchScore = n_inLig/n_ligand
    coverageScore = n_inP/n_pocket
    if((np.round(matchScore,2)>=matchScoreTh)and(np.round(coverageScore,2)>=coverageTh)):
        hit = True
    else:
        hit =False
    
    return hit,matchScore,coverageScore

def saveScores(top1,top3,top10,allHit,OS,VS,volume,nPockets,nAnalysed,progress):
    date = strftime("LOCAL TIME = %Y-%m-%d %H:%M:%S", localtime())
    of = open("rankingScore_fpocket.out","w")
    of.write("#"+date)
    of.write("\n# %d / %d Structures succesfully processed"%(progress[0],progress[1]))
    of.write("\n# Number of ligand-structure couples treated = %d"%nAnalysed)
    of.write("\n#Top1\tTop3\tTop10\tallHit\tOS\tVS\tAv volume\tnPockets")
    of.write("\n%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\t\t%.2f" 
    %(100*np.round(top1,4),100*np.round(top3,4),100*np.round(top10,4),100*np.round(allHit,4),100*np.round(OS,4),100*np.round(VS,4),np.round(volume,3),np.round(nPockets,2)))
    of.close()


def saveScoresCM(top1,top3,top10,volume,nAnalysed,progress):
    date = strftime("LOCAL TIME = %Y-%m-%d %H:%M:%S", localtime())
    of = open("CMrankingScore_fpocket.out","w")
    of.write("#"+date)
    of.write("\n# %d / %d Structures processed"%(progress[0],progress[1]))
    of.write("\n# Number of ligand-structure couples treated = %d"%nAnalysed)
    of.write("\n# Top1\tTop3\tTop10\tAv volume")
    of.write("\n%.2f\t%.2f\t%.2f\t%.3f" %(100*np.round(top1,4),100*np.round(top3,4),100*np.round(top10,4),np.round(volume,3)))
    of.close()

# ===========================================================
# ===========================


def main():
    
    if(onlyPeptides):
        print("CAREFUL: only peptides are analyzed")
        input('continue?')

    structures = getStructureIterator('ligandMap.txt')
    print(structures)
    n_structures = len(structures)
    print("N structures to analyse:", n_structures)
    input()


    ############ Main loop ################
    norm=0
    # check_norm = 0
    hitTop1 = 0
    hitTop3=0
    hitTop10=0
    allHit = 0
    OS_log=0
    VS_log=0
    volume_log=0
    failures=0
    failureList = []
    skipped = []
    nPockets = 0 # keep track of number of generated pockets
    counter=0
    for s in range(n_structures):
    
        proteinName = structures[s]['pdb']
        ligands_name = structures[s]['ligands']

        isPeptide = []
        for ln in ligands_name:
            if (bool(set(ln.split('_')) & amino)):
                isPeptide.append(True)
            else:
                isPeptide.append(False)
        # print("Loading ligand(s) coordinates")
        ligands = [(ln,isP) for ln,isP in zip(ligands_name,isPeptide)]


        try:
            proteinAtoms = get_proteinAtoms(readFolder+proteinName)
        except:
            # print("SKIPPING\n")
            print("FAILURE IN READING")
            sys.exit()
            failures+=1
            # continue

        # Preparing ligand coordindates
        ligands_coord=[]
        for ln in ligands:
            if((ln[1] == True) and excludePeptide):
                print("SKIPPING PEPTIDE "+ln[0])
                continue
            elif (onlyPeptides):
                if((ln[1] == True)):
                    print('Peptide, keeping..')
                else:
                    continue

            coord = readLigand(readFolder+ln[0],proteinCoord=proteinAtoms)

            if (len(coord)==0):
                print("Skipping "+ln[0]+" of structure "+proteinName+ " : No ligands heavy atoms within 5A from protein")
                continue
            ligands_coord.append({"name":ln[0],"coord":coord})
        if not ligands_coord:
            #skip current iteration: go to next structure
            if(not onlyPeptides):
                print("SKIPPING the structure: no ligands found")
                skipped.append(proteinName)
                failures+=1
            continue
        # Calling fpocket
        print(proteinName)
        print(ligands)
        counter +=1
        try:
            out = subprocess.check_output(['fpocket', '-f',proteinName+'.pdb'],cwd=runFolder)
        except subprocess.CalledProcessError as grepexc:                                                                                                   
            print ("error code", grepexc.returncode, grepexc.output)
            exit()
        # mv folder in results folder:
        subprocess.run(['cp','-r',runFolder+proteinName+'_out',outFolder])
        subprocess.run(['rm','-r',runFolder+proteinName+'_out'])
        # Gathering top 10 pockets
        # need to list all pockets contained in folder
        pList=[n for n in glob.glob(outFolder+proteinName+'_out/'+'pockets/*_atm.pdb')] #not ranked
        pList = sorted(pList,key=lambda x: int(re.sub('\D', '', x)))#ranked sorting

        infoFileName = outFolder+proteinName+'_out/'+proteinName+"_info.txt"
        infoFile = open(infoFileName,'r')
        infoLines = infoFile.readlines()
        
        
        top10=0
        
        volume =0 
        OS_kept = 0
        VS_kept =0
        numberEmptyPlaces = 0
        # for r,pname in enumerate(pList):
        ind=0
        r=0
        nPockets += len(pList)
        print('number of generated pockets = ',len(pList))
        ligandMap = np.zeros(len(ligands_coord))
        print('number ligands to check',len(ligandMap))
        hittenLigands=set()
        while (ind<len(pList)):
            pname = pList[ind]
            p_atoms=readPocket(pname)
            # print('empty places',r)
            # print('internal index', ind)
            gotHit = False 
            for k,l in enumerate(ligands_coord):
                if k in hittenLigands:
                    #already found ligand
                    print("skipping already hit ligand")
                    continue
                lc = l['coord']
                hit,OS,VS = matchScore(p_atoms,lc,coverageTh=VS_threshold,matchScoreTh=OS_threshold)
                if(hit):
                    hittenLigands.add(k)
                    gotHit=True
                    print(l['name'])
                    print('rank '+str(r)+" got hit")
                    print(OS,VS)
                    OS_kept += OS
                    VS_kept += VS
                    volume = +get_volume(infoLines,pnumber = r+1)
                    print('volume= ',volume)
                    allHit +=1
                    if(r<9):
                        top10+=1
                        hitTop10+=1
                        ligandMap[k] = 1
                        if(r<3):
                            hitTop3 += 1
                            if(r<1):
                                hitTop1 += 1
                    #break # out of ligands loop
            if gotHit:
                #one of the ligands is hit by the pocket
                if(top10==len(ligands_coord)):
                    print('breaking check, all ligand matched ')
                    #this avoid doule counting if more than 1 pocket matches (which happens for non severe scores)
                    break
            else:
                numberEmptyPlaces+=1
            r = numberEmptyPlaces # current ranking
            ind +=1 # actual index on pockets

        if(top10<len(ligands_coord)):
            for k,i in enumerate(ligandMap):
                # print(s,i)
                if(i==0):
                    #noHitMap.append([pqrName+': '+str(ligands[ln])])
                    failureList.append(proteinName+': '+ligands_coord[k]['name'])
                    print('NO TOP 10 MATCH FOUND FOR '+ligands_coord[k]['name'])

        norm+=len(ligands_coord)
        
       
        # print('hitTOP10 no norm=', hitTop10)
        OS_log += OS_kept
        VS_log += VS_kept
        volume_log += volume
        
    try:
        OS_log_av = OS_log/hitTop10 #average only over success
        VS_log_av = VS_log/hitTop10
        volume_av = volume_log/hitTop10
    except ZeroDivisionError:
        if(hitTop10==0):
            OS_log_av =0
            VS_log_av =0
            volume_av=0
        else:
            print("Problem!\n")
            sys.exit()
    # check_norm+=len(ligands_coord)
    print("n_strucutres:",counter)
    print("NORM (ligand-structures):",norm)

    hitTop1_av = hitTop1/norm
    hitTop3_av = hitTop3/norm
    hitTop10_av = hitTop10/norm
    av_allHit = allHit/norm

    saveScores(hitTop1_av,hitTop3_av,hitTop10_av,av_allHit,OS_log_av,VS_log_av,volume_av,nPockets/n_structures,nAnalysed=norm,progress=[s+1-failures,counter])
    

    fp = open("failure_list.txt",'w')
    fp.write(repr(skipped).replace("[","\n").replace("]","").replace(", ","\n"))
    fp.close()

    fp = open("noHit_list.txt",'w')
    # fp.write(repr(failureList).replace("[","\n").replace("]","").replace(", ","\n"))
    # fp.close()

    fp.write("#Total failures: %d/%d\n"%(len(failureList),norm))
    fp.write("#List of missed ligands\n")
    fp.write('#Thresholds: LC_Th=%.1f\t PC_Th=%.1f\n'%(np.round(OS_log_av,3)*100,np.round(VS_log_av,3)*100))
    for string in failureList:
        fp.write(repr(string).replace("[","\n").replace("]",""))
        fp.write('\n')
    fp.close()
    # print("Failures =", failures)
    print('number of not top 10 matching pockets:',len(failureList))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\nUser exit")
        sys.exit()


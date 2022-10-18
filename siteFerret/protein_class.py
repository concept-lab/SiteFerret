############### PROTEIN CLASS (contains loaded structure atoms and related computations) ####################

from siteFerret.functions import fetchRes
from siteFerret.functions import area,coplanar,compare_ids
from siteFerret.functions import Pdist_C,pcross,psum
import numpy as np



class Protein(object):
    def __init__(self):
        self.atoms = np.empty((0,4))
        self.resMap = None #Contains mapping between residue name and coordinates

    def load(self,filename,isPQR=True):
        # print("-- Loading atom coordinates --")
        if(isPQR):
            # print("here")
            self.resMap = fetchRes(filename) #also produces NS input file by default
            # global atoms
            for i in self.resMap:
                c = np.append(np.asarray(i['coord']),i['radius'])
                self.atoms = np.vstack([self.atoms,c])
        else:
            try:
                self.atoms = np.loadtxt(filename)
                # print("N atoms =",self.atoms.shape[0] )
            except:
                raise Exception
        # print("\n-- DONE --")
        # print("Number atoms: ",len(self.resMap))
        return 

    ############# MAIN FUNCTIONS ###################

    def matchCriterion(self,ligand_coordinates,pocket_atomList,distanceTh):
        """
        Ligand matching criterion based on the number of ligand atoms within a certain distance from any PROTEIN atom of the pocket.
        Returns two scores: 
        - matchScore : How many heavy atoms of the ligand are within given distance from a pocket's protein atom 
        divided by the total heavy atoms of the ligand
        - volumeScore : The ratio of these protein atoms with respect to all protein atoms of the pocket. 

        INPUT: heavy atoms of the ligand and atoms of the protein "facing" the pocket. 
        OPTIONARY: Distance criterion different from default
        """
        ligand_coordinates = np.array(ligand_coordinates)[:,0:3]

        resList = np.unique(np.array(pocket_atomList)) #unique list of tangent atoms to the pocket (1D, the triplets get flattened)
        pocket_tAtoms = self.atoms[resList,0:3]

        # print("\nTangent atoms:\n",pocket_tAtoms)
        # print("\nLigand atoms:\n",ligand_coordinates)

        n_ligand = ligand_coordinates.shape[0]
        n_pocket = pocket_tAtoms.shape[0]

        # print("Number of atoms ligand=%d\n Number of atoms pocket=%d" %(n_ligand,n_pocket))#CHECK

        d = np.empty(n_ligand*n_pocket)

        d,flag = Pdist_C(ligand_coordinates,pocket_tAtoms)
        if(flag!=0):
            raise ValueError("Something unexpected happened: ligand and pocket must be different")

        # print("Number of pair distances=",d.size) #CHECK
        # n_in=np.sum(d<=distanceTh)
        rowmap = np.ones(n_pocket,bool)*False
        n_inLig=0
        for i in range(n_ligand):
            rw = d[i*n_pocket:(i+1)*n_pocket]<=distanceTh #= row: distance relations between atom ligand and all pocket atoms"
            rowmap = np.logical_or(rw,rowmap)
            n_inLig+=np.any(rw) #at least one per row = one hit
        n_inP= np.sum(rowmap)#can be seen as a mask where a true entry is a hit, a pocket atom within the distance threshold

        # print(n_inLig,n_inP) #CHECK

        matchScore = n_inLig/n_ligand
        coverageScore = n_inP/n_pocket

        # i,j=getIndex(0,index,pocket_tAtoms.shape[0]) #first index refers to ligand atom; second index to pocket atom
        # can count unique indexes but I think there are more cleaver ways to do it..
        return matchScore,coverageScore



################## ACCESSORY FUNCTIONS (not important) ##################
    # def plotAtoms(self,ax):
    #     xa = self.atoms[:,0]
    #     ya = self.atoms[:,1]
    #     za = self.atoms[:,2]
    #     ra = self.atoms[:,3]
    #     for (xi,yi,zi,ri) in zip(xa,ya,za,ra):
    #         (xs,ys,zs) = drawSphere(xi,yi,zi,ri)
    #         ax.plot_wireframe(xs, ys, zs, linewidth = 0.1, color = 'black')



########### FUNCTIONS RELATED TO CLUSTERING PROCESS 

    def _get_plane(self,ind1,ind2 = None):
        """
            compute normals to the plane of the 2 triplets of atoms considered
            returns the normals and a reference atom

            IF one argument passed returns normal and all coordinates of trimming atoms 

            ----
            Parameters:
            2 list or np array, one for each triplet of atoms

        """
        #print("getting plane")
        if(ind2 is not None):

            a1_1 = self.atoms[ind1[0],0:3]
            a2_1 = self.atoms[ind1[1],0:3]
            a3_1 = self.atoms[ind1[2],0:3]

            a1_2 = self.atoms[ind2[0],0:3]
            a2_2 = self.atoms[ind2[1],0:3]
            a3_2 = self.atoms[ind2[2],0:3]

            #get normals 

            n1 = pcross((a2_1-a1_1),(a2_1-a3_1))#np.cross((a2_1-a1_1),(a2_1-a3_1))
            n2 = pcross((a2_2-a1_2),(a2_2-a3_2))#np.cross((a2_2-a1_2),(a2_2-a3_2))

            #normalize 
            n1 = n1/(np.sqrt(psum(n1*n1)))  #np.sum(n1**2)
            n2 = n2/(np.sqrt(psum(n2*n2)))  #np.sum(n2**2)
            
            return (n1,a1_1),(n2,a1_2)
        else:
            a1 = self.atoms[ind1[0],0:3]
            a2 = self.atoms[ind1[1],0:3]
            a3 = self.atoms[ind1[2],0:3]
            n =pcross((a2-a1),(a2-a3))#np.cross((a2-a1),(a2-a3))
            n = n/(np.sqrt(psum(n*n)))  #np.sum(n**2)
            ref_atoms = np.stack([a1,a2,a3])
            return n,ref_atoms

    def _inTriangle(self,n,ind,p):
        # COMMENTS 
        #Could we speculate about using these sub-areas to return also information on atoms closer to probe?
        # Maybe gives useful chemical clues?
        
        a1 = self.atoms[ind[0],0:3]
        a2 = self.atoms[ind[1],0:3]
        a3 = self.atoms[ind[2],0:3]

        #project everything on the plane
        refx = a2-a1
        refx = refx/(np.sqrt(psum(refx*refx)))  #np.sum(refx**2)
        refy = pcross(n,refx)#np.cross(n,refx)
        refy = refy/(np.sqrt(psum(refy*refy)))  #np.sum(refy**2)
        #refz = n of course
        a1p = np.array([np.dot(refx,a1),np.dot(refy,a1)])
        a2p = np.array([np.dot(refx,a2),np.dot(refy,a2)])
        a3p = np.array([np.dot(refx,a3),np.dot(refy,a3)])
        p_proj = p - np.dot(n,(p-a1))*n
        pp = np.array([np.dot(refx,p_proj),np.dot(refy,p_proj)])
        
        ######### 
        #Compare areas:
        #build triangles with permutations with p and compare them to the whole triangle area 
        A = area(a1p,a2p,a3p)

        A1 = area(a1p,a2p,pp)
        A2 = area(a2p,a3p,pp)
        A3 = area(a3p,a1p,pp)
        #print("areas:  ", A, A1+A2+A3)
        # if(abs(A-(A1+A2+A3))<=1.e-4):
        #     return True
        #relaxed version
        if(abs(A-(A1+A2+A3))<=1.e-4):
            return True
        else:
            return False


    def _inCircle(self,n,ind,p):
        a1 = self.atoms[ind[0],0:3]
        a2 = self.atoms[ind[1],0:3]
        a3 = self.atoms[ind[2],0:3]

        #project everything on the plane!!
        #fetch normal to the plane, then build reference frame on the plae
        refx = a2-a1
        refx = refx/(np.sqrt(psum(refx*refx))) #np.sum(refx**2)
        refy = pcross(n,refx)#np.cross(n,refx)
        refy = refy/(np.sqrt(psum(refy*refy))) #np.sum(refy**2)
        #refz = n of course
        a1p = np.array([np.dot(refx,a1),np.dot(refy,a1)])
        a2p = np.array([np.dot(refx,a2),np.dot(refy,a2)])
        a3p = np.array([np.dot(refx,a3),np.dot(refy,a3)])
        p_proj = p - np.dot(n,(p-a1))*n
        pp = np.array([np.dot(refx,p_proj),np.dot(refy,p_proj)])
        #or pp = (p-a1)-np.dot(n,(p-a1))*n
        ##########

        ax = a1p[0]-pp[0]
        ay = a1p[1]-pp[1]
        bx = a2p[0]-pp[0]
        by = a2p[1]-pp[1]
        cx = a3p[0]-pp[0]
        cy = a3p[1]-pp[1]

        det = np.sign((ax*ax+ ay*ay) * (bx*cy-cx*by) -(bx*bx + by*by) * (ax*cy-cx*ay) + (cx*cx + cy*cy) * (ax*by-bx*ay) )

        #check ordering of points (anticlockwise +1)
        s =(a2p[0]-a1p[0])*(a3p[1]-a1p[1]) - (a2p[1]-a1p[1]) * (a3p[0]-a1p[0])


        # if(s==0):
        #     print ("Collinear atoms.. PROBLEM")
        #     exit()
        return(np.sign(s)*det) #+1 inside, #-1  outside, 0 on the circumference



    def is_shift(self,probe1,ids1,probe2,ids2,getTag=False):
        """
        CAREFUL: this function is thought to compares probes which are close enough. 
        This check must be done externally so far.

        Returns True if probes are at the same side of both the planes relative to them AND share at least 2 identical atom ids.
        """


        plane1,plane2 = self._get_plane(ids1,ids2) #Even too much.. one of the two should be sufficient
        n1 = plane1[0]
        o1 = plane1[1]
        n2 = plane2[0]
        o2 = plane2[1]

        t1 = coplanar(n1,o1,probe1,probe2)
        t2 = coplanar(n2,o2,probe1,probe2)
        
        if(getTag):
            if(not(t1 and t2)):
                # print(t1,t2)
                return False,-1
            if(compare_ids(ids1,ids2)==3):
                #print(ids1,ids2)
                return True,3
            elif(compare_ids(ids1,ids2)==2):
                return True,4 #weak shift aggregation
            else:
                return False,-1
        else:
            if(not(t1 and t2)):
                # print(t1,t2)
                return False
            if(compare_ids(ids1,ids2)>=2):
                #print(ids1,ids2)
                return True
            else:
                return False 

    def is_btlnk(self,dist,probe1,ids1,probe2,ids2,tollerance,indulgent=False):
        """
        Returns 
        - True: if distances and geometrical properties correspond to a bottleneck
        - IDS of atoms defining the plane spanning the bottleneck entrance
        - Ref_probe coordinate, i.e. the probe defining the trimming plane (of the protein tangent atoms)
        ---
        Optional:
        Relax inTriangle requirement with only inCircle requirement
        --- 
        Note: 
        Probes can be of different radius, this is irrelevant.
        """

        plane1,plane2 = self._get_plane(ids1,ids2)
        n1 = plane1[0]
        o1 = plane1[1]
        n2 = plane2[0]
        o2 = plane2[1]
        

        pp = probe1-probe2
        ppn1 = abs(np.dot(pp,n1))
        ppn2 = abs(np.dot(pp,n2))

        conf1 = abs(ppn1-dist)
        conf2 = abs(ppn2-dist)


        test_side = False
        test_in= False


        if ((conf1<conf2)and(conf1<tollerance)):
            n=n1
            o=o1
            atomids = ids1

            ref_probe = probe1

            test_side = coplanar(n,o,probe1,probe2)

            if(test_side):
                return False,None,None

            if(indulgent):
                t1 = self._inCircle(n,atomids,probe1)
                t2 = self._inCircle(n,atomids,probe2)
            else:
                #Stronger requirement
                t1 = self._inTriangle(n,atomids,probe1)
                t2 = self._inTriangle(n,atomids,probe2)
            if(t1 and t2):
                test_in = True

            if(test_in):
                # print("confidence= ", conf)
                return True,atomids,ref_probe
        elif(conf2<tollerance):
            n = n2
            o = o2
            atomids = ids2
            ref_probe = probe2
            test_side = coplanar(n,o,probe1,probe2)
            
            if(test_side):
                return False,None,None

            if(indulgent):
                t1 = self._inCircle(n,atomids,probe1)
                t2 = self._inCircle(n,atomids,probe2)
            else:
                #Stronger requirement
                t1 = self._inTriangle(n,atomids,probe1)
                t2 = self._inTriangle(n,atomids,probe2)
            if(t1 and t2):
                test_in = True

            if(test_in):
                # print("confidence= ", conf)
                return True,atomids,ref_probe

        return False,None,None
    ###############################


# CLASSES/DATA CONTAINERS AND STRUCTURES

 #tag : tag of the node of last aggregation
 # -1 singleton or aggregation with singleton
 # aggregation between non-singleton clusters:
 #      0: lateral
 #      4: pyramidal
 # 2 bottleneck
 # 3 conical shift (only if at least one of the two clusters is a singleton)
 # 1 pocket (after post processing)


 #https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.ClusterNode.html#scipy.cluster.hierarchy.ClusterNode
 
 # preorder traversal function retrieves child nodes ids --> retrieves elements of a cluster
import numpy as np

from siteFerret.functions import naiveV,trueV,NSvolume,setup_NSInput,new_probe,OC,centerOfMass,extractOFF
from siteFerret.functions import Pdist_C

from  siteFerret import global_module
from siteFerret.global_module import Error
from siteFerret.protein_class import Protein

shift_Threshold = global_module.shift_Threshold


class ClusterNode(object):
    """
A tree node class for representing a cluster.
Leaf nodes correspond to original observations, while non-leaf nodes
correspond to non-singleton clusters.

----------
id : int
    The node id.
left : ClusterNode instance, optional
    The left child tree node.
right : ClusterNode instance, optional
    The right child tree node.
dist : float, optional
    Distance for this cluster in the linkage matrix.
count : int, optional
    The number of samples in this cluster.

-------------


tag
-1 singleton or aggregation to singleton
0 lateral aggregation with no special signature 
4 conical aggretation with no special signature (pyramid)
1 major candidate pocket 
2 simple bottleneck --> internal cavity candidate
3 conical shift

position = for leaves actual coordinates. For nodes meaningful probe position for tag event
t_atoms = IDS tangent atoms to the probe for leaves
          tangent atoms to bottleneck or to conical shift for non-leaf nodes
radius = 
    - Probe radius for leaves
    - None for nodes with no special tags. 
    - Probe defining the relevant plane for Bottlenecks
    - Larger probe for conical shift
NEW: coord = (position,radius) tuple
taglist: updated in gather_tag() traverse all tags and update signatures of pockets and cavities
"""
#Class variable
    wasAccurate = False #to remember if NS volume calculation is in light mode self.__class__.
    not_firstTime = False
    def __init__(self, id,coord,t_atoms, radius = None, left=None, right=None, tag=-1, count=1):
        if id < 0:
            raise ValueError('The id must be non-negative.')
        # if dist < 0:
        #     raise ValueError('The distance must be a non-negative integer.')
        if (left is None and right is not None) or \
            (left is not None and right is None):
            raise ValueError('Only full or proper binary trees are permitted.'
                                '  This node has one child.')
        if count < 1:
            raise ValueError('A cluster must contain at least one original '
                                'observation.')
        self.id = id
        self.coord = coord 
        self.t_atoms = t_atoms 
        self.r = radius
        self.left = left
        self.right = right
        
        self.tag = tag

        self._tagList = np.empty(0)
        self.index_btlnk = []
        self.index_Cshift = []

        self.vol_acc = None
        self.area_acc = None
        self.vol_Nacc = None
        self.area_Nacc = None

        self._idList = np.empty(0)
        self._shift = None
        #### 

        self._preorder = []
        self._coordinates = []
        self._atomList = []
        # self._wasDone = False
        self.protein = None


        if self.left is None:
            self.count = count
        else:
            self.count = left.count + right.count

    def load_protein(self,in_protein):
        if(self.protein==None):
            if(isinstance(in_protein,Protein)):
                self.protein = in_protein
                return 
            else:
                # try:
                raise TypeError ("Only protein class types can be loaded into a cluster node using load_protein method!")
            # except TypeError:
            #     info = "Only protein class types can be loaded into a cluster node using load_protein method!"
            #     return 1,info
        else:
            return #do nothing, already loaded

    def get_id(self):
        """
        Returns
        -------
        id : int
            The identifier of the target node.
        """
        return self.id

    def get_count(self):
        """
        The number of leaf nodes (original observations) belonging to
        the cluster node nd. If the target node is a leaf, 1 is
        returned.
        Returns
        -------
        get_count : int
            The number of leaf nodes below the target node.
        """
        return self.count

    def get_left(self):
        """
        Return a reference to the left child tree object.
        Returns
        -------
        left : ClusterNode
            The left child of the target node. If the node is a leaf,
            None is returned.
        """
        return self.left

    def get_right(self):
        """
        Return a reference to the right child tree object.
        Returns
        -------
        right : ClusterNode
            The left child of the target node. If the node is a leaf,
            None is returned.
        """
        return self.right

    def is_leaf(self):
        """
        Return True if the target node is a leaf. 
        Returns
        -------
        leafness : bool
            True if the target node is a leaf node.
        """
        return self.left is None
####################
#NOTE:
# The following "pre_order" type functions could all be wrapped in a single one, being the traversal identical.
# However, most of the time the algorithms uses the pre_order() function returning the ids, and much less the other types which are invoked
# on a very small subset of clusters which are recognised as pockets. So to save data I preferred to separate them.
    def getAtomsIds(self):
        """
        Returns ids of protein atoms tangent to the spheres contained in the cluster. 
        CAREFUL: these are ids not coordinates. 
        """
        if(self._atomList):
            return self._atomList

        n = self.count
        curNode = [None] * (2 * n)
        lvisited = set()
        rvisited = set()
        curNode[0] = self
        k = 0


        atomList = []

        while k >= 0:
            nd = curNode[k]
            ndid = nd.id
            if nd.is_leaf():
                atomList.append(nd.t_atoms)
                k = k - 1
            else:
                if ndid not in lvisited:
                    curNode[k + 1] = nd.left
                    lvisited.add(ndid)
                    k = k + 1
                elif ndid not in rvisited:
                    curNode[k + 1] = nd.right
                    rvisited.add(ndid)
                    k = k + 1
                # If we've visited the left and right of this non-leaf
                # node already, go up in the tree.
                else:
                    k = k - 1
        
        self._atomList = atomList
        return atomList

    def getCoordinates(self):
        """
        Returns x,y,z of each sphere contained in the cluster using a pre-order traversal.
        """
        if(self._coordinates):
            return self._coordinates

        n = self.count
        curNode = [None] * (2 * n)
        lvisited = set()
        rvisited = set()
        curNode[0] = self
        k = 0


        coordinates = []

        while k >= 0:
            nd = curNode[k]
            ndid = nd.id
            if nd.is_leaf():
                coordinates.append(np.append(nd.coord,nd.r))
                k = k - 1
            else:
                if ndid not in lvisited:
                    curNode[k + 1] = nd.left
                    lvisited.add(ndid)
                    k = k + 1
                elif ndid not in rvisited:
                    curNode[k + 1] = nd.right
                    rvisited.add(ndid)
                    k = k + 1
                # If we've visited the left and right of this non-leaf
                # node already, go up in the tree.
                else:
                    k = k - 1
        
        self._coordinates = coordinates
        return coordinates

    def pre_order(self):
        """
        Perform pre-order traversal without recursive function calls.
        Update attibutes of the node.
        Most important: update current node with relevant tag looking at child nodes.

        id_btlnk: lsit of IDS
        id_conShift: list of tuples (ID,persistence) where persistence counts how many times tag=3 found on same tangent atoms (conical shift)
        
        Returns a list of nodes ids corresponding to the leaf nodes of the tree as they appear from left to right.

        OBS. After creation each node is not supposed to be modified. Its attributes depend only on its child nodes
        """

        if self._preorder:
            #The traversal has been already performed, attributes have been already updated
            return self._preorder


        n = self.count
        curNode = [None] * (2 * n)
        lvisited = set()
        rvisited = set()
        curNode[0] = self
        k = 0

        
        ids = []

        ### used regarding tag retrieval and pocket identification
        # These are properties of non-singleton nodes

        tagList = []
        idList = []
        nodeAtomList = [] #relevant protein atoms for the event that generated the clustering
        

        ##

        while k >= 0:
            nd = curNode[k]
            ndid = nd.id
            if nd.is_leaf():
                ids.append(nd.id)
                k = k - 1
            else:
                if ndid not in lvisited:
                    curNode[k + 1] = nd.left
                    lvisited.add(ndid)
                    k = k + 1
                elif ndid not in rvisited:
                    curNode[k + 1] = nd.right
                    rvisited.add(ndid)
                    k = k + 1
                # If we've visited the left and right of this non-leaf
                # node already, go up in the tree.
                else:
                    #Here NON-leaf nodes
                    tagList.append(nd.tag)
                    idList.append((nd.id,nd._shift))
                    nodeAtomList.append(nd.t_atoms)
                    k = k - 1
        
        self._preorder = ids
        
        if((self.tag == 3)):
            tagList = np.asarray(tagList)
            idList = np.asarray(idList)
            #tag list was filled in reverse order from deep to top..
            tagList = tagList[::-1]
            idList = idList[::-1]
            nodeAtomList = nodeAtomList[::-1]
            c = (tagList == 3) | (tagList == 1)
            nodeAtomList = [np.sort(a) for (a, truth) in zip(nodeAtomList, c) if truth]
            shifts = (nodeAtomList == np.sort(self.t_atoms)).all(-1).sum() + 1

            if(shifts>=shift_Threshold):
                #Special tag only when it has structure. Useful for scoring and reward more channelled structures rather than isolated tag=3
                self.tag=1 
                tagList[0] = 1 #update tagList
                # if(n>=count_Threshold):
                self._shift = shifts
                idList[0][1] = shifts #update id list with persistence
                ib = np.where(tagList == 2) 
                index_btlnk = [i[0] for i in idList[ib]]
                index_Cshift = list(filter(lambda x: x[1]!= None, idList))# Pseudo-mouths

                self._tagList = tagList #taglist container filled only when useful, that is for pocket clusters
                #mouth lists fetched and data container updated only for successful pockets
                self._idList = idList[:,0]
                self.index_btlnk = index_btlnk
                self.index_Cshift = index_Cshift #Contains redundancies
        else:
            pass
        
        
        return ids


    def score(self,w,n_uniqueBtlnk = None):
        """
        EURISTIC SCORE
        Input : weight
        Returns a "tortuosity" score of the pocket (geometrical feature).

        - Pocket shift score = 1 (many contiguous are associated to deeper channels)
        - Bottleneck score weighted by w 
        - Might add tag=0 to weight more size BUT I loose sensitivity on important events.. 

        -When n_unique (= meant to be # botllenecks after filtering) is provided, it uses this to count btlnks events.

        Error raised if tagList was not built
        """
        if(self._tagList.size == 0):
            raise IndexError("A call to gather tag for this node should be implemented")
        if(n_uniqueBtlnk ==None):
            #s = np.sum((self._tagList==1) | (self._tagList==3) | w * (self._tagList == 2))
            s = np.sum((self._tagList==1) | w * (self._tagList == 2) | (self._tagList==0) | (self._tagList==3))
        else:
            # s = np.sum((self._tagList==1) | (self._tagList==3) ) + w * n_uniqueBtlnk
            s = np.sum((self._tagList==1) | (self._tagList==0)| (self._tagList==3)) + w * n_uniqueBtlnk
        #ramification = np.sum((self.tagList==3))
        return s


    def get_aggregations(self):
        if(self._tagList.size == 0):
            raise IndexError("A call to gather tag for this node should be implemented")
        idsAggregations_lat = [('lat',ids)for ids in self._idList[(self._tagList==0)]] #o: lateral, 1: "pyramidal"
        idsAggregations_pyr = [('pyr',ids)for ids in self._idList[(self._tagList==4)]]
        idsAggregations = idsAggregations_pyr+idsAggregations_lat #concatenate the 2 lists
        return idsAggregations


    def volume(self,consider_overlap = True):
        """
        ARVO volume calculator (anlytic method, slow but accurate)
        """
        err =Error()
        if(consider_overlap):
            if(self.vol_acc is None):
                coord = self.getCoordinates()
                try:
                    V,A=trueV(coord)
                except Exception as e:
                    err.value=1
                    err.info = e
            # self.wasAccurate = True #this is not overridden by class variable change, avoids reset if the calculation was already performed..
                if(V==0):
                    err.value=1
                    err.info = "Volume could not be computed neither with ARVO method.."
                    return self.vol_acc,self.area_acc,err
                self.vol_acc=V
                self.area_acc=A
            return self.vol_acc,self.area_acc,err
        else:
            coord = self.getCoordinates()
            return naiveV(coord),0,err
        

    #NEW NanoShaper based volume calculation..
    def buildTriang(self,triangName,isSkin=True,savePath = ".",accGridScale = 6):
        """
        Only to produce cluster triangulation.
        This is a reduced version of NSVolume() without all the checks of redundancies 
        (does not care if volume was already computed, or input file already created)
        """
        grid_selfInt=3
        maxProbes_selfInt=200
        err = Error()
        setup_NSInput(global_module.runFolder_path+global_module.conf,grid_scale=accGridScale,grid_selfInt=grid_selfInt,maxProbes_selfInt=maxProbes_selfInt,isSkin=isSkin,accTriang=True,pInput= True)
        new_probe(global_module.runFolder_path+global_module.conf,1.4) #use rp=1.4 to build pocket triangulation
        _V,_A,errline = NSvolume(self.getCoordinates(),triangName,savePath)
        if(errline):
            err.info = errline
            err.value = 1
        return err
    def NSVolume(self,accurate,triangName=False,savePath = ".",accGridScale = 6,isSkin=True):
        """
        Computes volume and optionally save triangualtion OFF file
        Pocket volume is a triangulation of several overlapping spheres. NS operative mode "skin" is more stable and is faster in these scenarios
        """
        defGridScale=4 #4 empirically the best to not incur into NS errors..
        grid_selfInt=3
        maxProbes_selfInt=200
        err = Error()
        errline=''
        #acc=false wasAccurate=True -> change NS input for fast triangulation and V calculation
        #acc = true wasAccurate=False-> change NS input for accurate triangulation
        #Other scenario skip changing input T T=T, F F = T 

        if(accurate):
            if (self.vol_acc is not None):
                # print("V accurate was already computed")
                return self.vol_acc,self.area_acc,err
            else:
                if ((not bool(accurate^self.wasAccurate))and self.not_firstTime): #XNOR + AND
                    # print('skipping change in NS input')
                    pass
                else:
                    # print("changing NS input")
                    self.__class__.wasAccurate = True
                    setup_NSInput(global_module.runFolder_path+global_module.conf,grid_scale=accGridScale,grid_selfInt=grid_selfInt,maxProbes_selfInt=maxProbes_selfInt,isSkin=isSkin,accTriang=True,pInput= True)
                    new_probe(global_module.runFolder_path+global_module.conf,1.4) #use rp=1.4 to build pocket triangulation
                # print("Computing accurate")
                V,A,errline = NSvolume(self.getCoordinates(),triangName,savePath)
                self.vol_acc=V
                self.area_acc = A
        else:
            if (self.vol_Nacc is not None):
                # print("V non accurate was already computed")
                return self.vol_Nacc,self.area_Nacc,err
            else:
                if ((not bool(accurate^self.wasAccurate))and self.not_firstTime): #XNOR + AND
                    # print('skipping change in NS input')
                    pass
                else:
                    # print("changing NS input")
                    self.__class__.wasAccurate = False
                    setup_NSInput(global_module.runFolder_path+global_module.conf,grid_scale=defGridScale,grid_selfInt=grid_selfInt,maxProbes_selfInt=maxProbes_selfInt,isSkin=isSkin,accTriang=False,pInput= True)
                    new_probe(global_module.runFolder_path+global_module.conf,1.4) #use rp=1.4 to build pocket triangulation
                # print("Computing non accurate")
                V,A,errline = NSvolume(self.getCoordinates(),triangName,savePath)
                self.vol_Nacc=V
                self.area_Nacc = A

        if(errline):
            # print("HERE!!")
            err.value=1
            err.info = errline
            # print(err.value)
        
        self.__class__.not_firstTime = True 
        return V,A,err

    #########

    def pArea_matchScore(self,ligand_coordinates,scoreTh=0.2):
        """
        Criterion based on pocket as aggregate of cluster sphere. Geometrically interesting but chemically weak.
        If matching found retuns 2 numbers representing overlap and size score.
        """
        coord = self.getCoordinates()
        success,overlap_score,ligand_volume = OC(ligand_coordinates,coord,scoreTh=scoreTh)
        if(success):
            print("Overlap score successful!")
            # print("Ligand volume=",ligand_volume)
            print("Overlap score=",overlap_score)
            pocket_volume=self.volume()
            if(pocket_volume>=ligand_volume):
                size_score =ligand_volume/pocket_volume
            else:
                size_score =pocket_volume/ligand_volume
            # score = 1-abs(overlap_score-size_score)
            print("Pocket volume:",pocket_volume)
            return overlap_score,size_score
        else:
            print("No sufficient overlap")
            print("Overlap score=",overlap_score)
            return overlap_score

    
    def matchScore(self,ligand_coordinates,distanceTh = 5, matchScoreTh=0.5, coverageTh=0.2):
        """
        Returns success if ligand identified according to criterion of matchCriterion() together with the 2 scores.
        NOTE: success is evaluated up to a rounding to the first decimal digit
        """
        pAtoms = self.getAtomsIds()
        if(self.protein==None):
            raise FileNotFoundError("Protein not loaded in cluster node: Cannot compute matching!")
        
        matchScore,coverageScore=self.protein.matchCriterion(ligand_coordinates,pAtoms,distanceTh)
        
        if((np.round(matchScore,2)>=matchScoreTh)and(np.round(coverageScore,2)>=coverageTh)):
            success = True
        else:
            success = False
        return success,matchScore,coverageScore

    def CMmatchScore(self,ligand_coordinates,distanceTh = 4):
        """
        Returns success if any ligand atom within 4A from pocket center of mass. The center of mass is extracted from the triangulation..
        """
        self.buildTriang(triangName='temp')
        verts,_faces = extractOFF('temp.off')
        CM = centerOfMass(verts)
        ligand_coordinates = np.array(ligand_coordinates)[:,0:3]
        d,_flag = Pdist_C(ligand_coordinates,CM.reshape((-1,1)).T)

        # print(d)
        d = np.round(d)
        # print(d)
    

        success = np.any(d<=distanceTh)

        # print(success)
        
        return success


    def getCentroid(self,weighted=False):
        '''Returns the center of mass of the cluster'''
        spheres = np.array(self.getCoordinates())
        coord = spheres[:,0:3]
        if(weighted):
            radius = spheres[:,3]
            return sum(radius[:,None]*coord)/(sum(radius))
        else:
            n=self.count
            return sum(coord)/n 

    

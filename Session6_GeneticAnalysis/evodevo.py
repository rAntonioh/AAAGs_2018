import numpy as np
from scipy import stats
import pandas as pd
import uuid
from abc import ABCMeta
from copy import deepcopy

### gene regulatory network class to iterate genetic interaction and expression dynamics ###

class GeneRegulatoryNetwork():
    """Class to make a gene regulatory network for use in evo-devo simulation."""

    __metaclass__ = ABCMeta

    def __init__(self, name):
        # founder properties
        self.uuid = str(uuid.uuid1())
        self.name = str(name)
        self.n = 0
        self.c = 0.
        self.a = 0

        # to run iterations
        self.bp = 0
        self.mxi = 0
        self.epsilon = 0.
        self.wc = 10  # value tau from Siegal & Bergman
        # self.wc = 3 # works well too but may exclude oscillators

        self.df = None
        
        # to hold values from iterations
        self.calc_var_buffer = []
        self.sbuffer = []
        self.dbuffer = []
        
        # for gene expression
        self.s0 = []
        self.s = []
        self.s_final = []
        self.ref_state = []
        
        #holds GRN
        self.w = []

        # to mutate the founder object
        self.delta = float
        self.mutated = False
        self.nom = 0  # number of mutations
        self.mut_index = []

        # quality of individual object values
        self.stabilized = False
        self.fitness = float
        self.distance = float
        self.path_at_stability = int
        self.is_accepted = bool
        
        # multicellular development
        self.corpus = object
        return
    
    def make_grn(self, n, c, a, mxi):
        """Set up basic GRN."""
        self.n = n
        self.c = c
        self.a = a
        self.mxi = mxi
        
        self.sbuffer = np.zeros((self.n, self.mxi))
        self.dbuffer = np.zeros((self.n, self.mxi))

        s0 = np.ndarray((self.n, 1))
        i = 0
        for _ in s0:
            s0[i] = (np.random.choice((0, 1)) * (np.random.rand(1) < .5))
            i+=1
        del i

        self.s0 = s0
        randweights = np.random.normal(size=(self.n, self.n))
        r = np.random.rand(self.n, self.n)
        self.w = np.multiply((r < self.c), randweights)
        del s0, randweights, r

    def __calculate_variance(self, bp, values):
        """Private method to calculate variance of gene expression during development."""
        bp = bp
        mxi = self.mxi
        wc = self.wc
        if bp < mxi:
            self.sbuffer[:, bp] = np.ndarray.flatten(values)
            result = 1.0
            if bp >= wc and bp < mxi:
                sp = bp - wc
                result = np.mean(np.var(self.sbuffer[:, sp:bp], axis=1)) / (4 * self.n)  # mean of the local variance
                del sp
                self.dbuffer[:, bp] = np.var(self.sbuffer[:, 0:bp], axis=1)  # distance measure from the start
            self.bp = bp + 1
            del bp, mxi, wc
            return result

    def __calculate_distance(self, ref, test):
        dist = np.vdot(np.subtract(test, ref), np.subtract(test, ref)) / (4*self.n)
        return dist 
    
    def sigmoid(self, a, x):
        np.seterr(all='ignore')
        r = np.round(1. / (1. + np.exp(-a * x)))
        del a, x
        return r

    def start_development(self, ss, epsilon):
        """Method to develop and grow"""
        # print('starting development')
        self.stabilized = False
        self.path_at_stability = 100
        self.s_final = []
        self.epsilon = epsilon
        self.calc_var_buffer = []
        self.sbuffer = np.zeros((self.n, self.mxi))
        self.dbuffer = np.zeros((self.n, self.mxi))

        s = ss
        w = self.w
        i = 0
        while not self.stabilized and i < self.mxi:
            s = self.sigmoid(self.a, np.dot(w, s))
            calc_var = self.__calculate_variance(i, s)
            self.calc_var_buffer.append(calc_var)
            if calc_var < self.epsilon:
                self.s_final = s
                self.stabilized = True
                self.path_at_stability = i
                self.distance = []
                self.distance = self.__calculate_distance(self.s0, self.s_final)
                # print('done developing!', self.path_at_stability)
            else:
                self.path_at_stability = i
                self.stabilized = False
                self.s_final = s
            i += 1
        return
    
    def mutate(self, delta, forced):
        """Method to mutate GRN, allows for new connections to be formed"""
        self.delta = np.divide(0.1,(self.c * (self.n * self.n)))
        if forced is True:
            self.delta = delta
        d = self.w
        u = np.random.uniform(0, 1, (self.n, self.n))
        mutate_mask = np.less(u,self.delta)
        self.nom = np.sum(mutate_mask)
        while self.nom == 0:
            del mutate_mask
            del d
            del u
            return
        else:
            normal_mask = u >= self.delta
            normal_w = np.multiply(normal_mask, d)
            mutant_w = np.multiply(mutate_mask, d)
            mutant_w = np.multiply(mutate_mask, np.random.normal(0, 1, (self.n, self.n)))
            rr, cc = np.nonzero(mutant_w) # find what is mutant
            self.mut_index = (rr,cc)
            # print(self.w[rr,cc]==0) # print if new connections
            new_w = normal_w + mutant_w
            self.w = new_w
            self.mutated = True
            del new_w
            del u
            del mutate_mask
            del mutant_w
            del normal_mask
            del normal_w
            del d
            return   

    def perturb(self):
        self.ref_state = self.s_final
        rows, cols = np.nonzero(self.w)
        ridx = np.random.choice(rows)
        cidx = np.random.choice(cols)
        self.w[ridx][cidx] = np.random.normal(0,1)
        return
    
    def multicellular_development(self):
        self.corpus.coordinated = []
        self.corpus.state_buffer = []
        for cell in self.corpus.cells:
            self.s_final = []
            self.sbuffer = []
            cell.state_buffer = []
            self.start_development(cell.initial_state, self.epsilon)
            cell.final_state = self.s_final
            cell.state_buffer = self.sbuffer
            cell.path_at_stability = self.path_at_stability
            cell.stabilized = self.stabilized
            # print('cell',cell.stabilized)
            if cell.path_at_stability in self.corpus.developmental_ranges:
                self.corpus.coordinated.append(True)
            else:
                self.corpus.coordinated.append(False)
        return

    def characterize(self):
        self.df = pd.Series({
                    'uuid': self.uuid,
                    'N':self.n,
                    'c':self.c,
                    'a':self.a,
                    'mxi':self.mxi,
                    'nom':self.nom,
                    'soma':self.corpus.cells[0].path_at_stability,
                    'gonad':self.corpus.cells[1].path_at_stability,
                    'path_at_stability': self.path_at_stability,
                    'stabilized': self.stabilized,
                    's0': self.s0.flatten(),
                    's_final': self.s_final.flatten(),
                    'dev_distance': self.distance,
                })
        return
        
    def __deepcopy__(self, memo):
        deepcopy_method = self.__deepcopy__
        self.__deepcopy__ = None
        cp = deepcopy(self, memo)
        self.__deepcopy__ = deepcopy_method
        cp.name = str(self.uuid) + '.copy'
        cp.bp = 0
        cp.stabilized = False
        return cp  
    

### organism class to govern GRN dynamics for simulation ###

class Organism(GeneRegulatoryNetwork):
    """Superclass of GeneRegulatoryNetwork to perform operations with"""
    pass
def generate_combinations(m, n, g, i, founder):
    """
        Function to generate all state vector combinations for a given founder.
        m = vector length, n = options [0.,1.], g = job index, i = chunk size of state space 
    
    """
    tmp = founder
    combinations = []
    vectors = []
    paths = []
    states = []
    attractors = []
    sbuff = []
    t = g + i if g + i < (len(n) ** m) else len(n) ** m
    for num in np.arange(g, t):
        combination = []
        while num != 0:
            combination.insert(0, n[num % len(n)])
            num = int(num/len(n))
        while len(combination) < m:
            combination.insert(0, n[0])
        if combination not in combinations:
            tmp.s0 = np.reshape(combination,(m,1))
            tmp.start_development(tmp.s0, founder.epsilon)
            combinations.append(combination)
            vectors.append(tmp.s0)
            paths.append(tmp.path_at_stability)
            states.append(tmp.s_final)
            attractors.append(tmp.sbuffer[:,1])
            sbuff.append(tmp.sbuffer)
        i -= 1
        if i == 0:
            break
    return vectors, np.asarray(states).squeeze(), paths, attractors, sbuff

def state_space(founder):
    options = [0.0, 1.0]
    chunk_size = len(options)**founder.n
    gob_id = np.arange(0,1)
    plas_buffer = []
    vector_buffer = []
    state_buffer = []
    for gob in gob_id:
        gob_idx = gob_id[gob] * chunk_size
        vectors, states, plas, attractors, sbuff  = generate_combinations(founder.n, options, gob_idx, chunk_size, founder)
        plas_buffer.append(plas)
        vector_buffer.append(vectors)
        state_buffer.append(states)
    plas_buffer = np.asarray(plas_buffer).squeeze()
    vector_buffer = np.asarray(vector_buffer).squeeze()
    state_buffer = np.asarray(state_buffer).squeeze()

    return plas_buffer, vector_buffer, state_buffer

def analyze_attractors(founder, vector_buffer, state_buffer):
    """Takes a founder, all vectors and all outputs as input and generates the attractor network in one iteration"""
    # sigmoid parameter
    a=100
    
    # sigmoid function
    def sigmoid(a, w, s):
        np.seterr(all='ignore')
        x = np.dot(w,s)
        r = np.round(1. / (1. + np.exp(-a * x)))
        del a, x
        return r

    # founder network
    w = founder.w
    
    ### starting and final vectors are supplied ###
    
    # vector buffer of starting vectors
    vb = vector_buffer
    
    # final state at stability buffer
    sb = state_buffer
    
    # next buffer for first iteration to find the attractor network mapping ###
    nb = []
    
    # compute and keep the next state from the network
    for v in vb: 
        nb.append(sigmoid(a,w,v))
    
    
    # create an edges list from location of final state vector found in vector_buffer    
    i=0
    sedge_list = []
    sidx = []
    for s in sb:
        j=0
        for v in vb:
            x = np.equal(s,v)
            if np.sum(x)==10:
                sedge_list.append((i,j))
                sidx.append(j)
            j+=1
        i+=1

    # finding the frequencies finds how many edges each state gets
    sfreqs = stats.itemfreq(sidx)


    # create an edges list and index from first vectors in vector_buffer to location of next vector found in vector_buffer
    i=0
    nedge_list = []
    nidx = []
    for n in nb:
        j=0
        for v in vb:
            x = np.equal(n,v)
            if np.sum(x)==10:
                nedge_list.append((i,j))
                nidx.append(j)
            j+=1
        i+=1

    ### pull out the item frequencies for next states identified in nb ###
    nfreqs = stats.itemfreq(nidx)
    
    # split out index of unique states in state_buffer (this is also the attractor index)
    x1 = np.asarray(list([row[0] for row in sfreqs]))
    # split out frequencies of unique states in state buffer
    y1 = np.asarray(list([row[1] for row in sfreqs]))

    # pull out to store the attractors and index together in tuple
    attractors = []
    for u in x1:
        attractors.append((u,state_buffer[u]))
 
    att_map = []
    i=0
    for att in attractors:
        for tta in attractors:
            diffs = np.subtract(att[1], tta[1])
            steps = np.vdot(diffs, diffs)
            att_map.append(steps)
            i+=1

    att_map = np.asarray(att_map).reshape((len(attractors),len(attractors)))
    
    # all the unique vectors present in nb
    x2 = np.asarray(list([row[0] for row in nfreqs]))
    y2 = np.asarray(list([row[1] for row in nfreqs]))
    
    return att_map,attractors,(x1,y1),(x2,y2)


# find state within a reference and return index
def find_state_index(states, reference):
    """takes inputs of state vectors and compares to reference, states,
    returns linear array the length of states."""
    res = []
    for s in states:
        s=np.asarray(s)
        i=0
        while i < 1024:
            t = np.asarray(reference[i])
            tf = np.equal(s,t)
            if np.sum(tf)==10:
                res.append(i)
            i+=1
    return np.asarray(res)

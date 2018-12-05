import random
import math
import matplotlib.pyplot as plt
from collections import Counter
import numpy
import time
import multiprocessing.dummy
from multiprocessing import Queue
import sys
import itertools

class Cell(object):
    """
    Representation of a HeLa cell with two properties (list of Bacterium) bound_bacteria and (bool) invaded. Additionally
    
    
    """
    def __init__(self):
        self.invaded = False
        self.bound_bacteria = []
        self.active_invaders = [] #Bacteria that invade by them selves
        self.passive_invaders = [] #Bacteria that coinvade with active invaders
        self.noninvaders = [] #Bound Bacteria that fail to invade
    
    def isInvaded(self):     
        """
        Returns True if the cell is invaded, False othervise
        """
        return self.invaded
           
    def getBoundBacteria(self):
        """
        Returns a list of Bactrium objects that have successfully bound the Cell.
        If no bacteria have bound an empty list is returned.
        """
        return self.bound_bacteria
    
    def getActiveInvaders(self):
        """
        Returns a list of Bactrium objects that have successfully invaded the Cell by themselves.
        If no bacteria have invaded an empty list is returned.
        """
        return self.active_invaders
    
    def getPassiveInvaders(self):
        """
        Returns a list of Bactrium objects that have successfully invaded the Cell through coinvasion.
        If no bacteria have coinvaded an empty list is returned.
        """
        return self.passive_invaders
    
    def getNonInvaders(self):
        """
        Returns a list of Bactrium objects that failed to invade the Cell.
        If no bacteria have bound and failed to invade empty list is returned.
        """
        return self.noninvaders
    
    def getBoundTags(self):
        """
        Returns a string consisting of the tags of the bound Bactrium objects.
        If no bacteria have bound an empty string is returned.
        """
        out = ''
        for b in self.bound_bacteria:
            out += b.getBarcode()
            
        return out
    
    def getInvadingTags(self, active=True, passive=True):
        """
        Returns a string consisting of the tags of invading Bactrium objects.
        Defaults to returning all invafing tags, both passive and active invaders
        If no bacteria have invaded an empty string is returned.
        """
        out = ''
        
        if active:
            for b in self.active_invaders:
                out += b.getBarcode()

        if passive:
            for b in self.passive_invaders:
                out += b.getBarcode()

        return out
    
    def getNumberOfBoundBacteria(self):
        """
        Returns the number (int) of bound Bactrium objects.
        """
        return len(self.bound_bacteria)
    
    def setInvaded(self):
        """
        Sets the state of the Cell to invaded
        """
        self.invaded = True
    
    def bindBacterium(self, bacterium):
        """
        Appends a Bacterium object to the list of bound bacteria,
        tests if the bacterium successfully invades, sets Cell as invaded if
        the invasion is the first successful invasion and appends the Bacterium to active invaders
        """
        self.bound_bacteria.append(bacterium)
        
        target_difficulty = random.random()
        
        if (bacterium.getP_inv() > target_difficulty) and not self.invaded:
            self.setInvaded()
            
        if (bacterium.getP_inv() > target_difficulty):
            self.active_invaders.append(bacterium)
        
        else:
            self.noninvaders.append(bacterium)
    
    def doCoinvasion(self):
        """
        Noninvading bacteria are added to passive invaders with a probability of p_coinvade
        """
        if self.invaded:
            for b in list(self.noninvaders): #iterating over a copy of the list to be able to modify it in the loop
                if b.getP_coinvade() > random.random():
                    self.noninvaders.remove(b)
                    self.passive_invaders.append(b)
        
        
            
class Bacterium(object):
    """
    Class to represent a bacterium binding and invading a Cell
    """
    def __init__(self, barcode, p_inv, p_bind, p_coinvade):
        
        self.barcode = barcode
        self.p_inv = p_inv
        self.p_bind = p_bind
        self.p_coinvade = p_coinvade
    
    def getBarcode(self):
        """
        Returns (str) barcode of bacterium
        """
        return self.barcode
    
    def getP_bind(self):
        """
        Returns (float) p_bind
        """
        return self.p_bind
    
    def getP_inv(self):
        """
        Returns (float) p_inv, the probability that the bacterium will invade the Cell.
    
        """
        return self.p_inv

    def getP_coinvade(self):
        """
        Returns (float) p_coinvade, the probability that the bacterium will invade the Cell if
        another bacterium has successfully invaded the cell.
    
        """
        return self.p_coinvade

class Experiment(object):
    """
    Represents running an invasion assay in one well including recovery, growth, and qPCR.
    Each Eperiment can be polled for absolute and relative amounts of tags.
    
    Experiments are created with the number of cells (int n_cells), 
    multiplicity of infection (float MOI), a dictionary of tag:(p_inv, p_bind) representing the invading strains,
    and (float) fraction_recovered, the fraction of invading bacteria that are counted in the end.
    (float) p_coinvade the probability that a non-invading bacterium binding a cell invades if the cell in invaded
    by another bacterium
    """
    
    def __init__(self, n_cells, MOI, innoculum_dict, fraction_recovered, doCoinvasion=False):
        
        self.n_cells = n_cells
        self.MOI = MOI
        self.innoculum_dict = innoculum_dict
        self.fraction_recovered = fraction_recovered 
        self.doCoinvasion = doCoinvasion
        
        self.cells = []
        self.innoculum = []
        self.bound_bacteria_counts = {}
        self.n_inv_cells = 0
        self.invaded_tags = ''
        self.recovered_tags = ''
        
        self._createCells()
        self._createInnoculumAndInfect()
        self._recoverTags()
    
    def _createCells(self):
        
        for i in range(self.n_cells):
            self.cells.append(Cell())
    
    def _createInnoculumAndInfect(self):
        
        #There can be unequal amounts of all strains, as determined by composition
        
        for tag in (self.innoculum_dict.keys()):
            bact_per_tag = int(self.n_cells * self.MOI * self.innoculum_dict[tag][3])
            p_inv = self.innoculum_dict[tag][0]
            p_bind = self.innoculum_dict[tag][1]
            p_coinv = self.innoculum_dict[tag][2]
            
            for i in range(bact_per_tag):
                b = Bacterium(tag, p_inv, p_bind, p_coinv)
                self.innoculum.append(b)
                
                if p_bind == 1.0:
                    random.choice(self.cells).bindBacterium(b)
                else:
                    if (p_bind > random.random()):
                        random.choice(self.cells).bindBacterium(b)
                
    
    def _recoverTags(self):
        
        for cell in self.cells:
            n_b = cell.getNumberOfBoundBacteria() 
            self.bound_bacteria_counts[n_b] = self.bound_bacteria_counts.get(n_b, 0) + 1
            
            if self.doCoinvasion:
                cell.doCoinvasion()
            
            if cell.isInvaded():
                self.n_inv_cells += 1
                self.invaded_tags += cell.getInvadingTags(active=True, passive=True)   
        
        if (self.fraction_recovered != 1):
            
            recovered_tags = random.sample(self.invaded_tags, int(len(self.invaded_tags)*self.fraction_recovered))
            self.recovered_tags = ''.join(recovered_tags)
        
        else:
            self.recovered_tags = self.invaded_tags
        
    def getNumberOfRecoveredTags(self):
        
        return len(list(self.getRecoveredTagCounts()))
    
    def getRecoveredTags(self):
        
        return self.recovered_tags
    
    def getRecoveredTagCounts(self):
        
        return Counter(self.getRecoveredTags())
    
    def getRecoveredTagFractions(self, precision=2):
        """
        precision is the rounding applied to the fraction, defaults to 2 significant digits
        Returns a dict containing the relative fractions of the innoculumtags in the recovered fraction
        """
        
        tagcounter = self.getRecoveredTagCounts()
        out = {}
        ntags = sum(tagcounter.values())
        
        for tag in self.innoculum_dict.keys():
            try:
                out[tag] = round(tagcounter[tag]/ntags,precision)
            except ZeroDivisionError:
                out[tag] = 0
                
        return out
    
    def getRecoveredTagFractionsNormToInnoculum(self):
        """
        precision is the rounding applied to the fraction, defaults to 2 significant digits
        Returns a dict containing the relative fractions of the innoculumtags in the recovered fraction
        """
        
        tagcounter = self.getRecoveredTagCounts()
        out = {}
        
        #ntags = sum(tagcounter.values())
        
        for tag in self.innoculum_dict.keys():
            bact_per_tag = int(self.n_cells * self.MOI * self.innoculum_dict[tag][3])
            
            try: 
                out[tag] = tagcounter[tag]/float(bact_per_tag)
            except ZeroDivisionError:
                out[tag] = 0
                
        return out




def createInnoculumDict(simpleInnoculumDict, sd_tag_composition):
    """
    simpleInnoculumDict: {'tag':(p_invade, p_bind, p_coinvade)}
    sd_tag_composition: (float) standard deviation of compositon

    It is possible for the tag_fraction to exceed 1

    returns: innoculumDict {'tag':(p_invade, p_bind, p_coinvade, tag_fraction)
    """
    out = {}
    mean_fraction = 1/len(simpleInnoculumDict.keys())
    
    for tag in simpleInnoculumDict.keys():
        fuzzy_fraction = max(0, random.gauss(mean_fraction, sd_tag_composition)) #Makes sure the fraction is not negative
        out[tag] = simpleInnoculumDict[tag] + (fuzzy_fraction,)
        
    return out


def worker_process(moi):
    
    simpleInnoculum = {'A':(p_invade, p_bind, p_coinvade), 
                      'B':(p_invade, p_bind, p_coinvade),
                      'C':(p_invade, p_bind, p_coinvade),
                      'D':(p_invade, p_bind, p_coinvade),
                      'E':(p_invade, p_bind, p_coinvade),
                      'F':(p_invade, p_bind, p_coinvade),
                      'G':(p_invade, p_bind, p_coinvade)
                      }
      
    innoculum = createInnoculumDict(simpleInnoculum, std_tag_composition)
    nr_cells = random.randint(125000, 175000)
    e = Experiment(n_cells=nr_cells, MOI=moi, innoculum_dict=innoculum, fraction_recovered=1, doCoinvasion=False)
    recovered_tags = e.getNumberOfRecoveredTags()
    return recovered_tags 

    
random.seed(666)        



std_tag_composition = 0.043
p_coinvade = 0
p_bind = 1.0
p_invade_invG = 1.018E-4
p_invade_wt = 1.146E-2

probs = [p_invade_wt, p_invade_invG]
MOIs = (200, 20, 2, 0.2, 0.02, 0.002)
nr_runs = (1000, 1000, 1000, 1000, 1000, 1000)

MOIs = (200, 20, 2, 0.2, 0.02, 0.002)
nr_runs = (100, 100, 1000, 1000, 1000, 1000)




if __name__== '__main__':
    
    t0 = time.time()
    f = open("invasion_simulation_results.csv", "w")
    header = "MOI,nr_tags_recovered,std_nr_tags,p_invade,nr_runs,time_to_process\n"
    f.write(header)
    for prob in probs:
        p_invade = prob
        with multiprocessing.dummy.Pool(6) as p:
            for i in range(len(MOIs)):
                print("Starting MOI:{}".format(MOIs[i]))
                t1 = time.time()
                arglist = [MOIs[i]]*nr_runs[i]
                results = p.map(worker_process, arglist)
                means = numpy.mean(results)
                stds = numpy.std(results)
                print("Parallel threads took {0:.1f} for MOI:{1} at pInv: {5} and {2} runs mean: {3} sd: {4}".format(time.time()-t1, MOIs[i], nr_runs[i], means, stds, p_invade))
                f.write(str(MOIs[i]) + "," +
                        str(means) +"," +
                        str(stds) + "," +
                        str(p_invade) + "," + 
                        str(nr_runs[i]) + "," +
                        str(time.time()-t1) + "\n")
        
    print("Parallel threads took de toto {0:.1f}".format(time.time()-t0))
    f.close()


    
#results = []
#while not q.empty():
#    results.append(q.get())

#print(results)


#print(results)



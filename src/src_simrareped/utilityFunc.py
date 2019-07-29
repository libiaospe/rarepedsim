#!/usr/bin/env python
# $File: srv_batch_simrareped.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the RarePedSim program
# Copyright (c) 2013-2015, Biao Li <libiaospe@gmail.com, biaol@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#
# Author: Biao Li
# purpose: utility functions for simRarePed program


import copy, random, sys, progressbar
import os, logging, math, time, warnings

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

import numpy as np
np.seterr(divide='raise', invalid='raise')
from scipy.stats import norm as gaussian

# This is necessary for cPickle to work correctly
import src_simrareped.gdata as gdata
sys.modules['gdata'] = gdata
from gdata import GFile
import itertools, pickle, collections

import simRarePed

try:
    from spower.calculator.manager import PowerData
    from spower.simulator.sampler import Sample
    if sys.version_info.major == 2:
        import spower.simulator.loci_py2 as L
        from spower.gsl_py2 import *
    else:
        import spower.simulator.loci_py3 as L
        from spower.gsl_py3 import *
except Exception:
    logging.warning("Fail to import spower")
    
import tarfile, shutil, glob, zipfile
try:
    import zlib
    ZIPMODE = zipfile.ZIP_DEFLATED 
except:
    ZIPMODE = zipfile.ZIP_STORED

import operator, bisect
from joblib import Parallel, delayed


class restoreStruct():
    '''
    obtain generation information for each individual in the family
    and restore family structure
    '''
    def __init__(self, familyID, familyDict, traitType="Qualitative", addMissingParents=False):
        '''
        if 'addMissingParents' is True: add missing parental individuals whose IDs begin with 1*, 2*, ...
        traitType - choose between 'Qualitative' and 'Quantitative'
        '''
        self.familyID = familyID
        self.familyDict = familyDict
        self.traitType = traitType
        self.missingnessCoding = ('na', 'n/a', 'none', 'unknown')
        self.maleCoding = ('1', 'm', 'male')
        self.femaleCoding = ('2', 'f', 'female')
        self.addMissingParents = addMissingParents
        # create a dict and use individual ids as keys
        famInfo = collections.OrderedDict({})
        # check individual IDs, no '0' allowed
        self.checkIndividualIDs(self.familyID, self.familyDict)
        # check father/mother IDs
        self.checkParentIDs(self.familyID, self.familyDict)
        # add missing founder individuals if self.addMissingParents == True
        if self.addMissingParents:
            self._addMissingParents() 
        # check sex (all fathers should be males, mothers should be females, male code to '1' and female to '2')
        self.checkSex(self.familyID, self.familyDict)
        # check and recode phenotypes (all unknown qualitative traits code to '0' and unknown quantitative traits to 'NA')
        self.checkPhenotypes()
        # retrieve family related info
        for idx, indID in enumerate(self.familyDict['individualID']):
            famInfo[indID] = {
                'indID': indID,
                'indIdx': idx,
                'sex': self.familyDict['sex'][idx],
                'trait': self.familyDict['trait'][idx],
                'fatherID': self.familyDict['fatherID'][idx],
                'motherID': self.familyDict['motherID'][idx],
                'fatherFatherID': self.getFatherFatherID(indID),
                'fatherMotherID': self.getFatherMotherID(indID),
                'motherFatherID': self.getMotherFatherID(indID),
                'motherMotherID': self.getMotherMotherID(indID),
                'childIDs': self.getChildIDs(indID),
                'spouseIDs': self.getSpouseIDs(indID),
                'ancestorIDs':[],
                'offspringIDs':[],
                'isFounder': self.checkIfFounder(indID),
                'ifGeno': self.familyDict['ifGeno'][idx]
            } 
        self.famInfo = famInfo
        # reconstruct info about which generation does each individual belong to
        # assuming the most distant common ancestors are at generation 1
        self.addGenInfo(indIDs = self.familyDict['individualID'])
        # add "ancestorIDs:[...]" and "offspringIDs":[...] info to each individual
        # to trace IDs of all ascending and descending relatives of such individual 
        self.addAncestorOffspringInfo()
        return
    
    
    def _addMissingParents(self):
        '''
        add missing founder individuals (with * appended to their indvidual ID)
        '''
        ## 1. single individual
        ## if only a single individual contained in the pedigree, add two virtual parents
        #if len(self.familyDict['individualID']) == 1:
        #    return
        # add zero-ID non-existent parents for occassional scenarios, which include
        indsToAdd = []
        # 2. single ind and sib-pairs
        if list(set(self.familyDict['fatherID'] + self.familyDict['motherID'])) == ['0']:
            numSibs = len(self.familyDict['individualID'])
            self.familyDict['fatherID'] = ['1*'] * numSibs
            self.familyDict['motherID'] = ['2*'] * numSibs
            indsToAdd.append([self.familyID, '1*', '0', '0', '1', '0', False])
            # last element 'False' of indsToAdd indicates that 'ifGeno' == 0 for missing inds
            indsToAdd.append([self.familyID, '2*', '0', '0', '2', '0', False])
            #for idx, ind in enumerate(self.familyDict['individualID']):
        # 3. for individuals contained in a pedigree that has missing parents (father or mother or both)
        # missing parents will be added and labeled by '*' with unknown phenotypes ('0') 
        else:
            flag = 1
            tmp = {}
            for idx, ind in enumerate(self.familyDict['individualID']):
                c = (self.familyDict['fatherID'][idx], self.familyDict['motherID'][idx])
                if c[0] == '0' and c[1] != '0':   # missing father
                    # check if the missing father has already been included
                    if c not in tmp.keys():   # if not, add ind
                        indsToAdd.append([self.familyID, str(flag)+'*', '0', '0', '1', 'NA' if self.traitType == 'Quantitative' else '0', False])
                        tmp.update({c:str(flag)+'*'})    
                        flag += 1
                    # update offspring father info
                    self.familyDict['fatherID'][idx] = tmp[c]
                elif c[0]!= '0' and c[1] == '0':   # missing mother
                    # check if the missing mother has already been included
                    if c not in tmp.keys():   # if not, add ind
                        indsToAdd.append([self.familyID, str(flag)+'*', '0', '0', '2', 'NA' if self.traitType == 'Quantitative' else '0', False])
                        tmp.update({c:str(flag)+'*'})
                        flag += 1
                    # update offspring mother info
                    self.familyDict['motherID'][idx] = tmp[c]
                elif c[0] == '0' and c[1] == '0' and (ind not in self.familyDict['fatherID']+self.familyDict['motherID']): # missing both parents
                    # add father
                    indsToAdd.append([self.familyID, str(flag)+'*', '0', '0', '1', 'NA' if self.traitType == 'Quantitative' else '0', False])
                    # add mother
                    indsToAdd.append([self.familyID, str(flag+1)+'*', '0', '0', '2', 'NA' if self.traitType == 'Quantitative' else '0', False])
                    # update offspring father & mother info
                    self.familyDict['fatherID'][idx] = str(flag)+'*'
                    self.familyDict['motherID'][idx] = str(flag+1)+'*'
                    flag += 2
                else:
                    continue
        # add missing inds
        for indInfo in indsToAdd:
            self.familyDict['individualID'].append(indInfo[1])
            self.familyDict['fatherID'].append(indInfo[2])
            self.familyDict['motherID'].append(indInfo[3])
            self.familyDict['sex'].append(indInfo[4])
            self.familyDict['trait'].append(indInfo[5])
            self.familyDict['ifGeno'].append(indInfo[6])
            #print indInfo
        return
    
    
    def checkIndividualIDs(self, familyID, familyDict):
        '''
        check individual IDs to be '0' free and unique
        '''
        if 0 in familyDict['individualID']:
            idx0 = familyDict['individualID'].index(0)
            raise ValueError("Individual ID cannot be 0 --> individual #%d in family %s" % (idx0+1, str(familyID)))
        # check uniqueness
        if len(set(familyDict['individualID'])) < len(familyDict['individualID']):
            raise ValueError("Duplicated individual IDs found in family %s" % familyID)
        return
    
    
    def checkParentIDs(self, familyID, familyDict):
        '''
        convert each non-zero but non-existent father/mother ID into '0'
        '''
        for idx, ind in enumerate(self.familyDict['individualID']):
            fID, mID = self.familyDict['fatherID'][idx], self.familyDict['motherID'][idx]
            if fID not in self.familyDict['individualID']:
                self.familyDict['fatherID'][idx] = '0'
            if mID not in self.familyDict['individualID']:
                self.familyDict['motherID'][idx] = '0'
        return 
    
    
    def checkPhenotypes(self):
        '''
        check phenotypes and recode missing values
        '''
        for idx, trait in enumerate(self.familyDict['trait']):
            if self.traitType in ('Mendelian', 'Qualitative'):
                if trait in ['0', '1', '2']:
                    continue
                elif trait.lower() in self.missingnessCoding:
                    self.familyDict['trait'][idx] = '0'
                else:
                    raise ValueError("Erroneous trait value '{}' for individual {} in family {}".format(trait, self.familyDict['individualID'][idx], self.familyID))
            else: # quantitative trait
                try:
                    float(trait)
                    continue
                except ValueError:
                    if trait.lower() in self.missingnessCoding:
                        self.familyDict['trait'][idx] = 'NA'
                    else:
                        raise ValueError("Erroneous trait value {} for individual {} in family {}".format(trait, self.familyDict['individualID'][idx], self.familyID))
        return
        
        
    def checkSex(self, familyID, familyDict):
        '''
        check sex for fathers and mothers 
        '''
        for (fatherID, motherID) in zip(familyDict['fatherID'], familyDict['motherID']):
            if fatherID in familyDict['individualID']:
                fatherIdx = familyDict['individualID'].index(fatherID)
                if str(familyDict['sex'][fatherIdx].lower()) not in  self.maleCoding:
                    raise ValueError("Individual %s in family %s is specified as father but with incorrect gender %s" % (str(fatherID), str(familyID), str(familyDict['sex'][fatherIdx])))
            if motherID in familyDict['individualID']:
                motherIdx = familyDict['individualID'].index(motherID)
                if str(familyDict['sex'][motherIdx].lower()) not in self.femaleCoding:
                    raise ValueError("Individual %s in family %s is specified as mother but with incorrect gender %s" % (str(motherID), str(familyID), str(familyDict['sex'][motherIdx])))
        return
    
    
    def _sortIndIDsByGen(self):
        indIDsByGen = {}
        for indID in self.familyDict['individualID']:
            if indIDsByGen.has_key(self.famInfo[indID]['gen']):
                indIDsByGen[self.famInfo[indID]['gen']].append(indID)
            else:
                indIDsByGen[self.famInfo[indID]['gen']] = [indID]
        indIDs_sorted = []
        gens = indIDsByGen.keys()
        gens.sort()
        for gen in gens:
            indIDs_sorted.extend(indIDsByGen[gen])
        return indIDs_sorted
    
    
    def addAncestorOffspringInfo(self):
        '''
        '''
        indIDs_sorted = self._sortIndIDsByGen()
        ## assign each ind's ancestors' IDs in sequential order of generational numbers (top down approach)
        for indID in indIDs_sorted:
            # append ind's parents to ind's ancestorIDs
            # and extend ind's partents' ancestorIDs to ind's ancestorIDs
            fatherID, motherID = self.famInfo[indID]['fatherID'], self.famInfo[indID]['motherID']
            if fatherID is not '0':
                self.famInfo[indID]['ancestorIDs'].append(fatherID)
                self.famInfo[indID]['ancestorIDs'].extend(self.famInfo[fatherID]['ancestorIDs'])
            if motherID is not '0':
                self.famInfo[indID]['ancestorIDs'].append(motherID)
                self.famInfo[indID]['ancestorIDs'].extend(self.famInfo[motherID]['ancestorIDs'])
            self.famInfo[indID]['ancestorIDs'] = list(set(self.famInfo[indID]['ancestorIDs']))
        ## assign each ind's offspringIDs according to ancestorIDs info
        for indID in indIDs_sorted:
            for ancInd in self.famInfo[indID]['ancestorIDs']:
                self.famInfo[ancInd]['offspringIDs'].append(indID)
        return
    
    
    def recodeGenInfo(self, indIDs):
        genInfo = [self.famInfo[i]['gen'] for i in indIDs]
        if min(genInfo) < 1:
            for i in indIDs:
                self.famInfo[i]['gen'] += (1-min(genInfo))
        return   
        
    
    def addGenInfo(self, indIDs):
        '''
        This function adds generational information to self.famInfo dict.
        Most ancestral founder individuals are given 'gen=1' while
        gen=2,3,... are assigned to offspring and married-in founder(s)/spouses(s)
        (it can handle the situation where a single family ID includes mixed pedigree structures <-->
        individuals belonging to same 'family ID' but are not necessarily all related to each other)
        '''
        # assign gen info to relatives of any individual (randomly selected) recursively
        self.indsAddedGenInfo = []
        ind = random.choice(indIDs)
        self.famInfo[ind]['gen'] = 1
        self.indsAddedGenInfo.append(ind)
        self._setGen(ind, gen=1)
        # recode gen info if min(ind['gen']) < 1
        self.recodeGenInfo(self.indsAddedGenInfo)
        # check if all inds have been assigned gen info
        try:
            assert [self.famInfo[i].has_key('gen') for i in self.familyDict['individualID']].count(True) == len(self.familyDict['individualID'])
        except Exception:
            probIDs = list(set(indIDs).symmetric_difference(set(self.indsAddedGenInfo)))
            self.addGenInfo(indIDs = probIDs)
        return
    
    
    def _funcGen(self, relativeID, relativeGen):
        if relativeID != '0' and relativeID not in self.indsAddedGenInfo:
            # assign generation info for indID's relative
            self.famInfo[relativeID]['gen'] = relativeGen
            self.indsAddedGenInfo.append(relativeID)
            # assign gen info for relative's relatives recursively
            self._setGen(indID=relativeID, gen=relativeGen)
        return
    
        
    def _setGen(self, indID, gen=1):
        '''
        set genenration info to direct relatives (parents, spouse and offspring) of a given ind recursively
        until all relatives of that ind have been successfully assigned gen info
        '''
        # father
        fatherID = self.famInfo[indID]['fatherID']
        self._funcGen(relativeID=fatherID, relativeGen=gen-1)
        # mother
        motherID = self.famInfo[indID]['motherID']
        self._funcGen(relativeID=motherID, relativeGen=gen-1)
        # spouses
        spouseIDs = self.famInfo[indID]['spouseIDs']
        [self._funcGen(relativeID=i, relativeGen=gen) for i in spouseIDs]
        # offspring
        childIDs = self.famInfo[indID]['childIDs']
        [self._funcGen(relativeID=i, relativeGen=gen+1) for i in childIDs]
        return    
        
    
    def idxByID(self, indID):
        return self.familyDict['individualID'].index(indID)
    
    
    def isMale(self, indID, maleCoding=['1', 'M', 'male', 'MALE']):
        sex = self.familyDict['sex'][self.idxByID(indID)]
        return True if sex in maleCoding else False
    
    
    def isFemale(self, indID, femaleCoding=['2', 'F', 'female', 'FEMALE']):
        sex = self.familyDict['sex'][self.idxByID(indID)]
        return True if sex in femaleCoding else False
            
    
    def _getIndices(self, element, myList):
        return [i for i, x in enumerate(myList) if x == element]
    
    
    def getFatherID(self, indID):
        if indID not in self.familyDict['individualID']:
            return '0'
        return self.familyDict['fatherID'][self.idxByID(indID)]
    
    
    def getMotherID(self, indID):
        if indID not in self.familyDict['individualID']:
            return '0'
        return self.familyDict['motherID'][self.idxByID(indID)]
    
    
    def checkIfFounder(self, indID):
        if (self.getFatherID(indID), self.getMotherID(indID)) == ('0', '0'):
            return True
        else:
            return False
    
    
    def getFatherFatherID(self, indID):
        fatherID = self.getFatherID(indID)
        return self.getFatherID(fatherID)
    
    
    def getFatherMotherID(self, indID):
        fatherID = self.getFatherID(indID)
        return self.getMotherID(fatherID)
    
    
    def getMotherFatherID(self, indID):
        motherID = self.getMotherID(indID)
        return self.getFatherID(motherID)
    
    
    def getMotherMotherID(self, indID):
        motherID = self.getMotherID(indID)
        return self.getMotherID(motherID)
    
    
    def getChildIDs(self, indID):
        # if male 
        if self.isMale(indID):
            offspringIdxes = self._getIndices(indID, self.familyDict['fatherID'])
        # if female
        elif self.isFemale(indID):
            offspringIdxes = self._getIndices(indID, self.familyDict['motherID'])
        # if unknown
        else:
            return []
        return [x for i,x in enumerate(self.familyDict['individualID']) if i in offspringIdxes]
        
    
    def getSpouseIDs(self, indID):
        '''
        For an individual it is possible to have multiple spouses
        '''
        childIDs = self.getChildIDs(indID)
        # return [] if individual does not have any offspring
        if len(childIDs) == 0:
            return []
        # if male
        if self.isMale(indID):
            tmp = list(set([self.getMotherID(x) for x in childIDs]))
        # if female
        elif self.isFemale(indID):
            tmp = list(set([self.getFatherID(x) for x in childIDs]))
        else:
            return []
        # remove '0' if any
        if '0' in tmp:
            tmp.remove('0')
        return tmp
    
    
def printDict(obj):
    try:
        type(obj) == type({})
        for key, value in zip(obj.keys(), obj.values()):
            print(key)
            print(value)
    except Exception, e:
        raise ValueError("Input object is not a dict")
    
    
def parseConfigFile(fileName):
    '''
    parse configuration file and return a dict of args and their values
    '''
    pars = {}
    lines = readFiles()._readlines(fileName, suffix='conf')
    for l in lines:
        if (not l.startswith('[')) and (not l.startswith('\n')):
            tmp = readFiles()._splitLine(l, sep='=')
            key, value = tmp[0], typeConvert(tmp[1])
            pars[key] = value
    return pars


def calRandSampleSizes(probs, N):
    inds = [0] * len(probs)
    cumuProbs = [sum(probs[:(i+1)]) for i in range(len(probs))]
    for n in range(N):
        rand = random.random()
        for i, p in enumerate(cumuProbs):
            if rand <= p:
                inds[i] += 1
                break
    return inds
        

def selectReps(currentPool, originalPool):
    '''
    Randomly select rep# (0,1,2,...) for each input file source, i, from ramaining reps given by currentPool
    and update currentPool to originalPool if any currentPool[i] == [].
    '''
    repsIdx = [None] * len(originalPool)
    for i,p in enumerate(currentPool):
        tmp = random.choice(p)
        repsIdx[i] = tmp
        p.remove(tmp)
        if len(p) == 0:
            currentPool[i] = copy.deepcopy(originalPool[i])
    return repsIdx, currentPool
    
    
def typeConvert(obj):
    '''
    convert obj in string to its correct python type
    '''
    try:
        return int(obj)
    except:
        try:
            return float(obj)
        except:
            try:
                return eval(obj)
            except:
                return obj
            

def wheelSpin(probs):
    '''
    return which idx is selected according to given probabilities
    '''
    choice = random.random() * sum(probs)
    for i, w in enumerate(probs):
        choice -= w
        if choice < 0:
            return i    
  
    
class readFiles():
    '''
    This class consists of functions that read multiple types of files
    '''
    def __init__(self, path=None, *args, **kwargs): 
        if path is not None:
            os.chdir(path)
    
    def _readlines(self, fileName, suffix):
        '''
        read file and return lines of strings
        '''
        # if filename ends with suffix remove suffix
        if fileName.endswith(suffix):
            fileName = fileName[:-(len(suffix)+1)]
        try:
            fi = open(fileName+'.'+suffix, 'r')
            lines = fi.readlines()
            fi.close()
            return lines
        except Exception:
            raise ValueError("Cannot find %s" % fileName+'.'+suffix)
            
    
    def _splitLine(self, line, sep=' '):
        '''
        split a single-line string by 'sep' and remove the end of line '\n'
        '''
        if sep == 'whitespace':
            return line.split()
        else:
            tmp = line.rstrip('\n')
            return tmp.split(sep)
        
        
    def sfs(self, sfsFileName, verbose, numGenes):
        '''
        read .sfs file and return a dict of dict objs, where each dict obj is a gene with info of 'maf', 'annotation' (selection coefficient in simulated data), whose key is 'genename' ("R*" in simulated data)
        and 'pos' information
        in the order of #gene_name, chr_name, pos, ref (optional), alt (optional), maf and annotation
        '''
        lines = self._readlines(sfsFileName, 'sfs')
        info = collections.OrderedDict({})
        if numGenes != -1:
            flag = 0
        if verbose >= 0:
            logging.info('Reading {} file'.format(os.path.basename(sfsFileName)))
            pbar = progressbar.ProgressBar(widgets=['', ' ', progressbar.Percentage(), ' ', progressbar.Bar('.'), ' ', progressbar.ETA(), ' '], maxval=len(lines)).start()
        geneNames = []
        for idx, l in enumerate(lines):
            tmp = self._splitLine(l, 'whitespace')
            if l.startswith('#'):
                continue
            if len(tmp) == 5:
                r, ch, pos, maf, sel = tmp
                if r in geneNames:
                    info[r]['maf'].append(float(maf))
                    info[r]['annotation'].append(float(sel))
                    info[r]['position'].append(pos)
                    info[r]['chr'].append(ch)
                    info[r]['ref'].append('.')
                    info[r]['alt'].append('.')
                else:
                    if numGenes != -1:
                        flag += 1
                        if flag > numGenes:
                            break
                    info[r] = {'maf':[float(maf)], 'annotation':[float(sel)], 'position':[pos], 'chr':[ch], 'ref':['.'], 'alt':['.']}
                    geneNames.append(r)
            if len(tmp) == 7:
                r, ch, pos, ref, alt, maf, sel = tmp
                if r in geneNames:
                    info[r]['maf'].append(float(maf))
                    info[r]['annotation'].append(float(sel))
                    info[r]['position'].append(pos)
                    info[r]['chr'].append(ch)
                    info[r]['ref'].append(ref)
                    info[r]['alt'].append(alt)
                else:
                    if numGenes != -1:
                        flag += 1
                        if flag > numGenes:
                            break
                    info[r] = {'maf':[float(maf)], 'annotation':[float(sel)], 'position':[pos], 'chr':[ch], 'ref':[ref], 'alt':[alt]}
                    geneNames.append(r)
            if verbose >= 0:
                pbar.update(idx)
        if verbose >= 0:
            pbar.finish()
        if numGenes != -1 and flag < numGenes:
            logging.warning("Input {} file does not have enough number of genes, {} < {}\n".format(os.path.basename(sfsFileName), flag, numGenes))        
        return info        
                    
    
    def config(self, configFileName):
        '''
        read *.conf file and return a dict obj with 'keys' as options and 'values' as values
        '''
        # line starts with '[' for section name
        lines = self._readlines(configFileName, 'conf')
        configInfo = {}
        for l in lines:
            l = l.rstrip('\n')
            if (l.startswith('[') and l.endswith(']')) or l is '':
                continue
            else:
                tmp = self._splitLine(l, '=')
                option, value = tmp[0].strip(), tmp[1].strip()
                configInfo[option] = convertString(value)
        # validate inputs
        # model: choose from LOGIT, PAR, LNR if traitType is Complex
        traitType = configInfo['trait_type']
        if traitType in ['Complex', 'C', 'c']:
            if configInfo['model'] not in ['LOGIT', 'PAR', 'LNR']:
                raise ValueError("Invalid disease model specified in configuration file, choose from LOGIT, PAR and LNR\n")
            #
            if configInfo['rare_only'] not in [True, False]:
                raise ValueError("Invalid option specified for 'rare_only' in configuration file, use 'True' or 'False'\n")
            if configInfo['model'] in ['LOGIT', 'PAR']:
                configInfo['p1'] = 0.5
        return configInfo
    
    
    def _recodeTraitInfo(self, trait):
        '''
        recode 'na, none, null, n/a, missing, unknown' to na, otherwise to numeric value
        '''
        try:
            if trait.lower() in ('na', 'null', 'none', 'missing', 'n/a', 'unknown'):
                return ('na')
            else:
                return trait
        except:
                return trait
    
    
    def ped(self, pedFileName, traitType='Qualitative'):
        '''
        read *.ped file and return a dict obj 'pedInfo{}' with keys as family ids.
        Each value of 'pedInfo{}' is a dict obj with lists of info,
        such as person ids, father ids, mother ids, gender and trait phenotype,
        for all individuals included in the corresponded family.
        Note! Only information contained in the first 7 columns will be parsed and useful.
        TraitType - choose between 'Qualitative' and 'Quantitative'
        Allowed labels for unknown phenotypes are: NA, N/A, None and Unknown. Unknown labels will be converted to '0' for qualitative trait type and 'NA' for quantitative trait, respectively.  
        '''
        lines = self._readlines(pedFileName, 'ped')
        pedInfo = collections.OrderedDict({})
        for idx, l in enumerate(lines): # skip lines of comments in ped file
            if l.startswith('#'):
                continue
            tmp = self._splitLine(l, 'whitespace')
            if len(tmp) == 7: # with ifGeno info in the 7th column
                familyID, individualID, fatherID, motherID, sex, trait, ifGeno = tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], self._recodeTraitInfo(tmp[5]), tmp[6]
                if pedInfo.has_key(familyID):
                    pedInfo[familyID]['individualID'].append(individualID)
                    pedInfo[familyID]['fatherID'].append(fatherID)
                    pedInfo[familyID]['motherID'].append(motherID)
                    pedInfo[familyID]['sex'].append(sex)
                    pedInfo[familyID]['trait'].append(trait)
                    pedInfo[familyID]['ifGeno'].append(False if ifGeno == '0' else True)
                else:
                    pedInfo[familyID] = {'individualID':[individualID],
                        'fatherID':[fatherID],
                        'motherID':[motherID],
                        'sex':[sex],
                        'trait':[trait],
                        'ifGeno':[False if ifGeno == '0' else True]}
            
            elif len(tmp) == 6: # without ifGeno info
                familyID, individualID, fatherID, motherID, sex, trait = tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], self._recodeTraitInfo(tmp[5])
                if pedInfo.has_key(familyID):
                    pedInfo[familyID]['individualID'].append(individualID)
                    pedInfo[familyID]['fatherID'].append(fatherID)
                    pedInfo[familyID]['motherID'].append(motherID)
                    pedInfo[familyID]['sex'].append(sex)
                    pedInfo[familyID]['trait'].append(trait)
                    pedInfo[familyID]['ifGeno'].append(None)
                else:
                    pedInfo[familyID] = {'individualID':[individualID],
                        'fatherID':[fatherID],
                        'motherID':[motherID],
                        'sex':[sex],
                        'trait':[trait],
                        'ifGeno':[None]}
            
            else: # raise error if any ind does not have enough fam info (7 column needed)
                logging.Error("Line {} in {} file need to have 6 or 7 values".format(idx, pedFileName))
                raise ValueError()
        return pedInfo
    
    
    def gdat(self, gdatFileName):
        '''
        read *.gdat file and return a list obj of file names compressed in *.gdat file
        '''
        if gdatFileName.endswith('.gdat'):
            return gdatItems(gdatFileName)
        else:
            return gdatItems(gdatFileName + '.gdat')


class SaveFiles:
    '''
    class with functions to save various types of files
    '''
    def __init__(self, ):
        pass
    
    @staticmethod
    def simulatedFamilyDataReps(fileLinks, data):
        
        for fi, rep in zip(fileLinks, data):
            for ind in rep:
                fi.write(' '.join([str(i) for i in ind]) + '\n')
    
    @staticmethod
    def simulatedData(fileLink, data):
        for ind in data:
            fileLink.write(' '.join([str(i) for i in ind]) + '\n')
        
    

def gdatItems(gdatFile):
    '''
    load gdat file and return list of data items
    '''
    with GFile.open(gdatFile) as g:
        names = g.getnames()
        data = [g.getdata(i) for i in names]
    # add key 'name' to data
    [i.update({'name':n}) for i,n in zip(data, names)]
    g.close()
    return data
    

def convertString(string):
    '''
    convert a string to its original python type
    '''
    # if bool/int/float/list/tuple/dict/None
    try:
        return eval(string)
    except: # if string
        return string
            

class OrganizePedInfo():
    '''
    re-organize pedInfo object (returned from readFiles().ped func) and sort all family samples by
    unique phenotypic pattern 
    '''
    def __init__(self, pedInfo):
        self.types = []
        self.famIDs = []
    
    
    def name(self, ):
        pass
    
        
class BalancedTernary:
    # Represented as a list of 0, 1 or -1s, with least significant digit first.
    '''
    This class has been adapted from Rosettacode.org
    URL: http://rosettacode.org/wiki/Balanced_ternary#Python
    to convert a ternary number to decimal one or vice versa.
    '''
 
    str2dig = {'+': 1, '-': -1, '0': 0} # immutable
    dig2str = {1: '+', -1: '-', 0: '0'} # immutable
    table = ((0, -1), (1, -1), (-1, 0), (0, 0), (1, 0), (-1, 1), (0, 1)) # immutable
 
    def __init__(self, inp):
        if isinstance(inp, str):
            self.digits = [BalancedTernary.str2dig[c] for c in reversed(inp)]
        elif isinstance(inp, int):
            self.digits = self._int2ternary(inp)
        elif isinstance(inp, BalancedTernary):
            self.digits = list(inp.digits)
        elif isinstance(inp, list):
            if all(d in (0, 1, -1) for d in inp):
                self.digits = list(inp)
            else:
                raise ValueError("BalancedTernary: Wrong input digits.")
        else:
            raise TypeError("BalancedTernary: Wrong constructor input.")
 
    @staticmethod
    def _int2ternary(n):
        if n == 0: return []
        if (n % 3) == 0: return [0] + BalancedTernary._int2ternary(n // 3)
        if (n % 3) == 1: return [1] + BalancedTernary._int2ternary(n // 3)
        if (n % 3) == 2: return [-1] + BalancedTernary._int2ternary((n + 1) // 3)
 
    def to_int(self):
        return reduce(lambda y,x: x + 3 * y, reversed(self.digits), 0)
 
    def __repr__(self):
        if not self.digits: return "0"
        return "".join(BalancedTernary.dig2str[d] for d in reversed(self.digits))
 
    @staticmethod
    def _neg(digs):
        return [-d for d in digs]
 
    def __neg__(self):
        return BalancedTernary(BalancedTernary._neg(self.digits))
 
    @staticmethod
    def _add(a, b, c=0):
        if not (a and b):
            if c == 0:
                return a or b
            else:
                return BalancedTernary._add([c], a or b)
        else:
            (d, c) = BalancedTernary.table[3 + (a[0] if a else 0) + (b[0] if b else 0) + c]
            res = BalancedTernary._add(a[1:], b[1:], c)
            # trim leading zeros
            if res or d != 0:
                return [d] + res
            else:
                return res
 
    def __add__(self, b):
        return BalancedTernary(BalancedTernary._add(self.digits, b.digits))
 
    def __sub__(self, b):
        return self + (-b)
 
    @staticmethod
    def _mul(a, b):
        if not (a and b):
            return []
        else:
            if   a[0] == -1: x = BalancedTernary._neg(b)
            elif a[0] ==  0: x = []
            elif a[0] ==  1: x = b
            else: assert False
            y = [0] + BalancedTernary._mul(a[1:], b)
            return BalancedTernary._add(x, y)
 
    def __mul__(self, b):
        return BalancedTernary(BalancedTernary._mul(self.digits, b.digits))
 
 
def main_BalancedTernary():
    a = BalancedTernary("+-0++0+")
    print "a:", a.to_int(), a
 
    b = BalancedTernary(-436)
    print "b:", b.to_int(), b
 
    c = BalancedTernary("+-++-")
    print "c:", c.to_int(), c
 
    r = a * (b - c)
    print "a * (b - c):", r.to_int(), r
 
 
def weighted_choice(weights, num):
    '''
    return random choices (indices) w.r.t. weights
    num - number of random choices needed
    '''
    totals = np.cumsum(weights)
    norm = totals[-1]
    throws = [np.random.rand()*norm for _ in xrange(num)]
    return [np.searchsorted(totals, throw) for throw in throws]
    

def compressSave(fileName, tempFolder, fm='bz2', numJobs=-1, verbose=1):
    '''
    compress files/dirs in tempFolder to fileName in given format and delete tempFolder
    format (fm) - zip, bz2, tar.gz
    numJobs - # of CPUs to run in parallel
    verbose - (0, 1); if 1 - show progress, otherwise - quiet  
    '''
    path = os.path.split(fileName)[0]
    if path != '':
        try:
            os.chdir(path)
        except:
            raise IOError("No such directory: %s" % path)
    geneNames = glob.glob(os.path.join(tempFolder, '*'))
    if fm == 'bz2':
        fiObj = tarfile.open(fileName+'.tar.bz2', 'w:bz2')
        os.chdir(tempFolder)
        # here it can parallel
        for name in geneNames:
        #tar.add(name)
            fiObj.add(os.path.basename(name))
    elif fm == 'gz':
        fiObj = tarfile.open(fileName+'.tar.gz', 'w:gz')
        os.chdir(tempFolder)
        for name in geneNames:
            fiObj.add(os.path.basename(name))
    elif fm == 'zip': 
        fiObj = zipfile.ZipFile(fileName+'.zip', "w", ZIPMODE)
        for name in geneNames:
            zfi = zipfile.ZipInfo(os.path.join(os.path.basename(name), ''))
            zfi.date_time = time.localtime()[:6]
            fiObj.writestr(zfi,'')
            repFiles = glob.glob(os.path.join(tempFolder, name, '*'))
            for rep in repFiles:
                fiObj.write(rep, os.path.join(os.path.basename(name), os.path.basename(rep)))
    else:
        pass
    fiObj.close()
    # remove tempFolder
    shutil.rmtree(tempFolder)
    return


def parallel_bz2(fiObj, name):
    fiObj.add(os.path.basename(name))
    return


def test_compress():
    pass


class JointProbMaps:
    '''
    This class calculates probabilities of join nuclear families (each represented by a ternaryProbabilityMap)
    '''
    def __init__(self, ):
        pass
    
    def func(self, ):
        pass
    
 
class TernaryProbabilityMap:
    '''
    Given total number of individuals in a nuclear
    family find likelihoods (represented by decimal
    numbers converted from ternary values -
    different combinations of possible genotypes of
    all family members) that contribute to the
    probability of each individual's genotype being
    "+" or "-" or "0".  
    '''
    def __init__(self, ):
        pass
        
    
    def createProbMap_Mendelian(self, hapVarFreq, aff, numReps, isFounder=[True, True], P_T=[0, 1, 1], nucIDs=[], familyID='', scale=5):
        '''
        Within the nuclear family return dict with {possible genotype : probability to occur}
        hapVarFreq - haplotype frequency (probabilities of different numbers of variants per haplotype)       
        aff - affection status (list of all individuals)
        For general case (allow for incomplete penetrance):
        P_T - probabilities of being unaffedted/affected
            given genotypes
            ({'+': a, '0': b, '-': c}),
            where a/b/c is probability of being affected given genotype '+'/'0'/'-';
            e.g. for dominant moi with full penetrance P_T = (Pr_- = 0, Pr_0 = 1, Pr_+ = 1)
            fully penetrant recessive mode P_T = (Pr_- = 0, Pr_0 = 0, Pr_+ = 1) 
        isFounder - list of two bool indicating if parent1 and/or parent2 are founders in the entire pedigree 
        '''
        numInds = len(aff)
        ternaryGenos = self._getTernaryGenosNucFam(numInds-2)
        probDict = collections.OrderedDict({})
        # use P_T to calculate P_T_X (can allow incomplete penetrance)
        a, b, c = P_T[2], P_T[1], P_T[0]
        # p - hap freq of var haps; q = hap freq of wild-type haps
        p, q = 1-hapVarFreq['0'], hapVarFreq['0']
        P_T_X = {'0': [1-b,b], '+': [1-a,a], '-': [1-c,c]}
        # 
        def P_X1(geno):
            '''
            Probability of 1st parent's genotype given that it is a founder
            '''
            if geno[0] == '+':
                return p**2
            elif geno[0] == '-':
                return q**2
            else:
                return 2*p*q
        #
        def P_X2(geno):
            '''
            Probability of 2nd parent's genotype given that it is a founder
            '''
            if geno[1] == '+':
                return p**2
            elif geno[1] == '-':
                return q**2
            else:
                return 2*p*q
        #
        # term A
        def termA(geno, pheno):
            prob = 1
            for gi, ti in zip(geno, pheno):
                if ti == 0:
                    continue
                else:
                    prob *= P_T_X[gi][ti-1]
            return prob
        #
        # term B
        def termB(geno):
            if geno[:2] in ['++', '+-', '-+', '--']:
                return 1
            elif geno[:2] in ['+0', '0+', '-0', '0-']:
                return (0.5)**(len(geno[2:]))
            else:
                return (1./3)**(len(geno[2:]))
        #
        # term C
        def termC(geno):
            prob = 1
            if isFounder[0]:
                prob *= P_X1(geno)
            if isFounder[1]:
                prob *= P_X2(geno)
            return prob
        ##
        for g in ternaryGenos:
            g_int = self.ternaryToInt(g)
            # calculate the probablity of observing genotype 'g' given affection status within nuclear family block
            prob = termA(g, aff) * termB(g) * termC(g)
            probDict[g_int] = prob
        try:
            assert sum(probDict.values()) > 0
        except:
            raise ZeroDivisionError("The input pedigree structure, mode of inheritance and/or penetrance incur Mendelian error for recorded phenotypes on individuals %s in family %s" % (';'.join(nucIDs), familyID))
        return probDict
    
    
    def createIndMap(self, numInds, updatedKeys=None):
        '''
        create individual-based map of possible ternary values (in int) for three possible genotypes, '0', '+' and '-'. 
        {ind1_ID:
            {'0':[list of ints that result in '0' for ind1's geno in ternary format],
             '+':[list of ints of '+' geno for ind1],
             '-':[list of ints of '-' geno for ind1]
            },
         ind2_ID:{...},
         ...,
         indN_ID:{...}
        }
        numInds - number of individuals in the nuclear family
        '''
        if updatedKeys is None:
            udpatedKeys = range(numInds)
        ternaryGenos = self._getTernaryGenosNucFam(numInds-2)
        indDict = collections.OrderedDict({})
        [indDict.update({ID:{'+':[], '-':[], '0':[]}}) for ID in updatedKeys]
        for g in ternaryGenos:
            g_int = self.ternaryToInt(g)
            for ID,i in zip(updatedKeys,g):
                indDict[ID][i].append(g_int)
        return indDict
 

    def createIndProbMap(self, indMap, probMap, updatedKeys=None):
        '''
        create individual-based probability map of three possible genotypes
        requires input indMap from self.createIndMap(...) and
        input probMap from self.createProbMap_Mendelian(...)
        '''
        if updatedKeys is None:
            updatedKeys = indMap.keys()
        tmp = {}
        for i, key in zip(indMap.keys(), updatedKeys):
            tmp.update({key:{'0':0, '+':0, '-':0}})
            for j in indMap[i].keys():
                p = sum([probMap[k] for k in indMap[i][j]])
                tmp[key][j] = p
        return tmp
    
    
    def nucConnectIndsProbMap(self, indMap, probMap, connectIndIDs, familyID):
        '''
        Only save genos that have non-zero probabilities with each geno as key and prob as value
        '''
        genoProbMap = {}
        if len(connectIndIDs) == 0:
            return genoProbMap
        keys = [''.join(i) for i in itertools.product('+-0', repeat=len(connectIndIDs))]
        for geno in keys:
            tmp = [indMap[indID][g] for indID, g in zip(connectIndIDs, geno)]
            intersectedTernaryVals = list(set(tmp[0]).intersection(*tmp[1:]))
            genoProb = sum([probMap[val] for val in intersectedTernaryVals])
            if genoProb != 0:
                genoProbMap[geno] = genoProb
        # Mendelian error occurs if genoProbMap is empty
        if len(genoProbMap) == 0:
            raise ValueError("Mendelian error occurs at individual(s) %s in family %s\n" % (','.join(connectIndIDs), str(familyID)))
        return genoProbMap
        
    
    @staticmethod
    def normalize(cutoff):
        def accumulate(iterable, func=operator.add):
            'Return running totals'
            # accumulate([1,2,3,4,5]) --> 1 3 6 10 15
            # accumulate([1,2,3,4,5], operator.mul) --> 1 2 6 24 120
            it = iter(iterable)
            total = next(it)
            yield total
            for element in it:
                total = func(total, element)
                yield total
        #        
        def func(probDict):
            s = sum(probDict.values())
            for key in probDict.keys():
                probDict[key] = probDict[key]/s
            if cutoff > 1: # remove items with cumu freqs < 1/cutoff * 100%
                sortedProb = sorted(probDict.iteritems(), key=operator.itemgetter(1))
                
                partSums = list(accumulate(sorted(probDict.values())))
                idx = bisect.bisect_right(partSums, 1./cutoff)
                for i in xrange(idx):
                    del probDict[sortedProb[i][0]]
        return func
    
   
    @classmethod
    def writeIndMap(numRange = range(3,11)):
        for n in numRange:
            with open('IndMaps/'+str(n)+'.pi', 'w') as fi:
                pickle.dump(self.createIndMap(n), fi)
        
        
    @staticmethod
    def loadIndMap(numInd):
        with open('IndMaps/'+str(numInd)+'.pi', 'r') as fi:
            return pickle.load(fi)


    @staticmethod
    def ternaryToInt(ternaryValue):
        a = BalancedTernary(ternaryValue)
        return a.to_int()
    
    
    def _getTernaryGenosNucFam(self, numOffspring):
        '''
        obtain all possible genotypes of nuclear family members
        and return in ternary format
        '''
        genos = []
        # parental geno = ++
        genos.append('+'*(numOffspring+2))
        # +-
        genos.append('+-' + '0'*numOffspring)
        # -+
        genos.append('-+' + '0'*numOffspring)
        # --
        genos.append('-'*(numOffspring+2))
        # +0
        genos.extend(['+0'+''.join(i) for i in list(itertools.product('+0', repeat=numOffspring))])
        # 0+
        genos.extend(['0+'+''.join(i) for i in list(itertools.product('+0', repeat=numOffspring))])
        # -0
        genos.extend(['-0'+''.join(i) for i in list(itertools.product('-0', repeat=numOffspring))])
        # 0-
        genos.extend(['0-'+''.join(i) for i in list(itertools.product('0-', repeat=numOffspring))])
        # 00
        genos.extend(['00'+''.join(i) for i in list(itertools.product('+-0', repeat=numOffspring))])
        return genos
   

def spowerArgsNData(configInfo, sfsInfo, geneName):
    '''
    parse phenotype-model configuration file and saved to 
    '''
    tmp = copy.deepcopy(sfsInfo[geneName])
    tmp['function_score'] = tmp['annotation']
    tmp['num_variants'] = len(tmp['maf'])
    del tmp['annotation']
    configInfoDict = configInfo
    return Bunch(**configInfoDict), tmp
    

class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        

def genotyping_artifact(powerData):
    '''
    return idices of sites that are missing; probablites of missing calls and error calls for detrimental, protective and neutral sites
    '''
    missingSites = []
    # missing_low_maf
    if powerData.args.missing_low_maf:
        missingSites = [idx for idx,maf in enumerate(powerData.data['maf']) if maf <= powerData.args.missing_low_maf]
    # missing_sites_[type]
    if powerData.args.missing_sites:
        missing_sites_d, missing_sites_p, missing_sites_n = powerData.args.missing_sites, powerData.args.missing_sites, powerData.args.missing_sites
    else:
        missing_sites_d = powerData.args.missing_sites_detrimental if powerData.args.missing_sites_detrimental else 0
        missing_sites_p = powerData.args.missing_sites_protective if powerData.args.missing_sites_protective else 0
        missing_sites_n = powerData.args.missing_sites_neutral if powerData.args.missing_sites_neutral else 0
    causality = powerData.data['direction']
    for idx, c in enumerate(causality):
        if (c == 'd' and random.random() < missing_sites_d) or \
           (c == 'p' and random.random() < missing_sites_p) or \
           (c == 'n' and random.random() < missing_sites_n):
            missingSites.append(idx)
    missingSites = list(set(missingSites))
    # missing_calls_[type]
    if powerData.args.missing_calls:
        missing_calls_d, missing_calls_p, missing_calls_n = powerData.args.missing_calls, powerData.args.missing_calls, powerData.args.missing_calls
    else:
        missing_calls_d = powerData.args.missing_calls_detrimental if powerData.args.missing_calls_detrimental else 0
        missing_calls_p = powerData.args.missing_calls_protective if powerData.args.missing_calls_protective else 0
        missing_calls_n = powerData.args.missing_calls_neutral if powerData.args.missing_calls_neutral else 0
    probMissingCalls = (missing_calls_d, missing_calls_p, missing_calls_n)
    # error_calls_[type]
    if powerData.args.error_calls:
        error_calls_d, error_calls_p, error_calls_n = powerData.args.error_calls, powerData.args.error_calls, powerData.args.error_calls
    else:
        error_calls_d = powerData.args.error_calls_detrimental if powerData.args.error_calls_detrimental else 0
        error_calls_p = powerData.args.error_calls_protective if powerData.args.error_calls_protective else 0
        error_calls_n = powerData.args.error_calls_neutral if powerData.args.error_calls_neutral else 0
    probErrorCalls = (error_calls_d, error_calls_p, error_calls_n)
    return missingSites, probMissingCalls, probErrorCalls


def applyGenotypeArtifact(geno, missingSites, probMissingCalls, probErrorCalls, causality):
    '''
    introduce genotyping artifact to a simulated genotype
    (assuming that marker genotypes/haplotypes are '1-2' coded, where '0' indicates missingness)
    '''
    for idx in range(len(geno)/2):
        # if site is missing:
        if idx in missingSites:
            geno[idx*2] = 0
            geno[idx*2+1] = 0
            continue
        # genotyping error calls
        c = causality[idx]
        pr = probErrorCalls[0] if c == 'd' else probErrorCalls[1] if c == 'p' else probErrorCalls[2]
        if pr > 0:
            if random.random() < pr:
                geno[idx*2] = 2 if geno[idx*2] == 1 else 1
            if random.random() < pr:
                geno[idx*2+1] = 2 if geno[idx*2+1] == 1 else 1
        # genotyping missing calls
        pr = probMissingCalls[0] if c == 'd' else probMissingCalls[1] if c == 'p' else probMissingCalls[2]
        if pr > 0:
            if random.random() < pr:
                geno[idx*2] = 0
            if random.random() < pr:
                geno[idx*2+1] = 0
    return geno


def generatePhenotype(sampleObj, powerdata, genotype=[]):
    '''
    generate and return an individual's phenotype given spower.Sample() instance object,
    disease model args in powerdata and specified genotype (= two given haplotypes if not [], otherwise, generate genotype on the fly)
    '''
    # get genotype
    if genotype == []:
        L.GenotypeGenerator().apply(sampleObj.data)
    else:
        sampleObj.set(haplotype1=genotype[0], haplotype2=genotype[1])
    model = powerdata.data['model']
    # 'LOGIT' model
    if model == 'LOGIT':
        L.ORModel(powerdata.data['odds_ratio']).apply(sampleObj.data)
        L.DiseaseEffectGenerator(powerdata.data['baseline_effect']).apply(sampleObj.data)
        L.DiseaseStatusGenerator().apply(sampleObj.data)
    # population attributable risk model
    elif model == 'PAR':
        L.PARModel(powerdata.data['par'], powerdata.data['PAR_variable']).apply(sampleObj.data)
        L.DiseaseEffectGenerator(powerdata.data['baseline_effect']).apply(sampleObj.data)
        L.DiseaseStatusGenerator().apply(sampleObj.data)
    # QT model
    elif model == 'LNR':
        L.MeanShiftModel(powerdata.data['mean_shift']).apply(sampleObj.data)
        L.QtEffectGenerator().apply(sampleObj.data)
        L.QtValueGenerator().apply(sampleObj.data)
    else:
        raise ValueError("Invalid input of phenotypic model, can only choose from (LOGIT, PAR, LNR), please modify the configuration file\n")
    return sampleObj.get('phenotype')


def getUpdatedMAF(sampleObj, powerdata, phenotype=2):
    '''
    return updated maf (for case samples if phenotype = 2 or control samples if phenotype = 1)
    '''
    # create a temporary sample object
    temp_sampleObj = sampleObj.clone()
    model = powerdata.data['model']
    if model == 'LOGIT':
        L.ORModel(powerdata.data['odds_ratio']).apply(temp_sampleObj.data)
        L.ORGFUpdater(phenotype, powerdata.data['baseline_effect']).apply(temp_sampleObj.data)
    elif model == 'PAR':
        L.PARModel(powerdata.data['par'], powerdata.data['PAR_variable']).apply(temp_sampleObj.data)
        L.PARGFUpdater(phenotype).apply(temp_sampleObj.data)
    else:
        raise ValueError("Invalid input of disease model, can only choose from (LOGIT, PAR)\n")
    return temp_sampleObj.get('maf')


def updateMAF(sampleObj, powerdata, phenotype=2):
    '''
    update maf on sampleObj (for case samples if phenotype = 2 or control samples if phenotype = 1)
    '''
    model = powerdata.data['model']
    if model == 'LOGIT':
        L.ORModel(powerdata.data['odds_ratio']).apply(sampleObj.data)
        L.ORGFUpdater(phenotype, powerdata.data['baseline_effect']).apply(sampleObj.data)
    elif model == 'PAR':
        L.PARModel(powerdata.data['par'], powerdata.data['PAR_variable']).apply(sampleObj.data)
        L.PARGFUpdater(phenotype).apply(sampleObj.data)
    else:
        raise ValueError("Invalid input of disease model, can only choose from (LOGIT, PAR)\n")
    return 


def initPowerdataObj(configInfo, sfsInfo, geneName):
    '''
    initialize a PowerData object for each replicate
    srvGeneRep - gene replicate # to be used in powerdata, randomly sample a gene replicate if None
    '''
    args, srvData = spowerArgsNData(configInfo, sfsInfo, geneName)
    powerdata = PowerData(args, [], srvData)
    powerdata.update_fixed()
    powerdata.update_random()
    return powerdata


def initSampleObj(powerDataObj):
    '''
    initialize a spower.Sample() object using information
    contained in a pre-processed 'PowerData' object
    '''
    sampleObj = Sample()
    sampleObj.seed(0, os.getpid())
    sampleObj.set(**powerDataObj.data)
    return sampleObj


def initSampleAndPowerdata(configFile, srvResFile, srvGeneRep=None):
    '''
    initialize a spower.Sample() object and a PowerData object
    srvGeneRep - gene replicate # to be used in powerdata, randomly sample a gene replicate if None
    '''
    args, srvData = spowerArgsNData(configFile, srvResFile, srvGeneRep)
    powerdata = PowerData(args, [], srvData)
    powerdata.update_fixed()
    powerdata.update_random()
    sampleObj = Sample()
    sampleObj.seed(0, os.getpid())
    sampleObj.set(**powerdata.data)
    return sampleObj, powerdata


def parseGeneInfo(sfsDict, causalVarIdx, causalVarMaf):
    '''
    parse input sfsDict and retrieve info of functional sites
    '''
    infoDict = {'position':[], 'maf':[], 'annotation':[], 'd_idx':[], 'd_maf':[]}
    infoDict['position'] = sfsDict['position']
    infoDict['maf'] = sfsDict['maf']
    infoDict['annotation'] = sfsDict['annotation']
    # d_idx: indices of causal variant sites (those that have sel > 0)
    infoDict['d_idx'] = causalVarIdx
    infoDict['nd_idx'] = [i for i in range(len(infoDict['annotation'])) if i not in causalVarIdx]
    # d_maf: mafs of causal variants
    infoDict['d_maf'] = causalVarMaf
    cumu_dMaf = [sum(infoDict['d_maf'][:i]) for i in range(1, len(infoDict['d_maf'])+1)]
    sum_dMaf = sum(infoDict['d_maf'])
    infoDict['cumuProbs_dMaf'] = [i/sum_dMaf for i in cumu_dMaf]
    return infoDict
    

class CreatePedStructure():
    '''
    This class includes functions to create a ped file
    from user-specified commonly used pedigree structures for complex trait data simulation,
    e.g. case-ctrl, trio, sibpair and nuclear family.
    '''
    def __init__(self, args, sfsDict, spowerData, missingSites, probMissingCalls, probErrorCalls, causality, hapCoding='1-2', nonProbandPheno=False):
        '''
        args - cmd input arguments
        hapCoding - haplotype coding scheme ('1-2' or '0-1')
        nonProbandPheno - bool, if or not to generate phenotypes for non-proband individuals
        '''
        self.args = args
        self.familyID = 1
        self.indID = 1
        self.sfsDict = sfsDict
        self.spowerData = spowerData
        self.missingSites = missingSites
        self.probMissingCalls = probMissingCalls
        self.probErrorCalls = probErrorCalls
        self.causality = causality
        self.hapCoding = hapCoding
        self.nonProbandPheno = nonProbandPheno
        # initiate a sample object
        self.sampleObj = initSampleObj(self.spowerData)
        # get updated MAF for cases and controls
        if self.spowerData.data['model'] in ['LOGIT', 'PAR']:
            # self.model -> type of phenotype model: 'b'-binary and 'q'-quantitative
            self.model = 'b'
            self.caseMaf = getUpdatedMAF(self.sampleObj, self.spowerData, phenotype=2)
            self.ctrlMaf = getUpdatedMAF(self.sampleObj, self.spowerData, phenotype=1)
        else:
            self.model = 'q'
        

    def getPedFile(self, fileName):
        fi = open(fileName, 'w')
        # case-control
        if self.args['casecontrol']:
            numCases, numCtrls = self.args['casecontrol'][0], self.args['casecontrol'][1]
            for _ in xrange(1, numCases+1):
                simRes = self.caseControl('case')
                SaveFiles.simulatedData(fi, simRes)
                self.familyID += 1
            for _ in xrange(1, numCtrls+1):
                simRes = self.caseControl('ctrl')
                SaveFiles.simulatedData(fi, simRes)
                self.familyID += 1
        # trio
        if self.args['trio']:
            for _ in xrange(1, self.args['trio']+1):
                simRes = self.trio()
                SaveFiles.simulatedData(fi, simRes)
                self.familyID += 1
        # sib-pair
        if self.args['sibpair']:
            for _ in xrange(1, self.args['sibpair']+1):
                simRes = self.sibpair()
                SaveFiles.simulatedData(fi, simRes)
                self.familyID += 1
        # nuclear
        if self.args['nuclear']:
            for _ in xrange(1, self.args['nuclear'][0]+1):
                simRes = self.nuclear(self.args['nuclear'][1:])
                SaveFiles.simulatedData(fi, simRes)
                self.familyID += 1
        fi.close()
        return 
    
    
    def caseControl(self, flag='case'):
        '''
        generate case control samples
        '''
        familyDict = {'individualID':[str(self.indID)],
                    'fatherID':['0'],
                    'motherID':['0'],
                    'sex':[str(random.randint(1,2))],
                    'trait':['2']}
        if flag == 'case':
            famInfo = simRarePed.simPed_linkage().simFamily(self.familyID, familyDict, self.sfsDict, maf_proband=self.caseMaf)
        else: # control'
            famInfo = simRarePed.simPed_linkage().simFamily(self.familyID, familyDict, self.sfsDict, maf_proband=self.ctrlMaf)
            famInfo[str(self.indID)]['trait'] = '1'
        self.indID += 1
        return self.outputInfo(self.familyID, famInfo)
    
    
    def trio(self):
        '''
        generate trio samples (with the child being proband)
        ''' 
        return self.nuclear([1,1])
    
    
    def sibpair(self):
        '''
        generate sibpair samples (with one sibling being proband)
        '''
        individualIDs = [str(self.indID), str(self.indID+1)]
        fatherIDs = ['0', '0']
        motherIDs = ['0', '0']
        sex = [str(random.randint(1,2)) for _ in range(2)]
        trait = ['2','0']
        familyDict = {'individualID':individualIDs,
                      'fatherID': fatherIDs,
                      'motherID': motherIDs,
                      'sex': sex,
                      'trait': trait}
        famInfo = famInfo = simRarePed.simPed_linkage().simFamily(self.familyID, familyDict, self.sfsDict, maf_proband=self.caseMaf, nonProbandPheno=self.nonProbandPheno, spowerData=self.spowerData, sampleObj=self.sampleObj)
        self.indID += 2
        return self.outputInfo(self.familyID, famInfo)
    
    
    def nuclear(self, numOffspringRange):
        '''
        generate nuclear family samples (with one child being proband)
        '''
        numOffspring = random.randint(numOffspringRange[0], numOffspringRange[1])
        individualIDs = [str(i) for i in range(self.indID, self.indID+numOffspring+2)]
        fatherIDs = ['0', '0'] + [str(self.indID)] * numOffspring
        motherIDs = ['0', '0'] + [str(self.indID+1)] * numOffspring
        sex = ['1','2'] + [str(random.randint(1,2)) for _ in range(numOffspring)]
        trait = ['0', '0', '2'] + ['0'] * (numOffspring-1)
        familyDict = {'individualID':individualIDs,
                      'fatherID': fatherIDs,
                      'motherID': motherIDs,
                      'sex': sex,
                      'trait': trait}        
        famInfo = simRarePed.simPed_linkage().simFamily(self.familyID, familyDict, self.sfsDict, maf_proband=self.caseMaf, nonProbandPheno=self.nonProbandPheno, spowerData=self.spowerData, sampleObj=self.sampleObj)
        self.indID += (numOffspring+2)
        return self.outputInfo(self.familyID, famInfo)
        
    
    def outputInfo(self, familyID, famInfo):
        '''
        convert simulated 'famInfo' into a list object
        of output information per individual in linkage format
        (func borrowed from simRarePed.simPed_forward.simulatedFamInfo)
        '''
        # note: remove individuals whose IDs end with '*' (missing parents)
        simFam = []
        for indID in famInfo.keys():
            if not indID.endswith('*'):
                pedFormatGeno = list(itertools.chain(*list(itertools.izip(famInfo[indID]['genotype'][0],famInfo[indID]['genotype'][1]))))
                if self.hapCoding == '1-2':
                    pedFormatGeno = [g+1 for g in pedFormatGeno]
                    # introduce genotyping artifact
                    pedFormatGeno = applyGenotypeArtifact(pedFormatGeno, self.missingSites, self.probMissingCalls, self.probErrorCalls, self.causality)
                simFam.append([familyID,
                               indID,
                               '0' if famInfo[indID]['fatherID'].endswith('*') else famInfo[indID]['fatherID'],
                               '0' if famInfo[indID]['motherID'].endswith('*') else famInfo[indID]['motherID'],
                               famInfo[indID]['sex'],
                               famInfo[indID]['trait'],
                              ] + pedFormatGeno
                    )
        return simFam


def writeVCF(vcfFileName, pedFileName, chrInfo, posInfo, refInfo, altInfo, varMafs):
    '''
    write ped file into vcf file
    '''
    fi_ped = open(pedFileName, 'r')
    lines = fi_ped.readlines()
    fi_ped.close()
    samples = [l.strip().split() for l in lines if l.strip() != '']
    numSamples = len(samples)
    #sample_names = ['{}-{}'.format(ind[0], ind[1]) for ind in samples]
    sample_names = ['{}'.format(ind[1]) for ind in samples]
    
    #print samples[0]
    #print len(samples[0])
    #sys.exit()
    varInfo = np.transpose(np.array(samples)[:,6:])
    numVars = len(varInfo)/2
    fi = open(vcfFileName, 'w')
    fi.write("##fileformat=VCFv4.0\n")
    fi.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(sample_names) + '\n')
    for idx in xrange(numVars):
        # skip wild-type sites
        if not ('2' in varInfo[2*idx] or '2' in varInfo[2*idx+1]):
            continue
        fi.write('\t'.join([str(chrInfo[idx]), str(posInfo[idx]), 'V'+str(idx+1), str(refInfo[idx]), str(altInfo[idx]), '.', 'PASS', 'EVSMAF='+str(varMafs[idx]), 'GT']) + '\t')
        for k, (i,j) in enumerate(zip(varInfo[2*idx], varInfo[2*idx+1])):
            if i == '1':
                x = '0'
            elif i == '0':
                x = '.'
            else:
                x = '1'
            if j == '1':
                y = '0'
            elif j == '0':
                y = '.'
            else:
                y = '1'
            if k < numSamples - 1:    
                fi.write(x+'|'+y+'\t')
            else:
                fi.write(x+'|'+y)
        fi.write('\n')
    fi.close()
    return


def calIndProbAffDict(penetrance=None, oddsRatio=None, baselineEffect=None, moi=None, max_vars=3):
    '''
    calculate probabilities of being affected given all possible genotypes of an individual and return in Dict obj
    For Mendelian trait use 'penetrance' - a list of 3 values; for Complex trait use 'oddsRatio' and 'baselineEffect'
    '''
    def _PrFromOR(r, k):
        return (k*r)/(1-k+k*r)
    if penetrance is not None: # a list of three values to denote (P(Aff|++), P(Aff|+-), P(Aff|--))
        return collections.OrderedDict({'00': penetrance[0], '01': penetrance[1], '11':penetrance[2]})
    elif oddsRatio is not None and baselineEffect is not None and moi is not None:
        indProbAffDict = collections.OrderedDict({})
        for geno in itertools.combinations_with_replacement(range(max_vars+1), 2):
            numVar = sum(geno)
            if moi == 'AAR': # Additive across region
                if numVar == 0:
                    odds = 1
                elif geno[0] == 0:
                    odds = oddsRatio
                else:
                    odds = 2 * oddsRatio - 1
            elif moi == 'MAR': # Multiplicative across region
                odds = oddsRatio ** (2-geno.count(0))
            elif moi == 'MAV': # Multiplicative across variant sites
                odds = oddsRatio ** numVar
            else:
                logging.error("WRONG mode of inheritance for Complex trait simulation, choose from 'D', 'R', 'DAR', 'RAR', 'AAR', 'MAR', 'MAV'")
                sys.exit(0)
            p = _PrFromOR(odds, baselineEffect)
            indProbAffDict[str(geno[0])+str(geno[1])] = p
        return indProbAffDict
    else:
        raise ValueError("Need specify penetrance or odds ratio\n")


def RetrieveCausalVarInfo(sfsInfo, configInfo, geneNames, minNumCausalVar, verbose):
    '''
    This func retrieves info about causal var sites, causal var mafs and calculates haplotype frequency for each gene
    '''
    causalVarSites, causalVarMafs, hapVarFreqs = [], [], []
    for geneName in geneNames:
        causalVarSite, max_vars = selectCausalVarSite(sfsInfo[geneName]['maf'], sfsInfo[geneName]['annotation'], configInfo['def_rare'], configInfo['rare_only'], configInfo['proportion_causal'], geneName=geneName, minNum=minNumCausalVar)
        causalVarMaf = [sfsInfo[geneName]['maf'][idx] for idx in causalVarSite]
        # calculate haplotype frequency
        if verbose >= 0:
            logging.info("Calculating gene-level causal haplotype frequency")
        hapVarFreq = calHapNumVarFreq(causalVarMaf, max_vars=max_vars)
        causalVarSites.append(causalVarSite)
        causalVarMafs.append(causalVarMaf)
        hapVarFreqs.append(hapVarFreqs)
    return causalVarSites, causalVarMafs, hapVarFreqs
        
    
def selectCausalVarSite(maf, sel, defRare, rareOnly, propCausal, geneName, minNum=1, defNeutral=None, defProtective=None):
    '''
    select causal "detrimental" variant sites for a gene region and return their corresponded indices
    (mainly used for pheno-to-geno simulation scheme)
    check if remained number of causal variant sites >= minNum
    # Note:
    If defNeutral and defProtective are not specified, any var site with non-zero sel score is treated as potential detrimental site. 
    '''
    if propCausal == None:
        propCausal = 1
    tmpCausalIdx = []
    for idx, s in enumerate(sel):
        if defNeutral and (defNeutral[0] <= s <= defNeutral[1]):
            continue
        if defProtective and (defProtective[0] <= s <= defProtective[1]):
            continue
        if s != 0:
            tmpCausalIdx.append(idx)
    #
    tmpCausalIdx = [idx for idx, s in enumerate(sel) if s != 0]
    rareCausalIdx = [idx for idx in tmpCausalIdx if maf[idx] <= defRare]
    # commonIdx = [idx for idx, f in enumerate(maf) if f >= defRare]
    causalIdx = rareCausalIdx if rareOnly else tmpCausalIdx
    random.shuffle(causalIdx)
    causalIdx = causalIdx[:int(round(propCausal*len(causalIdx)))]
    causalIdx.sort()
    if len(causalIdx) == 0:
        #FIXME!
        logging.warning('Gene {} does not have any causal variant site. '.format(geneName))
        #raise ValueError("Gene {} does not have any causal variant site. Check field info in input *.sfs file and/or 'proportion_causal' in input *.conf file\n".format(geneName))
    if len(causalIdx) < minNum:
        logging.warning("Gene {} results in having only {} causal detrimental variant sites, smaller than expected {}. Thus, maximum allowed number of variants per haplotype being {} will be used for gene {}".format(geneName, len(causalIdx), minNum, len(causalIdx), geneName))
        minNum = len(causalIdx)
    return causalIdx, minNum


def calHapNumVarFreq(maf, max_vars=1):
    '''
    calculate frequency of observing different numbers of causal variants on a haplotype (0,1,2,...)
    '''
    numVarFreq = {}
    cmaf = [1-m for m in maf]
    for num in range(0, max_vars):
        tmp = []
        for mIdx in itertools.combinations(range(len(maf)), num):
            cmIdx = list(set(range(len(maf))).difference(set(mIdx)))
            tmp.append([m for idx, m in enumerate(maf) if idx in mIdx] + [m for idx, m in enumerate(cmaf) if idx in cmIdx])
        numVarFreq[str(num)] = np.sum(np.prod(tmp, axis=1), axis=0)
    return numVarFreq


class ComplexTraitProbabilityMap():
    '''
    '''
    def __init__(self, ):
        pass
    
    def createProbMap(self, hapVarFreq, aff, isFounder, traitType, indProbAffDict=None, meanshift=None, debug=False, _precision=3, _famGenoFormat='tuple'):
        '''
        Within the nuclear family block return dict with {possible genotype : (relative) probability to occur}
        Genotype is represented by a string, where every other characters denote the genotype of each individual,
        (e.g. if there are 4 individuals contained in the nuclear family block, and genotype string is 01021202, then 01, 02, 12, 02 denote their genotypes, respectively)
        hapVarFreq - {0: p_h0, 1: p_h1, ..., n: p_hn, (n+1)+: q = 1-(p_h0+...+p_hn)}
        aff - affection status of nuclear family members (e.g. [0, 2, 1,1,2])
        isFounder - list of two bool indicating if parent1 and/or parent2 are founders in the entire pedigree
        indProbAffDict - (used for complex qualitative trait) probabilities of being affected given all possible genotypes for complex qualitative trait 
                  e.g. {'00': Pr(Aff|00), '01': Pr(Aff|01), ...}
        meanshift - (used for complex quantitative trait)
        traitType - choose between 'Qualitative' and 'Quantitative'
        '''
        # append 1-sum(hapVarFreq) to hapVarFreq
        hapFreq = copy.deepcopy(hapVarFreq)
        hapFreq[str(len(hapFreq.keys()))] = 1-sum(hapFreq.values())
        genos = [''.join(x) for x in list(itertools.combinations_with_replacement([str(i) for i in range(int(max(hapFreq))+1)], 2))]
        if traitType == 'Qualitative':
            assert len(indProbAffDict.keys()) == len(hapFreq)*(len(hapFreq)+1)/2
        else: # 'Quantitative'
            shiftVals = list(set([int(x[0])+int(x[1]) for x in genos]))
            rv = {}
            [rv.update({val:gaussian(val*meanshift,1)}) for val in shiftVals]
            # save and reuse already calculated pdfs
            pdfs = {}
            [pdfs.update({val:{}}) for val in shiftVals]
        # func to calculate probability of observed phenotypes given genotypes
        def calPhenoProb(geno, pheno):
            '''
            probability of observed phenotypes given genotypes of all individuals in the nuclear family block
            '''
            prob = 1
            if traitType == 'Qualitative':
                for gi, ti in zip(geno, pheno):
                    if ti == 0:
                        continue
                    elif ti == 1:
                        prob *= (1 - indProbAffDict[gi])
                    elif ti == 2:
                        prob *= indProbAffDict[gi]
                    else:
                        raise ValueError("Wrong input of affection status {}".format(aff))
            else: # quantitative trait, use ratio of pdf of shifited gaussian distribution to pdf of standard gaussian distribution
                for gi, ti in zip(geno, pheno):
                    if ti == 'NA' or gi == '00':
                        continue
                    else: # a quantitative value with >=1 variant
                        val = int(gi[0])+int(gi[1])
                        if val == 0:
                            continue
                        else:
                            try:
                                pv = pdfs[val][ti]
                            except KeyError:
                                pv = rv[val].pdf(ti)
                                pdfs[val][ti] = pv
                            try:
                                p0 = pdfs[0][ti]
                            except KeyError:
                                p0 = rv[0].pdf(ti)
                                pdfs[0][ti] = p0
                            prob *= pv / p0
            return prob
        #
        numKids = len(aff) - 2
        famGenoProbMap = collections.OrderedDict({})
        for fatherGeno in genos:
            hap1, hap2 = fatherGeno
            prFatherHap = hapFreq[hap1] * hapFreq[hap2] if isFounder[0] else 1
            for motherGeno in genos:
                hap3, hap4 = motherGeno
                prMotherHap = hapFreq[hap3] * hapFreq[hap4] if isFounder[1] else 1
                haps = [hap1, hap2, hap3, hap4]
                uniqueHaps = list(set(haps))
                uniqueHaps.sort()
                # if both parents are homozygotes (hap1==hap2 and hap3==hap4)
                if hap1 == hap2 and hap3 == hap4:
                    kidGeno = uniqueHaps[0]*2 if len(uniqueHaps) == 1 else uniqueHaps[0]+uniqueHaps[1]
                    famGeno = [fatherGeno, motherGeno] + [kidGeno] * numKids
                    prSeg = 1 # probability of random segregation of parental genotypes
                    prPheno = calPhenoProb(famGeno, aff)
                    famGenoProbMap[tuple(famGeno) if _famGenoFormat == 'tuple' else ''.join(famGeno)] = float('%.{}g'.format(_precision) % (prFatherHap * prMotherHap * prSeg * prPheno))
                # if father is homozygote or mother is homozygote
                elif (hap1 == hap2) or (hap3 == hap4):
                    if hap1 == hap2:
                        tmpGeno = [hap1+hap3 if hap1 <= hap3 else hap3+hap1, hap1+hap4 if hap1 <= hap4 else hap4+hap1]
                    else:
                        tmpGeno = [hap1+hap3 if hap1<= hap3 else hap3+hap1, hap2+hap3 if hap2 <= hap3 else hap3+hap2]
                    for kidsGenoIdx in itertools.product((0,1), repeat=numKids):
                        kidsGeno = [tmpGeno[idx] for idx in kidsGenoIdx]
                        famGeno = [fatherGeno, motherGeno] + kidsGeno
                        prSeg = (0.5) ** numKids
                        prPheno = calPhenoProb(famGeno, aff)
                        famGenoProbMap[tuple(famGeno) if _famGenoFormat == 'tuple' else ''.join(famGeno)] = float('%.{}g'.format(_precision) % (prFatherHap * prMotherHap * prSeg * prPheno))             
                # if both parents are heterozygotes (hap1!=hap2, hap3!=hap4)    
                else:
                    # if len(uniqueHaps) == 2
                    if len(uniqueHaps) == 2:
                        tmpGeno = [uniqueHaps[0]*2, uniqueHaps[0]+uniqueHaps[1], uniqueHaps[1]*2]
                        for kidsGenoIdx in itertools.product((0,1,2), repeat=numKids):
                            kidsGeno = [tmpGeno[idx] for idx in kidsGenoIdx]
                            famGeno = [fatherGeno, motherGeno] + kidsGeno
                            count1 = kidsGenoIdx.count(1)
                            prSeg = (0.5 ** count1) * (0.25 ** (numKids-count1))
                            prPheno = calPhenoProb(famGeno, aff)
                            famGenoProbMap[tuple(famGeno) if _famGenoFormat == 'tuple' else ''.join(famGeno)] = float('%.{}g'.format(_precision) % (prFatherHap * prMotherHap * prSeg * prPheno))
                    else:  # if len(uniqueHaps) == 3 or 4
                        if len(uniqueHaps) == 3:
                    
                            # find which hap occurred twice
                            tmp = collections.Counter(haps)
                            starHap = tmp.keys()[tmp.values().index(2)]
                            otherHaps = copy.deepcopy(uniqueHaps)
                            otherHaps.remove(starHap)
                            tmpGeno = [starHap+starHap,
                                starHap+otherHaps[0] if starHap <= otherHaps[0] else otherHaps[0]+starHap,
                                starHap+otherHaps[1] if starHap <= otherHaps[1] else otherHaps[1]+starHap,
                                otherHaps[0]+otherHaps[1]
                            ]
                        else: # if len(uniqueHaps) == 4
                            tmpGeno = [hap1+hap3 if hap1 <= hap3 else hap3+hap1,
                                    hap2+hap3 if hap2 <= hap3 else hap3+hap2,
                                    hap1+hap4 if hap1 <= hap4 else hap4+hap1,
                                    hap2+hap4 if hap2 <= hap4 else hap4+hap2
                            ]
                        for kidsGenoIdx in itertools.product((0,1,2,3), repeat=numKids):
                            kidsGeno = [tmpGeno[idx] for idx in kidsGenoIdx]
                            famGeno = [fatherGeno, motherGeno] + kidsGeno
                            prSeg = (0.25) ** numKids
                            prPheno = calPhenoProb(famGeno, aff)
                            famGenoProbMap[tuple(famGeno) if _famGenoFormat == 'tuple' else ''.join(famGeno)] = float('%.{}g'.format(_precision) % (prFatherHap * prMotherHap * prSeg * prPheno))
        return famGenoProbMap      
                    

    
    



    
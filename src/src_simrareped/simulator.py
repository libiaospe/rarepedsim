#!/usr/bin/env python
# $File: srv_batch_simrareped.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the RarePedSim program
# Copyright (c) 2013-2015, Biao Li <libiaospe@gmail.com, biaol@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#
# Author: Biao Li
# purpose: genotype and phenotype simulators

import random, os, sys, copy, tempfile, tarfile, shutil, glob, logging
from utilityFunc import calRandSampleSizes, readFiles, selectReps, wheelSpin
import numpy as np
try:
    from simuPOP.utils import ProgressBar
except Exception:
    logging.warning("Fail to import simuPOP")
from srv_batch_simrareped import testProgressBarGUI, bz2Save
# This is necessary for cPickle to work correctly
import gdata as gdata
sys.modules['gdata'] = gdata
from gdata import GData, GFile
# spower modules
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


class Phenotype():
    '''
    This class generated phenotypes for individuals
    contained in a family sample provided that
    genotypes are known
    '''
    def __init__(self, spowerdata, sampleObj):
        self.spowerdata = spowerdata
        self.sampleObj = sampleObj
        
    
    def setPheno(self):
        '''
        generate phenotypes given genotypes
        '''
        def operator(famInfo):
            for indID in famInfo.keys():
                pheno = generatePhenotype(self.sampleObj, self.spowerdata, famInfo[indID]['genotype'])
                famInfo[indID]['trait'] = str(int(pheno)) if self.spowerdata.data['model'] in ['LOGIT', 'PAR'] else str(pheno)
            return
        return operator
            


class Genotype():
    '''
    This class contains functional modules to generate genotype for individuals
    within a family structure for the following scenarios:
    <1> from *.ped file
    <2> using pre-defined commonly used pedigree structure
        (e.g. --sibpair, --trio, --nuclear, etc)
    sampled from forward-time simulation with specification of #gens, mating scheme, (using simuPOP).
    <3> simuPOP - specify #gens, mating scheme and #offspring
    Under any scenario, site frequency spectrum information (*.sfs) is required
    to generate genotype for founders
    '''
    def __init__(self, sfsDict, recRate
                # geneLength=5000, denovoMut=False, mutRate=1.8e-8
                 ):
        self.sfsDict = sfsDict
        self.recRate = recRate
        return
          
        
    def setGeno_withoutProband(self):
        '''
        Generate genotype for individuals within a family given by *.ped file without specified proband
        # FIXME! - current version does not include denovo  mutations
        Set genotype for all individuals from gen=1,2 ... to the most recent gen
        (assuming founders have already been taken care of)
        '''
        def operator(famInfo):
            # calculate gens{gen:[indIDs]}
            gens, indIDsByGen = self.indIDsByGen(famInfo)
            # set genotype for each individual by gen info (top-down)
            for gen in gens:
                [famInfo[indID].update({'genotype':self.indGeno_withoutProband(indID, self.sfsDict['maf'])(famInfo)}) for indID in indIDsByGen[gen]]
            return
        return operator
    
    
    def setGeno_withProband(self, probandID, maf_proband):
        '''
        Generate genotype for a family of individuals where a proband
        has been specified
        '''
        def operator(famInfo):
            # calculate gens{gen:[indIDs]}
            gens, indIDsByGen = self.indIDsByGen(famInfo)
            # locate and simulate proband's genotype first
            famInfo[probandID]['genotype'] = [self.generateHaplotype(maf_proband), self.generateHaplotype(maf_proband)]
            # assign two flags 'hapFather', 'hapMother' to indicate if each individual's two haplotypes from father and/or mother have been generated
            for indID in famInfo.keys():
                famInfo[indID]['hapFather'] = False
                famInfo[indID]['hapMother'] = False
            # First trace backward, assign each of proband's diploids to his/her parents, parents' parents, etc...
            self.assignIndGenoBackward(indID=probandID)(famInfo)
            # then top-down (forward-manner) simulate genotype for each individual
            for gen in gens:
                [self.indGeno_withProband(indID, self.sfsDict['maf'])(famInfo) for indID in indIDsByGen[gen]]
            return
        return operator
    
    
    def indGeno_withProband(self, indID, maf):
        '''
        '''
        def operator(famInfo):
            fatherID = famInfo[indID]['fatherID']
            motherID = famInfo[indID]['motherID']
            # if ind need generate two haplotypes
            if not famInfo[indID].has_key('genotype') or len(famInfo[indID]['genotype']) == 0:
                hap1 = self.generateHaplotype(maf) if fatherID is '0' else random.choice(famInfo[fatherID]['genotype'])
                hap2 = self.generateHaplotype(maf) if motherID is '0' else random.choice(famInfo[motherID]['genotype'])
                famInfo[indID]['genotype'] = [hap1, hap2]
            # if ind already has 1 haplotype, only need generate the other one
            elif len(famInfo[indID]['genotype']) == 1:
                if famInfo[indID]['hapFather'] and famInfo[indID]['hapMother']:
                    raise ValueError('ERROR! A haplotype for individual %s has to be inherited from either Father %s or Mother %s, cannot be both' % (indID, fatherID, motherID)) 
                # if existing hap comes from father
                if famInfo[indID]['hapFather']:
                    hap = self.generateHaplotype(maf) if motherID is '0' else random.choice(famInfo[motherID]['genotype'])
                # if exisiting hap comes from mother, vice versa
                elif famInfo[indID]['hapMother']:
                    hap = self.generateHaplotype(maf) if fatherID is '0' else random.choice(famInfo[fatherID]['genotype'])
                # if existing hap comes from population (indidivual is a founder)
                else:
                    hap = self.generateHaplotype(maf)
                famInfo[indID]['genotype'].append(hap)
            # if ind already has both haps         
            elif len(famInfo[indID]['genotype']) == 2:
                pass
            else:
                raise ValueError("ERROR! Check genotype for individual %s" % indID)
            # 
            if not famInfo[indID]['hapFather']:
                famInfo[indID]['hapFather'] = fatherID
            if not famInfo[indID]['hapMother']:
                famInfo[indID]['hapMother'] = motherID
            return
        return operator
                    

    def indIDsByGen(self, famInfo):
        '''
        calculate gen{gen:[indIDs]}
        '''
        gens = list(set([ind['gen'] for ind in famInfo.values()]))
        indIDsByGen = {}
        [indIDsByGen.update({gen:[]}) for gen in gens]
        [indIDsByGen[famInfo[indID]['gen']].append(indID) for indID in famInfo.keys()]
        return gens, indIDsByGen
       
       
    def indGeno_withoutProband(self, indID, maf):
        '''
        '''
        def operator(famInfo):
            fatherID = famInfo[indID]['fatherID']
            motherID = famInfo[indID]['motherID']
            hap1 = self.generateHaplotype(maf) if fatherID == '0' else self.genGametHap(famInfo[fatherID]['genotype'][0], famInfo[fatherID]['genotype'][1])
            hap2 = self.generateHaplotype(maf) if motherID == '0' else self.genGametHap(famInfo[motherID]['genotype'][0], famInfo[motherID]['genotype'][1])
            return [hap1, hap2]
        return operator
    
    
    def genGametHap(self, hap1, hap2):
        '''
        generate a gamet haplotype allowing for recombination
        '''
        if random.random() < self.recRate:
            spot = random.choice(range(len(hap1))[1:])
            hap3 = hap1[:spot] + hap2[spot:]
            hap4 = hap2[:spot] + hap1[spot:]
            return random.choice([hap1, hap2, hap3, hap4])
        else:
            return random.choice([hap1, hap2])
    
    
    def generateHaplotype(self, maf):
        '''
        generate haplotype from maf (list of minor allele frequencies)
        '''
        haplotype = []
        for i in maf:
            if random.random() <= i:
                haplotype.append(1)
            else:
                haplotype.append(0)
        return haplotype
  
    
    def insertHaplotype(self, indID, haplotype):
        '''
        Insert haplotype to indID['genotype']
        if indID does not have genotype yet, insert [haplotype]
        if len(genotype) == 1: append it as the 2nd haplotype
        if len(genotype) == 2: continue
        '''
        def operator(famInfo):
            if famInfo[indID].has_key('genotype'):
                currentGenotype = famInfo[indID]['genotype']
                if len(currentGenotype) in [0,1]:
                    famInfo[indID]['genotype'].append(haplotype)
            else:
                famInfo[indID]['genotype'] = [haplotype]
            return 
        return operator
    
    
    def assignIndGenoBackward(self, indID):
        '''
        Given genotype of indID (probandID), assign two haplotypes randomly to indID's parents,
        repeat this step for indID's parents' parents, and so on...
        '''
        def operator(famInfo):
            try:
                indGeno = famInfo[indID]['genotype']
            except:
                return
            fatherID, motherID = famInfo[indID]['fatherID'], famInfo[indID]['motherID']
            if len(indGeno) == 1:
                self._haploidToParents(indID, fatherID, motherID, indGeno[0])(famInfo)
            elif len(indGeno) == 2:
                self._diploidToParents(indID, fatherID, motherID, indGeno)(famInfo)
            else:
                raise ValueError("Genotype Error on individual %d!" % indID)
            # recursively run this func for ind's parents, ind's parents' parents, etc..
            if fatherID != '0' and famInfo[fatherID].has_key('genotype'):
                self.assignIndGenoBackward(fatherID)(famInfo)
            if motherID != '0' and famInfo[motherID].has_key('genotype'):
                self.assignIndGenoBackward(motherID)(famInfo)
        return operator
                                      
        
    def _haploidToParents(self, indID, fatherID, motherID, haplotype):
        '''
        Randomly assign offspring's haploid to either parent (father or mother)
        '''
        def operator(famInfo):
            # if inherited from father
            if random.random() >= 0.5:
                if fatherID != '0':
                    self.insertHaplotype(fatherID, haplotype)(famInfo)  
                    famInfo[indID]['hapFather'] = fatherID
            # if inherited from mother
            else:
                if motherID != '0':
                    self.insertHaplotype(motherID, haplotype)(famInfo)
                    famInfo[indID]['hapMother'] = motherID
            return
        return operator
    
    
    def _diploidToParents(self, indID, fatherID, motherID, diploid):
        '''
        Randomly assign offspring's genotype (diploid) to parents
        '''
        def operator(famInfo):
            idxes = [0,1]
            random.shuffle(idxes)
            hapFather, hapMother = diploid[idxes[0]], diploid[idxes[1]]
            if fatherID != '0':
                self.insertHaplotype(fatherID, hapFather)(famInfo)
                famInfo[indID]['hapFather'] = fatherID
            if motherID != '0':
                self.insertHaplotype(motherID, hapMother)(famInfo)
                famInfo[indID]['hapMother'] = motherID
            return
        return operator
        
        
        
class mixedPop():
    '''
    Generate population of mixed individuals to mimic population stratification.
    '''
    def __init__(self, inFileNames, probs, popSize, numReps, mode):
        '''
        mode : 'admix' (population admixture) or 'popstrat' (population stratification)
        '''
        self.inFileNames = inFileNames
        self.probs = probs
        self.popSize = popSize
        self.numReps = numReps
        self.mode = mode
    
    def get_data(self, output, variantPool):
        '''
        Return genotype pool and output sfs for mixed population pool
        '''
        dictGenos = {}
        # open output.sfs file to write sfs info
        sfsFile = open(output+'.sfs', 'w')
        # create temporary folder if variant pool is needed
        if variantPool:
            tempFolder = tempfile.mkdtemp()
        sfsInfo, gdatInfo, numRepsList = self.parseInputFiles(self.inFileNames)
        originalRepsPool, currentRepsPool = copy.deepcopy(numRepsList), copy.deepcopy(numRepsList)
        prog = ProgressBar("Creating %d replicates of the mixed population" % self.numReps, self.numReps, gui=testProgressBarGUI())
        for num in range(1, self.numReps+1):
            # get randomized idxes of reps from input pops
            repsIdxes, currentRepsPool = selectReps(currentRepsPool, originalRepsPool)
            if self.mode == 'popstrat':
                numInds = calRandSampleSizes(self.probs, self.popSize)
                # generate data obj (genotype, maf, sel and pos info) of mixed pop
                dataDict = self.generateDataDict_popstrat(numInds, repsIdxes, sfsInfo, gdatInfo, genoDataKey='rep'+str(num))
            elif self.mode == 'admix':
                dataDict = self.generateDataDict_admix(repsIdxes, sfsInfo, gdatInfo, genoDataKey='rep'+str(num))
            else:
                break
            # write sfs info to sfsFile
            self.writeSfsInfo(sfsFile, dataDict)
            # save variant pool if needed
            if variantPool:
                cwd = os.getcwd()
                os.chdir(tempFolder)
                obj = GData(data=dataDict, name=dataDict['name'])
                obj.compress()
                obj.sink('rep'+str(num))
                del obj
                os.chdir(cwd)
            # update progress bar
            prog.update()
        # close *.sfs file
        sfsFile.close()
        # create *.gdat if variant pool is needed
        if variantPool:
            bz2Save(output, tempFolder)
        return 
    
    
    def writeSfsInfo(self, outFile, dataDict):
        '''
        write maf, sel and pos info of one simulated replicate to sfs file
        '''
        rep = dataDict['name'].split('rep')[-1]
        print >> outFile, '#name maf annotation position'
        print >> outFile, '# Replicate #%s' % rep
        for maf, sel, pos in zip(dataDict['maf'], dataDict['annotation'], dataDict['position']):
            print >> outFile, ' '.join(['R'+rep, '%.8f' % maf, '%.8f' % sel, '%d' % pos,])
        return
        
    
    
    def parseInputFiles(self, inputFileNames):
        '''
        Return sfs and gdat info of each input file
        '''
        sfsInfo = {}
        gdatInfo = {}
        numRepsList = []
        for f in inputFileNames:
            # check if inFileName.gdat exists
            if os.path.isfile(f+'.gdat'):
                gdatInfo[f] = readFiles().gdat(f+'.gdat')
                numRepsList.append(range(len(gdatInfo[f])))
                sfsInfo[f] = None
            elif os.path.isfile(f+'.sfs'):
                sys.stdout.write("WARNING: Variant pool file %s not found, genotype will be restored by minor allele frequencies from %s\n" % (f+'.gdat', f+'.sfs'))
                sfsInfo[f] = readFiles().sfs(f+'.sfs')
                numRepsList.append(range(len(sfsInfo[f])))
                gdatInfo[f] = None
            else:
                raise NameError("File %s or %s not found\n" % (f+'.gdat', f+'.sfs'))
        return sfsInfo, gdatInfo, numRepsList
        
        
    def generateDataDict_popstrat(self, numInds, repsIdxes, sfsInfo, gdatInfo, genoDataKey):
        '''
        generate/sample num=numInds individuals from input sfs/gdat file as one replicate of the mixed population
        '''
        genoList = []
        posDict = {}
        posList, selList = [], []
        for idx, (n, inFileName, rep) in enumerate(zip(numInds, self.inFileNames, repsIdxes)):
            try: # use gdat info
                gdat = gdatInfo[inFileName][rep]
                samples = self.sampleInds(n, gdat[gdat['name']])
                # update dict of position information
                self.updatePosDict(idx, gdat['position'])(posDict)
                posList.append(gdat['position'])
                selList.append(gdat['annotation'])
            except Exception, e: # use sfs info
                sfs = sfsInfo[inFileName][rep]
                samples = np.array([genotype(sfs).generateHaplotype(sfs['maf']) for i in range(2*n)])
                # update dict of position information
                self.updatePosDict(idx, sfs['position'])(posDict)
                posList.append(sfs['position'])
                selList.append(sfs['annotation'])
            genoList.append(samples)
        # create dict of selection coefficients
        selDict = self.createSelDict(posDict, posList, selList)
        # create raw np.array of mixed genotype
        rawGenoArray = self.createRawGenotypeArray(posDict, numInds, genoList, posList)
        # trim raw genotype array by keeping only variant sites (remove cols that are all 0's)
        dataDict = self.trimGenotypeArray(rawGenoArray, genoDataKey, posDict, selDict)
        return dataDict
        
    
    def generateDataDict_admix(self, repsIdxes, sfsInfo, gdatInfo, genoDataKey):
        '''
        '''
        maf_or_hapArray = [None]*len(self.inFileNames)
        genoList = []
        posDict = {}
        posList, selList = [], []
        for idx, (inFileName, rep) in enumerate(zip(self.inFileNames, repsIdxes)):
            try: # use gdat info if available
                gdat = gdatInfo[inFileName][rep]
                tmp = gdatInfo[inFileName][rep][gdat['name']]
                maf_or_hapArray[idx] = np.array(tmp.todense())
                # update dict of position information
                self.updatePosDict(idx, gdat['position'])(posDict)
                posList.append(gdat['position'])
                selList.append(gdat['annotation'])
            except: # use sfs info if gdat data not available
                sfs = sfsInfo[inFileName][rep]
                maf_or_hapArray[idx] = sfs['maf']
                # update dict of position information
                self.updatePosDict(idx, sfs['position'])(posDict)
                posList.append(sfs['position'])
                selList.append(sfs['annotation'])
        # create dict of selection coefficients
        selDict = self.createSelDict(posDict, posList, selList)
        # create raw genotype array (np.array) of admixed population
        allPositions = posDict.keys()
        allPositions.sort()
        genoArray = np.zeros((2*self.popSize, len(allPositions)), dtype=np.int8)
        # create temp list of haplotype indexes of available gdat input files
        originalHapIdxes = [[range(len(i))] if type(i) == np.ndarray else None for i in maf_or_hapArray] 
        currentRemainingHapIdxes = copy.deepcopy(originalHapIdxes)
        for i in range(2*self.popSize):
            # determine from which population this haplo comes
            idx = wheelSpin(self.probs)
            # obtain a haplotype (either from gdat haplotype array or sfs-maf)
            # if haplotype if sampled from gdat pool without replacement
            if originalHapIdxes[idx] is not None:
                hapIdx, currentRemainingHapIdxes[idx] = selectReps(currentPool=currentRemainingHapIdxes[idx], originalPool=originalHapIdxes[idx])
                hap = maf_or_hapArray[idx][hapIdx[0], :] 
            # if haplotype if recovered from maf
            else: 
                hap = getHaplotypeFromMaf(maf_or_hapArray[idx])
            for s, var in enumerate(hap):
                if var == 0:
                    continue
                pos = posList[idx][s]
                col_idx = allPositions.index(pos)
                genoArray[i,col_idx] = var
        # trim raw genotype array by keeping only variant sites (remove cols that are all 0's)
        dataDict = self.trimGenotypeArray(genoArray, genoDataKey, posDict, selDict)
        return dataDict
   
    
    def getHaplotypeFromMaf(self, maf):
        '''
        recover a haplotype from maf list
        '''
        haplotype = []
        for i in maf:
            if random.random() <= i:
                haplotype.append(1)
            else:
                haplotype.append(0)
        return haplotype
            
    
    def trimGenotypeArray(self, rawGenoArray, genoDataKey, posDict, selDict):
        '''
        This function removes wild-type sites, updates posDict/selDict, calculates
        updated site-specific minor allele frequencies and
        returns a data dict obj of the mixed population 
        '''
        dataDict = {genoDataKey:[], 'maf':[], 'annotation':[], 'position':[], 'name': genoDataKey}
        allPositions = posDict.keys()
        allPositions.sort()
        varSiteIdxes = []
        def calMaf(tmp):
            return sum(tmp)/float(len(tmp))
        for idx, pos in enumerate(allPositions):
            maf = calMaf(rawGenoArray[:, idx])
            if maf != 0:
                varSiteIdxes.append(idx)
                if maf < 0.5:
                    dataDict['maf'].append(maf)
                else: # swap mutant and wild-type alleles
                    dataDict['maf'].append(1-maf)
                    for i,x in enumerate(rawGenoArray[:,idx]):
                        if x == 0:
                            rawGenoArray[i,idx] = 1
                        elif x == 1:
                            rawGenoArray[i,idx] = 0
                        else:
                            continue
                dataDict['position'].append(pos)
                dataDict['annotation'].append(selDict[pos])
        dataDict[genoDataKey] = rawGenoArray[:, varSiteIdxes] 
        return dataDict
        
    
    
    def createRawGenotypeArray(self, posDict, numInds, genoList, posList):
        '''
        create raw genotype array which contains geno information
        on all possible variant sites derived from
        all input files
        '''
        allPositions = posDict.keys()
        allPositions.sort()
        rangeList = [[0,2*numInds[0]]]
        for i,num in enumerate(numInds[1:]):
            startIdx = rangeList[i][-1]
            endIdx = 2*sum(numInds[:i+2])
            rangeList.append([startIdx, endIdx])
        genoArray = np.zeros((2*self.popSize, len(allPositions)), dtype=np.int8)
        for idx_col, pos in enumerate(allPositions):
            for file_idx in posDict[pos]:
                pos_idx = posList[file_idx].index(pos)
                colVars = genoList[file_idx][:, pos_idx]
                genoArray[range(rangeList[file_idx][0], rangeList[file_idx][1]), idx_col] = colVars
        return genoArray
        
    
    def createSelDict(self, posDict, posList, selList):
        '''
        create dict obj of selection coefficients for the mixed pop, with var positions as keys
        '''
        selDict = {}
        positions = posDict.keys()
        positions.sort()
        for p in positions:
            inFileIdxes = posDict[p]
            selIdxes = [posList[i].index(p) for i in inFileIdxes]
            tmp = [selList[i][j] for i,j in zip(inFileIdxes, selIdxes)]
            selDict[p] = sum(tmp)/len(tmp)
        return selDict
    

    
    def updatePosDict(self, inFileIndex, posInfoToUpdate):
        '''
        update dict obj of variant positions in the mixed pop
        '''
        def func(posDict):
            keys = posDict.keys()
            for pos in posInfoToUpdate:
                if pos in keys:
                    posDict[pos].append(inFileIndex)
                else:
                    posDict[pos] = [inFileIndex]
        return func
    
    
    def sampleInds(self, num, haploData):
        '''
        sample required number of individuals' genotype from a dataset (scipy.sparse matrix format or
        numpy.array format or list format) in gdat and return a numpy.array
        '''
        try:
            data = np.array(haploData.todense())
        except:
            try:
                assert type(haploData) in [np.ndarray, list] 
                data = haploData
            except:
                raise ValueError('Haplotype data are required to store in either numpy.array or scipy.sparse matrix or list format')
        popSize = len(data)/2
        numPops, numSamples = num/popSize, num%popSize
        popIdxes = range(popSize)
        random.shuffle(popIdxes)
        sampleIdxes = popIdxes[:numSamples]
        hapSampleIdxes = [i*2 for i in sampleIdxes] + [(i*2)+1 for i in sampleIdxes]
        hapSampleIdxes.sort()
        samples = data[hapSampleIdxes, ]
        for i in range(numPops):
            samples = np.vstack((samples, data))
        return samples
    


class indPhenotype():
    '''
    simulate an individual's phenotype given haplotype 1&2 and disease model
    This class also handles processing raw input parameters of disease model and
    transferring simulated data to required object
    '''
    def __init__(self, ):
        '''
        input pars are stored as a dict obj
        
        '''
        return
    
    

def spowerArgsNData(configFile, srvResFile, rep=None):
    '''
    parse phenotype-model configuration file and saved to 
    '''
    configInfoDict = readFiles().config(configFile)
    srvInfoDict = parseSrvInfo(srvResFile, rep=rep)
    return Bunch(**configInfoDict), srvInfoDict
    
    
def parseSrvInfo(simuResFile, rep=None):
    '''
    parse simulation result file(s) (*.sfs or *.gdat) and retrieve 'maf', 'position' and 'annotation' information for a single replicate (gene)
    '''
    # if filename ends with suffix remove suffix
    filePrefix = '.'.join(simuResFile.split('.')[:-1]) if '.' in simuResFile else simuResFile
    # use *.sfs if possible
    try:
        srvInfo = readFiles().sfs(filePrefix)
        resDict =  srvInfo[random.choice(range(len(srvInfo)))] if rep is None else srvInfo[rep-1]
        resDict['function_score'] = resDict['annotation']
        resDict['num_variants'] = len(resDict['maf'])
        del resDict['annotation']
        return resDict
    except:
        pass    
    # use *.gdat instead
    try:
        with GFile.open(filePrefix+'.gdat') as g:
            names = g.getnames()
            data = g.getdata(random.choice(names)) if rep is None else g.getdata('rep'+str(rep))
            return {'maf':data['maf'], 'function_score':data['annotation'], 'position':data['position'], 'num_variants':len(data['maf'])}
    except:
        raise ValueError("An error occurs possibly due to:\n1. Cannot find {0}.sfs file or {0}.gdat file.\n2. The gene database does not have enough gene replicates.\n".format(filePrefix))
    

class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
    

def generatePhenotype(sampleObj, powerdata, genotype=[]):
    '''
    generate and return an individual's phenotype given spower.Sample() instance object,
    disease model args in powerdata and 
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
        raise ValueError("Invalid input of disease model, can only choose from (LOGIT, PAR, LNR)\n")
    return sampleObj.get('phenotype')


def updatedMAF(sampleObj, powerdata, phenotype=2):
    '''
    return updated maf (for case samples if phenotype = 2 or control samples if phenotype = 1)
    '''
    model = powerdata.data['model']
    if model == 'LOGIT':
        L.ORGFUpdater(phenotype, powerdata.data['baseline_effect']).apply(sampleObj.data)
    elif model == 'PAR':
        L.PARGFUpdater(phenotype)
    else:
        raise ValueError("Invalid input of disease model, can only choose from (LOGIT, PAR)\n")
    return sampleObj.get('maf')


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

    
        
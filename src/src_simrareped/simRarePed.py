#!/usr/bin/env python
# $File: srv_batch_simrareped.py $
# $LastChangedDate: $
# $Rev: $
# This file is part of the RarePedSim program
# Copyright (c) 2013-2015, Biao Li <libiaospe@gmail.com, biaol@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#
# Author: Biao Li
# Purpose: Creation of a sample pedigree (a family of individuals within several
# generations) with assigned genotype, selection coefficient, recommbination rate
# and affection status for the linkage and association analysis of samples with rare # variants.
#

'''
In this script, the sample pedigree can be either simulated and drawn by
function 'drawPedSample(...)' or specified by a .ped file in MERLIN-ped-file
simuPOP-readable format. 
Genotype of founder individuals are assigned according to .sfs file (site-specific
information of minor allele frequencies and selection coefficients) outputted by
srv_batch_seqpower.py.
This script uses customized or randomized values of positions of recombination
hot spots and recombination rates together with founders' genotypic information
to set genotype for other family members in the pedigree.
A user-defined operator can be used to set phenotype for each individual
based on the genotype information.
The script can also depict the pedigree structure with affection status for
a qualitative trait.
The script returns the full sample pedigree as a list of instance objects of all
individuals with the completion of assigning genotype and phenotype for
each individual.
'''

try:
    import simuOpt
    simuOpt.setOptions(alleleType='long', optimized=True, quiet=True, version='1.0.5')
    import simuPOP as sim
    import simuPOP.sampling as simSample
    from simuPOP.utils import ProgressBar, migrIslandRates
    from simuPOP.sandbox import RevertFixedSites, revertFixedSites, MutSpaceSelector, MutSpaceMutator, MutSpaceRecombinator 
except Exception:
    import logging
    logging.warning("Fail to import simuPOP")
import os, sys, logging, math, time, random, warnings, copy
import utilityFunc, simulator, test_spower
import collections, itertools
from operator import itemgetter
import numpy as np
DEBUG = False

    

class SimPed_random():
    '''
    This class includes functions to create and evolve founder population in forward-time simulation
    from different sources, which can be sfs info, a previously saved simuPOP pop obj or
    ped info, etc. 
    '''
    def __init__(self):
        pass
    
    
    def sfsBased(self, sfsInfo):
        '''
        Create a simuPOP pop object from sfs information, which need to be a dict obj
        with 'maf', 'sel' and 'pos' information
        '''
        pass
    
    
    def _setIndGeno(self, ind, mafInfo, posInfo):
        '''
        set individual's genotype from recovering minor allele frequencies (mafInfo)
        and variant site position information (posInfo)
        '''
        geno = ind.genotype()
        numSites = len(geno)/2
        for maf, pos in zip(mafInfo, posInfo):
            p1, p2 = random.random(), random.random()
            if p1 <= maf:
                geno[pos-1] = 1
            if p2 <= maf:
                geno[pos-1+numSites] = 1
        ind.setGenotype(geno)
    
    
    def createPopAndSetGeno(self, popSize, geneLength, mafInfo, posInfo):
        '''
        create founder population and set genotype for each individual according to
        minor allele frequencies and position information (mafInfo and posInfo)
        '''
        # check gene length >= last position
        assert geneLength >= posInfo[-1]
        # create founder population with each nt site as a locus
        pop = sim.Population(size=[popSize], ploidy=2, loci=[geneLength])
        # set genotype for each individual
        [self._setIndGeno(ind, mafInfo, posInfo) for ind in pop.individuals()]
        return pop
            
    
    
    def evolve(self, pop, ancGens, mutRate, revertFixedSites, savePop='', verbose=0):
        ## add infoFields for sampling pedigrees
        pop.addInfoFields(['ind_id', 'father_id', 'mother_id'])
        pop.setAncestralDepth(ancGens)
        pop.evolve(
            # initiate sex and ind_ids
            initOps = [
                sim.InitSex(),
                sim.IdTagger()
            ],
            preOps = [
                # revert fixed sites if specified 
                sim.IfElse(revertFixedSites, ifOps=[sim.RevertFixedSites()]), 
            ],
            # monogamous mating and tag pedigrees
            matingScheme = sim.MonogamousMating(
                numOffspring = (sim.UNIFORM_DISTRIBUTION, 1, 4),
                ops = [
                    sim.MendelianGenoTransmitter(),
                    sim.IdTagger(),
                    sim.PedigreeTagger()
                ],
            ),
            # mutations
            postOps = [
                sim.SNPMutator(u=mutRate, v=0, output=''),
            ],
            finalOps = [
                sim.IfElse(len(savePop) > 0, ifOps=[sim.SavePopulation(savePop)]),
                # output statistics in verbose mode
                sim.IfElse(verbose > 0, ifOps=[
                    sim.Stat(popSize=True, numOfSegSites=range(geneLength), alleleFreq=range(geneLength)),
                    sim.PyEval(r'"%d %d %.3f %.3f %.3f %.3f %.3f\n" % (popSize, numOfSegSites, 1-alleleFreq[0][0], 1-alleleFreq[1][0], 1-alleleFreq[2][0], 1-alleleFreq[3][0], 1-alleleFreq[4][0])')])          
            ],
            gen =  ancGens
        )
        return pop
    
    
    def drawPedSample(self, pop, numSampler=2, ancGens=None):
        '''
        FIXME: refer to src_simrareped/sampling.py func plotPedigree and src_simrareped/srv_ped.py func plotPed 
        '''
        pass
    

    
class SimPed_linkage():
    '''
    This class includes functions to simulate genotype and phenotype for pedigrees given by linkage (*.ped) file
    '''
    def __init__(self, ):
        pass
    
    
    def simFamily(self, familyID, familyDict, sfsDict, phenoModel=None, maf_proband=None, nonProbandPheno=False, spowerData=None, sampleObj=None):
        '''
        simulate genotype & phenotype for a specified family structure,
        with genotype of founders recovered from sfs, genotype of other family members by segregation,
        phenotype by phenotypic model.
        If a proband is specified (disease model - 2 /affected, quantitative model - a numeric value), simulation begins with
        generating genotype for proband first.
        '''
        ### following funcs moved to utilityFunc.restoreStruct func, hence become deprecated here 
        ## check individual IDs, no '0' allowed
        #self.checkIndividualIDs(familyID, familyDict)
        ## check sex (all fathers should be male, mothers should be female)
        #self.checkSex(familyID, familyDict)
        # restore family structure (with generation, spouse, offspring info, etc.)
        obj = utilityFunc.restoreStruct(familyID, familyDict, addMissingParents=True)
        famInfo = obj.famInfo
        # locate 'proband' if there is any
        probandID = self.locateProband(familyID, familyDict, famInfo)
        # simulate genotype
        objGeno = simulator.Genotype(sfsDict)
        if probandID is None:
            objGeno.setGeno_withoutProband()(famInfo)
        else:
            # FIXME!! - calcluate/pass in maf_proband
            # use sfsDict['maf'] for test purpose now
            objGeno.setGeno_withProband(probandID, maf_proband=sfsDict['maf'])(famInfo)    
        # finish simulating genotype for each individual and save 'genotype' to dict famInfo
        # simulate phenotypes for nonproband individuals if needed
        if nonProbandPheno:
            simulator.Phenotype(spowerData, sampleObj).setPheno()(famInfo)    
            if probandID:
                famInfo[probandID]['trait'] = '2'
        # print genotype of each individual if DEBUG
        if DEBUG:
            for i in famInfo.values():
                if i.has_key('genotype'):
                    print i['indID'], i['genotype'], i['hapFather'], i['hapMother']
                
        return famInfo
    
    
    def locateProband(self, familyID, familyDict, famInfo):
        '''
        Locate 'proband' (whose disease status is 2 for binary traits
        or has a numeric value for quantitative traits.
        # Note! If more than 1 individuals are specified as diseased per family the
        current verison of simrareped will randomly choose one non-founder individual
        as proband if there exists at least one non-founder diseased individual, otherwise
        use a founder individual as proband. 
        Warning will be raised if there are more than 1 proband identified.
        '''
        # warning message
        warningMsg = ""
        # binary trait
        if set(familyDict['trait']) <= set(['0', '1', '2']):
            # count number of cases/probands included
            num = familyDict['trait'].count('2')
            if num >= 1: # choose only one case as proband
                probandIDs = [i for i in familyDict['individualID'] if famInfo[i]['trait']=='2']
                nonFounderProbandIDs = [i for i in probandIDs if (famInfo[i]['fatherID'], famInfo[i]['motherID']) != ('0','0')]
                if nonFounderProbandIDs != []: # if there exists non-founder proband(s)
                    return random.choice(nonFounderProbandIDs)
                else:
                    founderProbandIDs = [i for i in probandIDs if i not in nonFounderProbandIDs]
                    return random.choice(founderProbandIDs)
            else: # num == 0
                return None
    
    
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
    
        
    def checkSex(self, familyID, familyDict, maleCoding=['1', 'M', 'male', 'MALE'], femaleCoding=['2', 'F', 'female', 'FEMALE']):
        '''
        check sex for fathers and mothers 
        '''
        for (fatherID, motherID) in zip(familyDict['fatherID'], familyDict['motherID']):
            if fatherID in familyDict['individualID']:
                fatherIdx = familyDict['individualID'].index(fatherID)
                if str(familyDict['sex'][fatherIdx]) not in  maleCoding:
                    raise ValueError("Individual %s in family %s is specified as father but with incorrect gender %s" % (str(fatherID), str(familyID), str(familyDict['sex'][fatherIdx])))
            if motherID in familyDict['individualID']:
                motherIdx = familyDict['individualID'].index(motherID)
                if str(familyDict['sex'][motherIdx]) not in femaleCoding:
                    raise ValueError("Individual %s in family %s is specified as mother but with incorrect gender %s" % (str(motherID), str(familyID), str(familyDict['sex'][motherIdx])))
        return


class SimPed_forward:
    '''
    This class generates gene-based genotype for individuals in a known
    pedigree structure regardless of phenotypic information
    (equivalently to what "SimPed" does)
    and assign each individual's trait info according to generated genotype and phenotype model
    # Note: 11-17-14: allow recombination to occur at most once within the given gene
    '''
    def __init__(self, familyID, sfsDict, famInfo, missingSites, probMissingCalls, probErrorCalls, causality, recRate):
        # restore family structure (add missing parents)
        self.famInfo = famInfo
        self.familyID = familyID
        self.sfsDict = sfsDict
        self.missingSites = missingSites
        self.probMissingCalls = probMissingCalls
        self.probErrorCalls = probErrorCalls
        self.causality = causality
        self.recRate = recRate
    
    
    def generate(self, spowerdata=None, hapCoding='1-2', ascertain_func=None, genoOnly=False, model=None, ascertainment_qualitative=None, ascertainment_quantitative=None):
        '''
        Return generated family samples 
        '''
        while True:
            # generate genotypes
            objGeno = simulator.Genotype(self.sfsDict, self.recRate)
            objGeno.setGeno_withoutProband()(self.famInfo)
            # generate phenotypes (optional)
            if spowerdata:
                sampleObj = test_spower.initSampleObj(spowerdata)
                objPheno = simulator.Phenotype(spowerdata, sampleObj)
                objPheno.setPheno()(self.famInfo)
                # ascertainment (optional)
                if not callable(ascertain_func) or ascertain_func(self.famInfo, model=model, asc_qualitative=ascertainment_qualitative, asc_quantitative=ascertainment_quantitative):
                    break
                else:
                    continue
            else:
                for ind in self.famInfo.keys():
                    self.famInfo[ind]['trait'] = '0'
                break
        # successfully generated and ascertained pedigree sample
        return self.simulatedFamInfo(hapCoding, genoOnly)
                        
    
    def simulatedFamInfo(self, hapCoding, genoOnly):
        '''
        return linkage/ped format required format per individual from simulated family data
        if genoOnly is True: simFam only contains genotypes of each individual (without having info of first 6 cols of ped file)
        '''
        # note: remove individuals whose IDs end with '*' (missing parents)
        simFam = []
        keys = self.famInfo.keys()
        keys.sort()
        for indID in keys:
            if not indID.endswith('*'):
                pedFormatGeno = list(itertools.chain(*list(itertools.izip(self.famInfo[indID]['genotype'][0],self.famInfo[indID]['genotype'][1]))))
                if hapCoding == '1-2':
                    pedFormatGeno = [g+1 for g in pedFormatGeno]
                    # introduce genotyping artifact
                    pedFormatGeno = utilityFunc.applyGenotypeArtifact(pedFormatGeno, self.missingSites, self.probMissingCalls, self.probErrorCalls, self.causality)
                if genoOnly:
                    simFam.append(pedFormatGeno)
                else:
                    simFam.append([self.familyID,
                                   indID,
                                   '0' if self.famInfo[indID]['fatherID'].endswith('*') else self.famInfo[indID]['fatherID'],
                                   '0' if self.famInfo[indID]['motherID'].endswith('*') else self.famInfo[indID]['motherID'],
                                   self.famInfo[indID]['sex'],
                                   self.famInfo[indID]['trait'],
                                  ] + pedFormatGeno
                        )
        return simFam
    
    
        
class SimPed_backward:
    '''
    This class includes functions and operators to simulate
    genotype in a known pedigree structure with required phenotypic information (qualitative traits only) 
    in the "backward" manner.
    "Backward" means that the algorithm is rooted in numerical computation of Pr(X|T)
    instead of ascertaining desired phenotypic patterns from a large amount
    of randomly generated samples
    Note:! current version only supports simulation of Mendelian trait,
        mode of inheritance can be either 'D'ominant or 'R'ecessive
    '''
    def __init__(self, familyID, familyDict, numReps, scale, hapVarFreq, P_T=None, indProbAffDict=None, meanshift=None, traitType='Mendelian'):
        '''
        Requirements:
        1. for both single individual or any complex pedigree structure use backward manner (through P_T_X) if at least >= 1 individuals' phenotypes is known (otherwise, use func in class SimPed_forward) 
        2. all individuals are related, therefore, if there is missingness, missing individuals are only allowed to be founders
        P_T - penetrance for Mendelian trait
        indProbAffDict - used for complext qualitative trait, probabilities of being affected given all possible genotypes for complex case-ctrl trait e.g. {'00': Pr(Aff|00), '01': Pr(Aff|01), ...}
        meanshift - used for complex quantitative trait, mean shift of a rare detrimental variant
        traitType - choose from 'Mendelian', 'Qualitative' and 'Quantitative' (the latter two are for complex traits)
        '''
        self.numReps = numReps
        self.scale = scale
        self.hapVarFreq = hapVarFreq
        self.maxVar = int(max(hapVarFreq)) + 1
        self.P_T = P_T
        self.indProbAffDict = indProbAffDict
        self.meanshift = meanshift
        self.traitType = traitType
        # verify if traitType is matched with disease model
        if (self.traitType == 'Mendelian' and self.P_T == None) or (self.traitType == 'Qualitative' and self.indProbAffDict == None) or (self.traitType == 'Quantitative' and self.meanshift == None):
            raise ValueError("Check specified trait type and/or disease model")
        # check if ped only contains single individual
        # Note: 09/08/2015 - treat single individual as trio that is missing two parents
        #if len(familyDict['individualID']) == 1: # is single
        #    self.isSingle = True
        #    self.indPheno = int(familyDict['trait'][0])
        #    return
        if False:
            return
        else: # is sibpair or nuclear or large ped
            self.isSingle = False
            # restore family structure (add missing parents - missing parents' IDs end with '*', e.g. '1*', '2*', ...)
            obj = utilityFunc.restoreStruct(familyID, familyDict, traitType=self.traitType, addMissingParents=True)
            self.familyID = familyID
            self.familyDict = familyDict
            self.famInfo = obj.famInfo
            # pIDs - IDs of pairs of parents
            self.pIDs = self._getParentPairIDs()
            # connIDs - IDs of connecting individuals
            self.connIDs = self._getConnectionIndIDs()
            ## check if #pIDs - 1 == #connIDs
            # get detailed info of each nuclear family block
            self.nucInfo = self._createNuclearFamilyInfoDict()
            # sort list of dict objs of nucInfo by value of 'gen'
            self.nucInfo = sorted(self.nucInfo, key=itemgetter('gen'))
            # generate probability map of different genotypic configuration of connecting individuals
            # so that 'nucFamilyInfo' = self.nucInfo, 'connectIndsProbMap' = self.CIProbMap
            # for Mendelian trait
            if self.traitType == 'Mendelian':
                self.generateConnIndsProbMap()
            # for complex trait
            elif self.traitType in ('Qualitative', 'Quantitative'):
                self.generateConnIndsProbMap_complexTrait()
            else:
                raise ValueError("Incorrect trait type, choose from Mendelian, Qualitative and Quantitative")
            return
    
    
    @staticmethod
    def writeToFile(familyID, famInfo, genoDict, fileObj, missingSites, probMissingCalls, probErrorCalls, causality, hapCoding='1-2', configInfo=None, ifFillMissingPheno=False, causalSitesIdx=[], markerGeno=[], markerGeneOrder=2):
        '''
        convert simulated 'famInfo' into a list object
        of output information per individual in linkage format
        and write to file
        (func borrowed from utilityFunc.CreatePedStructure.outputInfo(...))
        markerGeno - given genotype of marker gene (for pairwise of genes in Mendelian setting)
        markerGeneOrder - 1 to be first and 2 to be second
        '''
        # incorporate genotype info into famInfo
        [famInfo[i].update({'genotype':genoDict[i]}) for i in genoDict.keys()]
        # note: remove individuals whose IDs end with '*' (missing parents)
        simFam = []
        keys = famInfo.keys()
        keys.sort()
        for idx, indID in enumerate(keys):
            if not indID.endswith('*'):
                # specify missing phenotype if 'ifFillMissingPheno' is True
                trait = famInfo[indID]['trait']
                
                if trait.lower() == 'na' and configInfo['model'] == 'LNR' and ifFillMissingPheno: # fill missing value of quantitative trait
                    hap1, hap2 = famInfo[indID]['genotype'][0], famInfo[indID]['genotype'][1]
                    cnt = SimPed_backward.numCausalVar(hap1, causalSitesIdx) + SimPed_backward.numCausalVar(hap2, causalSitesIdx)
                    trait = np.random.normal(cnt * configInfo['meanshift_rare_detrimental'], 1, 1)[0]
                    
                elif trait == '0' and ifFillMissingPheno:
                    # fill missing value of binary trait
                    hap1, hap2 = famInfo[indID]['genotype'][0], famInfo[indID]['genotype'][1]
                    # for mendelian trait
                    if configInfo['trait_type'].lower() in ['mendelian', 'm']:
                        penetr = SimPed_backward.ind_penetrance_mendelian(hap1, hap2, causalSitesIdx, configInfo['penetrance'], configInfo['moi'], configInfo['compound_hetero'])
                    # for complex trait
                    else:
                        penetr = SimPed_backward.ind_penetrance_complex(hap1, hap2, causalSitesIdx, configInfo['OR_rare_detrimental'], configInfo['baseline_effect'], configInfo['moi'])
                    
                    trait = '2' if random.random() < penetr else '1'
                    #print familyID, indID, num, prob, flag, trait
                
                else:
                    pass
                #     
                pedFormatGeno = list(itertools.chain(*list(itertools.izip(famInfo[indID]['genotype'][0],famInfo[indID]['genotype'][1]))))
                pedFormatMarkerGeno = markerGeno[idx] if len(markerGeno) == len(keys) else []
                if hapCoding == '1-2':
                    pedFormatGeno = [g+1 for g in pedFormatGeno]
                    # introduce missing sites, genotype missing calls and genotype error calls (genotyping artifact)
                    pedFormatGeno = utilityFunc.applyGenotypeArtifact(pedFormatGeno, missingSites, probMissingCalls, probErrorCalls, causality)
                geno = pedFormatGeno + pedFormatMarkerGeno if markerGeneOrder == 2 else pedFormatMarkerGeno + pedFormatGeno
                if famInfo[indID]['ifGeno'] == False:
                    geno = [0] * len(geno)
                simFam.append([familyID,
                               indID,
                               '0' if famInfo[indID]['fatherID'].endswith('*') else famInfo[indID]['fatherID'],
                               '0' if famInfo[indID]['motherID'].endswith('*') else famInfo[indID]['motherID'],
                               famInfo[indID]['sex'],
                               #famInfo[indID]['trait']
                               trait
                              ] + geno 
                    )
            #if famInfo[indID]['ifGeno'] is False:
            #    print geno
            #    print len(geno)
            #    print type(geno[0]), type(geno[-1])
            #    sys.exit()
            #    
        for ind in simFam:
            fileObj.write(' '.join([str(i) for i in ind]) + '\n')
        return
    
    
    @staticmethod
    def _writeToLineageFile(familyID, famInfo, lineageDict, fileObj):
        info = []
        keys = famInfo.keys()
        keys.sort()
        for idx, indID in enumerate(keys):
            if not indID.endswith('*'):
                info.append([familyID,
                             indID,
                             '0' if famInfo[indID]['fatherID'].endswith('*') else famInfo[indID]['fatherID'],
                             '0' if famInfo[indID]['motherID'].endswith('*') else famInfo[indID]['motherID'],
                             famInfo[indID]['sex'],
                             famInfo[indID]['trait'],
                             lineageDict[indID][0]+1,
                             lineageDict[indID][1]+1
                    ])
        for ind in info:
            fileObj.write(' '.join([str(i) for i in ind]) + '\n')
        return
    
    
    @staticmethod
    def createNewHaplotypes(ternaryGeno, moi, geneInfo, dVarIdx, ifCompoundHeterozygote=False):
        '''
        create a pair of haplotypes for a founder individual
        moi - choose between 'D' and 'R'
        ifCompoundHeterozygote - if 'True', compound dominant or compound recessive
        '''
        # if geno is '++'
        if ternaryGeno == '+': # two causal haps
            if ifCompoundHeterozygote: # compound heterozygote
                # causal vars need to be randomly redrawn
                dVarIdx1 = geneInfo['d_idx'][utilityFunc.wheelSpin(geneInfo['d_maf'])]
                dVarIdx2 = geneInfo['d_idx'][utilityFunc.wheelSpin(geneInfo['d_maf'])]
                hap1 = SimPed_backward.generateHaplotype(geneInfo, includeVars=[dVarIdx1])
                hap2 = SimPed_backward.generateHaplotype(geneInfo, includeVars=[dVarIdx2])
            else: # without compound heterozygote, need two causal haps with given causal variant
                hap1 = SimPed_backward.generateHaplotype(geneInfo, includeVars=[dVarIdx])
                hap2 = SimPed_backward.generateHaplotype(geneInfo, includeVars=[dVarIdx])
        # if geno is '+-'     
        elif ternaryGeno == '0': # one causal and one wild-type hap
            if ifCompoundHeterozygote: # compound heterozygote
                # randomly redraw a causal variant
                dVarIdx = geneInfo['d_idx'][utilityFunc.wheelSpin(geneInfo['d_maf'])]
            hap1 = SimPed_backward.generateHaplotype(geneInfo, includeVars=[dVarIdx])
            hap2 = SimPed_backward.generateHaplotype(geneInfo, excludeVars=geneInfo['d_idx'])
        else: # ternaryGeno = '--', need two wild-type haps
            hap1 = SimPed_backward.generateHaplotype(geneInfo, excludeVars=geneInfo['d_idx'])
            hap2 = SimPed_backward.generateHaplotype(geneInfo, excludeVars=geneInfo['d_idx'])
        return [hap1, hap2]
    
    
    @staticmethod
    def createNewHaplotypes_complexTrait(indGeno, geneInfo):
        '''
        For complex trait create a pair of haplotypes for a founder individual
        '''
        numVarHap1, numVarHap2 = int(indGeno[0]), int(indGeno[1])
        # create hap1
        hap1 = SimPed_backward.generateHaplotype_complexTrait(numVarHap1, geneInfo)
        # create hap2
        hap2 = SimPed_backward.generateHaplotype_complexTrait(numVarHap2, geneInfo)
        return [hap1, hap2]
        
    
    @staticmethod
    def generateHaplotype(geneInfo, includeVars=[], excludeVars=[]):
        '''
        generate a haplotype according to gene info with specified variants to include and/or exclude
        '''
        hap = [0] * len(geneInfo['maf'])
        tmpIndices = list(set(range(len(geneInfo['maf']))).difference(includeVars).difference(excludeVars))
        for idx in tmpIndices:
            if random.random() < geneInfo['maf'][idx]:
                hap[idx] = 1
        # add causal variantss on given sites
        for idx in includeVars:
            hap[idx] = 1
        return hap
    
    
    @staticmethod
    def generateHaplotype_complexTrait(numCausalVar, geneInfo):
        '''
        For complex trait generate a haplotype according to gene info with specified variants to include and/or exclude
        '''
        varIdxIncluded = []
        weights = copy.deepcopy(geneInfo['d_maf'])
        idxes = copy.deepcopy(geneInfo['d_idx'])
        for num in range(numCausalVar):
            tmpIdx = utilityFunc.weighted_choice(weights, 1)[0]
            varIdxIncluded.append(idxes.pop(tmpIdx))
            weights.pop(tmpIdx)
            
            
        if len(set(varIdxIncluded)) != numCausalVar:
            print varIdxIncluded, numCausalVar
            sys.exit()
        
        varIdxExcluded = list(set(geneInfo['d_idx']).difference(varIdxIncluded))
        return SimPed_backward.generateHaplotype(geneInfo, includeVars=varIdxIncluded, excludeVars=varIdxExcluded)
    
    
    @staticmethod
    def pairOfHaplotypes(ternaryGeno, moi=None, geneInfo=None, dVarIdx=None, ifCompoundHeterozygote=False, fatherTernaryGeno=None, fatherHaps=None, motherTernaryGeno=None, motherHaps=None):
        '''
        create an individual's haplotypes given ternary genotype ('+', '0' or '-') 
        if '+': two causal haps; '-': two wildtype haps; '0': first hap being causal and second wildtype
        '''
        if fatherTernaryGeno == None: # is founder, need create new haplotypes
            return SimPed_backward.createNewHaplotypes(ternaryGeno, moi, geneInfo, dVarIdx, ifCompoundHeterozygote)
        else:  # is non-founder, need two haps segregated from parental genotypes according to ternary genotypic origins of trio
            parTernGenos = fatherTernaryGeno + motherTernaryGeno
            # parental genos = ++ or +-
            if parTernGenos in ['++', '+-']:
                return (random.choice(fatherHaps), random.choice(motherHaps))
            # -+ or --
            elif parTernGenos in ['-+', '--']:
                return (random.choice(motherHaps), random.choice(fatherHaps))
            # +0
            elif parTernGenos == '+0':
                idx = 0 if ternaryGeno == '+' else 1
                return (random.choice(fatherHaps), motherHaps[idx])
            # 0+
            elif parTernGenos == '0+':
                idx = 0 if ternaryGeno == '+' else 1
                return (fatherHaps[idx], random.choice(motherHaps))
            # -0
            elif parTernGenos == '-0':
                idx = 0 if ternaryGeno == '0' else 1
                return (random.choice(fatherHaps), motherHaps[idx])
            # 0-
            elif parTernGenos == '0-':
                idx = 0 if ternaryGeno == '0' else 1
                return (fatherHaps[idx], random.choice(motherHaps))
            # 00
            else:
                if ternaryGeno == '+':
                    return (fatherHaps[0], motherHaps[0])
                elif ternaryGeno == '-':
                    return (fatherHaps[1], motherHaps[1])
                else: # 0
                    if random.random() < 0.5:
                        return (fatherHaps[0], motherHaps[1])
                    else:
                        return (motherHaps[0], fatherHaps[1])
    
    
    @staticmethod
    def pairOfHaplotypes_complexTrait(indGeno, geneInfo=None, fatherGeno=None, fatherHaps=None, motherGeno=None, motherHaps=None, recRate=0):
        '''
        For complex trait create an individual's haplotypes given its '0-1-2...' coded genotype
        '''
        if fatherGeno == None: # is founder, need create new haplotypes
            return SimPed_backward.createNewHaplotypes_complexTrait(indGeno, geneInfo)
        else: # is non-founder, need two haps from parental genotypes
            ## first determine if any paternal or maternal gamet is recombined
            ifRec = [False, False]
            if recRate:
                if random.random() < recRate:
                    ifRec[0] = True
                if random.random() < recRate:
                    ifRec[1] = True
            # if no recombination
            if True not in ifRec:
                parHaps = fatherHaps + motherHaps
                parGenos = fatherGeno + motherGeno
                # can have 4 different kid genos from parental haps
                validKidHapIdx = []
                for idx1, idx2 in zip([0,0,1,1], [2,3,2,3]):
                    tmp, hapIdx = (parGenos[idx1]+parGenos[idx2], (idx1, idx2)) if parGenos[idx1] <= parGenos[idx2] else (parGenos[idx2]+parGenos[idx1], (idx2, idx1))
                    if tmp == indGeno:
                        validKidHapIdx.append(hapIdx)
                # randomly choose a possible kid geno
                kidHapIdx = random.choice(validKidHapIdx)
                return [parHaps[kidHapIdx[0]], parHaps[kidHapIdx[1]]], kidHapIdx
            # if either paternal gamet is recombined or maternal or both
            else:
                while True:
                    kidHap1 = SimPed_backward.genGametHap(fatherHaps[0], fatherHaps[1], 1 if ifRec[0] == True else 0)
                    kidHap2 = SimPed_backward.genGametHap(motherHaps[0], motherHaps[1], 1 if ifRec[1] == True else 0)
                    tmp = SimPed_backward.validateHaps_complexTrait(indGeno, kidHap1, kidHap2, geneInfo['d_idx'])
                    if tmp:
                        return tmp
    
    
    @staticmethod
    def validateHaps_complexTrait(indGeno, hap1, hap2, causalSitesIdx):
        indHapCoding = (int(indGeno[0]), int(indGeno[1]))
        hapCausalVarCounts = (SimPed_backward.numCausalVar(hap1, causalSitesIdx), SimPed_backward.numCausalVar(hap2, causalSitesIdx))
        if indHapCoding == hapCausalVarCounts:
            return [hap1, hap2]
        elif indHapCoding[0] == hapCausalVarCounts[1] and indHapCoding[1] == hapCausalVarCounts[0]:
            return [hap2, hap1]
        else:
            return False
    
    
    @staticmethod
    def generateGenotype(nucFamilyInfo, indGenoOrigins, sfsDict, causalVarIdx, causalVarMaf, moi, recRate, ifCompoundHeterozygote=False, map_file=None):
        '''
        generate genotypes of all family members given a randomly sampled configuration of genotypic origin
        NEW! 11-05-14: allow for at most 1 recombination to occur per mating event
        '''
        # use recRate in map_file if specified
        if map_file:
            recRate = SimPed_backward._read_map_file(map_file)
        # first parse gene info (from sfs dict object)
        geneInfo = utilityFunc.parseGeneInfo(sfsDict, causalVarIdx, causalVarMaf)
        # randomly choose a disease causal variant according to maf weights and obtain its index if there exists at least one causal var site, otherwise return None
        dVarIdx = None if geneInfo['d_idx'] == [] else geneInfo['d_idx'][utilityFunc.wheelSpin(geneInfo['d_maf'])]
        # iterate over all nuclear family blocks contained in the pedigree, restore genotypes of all individuals and save into genoDict
        genoDict = {}
        while len(nucFamilyInfo) > 0:
            tmpList = []
            for nuc in nucFamilyInfo:
                #print nuc['IDs']
                # restore parental genotypes in the nuc family
                fID, mID = nuc['IDs'][0], nuc['IDs'][1]
                if (not nuc['isFounder'][0] and not genoDict.has_key(fID)) or (not nuc['isFounder'][1] and not genoDict.has_key(mID)): # when there exists loops in pedigree structure
                    tmpList.append(nuc)
                    continue
                if nuc['isFounder'][0] and not genoDict.has_key(fID):   # if father is founder and not been assigned geno yet, create a pair of new haplotypes
                    genoDict[fID] = SimPed_backward.pairOfHaplotypes(indGenoOrigins[fID], moi, geneInfo, dVarIdx, ifCompoundHeterozygote)
                if nuc['isFounder'][1] and not genoDict.has_key(mID):   # if mother is founder and not been assigned geno yet, create a pair of new haps
                    genoDict[mID] = SimPed_backward.pairOfHaplotypes(indGenoOrigins[mID], moi, geneInfo, dVarIdx, ifCompoundHeterozygote)
                for oID in nuc['IDs'][2:]: #  for each offspring individual
                    if not genoDict.has_key(oID):
                        # if without recombination can also use the following commented out function; if with recombination need to use func 'genOffspringHaplotypes' 
                        genoDict[oID] = SimPed_backward.genOffspringHaplotypes(indGenoOrigins[oID], moi, geneInfo, genoDict[fID], genoDict[mID], dVarIdx, recRate, ifCompoundHeterozygote)
            nucFamilyInfo = tmpList
        return genoDict  
 
    
    @staticmethod
    def _read_map_file(fi):
        '''
        read map file and convert physical distance in cm to recombination rate in probability
        '''
        tmp = np.genfromtxt(fi)
        centi_morgans = tmp[:, 2].tolist()
        rec_rates = [0.01 * (centi_morgans[idx+1] - centi_morgans[idx]) for idx in range(len(centi_morgans)-1)]
        return rec_rates

    
    @staticmethod
    def generateGenotype_complexTrait(nucFamilyInfo, indGenoOrigins, sfsDict, causalVarIdx, causalVarMaf, recRate):
        '''
        For complex trait generate genotypes of all family members
        given a randomly sampled configuratioin of genotypic origin - 'indGenoOrigins'.
        Allow for at most 1 recombination to occur during mating
        This function also tracks lineage of haplotype flow. Eg. for each non-founder, if its father's haplotypes are indexed as 1,2 and mother's as 3,4; keep track of which two haps (by indices) are inherited by the offspring
        '''
        # first pare geneInfo (from sfs dict object)
        geneInfo = utilityFunc.parseGeneInfo(sfsDict, causalVarIdx, causalVarMaf)
        # iterate over all nuclear family blocks contained in the pedigree, restore genotypes of all individuals and save into genoDict
        genoDict = {}
        lineage = {}
        while len(nucFamilyInfo) > 0:
            tmpList = []
            for nuc in nucFamilyInfo:
                fID, mID = nuc['IDs'][0], nuc['IDs'][1]
                if (not nuc['isFounder'][0] and not genoDict.has_key(fID)) or (not nuc['isFounder'][1] and not genoDict.has_key(mID)): # when there exists loops in pedigree structure
                    tmpList.append(nuc)
                    continue
                if nuc['isFounder'][0] and not genoDict.has_key(fID):
                    # if father is founder and not been assigned geno yet, create a pair of new haplotypes
                    genoDict[fID] = SimPed_backward.pairOfHaplotypes_complexTrait(indGenoOrigins[fID], geneInfo)
                    lineage[fID] = (-1,-1)
                if nuc['isFounder'][1] and not genoDict.has_key(mID):
                    # if mother is founder and has not been assigned geno yet, create a pair of new haps
                    genoDict[mID] = SimPed_backward.pairOfHaplotypes_complexTrait(indGenoOrigins[mID], geneInfo)
                    lineage[mID] = (-1,-1)
                for oID in nuc['IDs'][2:]:
                    # for each offspring 
                    if not genoDict.has_key(oID):
                        ## FIXME! (no recombination for now, refer to func genOffspringHaplotypes to add recombination later)
                        genoDict[oID], lineage[oID] = SimPed_backward.pairOfHaplotypes_complexTrait(indGenoOrigins[oID], geneInfo, fatherGeno=indGenoOrigins[fID], fatherHaps=genoDict[fID], motherGeno=indGenoOrigins[mID], motherHaps=genoDict[mID], recRate=recRate)
            nucFamilyInfo = tmpList
        return genoDict, lineage
                
    
    @staticmethod
    def genGametHap(hap1, hap2, recRate):
        '''
        generate a gamet haplotype allowing for recombination
        '''
        if type(recRate) == list: # use map file with site-specific recombination rate (args.map_file for LH-NPL project)
            recSites = [] if random.random() < 0.5 else [idx+1 for idx, p in enumerate(recRate) if random.random() < p]
            if recSites == []: # nonrecombinant gamet
                return random.choice((hap1, hap2))
            else: # recombinant
                recSites = [0] + recSites + [len(hap1)]
                #
                #print recSites
                #
                flag, gamet1, gamet2 = 1, [], []
                for idx in range(len(recSites)-1):
                    if flag == 1:
                        gamet1 += hap1[recSites[idx] : recSites[idx+1]]
                        gamet2 += hap2[recSites[idx] : recSites[idx+1]]
                        flag = 2
                    else:
                        gamet1 += hap2[recSites[idx] : recSites[idx+1]]
                        gamet2 += hap1[recSites[idx] : recSites[idx+1]]
                        flag = 1
                return random.choice((gamet1, gamet2))
        # use args.rec_rate
        else:    
            if random.random() < recRate:
                spot = random.choice(range(len(hap1))[1:])
                hap3 = hap1[:spot] + hap2[spot:]
                hap4 = hap2[:spot] + hap1[spot:]
                return random.choice([hap3, hap4])
            else:
                return random.choice([hap1, hap2])
    
    
    @staticmethod
    def genOffspringHaplotypes(ternaryGeno, moi, geneInfo, fatherHaps, motherHaps, dVarIdx, recRate, ifCompoundHeterozygote):
        '''
        generate offspring haplotypes and allow for recombination (at most once)
        '''
        mapDict = {'+': 2, '0': 1, '-': 0}
        while True:
            offHap1 = SimPed_backward.genGametHap(fatherHaps[0], fatherHaps[1], recRate)
            offHap2 = SimPed_backward.genGametHap(motherHaps[0], motherHaps[1], recRate)
            causality = [SimPed_backward.ifHapCausal(offHap1, [dVarIdx]), SimPed_backward.ifHapCausal(offHap2, [dVarIdx])] if ifCompoundHeterozygote == False  else [SimPed_backward.ifHapCausal(offHap1, geneInfo['d_idx']), SimPed_backward.ifHapCausal(offHap2, geneInfo['d_idx'])]
            if mapDict[ternaryGeno] == causality.count(True):
                if causality[0] >= causality[1]:
                    return (offHap1, offHap2)
                else:
                    return (offHap2, offHap1)
                
                        
    @staticmethod
    def ifHapCausal(hap, causalSitesIdx):
        if causalSitesIdx[0] == None:
            return False
        for idx in causalSitesIdx:
            if hap[idx] == 1:
                return True
        return False
    
    
    @staticmethod
    def numCausalVar(hap, causalSitesIdx):
        count = 0
        for idx in causalSitesIdx:
            if hap[idx] == 1:
                count += 1
        return count
        
    @staticmethod
    def ind_penetrance_mendelian(hap1, hap2, causalSitesIdx, penetrance, moi='D', compound_hetero=False):
        if not compound_hetero: # without compound heterogeneity
            flag = 0
            for idx in causalSitesIdx:
                if hap1[idx] == hap2[idx] == 0:
                    continue
                elif hap1[idx] == hap2[idx] == 1:
                    flag = 2
                    break
                else:
                    flag = 1
            return penetrance[flag]
        else:   # with compound heterogeneity
            numVarHap1, numVarHap2 = SimPed_backward.numCausalVar(hap1, causalSitesIdx), SimPed_backward.numCausalVar(hap2, causalSitesIdx)
            return penetrance[-1 - (numVarHap1, numVarHap2).count(0)]
    
    @staticmethod
    def ind_penetrance_complex(hap1, hap2, causalSitesIdx, oddsRatio, baseline, moi):
        def _PrFromOR(r, k):
            return (k*r)/(1-k+k*r)
        
        if moi in ['D', 'R', 'DAR', 'RAR']:
            pr = _PrFromOR(oddsRatio, baseline)
            penetrance = (baseline, baseline if moi.startswith('R') else pr, pr)
            compound_hetero = False if len(moi) == 1 else True
            return SimPed_backward.ind_penetrance_mendelian(hap1, hap2, causalSitesIdx, penetrance, moi[0], compound_hetero)
        else:
            numVarHap1, numVarHap2 = SimPed_backward.numCausalVar(hap1, causalSitesIdx), SimPed_backward.numCausalVar(hap2, causalSitesIdx)
            count0 = (numVarHap1, numVarHap2).count(0)
            if moi == 'AAR':
                tmp = (1, oddsRatio, 2*oddsRatio-1)
                OR = tmp[-1 - count0]
            elif moi == 'MAR':
                tmp = (1, oddsRatio, oddsRatio**2)
                OR = tmp[-1 - count0]
            else: # moi == 'MAV'
                OR = oddsRatio ** (numVarHap1 + numVarHap2)
            
            return _PrFromOR(OR, baseline)
            
                
        
    
    
        
    #@staticmethod
    #def ifHapCausal_recessive(hap, dVarIdx):
    #    return True if hap[dVarIdx] == 1 else False
    
    
    def _getParentPairIDs(self):
        '''
        Return IDs of parental individuals in pairs.
        '''
        parentPairs = []
        for indID, ind in zip(self.famInfo.keys(), self.famInfo.values()):
            pIDs = [ind['fatherID'], ind['motherID']]
            tmp = pIDs.count('0')
            if tmp == 1:
                raise ValueError("Missing parent not allowed --> check individual %s in family %s or rerun command by switching on 'addMissingParents' in class restoreStruct(...) in utilityFunc.py" % (indID, self.familyID))
            elif tmp == 0 and pIDs not in parentPairs:
                parentPairs.append(pIDs)
            else:
                continue
        return parentPairs
    

    def _getConnectionIndIDs(self):
        '''
        obtain IDs of connecting individuals (those through whom nuclear family blocks are connected)
        '''
        tmp = []
        for i in self.pIDs:
            for j in i:
                if (self.famInfo[j]['fatherID'] != '0' and self.famInfo[j]['motherID'] != 0) or (len(self.famInfo[j]['spouseIDs']) > 1):
                    tmp.append(j)
        connIndIDs = list(set(tmp))
        return connIndIDs
    
    
    def _createNuclearFamilyInfoDict(self):
        '''
        For each nuclear family block retrieve the following information:
            "IDs": [father, mother, offspring, ...];
            "isFounder": [True/False, True/False];
            "gen": [parental generation]
            "pheno": [1,2,...];
            "cIndIDs: [c1ID, c2ID,...]; (IDs of connecting individuals)
            "cIndIdxInNuc": [c1IdxInNuc, c2IdxInNuc, ...];
            "cIndIdxInAll": [c1IdxInAll, c2IdxInAll, ...] (indices of connection individuals in self.connIDs);
        '''
        def _convertTraitValue(val):
            try:
                return int(val)
            except:
                try:
                    return float(val)
                except:
                    return val
        nucInfo = []
        for pID in self.pIDs:
            infoDict = {"IDs":[], "isFounder":[False, False], "pheno":[], "cIndIDs":[], "cIndIdxInNuc":[], "cIndIdxInAll":[]}
            fID, mID = pID[0], pID[1]
            childIDs = list(set(self.famInfo[fID]['childIDs']).intersection(self.famInfo[mID]['childIDs']))
            #childIDs.sort()
            infoDict["IDs"] = [fID, mID] + childIDs
            infoDict["isFounder"][0] = True if self.famInfo[fID]['fatherID'] == '0' and self.famInfo[fID]['motherID'] == '0' else False
            infoDict["isFounder"][1] = True if self.famInfo[mID]['fatherID'] == '0' and self.famInfo[mID]['motherID'] == '0' else False
            infoDict["pheno"] = [_convertTraitValue(self.famInfo[iID]['trait']) for iID in infoDict["IDs"]]
            infoDict["cIndIDs"] = [iID for iID in infoDict["IDs"] if iID in self.connIDs]
            infoDict["cIndIdxInNuc"] = [infoDict["IDs"].index(iID) for iID in infoDict["cIndIDs"]]
            infoDict["cIndIdxInAll"] = [self.connIDs.index(iID) for iID in infoDict["cIndIDs"]]
            infoDict['gen'] = min(self.famInfo[fID]['gen'], self.famInfo[mID]['gen'])
            #
            nucInfo.append(infoDict)
        return nucInfo
    
    
    def generateConnIndsProbMap_complexTrait(self, _precision=3):
        '''
        '''
        # calculate geno-probability map per nuclear family block within the given pedigree
        ins = utilityFunc.ComplexTraitProbabilityMap()
        # create family-wide CIProbMap
        CIProbMap = collections.OrderedDict({})
        if len(self.connIDs) > 0:
            genos = [''.join(x) for x in list(itertools.combinations_with_replacement([str(i) for i in range(self.maxVar+1)], 2))]
            [CIProbMap.update({i: 1}) for i in list(itertools.product(genos, repeat=len(self.connIDs)))]
        for idx, nuc in enumerate(self.nucInfo):
            # create IndProb map for each nuclear family block
            probMap = ins.createProbMap(self.hapVarFreq, aff=nuc['pheno'], isFounder=nuc['isFounder'], traitType=self.traitType, indProbAffDict=self.indProbAffDict, meanshift=self.meanshift)
            # add 'probMap' dict to nucInfo
            nuc['probMap'] = probMap
            # update CIProbMap using probMap info contained in the idx'th nuc family block
            if len(nuc['cIndIdxInNuc']) > 0:
                # first compute indProbMap of cInd(s)
                nucCIProbMap = {}
                [nucCIProbMap.update({i: 0}) for i in list(itertools.product(genos, repeat=len(nuc['cIndIdxInNuc'])))]
                for g in probMap:
                    cIndGeno = tuple([g[i] for i in nuc['cIndIdxInNuc']])
                    nucCIProbMap[cIndGeno] += probMap[g]
                    nucCIProbMap[cIndGeno] = float(format(nucCIProbMap[cIndGeno], '.{}g'.format(_precision)))
                # normalize dict values
                self._normalizeDictValues(nucCIProbMap, _precision)
                # then update family-wide CIProbMap using nucCIProbMap
                for k in nucCIProbMap:
                    for g in CIProbMap:
                        if tuple([g[i] for i in nuc['cIndIdxInAll']]) == k:
                            CIProbMap[g] *= nucCIProbMap[k]
        # normalize CIProbMap values
        self._normalizeDictValues(CIProbMap, _precision)
        
        if self.scale != -1:
            utilityFunc.TernaryProbabilityMap.normalize(self.numReps*self.scale)(CIProbMap)
        self.CIProbMap = collections.OrderedDict(CIProbMap)
        
        
    def _normalizeDictValues(self, dictObj, _precision=3):
        s = sum(dictObj.values())
        for k in dictObj.keys():
            dictObj[k] = float(format(dictObj[k]/s, '.{}g'.format(_precision)))
        return
    
    
    def generateConnIndsProbMap(self):
        '''
        generate possible configuration of genotypic origins for the given pedigree structure
        '''
        # calculate geno-probability map per nuclear family block within the given pedigree
        ins = utilityFunc.TernaryProbabilityMap()
        CIProbMap = {}
        if len(self.connIDs) > 0:
            ## FIXME: Need optimize the following code when
            ## number of connection individuals gets large (>15)
            [CIProbMap.update({''.join(i):1}) for i in list(itertools.product('+-0', repeat=len(self.connIDs)))]
        # update self.nucInfo by adding prob of genos of connecting inds in each nuclear family
        for idx, nuc in enumerate(self.nucInfo):
            # create IndProb map for each nuclear family block
            probMap = ins.createProbMap_Mendelian(self.hapVarFreq, aff=nuc['pheno'], numReps=self.numReps, isFounder=nuc['isFounder'], P_T=self.P_T, nucIDs=nuc['IDs'], familyID=self.familyID, scale=self.scale)
            indMap = ins.createIndMap(len(nuc['IDs']), updatedKeys=nuc['IDs'])
            # remove 0-prob geno entries in both 'probMap' and 'indMap'
            for key, value in zip(probMap.keys(), probMap.values()):
                if value == 0:
                    del probMap[key]
            pmKeys = probMap.keys()
            for indID in indMap.keys():
                for g in '+-0':
                    indMap[indID][g] = list(set(indMap[indID][g]).intersection(pmKeys))
            # add updated 'indMap' and 'probMap' to nucInfo
            nuc['indMap'] = indMap
            nuc['probMap'] = probMap
            # create individual based probability map
            # indProbMap = ins.createIndProbMap(indMap, probMap)
            if len(nuc['cIndIDs']) > 0:
                self.nucConnIndsProbMap = ins.nucConnectIndsProbMap(indMap, probMap, nuc['cIndIDs'], familyID=self.familyID)
                ### WARNING!NOTE!: nucConnInds can be none (single nuclear family pedigree or multiple disjoint nuc fams but all individuals share same family ID)
                # update CIProbMap
                self._updateCIProbMap(nuc['cIndIdxInAll'])(CIProbMap)
        if self.scale != -1:
            ins.normalize(self.numReps*self.scale)(CIProbMap)
        self.CIProbMap = collections.OrderedDict(CIProbMap)
        #print self.connIDs
        return
    
    
    @staticmethod
    def restoreGenotypeOrigin(nucFamilyInfo, connectIndsProbMap):
        '''
        "For a Mendelian trait" restore genotypic origin of all individuals given genotype probability map of
        connecting individuals and info of all nuclear (joint) families
        '''
        def _intToTernary(integer, numInds):
            t = str(utilityFunc.BalancedTernary(integer))
            return '0'*(numInds-len(t)) + t
        # 
        indGenos = {}
        # Randomly choose genotypes of connecting individuals according to their probability map 'connectIndsProbMap'
        if len(connectIndsProbMap) == 0:
            connIndsGeno = ''
        else:
            idx = utilityFunc.wheelSpin(connectIndsProbMap.values())
            connIndsGeno = connectIndsProbMap.keys()[idx]
        #
        for nuc in nucFamilyInfo:
            # w/o connecting individual(s)
            if len(nuc['cIndIDs']) == 0:
                intGenoKeys = nuc['probMap'].keys()
            # just one connecting individual
            elif len(nuc['cIndIDs']) == 1:
                intGenoKeys = nuc['indMap'][nuc['cIndIDs'][0]][connIndsGeno[nuc['cIndIdxInAll'][0]]]
            # more than one connecting individual(s)
            else:
                nucConnGenos = [connIndsGeno[idx] for idx in nuc['cIndIdxInAll']]
                tmp = [nuc['indMap'][indID][g] for indID, g in zip(nuc['cIndIDs'], nucConnGenos)]
                intGenoKeys = list(set(tmp[0]).intersection(*tmp[1:]))
            intGenoProbs = [nuc['probMap'][key] for key in intGenoKeys]
            nucTernaryGeno = _intToTernary(intGenoKeys[utilityFunc.wheelSpin(intGenoProbs)], numInds=len(nuc['IDs']))
            # update dict 'indGenos'
            [indGenos.update({indID:g}) for g, indID in zip(nucTernaryGeno, nuc['IDs'])]
        return indGenos  
            
    
    @staticmethod
    def restoreGenotypeOrigin_complexTrait(nucFamilyInfo, connectIndsProbMap):
        '''
        For a Complex trait resotre genotypic origins of all individuals
        '''
        indGenos = {}
        # randomly choose genotypes connecting inds
        if len(connectIndsProbMap) == 0:
            connIndsGeno = ''
        else:
            idx = utilityFunc.wheelSpin(connectIndsProbMap.values())
            connIndsGeno = connectIndsProbMap.keys()[idx]
        def _checkIfMatch(g, i, k):
            return True if k[i] == g else False
        for nuc in nucFamilyInfo:
            indGenoKeys = []
            indGenoProbs = []
            # no connecting individual
            if len(nuc['cIndIDs']) == 0:
                indGenoKeys = nuc['probMap'].keys()
                indGenoProbs = nuc['probMap'].values()
            #
            else: # choose all possible genotypes of nuc fam that correspond to connInds genos
                connIndsGenoInNuc = [connIndsGeno[idx] for idx in nuc['cIndIdxInAll']]
                for key, prob in zip(nuc['probMap'].keys(), nuc['probMap'].values()):
                    if not False in [_checkIfMatch(geno, idx, key) for geno, idx in zip(connIndsGenoInNuc, nuc['cIndIdxInNuc'])]:
                        indGenoKeys.append(key)
                        indGenoProbs.append(prob)
            nucGeno = indGenoKeys[utilityFunc.wheelSpin(indGenoProbs)]
            # update dict 'indGenos'
            [indGenos.update({indID:g}) for g, indID in zip(nucGeno, nuc['IDs'])]
        return indGenos
    
    
    def _updateCIProbMap(self, nucConnIndsIdx):
        '''
        update probability map in '+-0' format for all connection individuals
        '''
        availKeys = self.nucConnIndsProbMap.keys()
        def operator(probMap):
            genoToDelete = []
            for geno in probMap.keys():
                nucGeno = ''.join([geno[i] for i in nucConnIndsIdx])
                if nucGeno in availKeys:
                    probMap[geno] *= self.nucConnIndsProbMap[nucGeno]
                else:
                    genoToDelete.append(geno)
            for geno in genoToDelete:
                del probMap[geno]   
            return
        return operator
    
    
    def _genoMapSingleInd(self):
        genoMap = collections.OrderedDict({})
        [genoMap.update({i:0}) for i in ['+', '0', '-']]
        if self.indPheno == 0:
            genoMap['+'] = self.k ** 2
            genoMap['-'] = (1-self.k) ** 2
            genoMap['0'] = 2 * self.k * (1-self.k)
        else:
            p_T = self.k if self.indPheno == 2 else 1-self.k
            genoMap['+'] = self.P_T_X['+'][self.indPheno-1]*(self.k**2)/p_T
            genoMap['0'] = self.P_T_X['0'][self.indPheno-1]*(2*self.k*(1-self.k))/p_T
            genoMap['-'] = self.P_T_X['-'][self.indPheno-1]*((1-self.k)**2)/p_T
        return genoMap
            
            


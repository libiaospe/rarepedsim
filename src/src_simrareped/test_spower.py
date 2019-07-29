#!/usr/bin/env python
# $File: test_spower.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the RarePedSim program
# Copyright (c) 2013-2015, Biao Li <libiaospe@gmail.com, biaol@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#
# Author: Biao Li
# purpose: calculate sample individual's phenotype based on genotype and disease model

import sys, random, os, logging
from utilityFunc import readFiles
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
from gdata import GFile
import collections

DEBUG = True

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
        
        ## FIXME!
        resDict['pool'] = None
        resDict['missing'] = None
        
        
        return resDict
    except Exception, e:
        if DEBUG: print e
        else: pass    
    # use *.gdat instead
    try:
        with GFile.open(filePrefix+'.gdat') as g:
            names = g.getnames()
            data = g.getdata(random.choice(names)) if rep is None else g.getdata('rep'+str(rep))
            return {'maf':data['maf'], 'function_score':data['annotation'], 'position':data['position'], 'num_variants':len(data['maf'])}
    except Exception, e:
        if DEBUG: print e
        else:
            raise ValueError("An error occurs possibly due to:\n1. Cannot find {0}.sfs file or {0}.gdat file.\n2. The gene database does not have enough gene replicates.\n".format(filePrefix))
    

class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
    

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


def initPowerdataObj(configFile, srvResFile, srvGeneRep=None):
    '''
    initialize a PowerData object for each replicate
    srvGeneRep - gene replicate # to be used in powerdata, randomly sample a gene replicate if None
    '''
    args, srvData = spowerArgsNData(configFile, srvResFile, srvGeneRep)
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


def calProb_T_X(powerdata, n=10000):
    '''
    calculate probablity of being unaffected/affected given haplotype origins ('+', '0', '-')
    n - sample size to create for calculation of P_T_X
    '''
    # initialize a sample object first
    sampleObj = initSampleObj(powerdata)
    # update case and control MAFs
    caseMAF = getUpdatedMAF(sampleObj, powerdata, phenotype=2)
    ctrlMAF = getUpdatedMAF(sampleObj, powerdata, phenotype=1)
    phenoPool = {'+':[], '0':[], '-':[]}
    def _genHap(maf):
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
    phenoPool = {
        '+': [generatePhenotype(sampleObj, powerdata, genotype=[_genHap(caseMAF), _genHap(caseMAF)]) for _ in xrange(n)],
        '0': [generatePhenotype(sampleObj, powerdata, genotype=[_genHap(caseMAF), _genHap(ctrlMAF)]) for _ in xrange(n)],
        '-': [generatePhenotype(sampleObj, powerdata, genotype=[_genHap(ctrlMAF), _genHap(ctrlMAF)]) for _ in xrange(n)]
    }
    p = [(phenoPool['+'].count(2.0))/float(n), (phenoPool['0'].count(2.0))/float(n), (phenoPool['-'].count(2.0))/float(n)]
    Prob_T_X = {'+':[1-p[0], p[0]], '0':[1-p[1], p[1]], '-':[1-p[2], p[2]]}
    return Prob_T_X


def temp(powerdata, n=100000):
    sampleObj = initSampleObj(powerdata)
    # update case and control MAFs
    caseMAF = getUpdatedMAF(sampleObj, powerdata, phenotype=2)
    ctrlMAF = getUpdatedMAF(sampleObj, powerdata, phenotype=1)
    phenoPool = {'+':[], '0':[], '-':[]}
    def _genHap(maf):
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
    #
    hapPool = {
        '+': [_genHap(caseMAF).count(1) for _ in xrange(n)],
        '-': [_genHap(ctrlMAF).count(1) for _ in xrange(n)]
    }
    counts = {'+': dict(collections.Counter(hapPool['+'])),
              '-': dict(collections.Counter(hapPool['-']))
    }
    
    freqs = {'+': [float(i)/n for i in dict(collections.Counter(hapPool['+'])).values()],
              '-': [float(i)/n for i in dict(collections.Counter(hapPool['-'])).values()]
    }
    
    caseHapPool = []
    for i in xrange(1000000):
        hap1 = _genHap(powerdata.data['maf'])
        hap2 = _genHap(powerdata.data['maf'])
        if generatePhenotype(sampleObj, powerdata, genotype=[hap1, hap2]) == 2:
            caseHapPool.append(hap1)
            caseHapPool.append(hap2)
    
    
    print freqs

    tmp = [i.count(1) for i in caseHapPool]
    caseHapsVarCount = dict(collections.Counter(tmp))
    print caseHapsVarCount
    print sum(tmp)/float(len(tmp))
    print [float(i)/sum(caseHapsVarCount.values()) for i in caseHapsVarCount.values()]
    
    
    return counts
    
    


        
if __name__ == '__main__':
    configFileName = '/Users/biao/MySVN/simRarePed/temp/phenotype-model.conf'
    srvResFileName = '/Users/biao/MySVN/simRarePed/temp/Kryukov2009European1800Rep500Cp1.sfs'
    
    srvResFileName = '/Users/biao/MySVN/simRarePed/temp/Boyko2008African1800.sfs'
    
    powerdata = initPowerdataObj(configFileName, srvResFileName, 1)
    sampleObj = initSampleObj(powerdata)

    #sampleObj, powerdata = initSampleAndPowerdata(configFileName, srvResFileName, 2)
    #
    #print sampleObj.get('maf')[:10]
    
    print temp(powerdata, 100000)
    sys.exit(0)
    
    print calProb_T_X(sampleObj, powerdata, 10000)
    
    sys.exit()
    
    phenos = [generatePhenotype(sampleObj, powerdata) for _ in xrange(5000)]
    print phenos.count(1.0)/5000.
    
    temp_sampleObj = sampleObj.clone()
    updateMAF(temp_sampleObj, powerdata, 1)
    phenos = [generatePhenotype(temp_sampleObj, powerdata) for _ in xrange(5000)]
    print (phenos.count(1.0)/5000.)
    
    temp_sampleObj = sampleObj.clone()
    updateMAF(temp_sampleObj, powerdata, 2)
    phenos = [generatePhenotype(temp_sampleObj, powerdata) for _ in xrange(5000)]
    print phenos.count(1.0)/5000.
    
    
    
    
    sys.exit(0)
    
    
    
    #print powerdata.data
    #dir(powerdata.data)
    #attrs = dir(powerdata.args)
    #for a in attrs:
    #    print a, getattr(powerdata.args, a)
    #
    #### generate disease by ORs
    #s = Sample()
    #s.seed(0, os.getpid())
    #s.set(**powerdata.data)
    #
    #tmp = []
    #for i in range(10000):
    #    L.GenotypeGenerator().apply(s.data)
    #    L.ORModel(powerdata.data['odds_ratio']).apply(s.data)
    #    L.DiseaseEffectGenerator(powerdata.data['baseline_effect']).apply(s.data)
    #    L.DiseaseStatusGenerator().apply(s.data)
    #    tmp.append(s.get('phenotype'))
    #
    #sys.exit(0)
    
    #L.GenotypeGenerator().apply(s.data)
    #L.ORModel(powerdata.data['odds_ratio']).apply(s.data)
    #L.DiseaseEffectGenerator(powerdata.data['baseline_effect']).apply(s.data)
    #L.DiseaseStatusGenerator().apply(s.data)
    #print powerdata.data['maf']
    #print s.get('haplotype1'), s.get('haplotype2')
    #print s.get('effect'),
    #print s.get('phenotype')
    #
    #
    ##### generate disease by PAR
    #s = Sample()
    #s.seed(0)
    #s.set(**powerdata.data)
    #L.GenotypeGenerator().apply(s.data)
    #L.PARModel(powerdata.data['par'], powerdata.data['PAR_variable']).apply(s.data)#L.ORModel(powerdata.data['odds_ratio']).apply(s.data)
    #L.DiseaseEffectGenerator(powerdata.data['baseline_effect']).apply(s.data)
    #L.DiseaseStatusGenerator().apply(s.data)
    #print powerdata.data['maf']
    #print s.get('haplotype1'), s.get('haplotype2')
    #print s.get('effect'),
    #print s.get('phenotype')
    
    
    #### generate qt
    s = Sample()
    s.seed(0)
    
    print dir(s)
    
    s.set(**powerdata.data)
    #L.GenotypeGenerator().apply(s.data)
    s.set(haplotype1=[0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0], haplotype2=[1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1])
    L.MeanShiftModel(powerdata.data['mean_shift']).apply(s.data)
    L.QtEffectGenerator().apply(s.data)
    L.QtValueGenerator().apply(s.data)
    print s.get('haplotype1'), s.get('haplotype2')
    print s.get('effect')
    print s.get('phenotype')
    
    
    
    
  
    
    
    
#!/usr/bin/env python
# $File: parallel.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the RarePedSim program
# Copyright (c) 2013-2015, Biao Li <libiaospe@gmail.com, biaol@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#
# Author: Biao Li
# purpose: wrapper functions to run simulation in parallel

from simRarePed import SimPed_backward, SimPed_forward
from utilityFunc import restoreStruct, SaveFiles

def parallel_restoreGenoOrigin(familyID, familyDict, spowerData, moi, numReps, scale, P_T=None):
    ## FIXME:! add progressbar to show status when SimPed_backward is computing CIProbMap
    ins = SimPed_backward(familyID, familyDict, spowerData, numReps, scale, moi, P_T)
    return ins.famInfo, ins.nucInfo, ins.CIProbMap


def parallel_generateGeno(rep, pedInfo, sfsDict, famInfoList, nucInfoList, CIProbMapList, moi, missingSites, probMissingCalls, probErrorCalls, causality, hapCoding='1-2'):
    fileName = 'rep_'+str(rep)+'.ped'
    fi = open(fileName, 'w')
    for idx, familyID in enumerate(pedInfo.keys()):
        # generate genotype origins (++ or -- or +-/0)
        indGenoOrigins = SimPed_backward.restoreGenotypeOrigin(nucInfoList[idx], CIProbMapList[idx])
        # restore actual genotypes
        genoDict = SimPed_backward.generateGenotype(nucInfoList[idx], indGenoOrigins, sfsDict, moi=moi)
        # output simu data in Linkage format
        SimPed_backward.writeToFile(familyID, famInfoList[idx], genoDict, fi, missingSites, probMissingCalls, probErrorCalls, causality, hapCoding=hapCoding)
    fi.close()
    return


def parallel_restoreStruct(familyID, familyDict):
    return restoreStruct(familyID, familyDict, addMissingParents=True).famInfo


def parallel_simForward(rep, pedInfo, sfsDict, famInfoList, spowerData, missingSites, probMissingCalls, probErrorCalls, causality, hapCoding='1-2'):
    fileName = 'rep_'+str(rep)+'.ped'
    fi = open(fileName, 'w')
    for familyID, famInfo in zip(pedInfo.keys(), famInfoList):
        simRes = SimPed_forward(familyID, sfsDict, famInfo, missingSites, probMissingCalls, probErrorCalls, causality).generate(spowerData, hapCoding)
        SaveFiles.simulatedData(fi, simRes)
    fi.close()
    return


def parallel_restoreGenotypeOrigins_Mendelian(familyID, familyDict, numReps, scalar, hapVarFreq, penetrance):
    '''
    '''
    ins = SimPed_backward(familyID, familyDict, numReps, scalar, hapVarFreq, P_T=penetrance)
    return ins.famInfo, ins.nucInfo, ins.CIProbMap


def parallel_restoreGenotypeOrigins_Complex(familyID, familyDict, configInfo, args, hapVarFreq, indProbAffDict=None):
    '''
    '''
    if configInfo['model'] == 'LOGIT':
        ins = SimPed_backward(familyID, familyDict, args.num_reps, args.scalar, hapVarFreq, indProbAffDict=indProbAffDict, traitType='Qualitative')
    elif configInfo['model'] == 'LNR':
        ins = SimPed_backward(familyID, familyDict, args.num_reps, args.scalar, hapVarFreq, traitType='Quantitative', meanshift=configInfo['meanshift_rare_detrimental'])
    else:
        raise ValueError("Can only choose model between 'LOGIT' and 'LNR'")
    return ins.famInfo, ins.nucInfo, ins.CIProbMap



def parallel_generateGenotype_Mendelian_single(rep, geneIdx, geneName, pedInfo, famInfoList, nucInfoList, CIProbMapList, sfsInfo, causalVarSites, causalVarMafs, configInfo, recRate, hapCoding='1-2'):
    '''
    '''
    fileName = 'rep{}.ped'.format(rep)
    fi = open(fileName, 'w')
    for famIdx, familyID in enumerate(pedInfo.keys()):
        # sample genotype origins for all individuals contained in the 'famIdx'th family
        indGenoOrigins = SimPed_backward.restoreGenotypeOrigin(nucInfoList[famIdx], CIProbMapList[famIdx])
        # restore actual genotypes
        genoDict = SimPed_backward.generateGenotype(nucInfoList[famIdx], indGenoOrigins, sfsInfo[geneName], causalVarSites[geneIdx], causalVarMafs[geneIdx], configInfo['moi'], recRate, configInfo['compound_hetero'])
        # write simulated replicate to file
        if configInfo['missing_low_maf'] == None:
            configInfo['missing_low_maf'] = 0
        if configInfo['missing_sites'] == None:
            configInfo['missing_sites'] = 0
        if configInfo['missing_calls'] == None:
            configInfo['missing_calls'] = 0
        if configInfo['error_calls'] == None:
            configInfo['error_calls'] = 0
        lowMafSitesIdx = [idx for idx, maf in enumerate(sfsInfo[geneName]['maf']) if maf < configInfo['missing_low_maf']]
        tmpMissingSitesIdx = [idx for idx in range(len(sfsInfo[geneName]['maf'])) if random.random() < configInfo['missing_sites']]
        _missingSites = list(set.union(*[set(lowMafSitesIdx), set(tmpMissingSitesIdx)]))
        _probMissingCalls = [configInfo['missing_calls']] * 3
        _probErrorCalls = [configInfo['error_calls']] * 3
        _causality = ['d'] * len(sfsInfo[geneName]['maf'])
        SimPed_backward.writeToFile(familyID, famInfoList[famIdx], genoDict, fi, _missingSites, _probMissingCalls, _probErrorCalls, _causality, hapCoding=hapCoding)
    fi.close()

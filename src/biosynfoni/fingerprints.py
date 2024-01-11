# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 12:20:39 2023
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: fingerprints                 ||
created: 2023-09-28 12:20           ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
____________________________________
 
||||||||||||  ()()()  |||||||||||||||

description: fingerprints and related functionality.
"""
import sys
import numpy as np
from numpy import dot
from numpy.linalg import norm

from rdkit import Chem, DataStructs
from rdkit.DataManip import Metric
from rdkit.Chem import AllChem

# my imports
from biosynfoni.def_biosynfoni import defaultVersion
from biosynfoni.concerto_fp import Biosynfoni


# circular fingerprint -------------------------------------------------
def morgan_getter(
    mol: Chem.Mol, useChirality: bool = False, radius: int = 2, nBits: int = 2048
) -> np.array:
    """returns explicit bit vector fp"""
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=useChirality, radius=radius, nBits=nBits
    )
    return np.array(fingerprint)


# topology path fingerprint ---------------------------------------------
def rdk_fp_getter(mol: Chem.Mol, nbits: str = 2048) -> np.array:
    """returns explicit bit vector fp"""

    fingerprint = Chem.RDKFingerprint(mol, fpSize=nbits)
    return np.array(fingerprint)


# substructure key fingerprints ------------------------------------------
def maccs_getter(mol: Chem.Mol) -> np.array:
    """returns explicit bit vector fp"""

    fingerprint = AllChem.GetMACCSKeysFingerprint(mol)
    return np.array(fingerprint)


def biosynfoni_getter(
    mol: Chem.Mol, version: str = defaultVersion, *args, **kwargs
) -> np.array:
    """returns counted fingerprint list"""
    counted_fingerprint = Biosynfoni(
        mol, version_name=version, *args, **kwargs
    ).fingerprint
    return np.array(counted_fingerprint)


def maccsynfoni_getter(
    mol: Chem.Mol, version: str = defaultVersion, *args, **kwargs
) -> np.array:
    """returns counted fingerprint list"""
    counted_fingerprint = biosynfoni_getter(mol, version=version, *args, **kwargs)
    maccs = maccs_getter(mol)
    return np.concatenate((counted_fingerprint, maccs))


def bino_maccs_getter(
    mol: Chem.Mol, version: str = defaultVersion, *args, **kwargs
) -> np.array:
    binosynfoni = binosynfoni_getter(mol, version=version, *args, **kwargs)
    maccs = maccs_getter(mol)
    return np.concatenate((binosynfoni, maccs))


def binosynfoni_getter(
    mol: Chem.Mol, version: str = defaultVersion, *args, **kwargs
) -> np.array:
    """returns explicit bit vector"""
    counted = biosynfoni_getter(mol, version=version, *args, **kwargs)
    binary = counted.copy()
    binary[np.where(counted > 0)] = 1
    # binosynfoni = list_to_bitvect(binary)
    return binary


# ============================= distance, similarity =========================
# ----------------------------- distance -------------------------------------


def bitvect_to_tanimoto(expl_bitvectors: list[DataStructs.ExplicitBitVect]) -> float:
    array = Metric.rdMetricMatrixCalc.GetTanimotoDistMat(expl_bitvectors)
    return array.tolist()[0]


def bitvect_to_euclidean(expl_bitvectors: list[DataStructs.ExplicitBitVect]) -> float:
    array = Metric.rdMetricMatrixCalc.GetEuclideanDistMat(expl_bitvectors)
    return array.tolist()[0]


def bitvect_to_manhattan(expl_bitvectors: list[DataStructs.ExplicitBitVect]) -> float:
    array = Metric.rdMetricMatrixCalc.GetManhattanDistMat(expl_bitvectors)
    return array.tolist()[0]


def bitvect_to_cosine(expl_bitvectors: list[DataStructs.ExplicitBitVect]) -> float:
    array = Metric.rdMetricMatrixCalc.GetCosineDistMat(expl_bitvectors)
    return array.tolist()[0]


# ----------------------------- similarity -----------------------------------
def counted_tanimoto_sim(fp1: np.array, fp2: np.array) -> float:
    nom = sum(np.minimum(fp1, fp2))  # overlap
    denom = float(sum(np.maximum(fp1, fp2)))  # all bits that are 'on'
    if denom == 0:
        return -1  # fingerprint does not describe either molecule
    else:
        return nom / denom


def countanimoto(pair: list) -> float:
    return counted_tanimoto_sim(np.array(pair[0]), np.array(pair[1]))


def cosine_sim(pair: list[np.array]) -> float:
    cos_sim = dot(np.array(pair[0]), np.array(pair[1])) / (
        norm(np.array(pair[0])) * norm(np.array(pair[1]))
    )
    return cos_sim


# ============================= formatting ===================================


def list_to_bitvect(lst: list):
    """returns DataStructs.cDataStructs.ExplicitBitVect type"""
    binarystring = "".join([str(x) for x in lst])
    return DataStructs.cDataStructs.CreateFromBitString(binarystring)


def array_to_bitvect(numpy_array):
    return list_to_bitvect(numpy_array.tolist())

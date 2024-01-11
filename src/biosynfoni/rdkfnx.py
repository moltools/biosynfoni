# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:59:28 2023
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: RDKit functions              ||
created: 2023-09-13                 ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
____________________________________
 
||||||||||||  ()()()  |||||||||||||||

contains general reusable functionalities depending on RDKit packages,
mainly supplier-handling
"""

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")  # for muting warnings

from biosynfoni.subkeys import (
    substructureSmarts,
    fpVersions,
    defaultVersion,
    get_smarts,
    get_names,
    get_values,
)
from biosynfoni.inoutput import outfile_namer
from biosynfoni.moldrawing import drawfp, pathway_colours


def supplier(sdf_file: str) -> Chem.SDMolSupplier:
    """returns a mol supplier from an sdf file"""
    suppl = Chem.SDMolSupplier(sdf_file)
    return suppl


def sdf_writer(outfile: str) -> None:
    """returns an sdf writer object to write to outfile"""
    return Chem.SDWriter(outfile)


# used to be get_subsset
def get_subs_set(
    fp_version_name: str,
    subs_set: dict[dict] = substructureSmarts,
    fp_versions: dict[list] = fpVersions,
) -> list[Chem.Mol]:
    """gives list of rdkit.Chem.Mols of substructures of choice
    input:   fp_version_name (str) -- name of the version
             subs_smarts (dict)
                         (key) substructure names (e.g. 'fp1')
                         (val) (dict)
                            'smarts': substructure RDK molfiles (from SMARTS)
             fp_versions (dict)
                         (key) version name (e.g. fps_full_2)
                         (val) (list) substructure names (e.g. 'd_isoprene')

    output: (list) rdkit.Chem.rdchem.Mol files for substructure keys
    """

    if not fp_version_name:
        raise "No version name provided to select substructure set"

    substructures_smarts = []

    chosen_version = fp_versions[fp_version_name]
    for substructure_key in chosen_version:
        substructure_properties = subs_set[substructure_key]
        smartsmol = Chem.MolFromSmarts(substructure_properties["smarts"])
        if smartsmol is None:
            print(f"{substructure_key} could not be converted to mol: skipped")
            continue
        substructures_smarts.append(smartsmol)

    return substructures_smarts


class BiosynfoniVersion:
    def __init__(self, fp_version: str = defaultVersion):
        if fp_version.lower() == "default":
            fp_version = defaultVersion
        self.fp_version = fp_version
        self.substructures = get_subs_set(fp_version)
        self.smarts = get_smarts(fp_version)
        self.subs_ids = fpVersions[fp_version]
        self.subs_names = get_names(fp_version)
        self.subs_pathways = get_values("pathway", fp_version)
        self.subs_colors = None

    def save_smarts(self):
        outfilename = outfile_namer("version", self.fp_version)
        with open(f"{outfilename}.smarts", "w") as f:
            for smart in self.smarts:
                f.write(f"{smart[0]}\t{smart[1]}\n")
        return None

    # def save_svg(self, window_size=(1000, 1000)):
    #     outfilename = outfile_namer("version", self.fp_version)
    #     svg_text = Draw.MolsToGridImage(
    #         self.substructures,
    #         molsPerRow=4,
    #         subImgSize=window_size,
    #         legends=self.subs_names,
    #     )
    #     with open(f"{outfilename}.svg", "w") as f:
    #         f.write(svg_text)
    #     return None
    def save_svg(self, window_size=(1000, 1000)):
        outfilename = outfile_namer("version", self.fp_version)
        svg_text = drawfp(
            subs_set=self.substructures, subs_ids=self.subs_ids, window_size=window_size
        )
        print(svg_text)
        with open(f"{outfilename}.svg", "w") as f:
            f.write(svg_text)
        return None

    def _get_sub_color(self, ind: int):
        pathway = self.subs_pathways[ind]
        if pathway:
            first_pathway = pathway[0]
            color = pathway_colours[first_pathway]
        else:
            color = (0.75, 0.75, 0.75)  # grey
        return color

    def set_subs_colors(self):
        colors = []
        for i in range(len(self.subs_pathways)):
            colors.append(self._get_sub_color(i))
        self.subs_colors = colors
        return None

    def get_subs_colors(self):
        if self.subs_colors is None:
            self.set_subs_colors()
        return self.subs_colors


# def save_version(
#     fp_version: str, window_size=(1000, 1000), extra_text: str = ""
# ) -> None:
#     """depracated, use BiosynfoniVersion instead"""
#     outfilename = outfile_namer("version", f"{fp_version}_{extra_text}")

#     # svg_text = fm.drawfp(fp_version, window_size=window_size)
#     # with open(f'{outfilename}.svg', 'w') as f:
#     #    f.write(svg_text)

#     smarts = get_smarts(fp_version)
#     with open(f"{outfilename}.smarts", "w") as f:
#         for smart in smarts:
#             f.write(f"{smart[0]}\t{smart[1]}\n")
#     return None


def save_version(version: str, window_size=(1000, 1000)) -> None:
    saver = BiosynfoniVersion(version)
    saver.save_smarts()
    saver.save_svg(window_size=window_size)
    return None


# ===========================individual molecule functions=======================================


def smiles_to_mol(smiles: str) -> Chem.Mol:
    """returns rdkit.Chem.rdchem.Mol from smiles"""
    return Chem.MolFromSmiles(smiles)


def inchi_to_mol(inchi: str) -> Chem.Mol:
    """returns rdkit.Chem.rdchem.Mol from inchi"""
    return Chem.MolFromInchi(inchi)


def nonh_atomcount(mol: Chem.Mol) -> int:
    """returns number of non-hydrogen atoms in mol. GetNumAtoms on its own acts weird, so this is a workaround"""
    nonh_mol = Chem.rdmolops.RemoveHs(mol, implicitOnly=False, sanitize=True)
    return nonh_mol.GetNumAtoms()


def get_sub_matches(mol: Chem.Mol, substructure: Chem.Mol) -> list[list[int]]:
    return [list(matchatoms) for matchatoms in mol.GetSubstructMatches(substructure)]

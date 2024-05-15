# -*- coding: utf-8 -*-

"""Script for parsing MetaCyc database in mol-mol pairs and metabolic distance."""

import argparse
import logging
import os
import typing as ty

from rdkit import Chem, RDLogger
from tqdm import tqdm


# Turn any RDKit message off
RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")


def smiles_is_valid(smiles: str) -> bool:
    """Check if SMILES string is valid.
    
    :param smiles: SMILES string.
    :type smiles: str
    :return: True if valid, False otherwise.
    :rtype: bool
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return True
    except:
        pass

    return False


def inchi_is_valid(inchi: str) -> bool:
    """Check if InChI string is valid.
    
    :param inchi: InChI string.
    :type inchi: str
    :return: True if valid, False otherwise.
    :rtype: bool
    """
    try:
        mol = Chem.MolFromInchi(inchi)
        if mol is not None:
            return True
    except:
        pass

    return False


class CompoundRecord:
    """Compound record."""

    def __init__(
        self,
        unique_name: ty.Optional[str] = None,
        smiles: ty.Optional[str] = None,
        inchi: ty.Optional[str] = None
    ) -> None:
        """Initialize CompoundRecord.
        
        :param unique_name: Unique name.
        :type unique_name: ty.Optional[str]
        :param smiles: SMILES string.
        :type smiles: ty.Optional[str]
        :param inchi: InChI string.
        :type inchi: ty.Optional[str]
        """
        self.unique_name = unique_name
        self.smiles = smiles
        self.inchi = inchi

    def set_unique_name(self, unique_name: str) -> None:
        """Set unique name.
        
        :param unique_name: Unique name.
        :type unique_name: str
        """
        self.unique_name = unique_name

    def set_smiles(self, smiles: str) -> None:
        """Set SMILES string.
        
        :param smiles: SMILES string.
        :type smiles: str
        """
        if smiles_is_valid(smiles):
            self.smiles = smiles

    def set_inchi(self, inchi: str) -> None:
        """Set InChI string.
        
        :param inchi: InChI string.
        :type inchi: str
        """
        if inchi_is_valid(inchi):
            self.inchi = inchi

    def is_valid(self) -> bool:
        """Check if record is valid.
        
        :return: True if valid, False otherwise.
        :rtype: bool
        """
        if (
            self.unique_name is not None
            and (self.smiles is not None or self.inchi is not None)
        ):
            return True
        
        return False
    
    def auto_complete(self) -> None:
        """Auto complete record."""
        if self.smiles is None and self.inchi is None:
            return

        if self.smiles is None:
            mol = Chem.MolFromInchi(self.inchi)
            self.smiles = Chem.MolToSmiles(mol)

        if self.inchi is None:
            mol = Chem.MolFromSmiles(self.smiles)
            self.inchi = Chem.MolToInchi(mol)

    def is_complete(self) -> bool:
        """Check if record is complete.
        
        :return: True if complete, False otherwise.
        :rtype: bool
        """
        if (
            self.unique_name is not None
            and self.smiles is not None 
            and self.inchi is not None
        ):
            return True

        return False

    def to_tsv(self) -> str:
        """Convert record to TSV.
        
        :return: Record as TSV.
        :rtype: str
        """
        return f"{self.unique_name}\t{self.smiles}\t{self.inchi}\n"


def cli() -> argparse.Namespace:
    """Command line interface.
    
    :return: Parsed arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Parse MetaCyc database.")
    parser.add_argument("--data", type=str, required=True, help="Input MetaCyc database data folder.")
    parser.add_argument("--output", type=str, required=False, default=".", help="Output folder (default: current folder).")
    parser.add_argument("--log-file", type=str, required=False, default=None, help="Log file (default: None).")
    parser.add_argument("--log-level", type=str, required=False, default="INFO", help="Logging level (default: INFO).")
    return parser.parse_args()


def setup_logger(log_file: str, log_level: str) -> None:
    """Setup logger.
    
    :param log_file: Log file.
    :type log_file: str
    :param log_level: Logging level.
    :type log_level: str
    """
    # Remove previous log file if exists.
    if os.path.exists(log_file):
        os.remove(log_file)

    logging.basicConfig(
        filename=log_file,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        level=getattr(logging, log_level.upper())
    )


def validate_path(path: str) -> None:
    """Check if path exists.
    
    :param path: Path to check.
    :type path: str
    """
    logger = logging.getLogger(__name__)

    if not os.path.exists(path):
        logger.error(f"Path '{path}' does not exist.")
        raise FileNotFoundError(f"Path '{path}' does not exist.")
    else:
        logger.debug(f"Path '{path}' exists.")


def parse_compound_records(compounds_path: str, out_dir: str) -> None:
    """Parse compounds.
    
    :param compounds_path: Path to compounds file.
    :type compounds_path: str
    :param out_dir: Output directory.
    :type out_dir: str
    """
    logger = logging.getLogger(__name__)

    compounds_path_out = os.path.join(out_dir, "compounds.tsv")
    compounds_file_out = open(compounds_path_out, "w", encoding="utf-8")
    compounds_file_out.write("name\tSMILES\tInChI\n")

    valid_compound_count = 0
    record = CompoundRecord()

    with open(compounds_path, "r", encoding="utf-8", errors="ignore") as compounds_file:
        for line in tqdm(compounds_file, desc="Parsing compounds"):
            line = line.strip()
            
            if line.startswith("#"):
                continue

            if line.startswith("UNIQUE-ID"):
                
                if record.is_valid():
                    record.auto_complete()
                    valid_compound_count += 1
                    if record.is_complete():
                        compounds_file_out.write(record.to_tsv())

                unique_name = line.split(" - ")[1]
                record = CompoundRecord()
                record.set_unique_name(unique_name)

            if line.startswith("SMILES"):
                smiles = line.split(" - ")[1]
                record.set_smiles(smiles)
            
            if line.startswith("INCHI"):
                inchi = line.split(" - ")[1]
                record.set_inchi(inchi)

    logger.debug(f"Number of valid compound records: {valid_compound_count}")


def parse_compound_links(compound_links_path: str, out_dir: str) -> None:
    """Parse compound links.
    
    :param compound_links_path: Path to compound links file.
    :type compound_links_path: str
    :param out_dir: Output directory.
    :type out_dir: str
    """
    logger = logging.getLogger(__name__)

    compound_links_path_out = os.path.join(out_dir, "compound-links.tsv")
    compound_links_file_out = open(compound_links_path_out, "w", encoding="utf-8")
    compound_links_file_out.write("source\ttarget\tmetabolic_distance\n")

    valid_link_count = 0

    with open(compound_links_path, "r", encoding="utf-8", errors="ignore") as compound_links_file:
        for line in tqdm(compound_links_file, desc="Parsing compound links"):
            line = line.strip()
            
            if line.startswith("#"):
                continue
            
            logger.error(f"Unable to parse line: {line}")


def main() -> None:
    """Run main function."""
    args = cli()

    # Setup logger.
    setup_logger(args.log_file, args.log_level)
    logger = logging.getLogger(__name__)

    # Define data paths.
    compounds_path = args.data + "/compounds.dat"
    compound_links_path = args.data + ""

    # Check if paths exist.
    validate_path(compounds_path)
    validate_path(compound_links_path)

    # Parse compounds.
    parse_compound_records(compounds_path, args.output)

    # Parse compound links.
    parse_compound_links(compound_links_path, args.output)
    

if __name__ == "__main__":
    main()

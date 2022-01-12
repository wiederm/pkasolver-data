import argparse
import gzip
from rdkit import Chem
import tqdm
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize


def main():
    """
    takes sdf file of initial training molecules and
    sdf file of training molecules (both optionally gzipped)
    and returns only those initial training molecules
    not contained in the training molecules file as sdf file.
    """
    RDLogger.DisableLog("rdApp.*")
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input filename, type: .sdf.gz or .sdf")
    parser.add_argument("--filter", help="filter filename, type: .sdf or .sdf.gz")
    parser.add_argument("--output", help="output filename, type: .sdf.gz or .sdf")
    args = parser.parse_args()
    input_zipped = False
    print("inputfile:", args.input)
    print("outputfile:", args.output)
    ini_list = []
    smi_list = []
    # start with generating INCHI and SMILES for mols in the filter set
    for i in args.filter.split(","):
        # test if it's gzipped
        with gzip.open(i, "r") as fh:
            try:
                fh.read(1)
                input_zipped = True
            except gzip.BadGzipFile:
                input_zipped = False

        if input_zipped:
            with gzip.open(i, "r") as fh:
                suppl = Chem.ForwardSDMolSupplier(fh, removeHs=True)
                ini_list.extend(ini_filter(suppl))
            with gzip.open(i, "r") as fh:
                suppl = Chem.ForwardSDMolSupplier(fh, removeHs=True)
                smi_list.extend(smi_filter(suppl))
        else:
            with open(i, "rb") as fh:
                suppl = Chem.ForwardSDMolSupplier(fh, removeHs=True)
                ini_list.extend(ini_filter(suppl))
            with open(i, "rb") as fh:
                suppl = Chem.ForwardSDMolSupplier(fh, removeHs=True)
                smi_list.extend(smi_filter(suppl))

    print(f"{len(ini_list)} inchi test molecules found")
    print(f"{len(smi_list)} smi test molecules found")
    # test if it's gzipped
    with gzip.open(args.input, "r") as fh:
        try:
            fh.read(1)
            input_zipped = True
        except gzip.BadGzipFile:
            input_zipped = False

    if input_zipped:
        with gzip.open(args.input, "r") as fh:
            suppl = Chem.ForwardSDMolSupplier(fh, removeHs=True)
            processing(suppl, args, ini_list, smi_list)
    else:
        with open(args.input, "rb") as fh:
            suppl = Chem.ForwardSDMolSupplier(fh, removeHs=True)
            processing(suppl, args, ini_list, smi_list)


def smi_filter(suppl):
    un = rdMolStandardize.Uncharger()
    smi_list = []
    # SMILES are generated
    for mol in suppl:
        # the mol is neutralized
        mol = un.uncharge(mol)
        smi_list.append(Chem.MolToSmiles(mol))
    return smi_list


def ini_filter(suppl):
    un = rdMolStandardize.Uncharger()
    ini_list = []
    # InCHIs are generated
    for mol in suppl:
        # the mol is neutralized
        mol = un.uncharge(mol)
        ini_list.append(Chem.inchi.MolToInchi(mol))
    return ini_list


def processing(suppl, args, exclude_ini_list, exclude_smi_list):
    dup = 0
    skipped = 0
    written = 0
    # iterate through dataset for which molecules are filtered
    un = rdMolStandardize.Uncharger()
    with gzip.open(args.output, "wt+") as sdf_zip:
        with Chem.SDWriter(sdf_zip) as writer:
            for idx, mol in enumerate(tqdm.tqdm(suppl)):
                if mol:
                    # uncharge
                    mol_uncharged = un.uncharge(mol)
                    smiles = Chem.MolToSmiles(mol_uncharged)
                    try:
                        inchi = Chem.inchi.MolToInchi(mol_uncharged)
                    except Chem.rdchem.KekulizeException:
                        print(smiles)
                    # test if either an inchi or a smiles are in the exclude lists
                    if inchi in exclude_ini_list or smiles in exclude_smi_list:
                        dup += 1
                    else:
                        # if not write mol to filtered data set
                        written += 1
                        writer.write(mol)

                else:
                    skipped += 1

    print(f"{dup} duplicate molecules found and discarted")
    print(f"{skipped} molecules skipped")
    print(f"{written} molecules")


if __name__ == "__main__":
    main()

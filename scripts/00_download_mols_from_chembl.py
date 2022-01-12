from chembl_webresource_client.new_client import new_client
from tqdm import tqdm
import gzip
import argparse


def main():
    """
    Filters the molecules of the Chembl database
    by the specified criteria
    (e.g. max number of rule of five violation = 1)
    and save them in a gzipped sdf file.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="output filename, type: .sdf.gz")
    args = parser.parse_args()

    print("outputfile:", args.output)

    molecule = new_client.molecule

    # Filters for chembl query are set here
    mols = molecule.filter(molecule_type="Small molecule").filter(
        molecule_properties__num_ro5_violations=1
    )
    print(len(mols))

    with gzip.open(args.output, "wb+") as output:
        for mol in tqdm(mols):
            if mol["molecule_structures"]:
                output.write(mol["molecule_structures"]["molfile"].encode())
                output.write(b"$$$$\n")


if __name__ == "__main__":
    main()

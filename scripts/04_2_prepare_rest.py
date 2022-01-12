from rdkit import Chem
from pkasolver.data import iterate_over_acids, iterate_over_bases
import argparse
import gzip
from molvs import Standardizer
from copy import deepcopy
import pickle

s = Standardizer()

PH = 7.4


def main():
    """
    takes sdf of molecule set containing pka data and returns a pkl file of the data.
    ?Precise purpose?
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input filename, type: .sdf.gz or .sdf")
    parser.add_argument("--output", help="output filename, type: .pkl")
    args = parser.parse_args()
    input_zipped = False
    print(f"pH splitting used: {PH}")
    print("inputfile:", args.input)
    print("outputfile:", args.output)

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
            processing(suppl, args)
    else:
        with open(args.input, "rb") as fh:
            suppl = Chem.ForwardSDMolSupplier(fh, removeHs=True)
            processing(suppl, args)


def processing(suppl, args):
    GLOBAL_COUNTER = 0
    nr_of_skipped_mols = 0
    all_protonation_states_enumerated = dict()
    for nr_of_mols, mol in enumerate(suppl):

        if not mol:
            continue

        skipping_bases = 0
        skipping_acids = 0

        try:
            props = mol.GetPropsAsDict()
        except AttributeError as e:
            # this mol has no pka value
            nr_of_skipped_mols += 1
            print(e)
            continue

        pkas = []

        pkas.append(
            {
                "pka_value": float(props[f"pKa"]),
                "atom_idx": int(props[f"marvin_atom"]),
                "chembl_id": f"mol{nr_of_mols}",
            }
        )

        # calculate number of acidic and basic pka values
        nr_of_acids = sum(
            pka["pka_value"] <= PH and pka["pka_value"] > 0.5 for pka in pkas
        )
        nr_of_bases = sum(
            pka["pka_value"] > PH and pka["pka_value"] < 13.5 for pka in pkas
        )
        assert nr_of_acids + nr_of_bases <= len(pkas)

        acidic_mols_properties = [
            mol_pka
            for mol_pka in pkas
            if mol_pka["pka_value"] <= PH and mol_pka["pka_value"] > 0.5
        ]
        basic_mols_properties = [
            mol_pka
            for mol_pka in pkas
            if mol_pka["pka_value"] > PH and mol_pka["pka_value"] < 13.5
        ]

        if len(acidic_mols_properties) != nr_of_acids:
            raise RuntimeError(f"{acidic_mols_properties=}, {nr_of_acids=}")
        if len(basic_mols_properties) != nr_of_bases:
            raise RuntimeError(f"{basic_mols_properties=}, {nr_of_bases=}")

        # clear porps
        for prop in props.keys():
            mol.ClearProp(prop)

        # save values
        pka_list = []
        smiles_list = []
        counter_list = []

        # add mol at pH=PH
        mol_at_ph7 = mol
        print(Chem.MolToSmiles(mol_at_ph7))
        acidic_mols = []
        partner_mol = deepcopy(mol_at_ph7)
        (
            acidic_mols,
            nr_of_skipped_mols,
            GLOBAL_COUNTER,
            skipping_acids,
        ) = iterate_over_acids(
            acidic_mols_properties,
            nr_of_mols,
            partner_mol,
            nr_of_skipped_mols,
            pka_list,
            GLOBAL_COUNTER,
            PH,
            counter_list,
            smiles_list,
        )

        # same workflow for basic mols
        basic_mols = []
        partner_mol = deepcopy(mol_at_ph7)
        (
            basic_mols,
            nr_of_skipped_mols,
            GLOBAL_COUNTER,
            skipping_bases,
        ) = iterate_over_bases(
            basic_mols_properties,
            nr_of_mols,
            partner_mol,
            nr_of_skipped_mols,
            pka_list,
            GLOBAL_COUNTER,
            PH,
            counter_list,
            smiles_list,
        )

        # combine basic and acidic mols, skip neutral mol for acids
        combined_mols = acidic_mols + basic_mols
        if (
            len(combined_mols)
            != len(acidic_mols_properties)
            - skipping_acids
            + len(basic_mols_properties)
            - skipping_bases
        ):
            raise RuntimeError(
                combined_mols,
                acidic_mols_properties,
                skipping_acids,
                basic_mols_properties,
                skipping_bases,
            )

        if len(combined_mols) != 0:
            chembl_id = combined_mols[0][0].GetProp("CHEMBL_ID")
            print(f"CHEMBL_ID: {chembl_id}")
            for mols in combined_mols:
                if mols[0].GetProp("pKa") != mols[1].GetProp("pKa"):
                    raise AssertionError(mol[0].GetProp("pKa"), mol[1].GetProp("pKa"))

                mol1, mol2 = mols
                pka = mol1.GetProp("pKa")
                counter = mol1.GetProp("INTERNAL_ID")
                print(
                    f"{counter=}, {pka=}, {mol1.GetProp('mol-smiles')}, prot, {mol1.GetProp('epik_atom')}"
                )
                pka = mol2.GetProp("pKa")
                counter = mol2.GetProp("INTERNAL_ID")
                print(
                    f"{counter=}, {pka=}, {mol2.GetProp('mol-smiles')}, deprot, {mol1.GetProp('epik_atom')}"
                )
            print(pka_list)
            if chembl_id in all_protonation_states_enumerated.keys():
                raise RuntimeError("Repeated chembl id!")

            all_protonation_states_enumerated[chembl_id] = {
                "mols": combined_mols,
                "pKa_list": pka_list,
                "smiles_list": smiles_list,
                "counter_list": counter_list,
            }

    print(f"finished splitting {nr_of_mols} molecules")
    print(f"skipped mols: {nr_of_skipped_mols}")
    pickle.dump(all_protonation_states_enumerated, open(args.output, "wb+"))


if __name__ == "__main__":
    main()

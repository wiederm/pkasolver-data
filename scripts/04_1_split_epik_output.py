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
    takes sdf file with molcules containing epik pka predictions in their properties
    and outputs a pkl file in which pairs of molecules are deposited
    that describe the protonated and deprotonated species for each pka value.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input filename, type: .sdf.gz or .sdf")
    parser.add_argument("--output", help="output filename, type: .pkl")
    args = parser.parse_args()
    input_zipped = False
    print(f"pH splitting used: {PH}")
    print("inputfile:", args.input)
    print("outputfile:", args.output)

    #  test if it's gzipped
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

    # iterating through mols
    for nr_of_mols, mol in enumerate(suppl):
        # skit if mol can not be read
        if not mol:
            continue

        skipping_bases = 0
        skipping_acids = 0

        # test if mol has pka values
        try:
            props = mol.GetPropsAsDict()
        except AttributeError as e:
            # this mol has no pka value
            nr_of_skipped_mols += 1
            print(e)
            continue

        # count  number of pka states that epik predicted
        nr_of_protonation_states = len([s for s in props.keys() if "r_epik_pKa" in s])

        # for each protonation state extract the pka value, the atom idx and the chembl id
        properties_for_each_protonation_state = []
        for i in range(nr_of_protonation_states):
            properties_for_each_protonation_state.append(
                {
                    "pka_value": float(props[f"r_epik_pKa_{i+1}"]),
                    "atom_idx": int(props[f"i_epik_pKa_atom_{i+1}"]) - 1,
                    "chembl_id": props[f"chembl_id"],
                }
            )

        # eventhough we restricted epik predictions within a pH range of 0 to 14 there were some
        # additional pka values predicted. We introduce here a cutoff for these extrem pKa values
        upper_pka_limit = 16
        lower_pka_limit = -2

        # calculate number of acidic and basic pka values
        nr_of_acids = sum(
            pka["pka_value"] <= PH and pka["pka_value"] > lower_pka_limit
            for pka in properties_for_each_protonation_state
        )
        nr_of_bases = sum(
            pka["pka_value"] > PH and pka["pka_value"] < upper_pka_limit
            for pka in properties_for_each_protonation_state
        )
        # make sure that the sum of acid and bases equals to the number of extracted pka values
        assert nr_of_acids + nr_of_bases <= len(properties_for_each_protonation_state)
        # split properties_for_each_protonation_state into acids and bases
        acidic_mols_properties = [
            mol_pka
            for mol_pka in properties_for_each_protonation_state
            if mol_pka["pka_value"] <= PH and mol_pka["pka_value"] > lower_pka_limit
        ]
        basic_mols_properties = [
            mol_pka
            for mol_pka in properties_for_each_protonation_state
            if mol_pka["pka_value"] > PH and mol_pka["pka_value"] < upper_pka_limit
        ]
        # double check
        if len(acidic_mols_properties) != nr_of_acids:
            raise RuntimeError(f"{acidic_mols_properties=}, {nr_of_acids=}")
        if len(basic_mols_properties) != nr_of_bases:
            raise RuntimeError(f"{basic_mols_properties=}, {nr_of_bases=}")

        # clear porps for the mol at pH 7.4
        for prop in props.keys():
            mol.ClearProp(prop)

        # prepare lists in which we save the pka values, smiles and atom_idxs
        pka_list = []
        smiles_list = []
        counter_list = []

        # add mol at pH=7.4
        mol_at_ph7 = mol

        # generate states for acids and save them in acidic_mols list
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

        # generate states for bases and save them in acidic_mols list
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
        # make sure that the number of acids and bases make sense
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

        # if protonation states were present
        if len(combined_mols) != 0:
            # extract chembl id
            chembl_id = combined_mols[0][0].GetProp("CHEMBL_ID")
            print(f"CHEMBL_ID: {chembl_id}")
            # iterate over protonation states
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
    # save everything
    pickle.dump(all_protonation_states_enumerated, open(args.output, "wb+"))


if __name__ == "__main__":
    main()

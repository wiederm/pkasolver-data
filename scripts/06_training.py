import argparse
import os
import pickle
import random

import torch
from pkasolver.constants import DEVICE
from pkasolver.data import calculate_nr_of_features
from pkasolver.ml import dataset_to_dataloader
from pkasolver.ml_architecture import GINPairV1, gcn_full_training
from sklearn.model_selection import train_test_split

# all used node features
node_feat_list = [
    "element",
    "formal_charge",
    "hybridization",
    "total_num_Hs",
    "aromatic_tag",
    "total_valence",
    "total_degree",
    "is_in_ring",
    "reaction_center",
    "smarts",
]

# all possible edge features (none are used)
edge_feat_list = ["bond_type", "is_conjugated", "rotatable"]

num_node_features = calculate_nr_of_features(node_feat_list)
num_edge_features = calculate_nr_of_features(edge_feat_list)


def main():
    """
    takes training set as pkl file and trains new model or retrains existing one.

    Parameters:
    --input: set of training molecules as pyg graphs (pkl)
    --model_name: name of model architecture used
    --path: path for saving model or containing model for retraining

    optional parameters:
    --val: set of validation molecules as pyg graphs (pkl)
    --epochs: set number of training epochs (default == 1000)
    --reg: optional regularization training set (pkl)
    -r: flag for retraining model at path give by --path
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="training set path, type: .pkl")
    parser.add_argument("--path", help="training directory path")
    parser.add_argument(
        "--epochs",
        nargs="?",
        default="1000",
        help="set number of epochs (default=1000)",
    )
    parser.add_argument(
        "--reg", nargs="?", default="", help="regularization set filename"
    )
    parser.add_argument("-r", action="store_true", help="retraining run")
    args = parser.parse_args()

    if args.r:
        BATCH_SIZE = 64
    else:
        BATCH_SIZE = 512

    ######################################
    LEARNING_RATE = 0.001
    hidden_channels = 96
    NUM_EPOCHS = int(args.epochs)
    print(f"Number of epochs set to {NUM_EPOCHS}")
    model_name, model_class = "GINPairV1", GINPairV1
    # where to save training progress
    print(f"Used model: {model_name}")
    # call out that this is a fine tuning run
    if args.r:
        print("THIS IS A FINE TUNING RUN")
    # check if there is a regularization dataset present
    if args.reg:
        print(f"load regularization dataset from: {args.reg}")
    ######################################

    # generate output directory
    print(f"Write models and training progress to: {args.path}")
    os.makedirs(args.path, exist_ok=True)

    # read training set
    with open(args.input, "rb") as f:
        train_dataset = pickle.load(f)

    # reload randint if present
    if os.path.isfile(f"{args.path}/randint.pkl"):
        rs = pickle.load(open(f"{args.path}/randint.pkl", "rb"))
        print(f"Loading randing: {rs}")
    else:
        rs = random.randint(
            0, 1_000_000
        )  # save random_state to reproduce splitting if needed!
        print(rs)
        with open(f"{args.path}/randint.pkl", "wb+") as f:
            pickle.dump(rs, f)

    # randomly split training set in training and validation set
    print(f"load training dataset from: {args.input}")
    print(f"random 90:10 split is used to generate validation set.")

    # split dataset
    train_dataset, validation_dataset = train_test_split(
        train_dataset, test_size=0.1, shuffle=True, random_state=rs
    )

    train_loader = dataset_to_dataloader(train_dataset, BATCH_SIZE, shuffle=True)
    val_loader = dataset_to_dataloader(validation_dataset, BATCH_SIZE, shuffle=True)

    # if retraining
    if args.r:
        # do we have a regularization dataset
        with open(args.reg, "rb") as f:
            reg_dataset = pickle.load(f)
        reg_loader = dataset_to_dataloader(reg_dataset, 1024, shuffle=True)
    else:
        reg_loader = None

    model = model_class(
        num_node_features, num_edge_features, hidden_channels=hidden_channels
    )

    # only load model when in retraining mode, otherwise generate new one
    if args.r:
        checkpoint = torch.load(f"{args.path}/pretrained_best_model.pt")
        model.load_state_dict(checkpoint["model_state_dict"])
        prefix = "fine_tuned_"
        print("Attention: RELOADING model and extracting all layers")
        optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE,)

    else:
        prefix = "pretrained_"
        optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE,)

    # put model in training mode
    model.train()
    print(
        "Number of parameters: ",
        sum(p.numel() for p in model.parameters() if p.requires_grad == True),
    )
    print(f'Training {model_name} at epoch {model.checkpoint["epoch"]} ...')
    print(f"LR: {LEARNING_RATE}")
    print(f"Batch-size: {BATCH_SIZE}")
    print(f"Training on {DEVICE}.")
    print(f"Saving models to: {args.path}")
    gcn_full_training(
        model.to(device=DEVICE),
        train_loader,
        val_loader,
        optimizer,
        NUM_EPOCHS=NUM_EPOCHS,
        path=args.path,
        prefix=prefix,
        reg_loader=reg_loader,
    )


if __name__ == "__main__":
    main()

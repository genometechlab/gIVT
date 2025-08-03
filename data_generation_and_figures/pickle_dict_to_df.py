import sys
import os
import pandas as pd
import pickle

def update_dict(old_dict, new_dict):
    for kmer in new_dict:
        if kmer not in old_dict:
            old_dict[kmer] = {}
        for mod in new_dict[kmer]:
            if mod not in old_dict[kmer]:
                old_dict[kmer][mod] = {}
                for i in range(30):
                    old_dict[kmer][mod][round(0.7 + i / 100, 4)] = {"mod":0, "canon":0}
            for threshold in old_dict[kmer][mod]:
                old_dict[kmer][mod][threshold]["mod"] += new_dict[kmer][mod][threshold]["mod"]
                old_dict[kmer][mod][threshold]["canon"] += new_dict[kmer][mod][threshold]["canon"]
    return old_dict
    
def mer_dict_to_df(mer_dict):
    total = 0
    for mer in mer_dict:
        for mod in mer_dict[mer]:
            total += 1
    false_positive = {round(0.7 + i / 100, 4):[0] * total for i in range(30)}
    mods = [0] * total
    mer_list = [0] * total
    support = [0] * total

    i = 0
    for mer in mer_dict:
        for mod in mer_dict[mer]:
            for j in range(30):
                false_positive[round(0.7 + j/100, 4)][i] = mer_dict[mer][mod][round(0.7 + j/100, 4)]["mod"] / (mer_dict[mer][mod][round(0.7 + j/100, 4)]["canon"] + mer_dict[mer][mod][round(0.7 + j/100, 4)]["mod"])
            mods[i] = mod
            mer_list[i] = mer
            support[i] = mer_dict[mer][mod][round(0.7 + j/100, 4)]["canon"] + mer_dict[mer][mod][round(0.7 + j/100, 4)]["mod"]
            i += 1

    mer_dict_df = pd.DataFrame()
    mer_dict_df["mer"] = mer_list
    mer_dict_df["mod"] = mods
    for i in range(30):
        mer_dict_df[f"false_positive_{round(0.7 + i/100, 4)}_threshold"] = false_positive[round(0.7 + i/100, 4)]
    mer_dict_df["support"] = support
    mer_dict_df.sort_values(by="false_positive_0.7_threshold", ascending=False, inplace=True)

    return mer_dict_df


def main(pickle_dir, outpath_ref):
    merged_dict_ref = {}
    
    for file in os.listdir(pickle_dir):
        print(file)
        with open(os.path.join(pickle_dir, file), 'rb') as infile:
            sub_dict_ref = pickle.load(infile)

            merged_dict_ref = update_dict(merged_dict_ref, sub_dict_ref)

    merged_ref_df = mer_dict_to_df(merged_dict_ref)
    merged_ref_df.to_csv(outpath_ref, sep="\t", index=False)
                         

if __name__ == "__main__":
    main(sys.argv[1],
         sys.argv[2])
    
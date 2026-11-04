"""Verbatim logic of the published classification notebook ("desicion tree and plots.ipynb").

Only two changes were made so it can run as a function: the working directory
is a parameter, and the positional ``Series[i]`` look-ups were replaced by
``.iloc[i]`` (same result; the original relied on a deprecated pandas
behaviour). Used as the ground truth for the ``paper`` preset.
"""

import os
from glob import glob

import pandas as pd


def classify_like_notebook(folder):
    probbCutoff = 0.05
    fcCutoff = 0
    cases = {
        "stress-ready": {"1": [], "2": [], "3": [], "4": []},
        "shared-response": {"1": [], "2": []},
        "unique-response": {"1": [], "2": [], "3": [], "4": []},
        "opposite-response": {"1": [], "2": []},
        "no-response": {"1": [], "2": [], "3": []},
    }
    files = [os.path.basename(f) for f in glob(os.path.join(folder, "*.csv")) if "ortho" not in f]
    data = {}
    for file in files:
        dfAux = pd.read_csv(os.path.join(folder, file))
        data[file] = {}
        for i in range(len(dfAux["Unnamed: 0"])):
            og = dfAux["Row.names"][i]
            de = 0
            if dfAux["padj"][i] < probbCutoff:
                de = 1
            data[file][og] = de
    dfDE = pd.DataFrame(data)
    data = {}
    for file in files:
        dfAux = pd.read_csv(os.path.join(folder, file))
        data[file] = {}
        for i in range(len(dfAux["Unnamed: 0"])):
            og = dfAux["Row.names"][i]
            fc = 0
            if round(dfAux["log2FoldChange"][i], 2) > fcCutoff:
                fc = 1
            elif round(dfAux["log2FoldChange"][i], 2) < fcCutoff * -1:
                fc = -1
            data[file][og] = fc
    dfFC = pd.DataFrame(data)

    for i, og in enumerate(dfDE.index):
        DE_status_sly = dfDE["Sly.drought_Sly.control.csv"].iloc[i]
        DE_status_spenn = dfDE["Spenn.drought_Spenn.control.csv"].iloc[i]
        DE_status_spenn_control_sly_control = dfDE["Spenn.control_Sly.control.csv"].iloc[i]
        DE_status_spenn_drought_sly_control = dfDE["Spenn.drought_Sly.control.csv"].iloc[i]
        DE_status_sly_drought_spenn_control = dfDE["Sly.drought_Spenn.control.csv"].iloc[i]
        FC_status_sly = dfFC["Sly.drought_Sly.control.csv"].iloc[i]
        FC_status_spenn = dfFC["Spenn.drought_Spenn.control.csv"].iloc[i]
        FC_status_spenn_control_sly_control = dfFC["Spenn.control_Sly.control.csv"].iloc[i]

        if DE_status_sly == 0 and DE_status_spenn == 1:
            if FC_status_spenn == 1:
                if DE_status_spenn_drought_sly_control == 0:
                    cases["stress-ready"]["1"].append(og)
                else:
                    cases["unique-response"]["3"].append(og)
            else:
                if DE_status_spenn_drought_sly_control == 0:
                    cases["stress-ready"]["2"].append(og)
                else:
                    cases["unique-response"]["4"].append(og)
        elif DE_status_sly == 1 and DE_status_spenn == 0:
            if FC_status_sly == 1:
                if DE_status_sly_drought_spenn_control == 0:
                    cases["stress-ready"]["3"].append(og)
                else:
                    cases["unique-response"]["1"].append(og)
            else:
                if DE_status_sly_drought_spenn_control == 0:
                    cases["stress-ready"]["4"].append(og)
                else:
                    cases["unique-response"]["2"].append(og)
        elif DE_status_sly == 1 and DE_status_spenn == 1:
            if FC_status_sly == 1 and FC_status_spenn == 1:
                cases["shared-response"]["1"].append(og)
            elif FC_status_sly == -1 and FC_status_spenn == -1:
                cases["shared-response"]["2"].append(og)
            elif FC_status_sly == 1 and FC_status_spenn == -1:
                cases["opposite-response"]["1"].append(og)
            else:
                cases["opposite-response"]["2"].append(og)
        else:
            if DE_status_spenn_control_sly_control == 1 and FC_status_spenn_control_sly_control == 1:
                cases["no-response"]["1"].append(og)
            elif DE_status_spenn_control_sly_control == 1 and FC_status_spenn_control_sly_control == -1:
                cases["no-response"]["2"].append(og)
            else:
                cases["no-response"]["3"].append(og)

    rename = {"stress-ready": "stress_ready", "shared-response": "shared", "unique-response": "unique",
              "opposite-response": "opposite", "no-response": "no_response"}
    labels = {}
    for model, d in cases.items():
        for case, ogs in d.items():
            for og in ogs:
                labels[og] = f"{rename[model]}_{case}"
    return pd.Series(labels)


PAPER_FILES = {
    "a_response": "Sly.drought_Sly.control.csv",
    "b_response": "Spenn.drought_Spenn.control.csv",
    "basal": "Spenn.control_Sly.control.csv",
    "b_stress_vs_a_control": "Spenn.drought_Sly.control.csv",
    "a_stress_vs_b_control": "Sly.drought_Spenn.control.csv",
}

"""Loaders. One function per shipped table; every one returns a tidy frame.

Nothing here computes anything. Analysis lives in `diversity`, `jacobian` and
`stats`; drawing lives in `figures/`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import CLONE_RELABEL, COLONISED, DATA, group_of


def _taxa_suffix(mouse: str) -> str:
    return "_16S_family" if mouse.startswith("c_") else "_16S_taxa"


def load_family(mouse: str) -> pd.DataFrame:
    """Family-level relative abundance. Columns: Family, Time, abundance."""
    df = pd.read_csv(DATA / "16s" / "family" / f"{mouse}_family.csv")
    df = df.rename(columns={"Abundance.family": "abundance"})
    df["abundance"] = df["abundance"].fillna(0.0)
    return df[["Family", "Time", "abundance"]]


def load_diversity_16s() -> pd.DataFrame:
    """Hill numbers for every sample. Columns: mouse, index, q0, q1, qinf."""
    df = pd.read_csv(DATA / "16s" / "diversity.csv")
    df = df.rename(columns={"Sample": "mouse", "Time": "index",
                            "q_0": "q0", "q_1": "q1", "q_inf": "qinf"})
    return df[["mouse", "index", "q0", "q1", "qinf"]]


def load_diversity_barcode() -> pd.DataFrame:
    """Barcode lineage Hill numbers. Columns: mouse, index, q0, q1, qinf."""
    df = pd.read_csv(DATA / "barcodes" / "diversity.csv")
    df = df.rename(columns={"Sample": "mouse", "Time": "index",
                            "q_0": "q0", "q_1": "q1", "q_inf": "qinf"})
    return df[["mouse", "index", "q0", "q1", "qinf"]]


def load_taxa_matrix(mouse: str) -> tuple[np.ndarray, list[str], np.ndarray]:
    """The curated 16S taxon set that feeds the Jacobian.

    Returns (index, taxa names, matrix) with matrix shaped (n_index, n_taxa).
    This is a different table from `load_family`: it holds only the taxa that
    survived filtering into the interaction model, and Hill diversity computed
    from it differs slightly from the diversity table for that reason.
    """
    path = DATA / "16s" / "for_jacobian" / f"{mouse}{_taxa_suffix(mouse)}.csv"
    df = pd.read_csv(path)
    taxa = [c for c in df.columns if c != "Time"]
    return df["Time"].to_numpy(float), taxa, df[taxa].to_numpy(float)


def load_clone_loess(mouse: str, relabel: bool = True) -> pd.DataFrame:
    """LOESS consensus per clonal cluster, on the dense grid, log10 scale.

    Columns: cluster (int), clone ('C1'...), index, log10_freq.
    `relabel` applies the m2/m8 C1<->C2 swap described in config.CLONE_RELABEL.
    """
    path = DATA / "barcodes" / "clusters" / f"{mouse}_loess_clusters.csv"
    if not path.exists():
        # Control mice were never gavaged, so they have no barcode dimension and
        # their Jacobian is 16S-only. Callers branch on the empty frame.
        return pd.DataFrame(columns=["cluster", "clone", "index", "log10_freq"])
    df = pd.read_csv(path)
    df = df.rename(columns={"time": "index", "loess_value": "log10_freq"})
    df["cluster"] = df["cluster"].astype(int)
    if relabel and mouse in CLONE_RELABEL:
        df["cluster"] = df["cluster"].map(lambda c: CLONE_RELABEL[mouse].get(c, c))
    df["clone"] = "C" + df["cluster"].astype(str)
    return df[["cluster", "clone", "index", "log10_freq"]].sort_values(["cluster", "index"])


def load_cfu() -> pd.DataFrame:
    """E. coli CFU per gram. Columns: mouse, hours, cfu (zeros dropped)."""
    out = []
    for fname, mice in (("cfu_m1-m4.csv", ["m1", "m2", "m3", "m4"]),
                        ("cfu_m5-m8.csv", ["m5", "m6", "m7", "m8"])):
        df = pd.read_csv(DATA / "cfu" / fname)
        long = df.melt(id_vars="Time", value_vars=[m for m in mice if m in df],
                       var_name="mouse", value_name="cfu")
        out.append(long.rename(columns={"Time": "hours"}))
    cfu = pd.concat(out, ignore_index=True)
    return cfu[cfu["cfu"] > 0].reset_index(drop=True)


def load_jacobian_timeseries(mouse: str) -> pd.DataFrame:
    """Long-form stored Jacobian. Columns: index, window, target, driver, strength.

    strength(i, j) = cov(dz_i/dt, z_j) = J[i <- j]; `target` responds, `driver` acts.
    """
    df = pd.read_csv(DATA / "jacobian" / f"{mouse}_jacobian_timeseries.csv")
    return df.rename(columns={"time": "index", "effector_j": "driver",
                              "target_i": "target"})


def load_jacobian_matrices(mouse: str) -> list[tuple[float, np.ndarray, list[str]]]:
    """Square Jacobians per window, for eigenvalue work.

    Returns [(index, J, species)], J[i, j] with rows = targets, columns = drivers.
    """
    df = pd.read_csv(DATA / "jacobian" / f"{mouse}_jacobian_matrices_by_time.csv")
    species = list(df["target_i"].unique())
    out = []
    for idx, grp in df.groupby("time"):
        block = grp.set_index("target_i").reindex(species)
        out.append((float(idx), block[species].to_numpy(float), species))
    return sorted(out, key=lambda r: r[0])


def evaluation_indices(mouse: str) -> np.ndarray:
    """The evaluation grid, taken from the stored Jacobian timeseries.

    17 values for m1-m4 (2..17 then 18.9), 16 for m5, 15 for m6-m8. This is not
    the list of sampling slots: slot 1 has no derivative, and the last two slots
    are represented by a single dense-grid endpoint.
    """
    df = pd.read_csv(DATA / "jacobian" / f"{mouse}_jacobian_timeseries.csv",
                     usecols=["time"])
    return np.sort(df["time"].unique().astype(float))


def load_sbd_distance(mouse: str) -> pd.DataFrame:
    """Pairwise shape-based distances between clone and family trajectories."""
    return pd.read_csv(DATA / "coclustering" / f"{mouse}_sbd_distance.csv")


def load_barcode_trajectories(mice=None) -> pd.DataFrame:
    """Every barcode at every timepoint. Columns: mouse, index, ID, freq, hex.

    ~10.4 M rows, 36 MB compressed. Only Figure 1C-D needs it; load once and
    drop it afterwards.
    """
    df = pd.read_csv(DATA / "barcodes" / "trajectories.csv.gz")
    df = df.rename(columns={"Sample": "mouse", "Time": "index",
                            "Freq": "freq", "hex_line": "hex"})
    if mice is not None:
        df = df[df["mouse"].isin(mice)]
    return df


def load_rpse_alignment() -> dict[str, str]:
    """Reference and isolate RpsE protein sequences."""
    seqs, name = {}, None
    for line in (DATA / "genomics" / "rpsE_alignment.fasta").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip()
    return seqs


def load_card_hits() -> pd.DataFrame:
    """CARD/RGI strict hits for the isolate genome."""
    return pd.read_csv(DATA / "genomics" / "card_rgi_hits.csv")


__all__ = [n for n in dir() if n.startswith("load_") or n == "evaluation_indices"]

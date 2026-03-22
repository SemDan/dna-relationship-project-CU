from itertools import combinations

import pandas as pd
from tqdm.auto import tqdm

from src.similarity import similarity_score


def build_similarity_matrix(sequences_by_species: dict[str, str], use_progress: bool = True) -> pd.DataFrame:
    """
    Строит симметричную матрицу сходства между организмами.

    В ячейке [i][j] лежит нормированное сходство:
    similarity = LCS(seq_i, seq_j) / min(len(seq_i), len(seq_j))

    Parameters
    ----------
    sequences_by_species : dict[str, str]
        Словарь вида {SpeciesName: DNASequence}
    use_progress : bool
        Если True, показывает progress bar.
    """
    species_names = sorted(sequences_by_species.keys())

    similarity_df = pd.DataFrame(
        0.0,
        index=species_names,
        columns=species_names,
    )

    for species_name in species_names:
        similarity_df.loc[species_name, species_name] = 1.0

    species_pairs = list(combinations(species_names, 2))

    iterator = species_pairs
    if use_progress:
        iterator = tqdm(species_pairs, desc="Building similarity matrix")

    for left_name, right_name in iterator:
        left_sequence = sequences_by_species[left_name]
        right_sequence = sequences_by_species[right_name]

        score = similarity_score(left_sequence, right_sequence)

        similarity_df.loc[left_name, right_name] = score
        similarity_df.loc[right_name, left_name] = score

    return similarity_df


def build_rank_matrix(similarity_df: pd.DataFrame) -> pd.DataFrame:
    """
    Строит матрицу рангов родственности.

    Правила:
    - диагональ всегда = 1
    - остальные организмы ранжируются по убыванию сходства
    - меньший ранг = большее сходство
    """
    species_names = list(similarity_df.index)

    rank_df = pd.DataFrame(
        0,
        index=species_names,
        columns=species_names,
        dtype=int,
    )

    for species_name in species_names:
        # берем сходства, кроме самого себя
        row = similarity_df.loc[species_name].drop(index=species_name)

        # сортируем по убыванию сходства
        sorted_species = row.sort_values(ascending=False).index.tolist()

        # сам с собой — всегда 1
        rank_df.loc[species_name, species_name] = 1

        # остальные — начиная с 2
        for rank, other_species in enumerate(sorted_species, start=2):
            rank_df.loc[species_name, other_species] = rank

    return rank_df

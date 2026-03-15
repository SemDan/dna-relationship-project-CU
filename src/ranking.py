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
    По матрице сходства строит матрицу рангов родственности.

    Для каждой строки:
    - 1 означает максимальную близость
    - большее число означает меньшую близость

    Диагональ всегда равна 1.
    """
    species_names = list(similarity_df.index)

    rank_df = pd.DataFrame(
        0,
        index=species_names,
        columns=species_names,
        dtype=int,
    )

    for species_name in species_names:
        row_scores = similarity_df.loc[species_name].sort_values(ascending=False)

        for rank_value, other_species_name in enumerate(row_scores.index, start=1):
            rank_df.loc[species_name, other_species_name] = rank_value

    return rank_df

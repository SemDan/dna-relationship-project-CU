def lcs_length(seq1: str, seq2: str) -> int:
    """
    Вычисляет длину наибольшей общей подпоследовательности (LCS)
    между двумя последовательностями.

    Используется оптимизация по памяти:
    хранится только две строки динамического программирования.

    Time complexity: O(n * m)
    Memory complexity: O(min(n, m))
    """

    # чтобы использовать меньше памяти
    if len(seq1) < len(seq2):
        seq1, seq2 = seq2, seq1

    previous = [0] * (len(seq2) + 1)
    current = [0] * (len(seq2) + 1)

    for i in range(1, len(seq1) + 1):

        for j in range(1, len(seq2) + 1):

            if seq1[i - 1] == seq2[j - 1]:
                current[j] = previous[j - 1] + 1
            else:
                current[j] = max(previous[j], current[j - 1])

        previous, current = current, previous

    return previous[len(seq2)]

def similarity_score(seq1: str, seq2: str) -> float:
    """
    Нормированное сходство между последовательностями.
    """

    lcs = lcs_length(seq1, seq2)

    return lcs / min(len(seq1), len(seq2))

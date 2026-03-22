import matplotlib.pyplot as plt
import seaborn as sns

def plot_similarity_heatmap(similarity_df, title="Матрица сходства организмов по ДНК"):
    """
    Построение heatmap для матрицы сходства.

    Args:
        similarity_df (pd.DataFrame): Матрица сходства.
        title (str): Заголовок графика.
    """
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        similarity_df,
        annot=True,
        fmt=".2f",
        cmap="viridis",
        square=True,
        linewidths=0.5,
    )

    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()
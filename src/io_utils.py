from pathlib import Path
import re


def extract_species_name(header_line: str) -> str:
    """
    Извлекает имя организма из строки заголовка FASTA/.fna.

    Пример заголовка:
    >NC_037345.1:c55263944-55260179 RPL18 [organism=Bos taurus] [GeneID=509163] [chromosome=18]

    Возвращает:
    Bos taurus

    Если organism=... не найден, возвращает заголовок без символа '>'.
    """
    clean_header = header_line.strip()

    if clean_header.startswith(">"):
        clean_header = clean_header[1:]

    match = re.search(r"\[organism=([^\]]+)\]", clean_header)
    if match is not None:
        return match.group(1).strip()

    return clean_header


def read_fna_records(file_path: str | Path) -> list[tuple[str, str, str]]:
    """
    Читает все FASTA-записи из .fna файла.

    Возвращает список кортежей:
    [
        (species_name, header_line, sequence),
        ...
    ]
    """
    path = Path(file_path)

    with path.open("r", encoding="utf-8") as file:
        lines = [line.strip() for line in file if line.strip()]

    if not lines:
        raise ValueError(f"Файл пустой: {path}")

    records = []
    current_header = None
    current_sequence_parts = []

    for line in lines:
        if line.startswith(">"):
            if current_header is not None:
                sequence = "".join(current_sequence_parts)
                if not sequence:
                    raise ValueError(
                        f"В файле {path} найдена запись без последовательности: {current_header}"
                    )

                species_name = extract_species_name(current_header)
                records.append((species_name, current_header, sequence))

            current_header = line
            current_sequence_parts = []
        else:
            if current_header is None:
                raise ValueError(
                    f"Некорректный формат файла {path}: "
                    f"последовательность встретилась раньше заголовка."
                )

            cleaned_line = re.sub(r"[^A-Za-z]", "", line).upper()
            current_sequence_parts.append(cleaned_line)

    if current_header is not None:
        sequence = "".join(current_sequence_parts)
        if not sequence:
            raise ValueError(
                f"В файле {path} найдена запись без последовательности: {current_header}"
            )

        species_name = extract_species_name(current_header)
        records.append((species_name, current_header, sequence))

    if not records:
        raise ValueError(f"В файле {path} не найдено ни одной FASTA-записи")

    return records


def read_fna_file(file_path: str | Path) -> tuple[str, str]:
    """
    Читает .fna файл и возвращает одну 대표- последовательность для организма.

    Если в файле несколько FASTA-записей, выбирается самая длинная.
    Возвращает:
    - имя организма
    - последовательность
    """
    records = read_fna_records(file_path)

    longest_record = max(records, key=lambda record: len(record[2]))
    species_name, _, sequence = longest_record

    return species_name, sequence


def load_sequences_from_folder(folder_path: str | Path) -> dict[str, str]:
    """
    Читает все .fna файлы из папки и возвращает словарь:
    {
        "Bos taurus": "ATGC...",
        ...
    }

    Если имена организмов повторяются, выбрасывает ошибку.
    """
    folder = Path(folder_path)

    if not folder.exists():
        raise FileNotFoundError(f"Папка не найдена: {folder}")

    if not folder.is_dir():
        raise NotADirectoryError(f"Это не папка: {folder}")

    fna_files = sorted(folder.glob("*.fna"))

    if not fna_files:
        raise ValueError(f"В папке {folder} нет файлов .fna")

    sequences_by_species = {}

    for file_path in fna_files:
        species_name, sequence = read_fna_file(file_path)

        if species_name in sequences_by_species:
            raise ValueError(
                f"Повторяющееся имя организма '{species_name}'. "
                f"Проверь заголовки файлов."
            )

        sequences_by_species[species_name] = sequence

    return sequences_by_species


def build_sequences_info(sequences_by_species: dict[str, str]) -> list[dict[str, int | str]]:
    """
    Возвращает список словарей с краткой информацией по считанным последовательностям.
    """
    info_rows = []

    for species_name, sequence in sequences_by_species.items():
        info_rows.append(
            {
                "Species": species_name,
                "SequenceLength": len(sequence),
            }
        )

    info_rows.sort(key=lambda row: row["Species"])
    return info_rows
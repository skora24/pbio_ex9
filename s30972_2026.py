import random
import csv
import matplotlib.pyplot as plt

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """
    Bartosz Skóra S30972, 07.05.2026r
    Pobiera liczbę od użytkownika i sprawdza, czy mieści się w zakresie.
    Jeśli użytkownik wpisze tekst lub liczbę spoza zakresu, funkcja zapyta ponownie.
    """
    while True:
        try:
            value = int(input(prompt))
            if min_val <= value <= max_val:
                return value
            print(f"Błąd: Wartość musi być liczbą z przedziału [{min_val}, {max_val}].")
        except ValueError:
            print(f"Błąd: To nie jest liczba całkowita. Spróbuj ponownie.")

def get_nucleotide_weights() -> list:
    """
    Pyta użytkownika o procentowy udział każdego nukleotydu.
    Sprawdza, czy suma udziałów wynosi dokładnie 100%.
    """
    while True:
        print("\n--- Konfiguracja składu nukleotydowego ---")
        a = validate_positive_int("Udział A (%): ", 0, 100)
        c = validate_positive_int("Udział C (%): ", 0, 100)
        g = validate_positive_int("Udział G (%): ", 0, 100)
        t = validate_positive_int("Udział T (%): ", 0, 100)

        if a + c + g + t == 100:
            return [a, c, g, t]
        print(f"Błąd: Suma udziałów wynosi {a+c+g+t}%. Musi wynosić dokładnie 100%!")

def generate_sequence(length: int, weights: list = None) -> str:
    """
    Tworzy ciąg liter A, C, G, T.
    Używa wag (weights), aby kontrolować prawdopodobieństwo wylosowania danej litery.
    """
    nukleotydy = ['A', 'C', 'G', 'T']
    return "".join(random.choices(nukleotydy, weights=weights, k=length))

def calculate_stats(sequence: str) -> dict:
    """
    Liczy ile jest jakich liter i oblicza ich procentowy udział.
    GC-content to suma udziałów G oraz C.
    """
    length = len(sequence)
    if length == 0:
        return {"A": 0, "C": 0, "G": 0, "T": 0, "GC": 0}

    stats = {n: (sequence.count(n) / length) * 100 for n in "ACGT"}
    stats["GC"] = stats["G"] + stats["C"]
    return stats

def sliding_window_gc(sequence: str, window_size: int) -> list:
    """
    Przesuwa 'okno' po sekwencji i liczy GC dla każdego fragmentu.
    Zwraca listę par: (pozycja, procent_GC).
    """
    results = []
    # Przesuwamy się co 1 znak
    for i in range(len(sequence) - window_size + 1):
        fragment = sequence[i : i + window_size]
        gc = ((fragment.count('G') + fragment.count('C')) / window_size) * 100
        results.append((i, gc))
    return results

def find_orfs(sequence: str, min_len: int) -> list:
    """
    Szuka fragmentów zaczynających się od ATG i kończących na TAA, TAG lub TGA.
    Sprawdza tylko fragmenty o długości co najmniej min_len.
    """
    stops = ["TAA", "TAG", "TGA"]
    found = []
    for i in range(len(sequence) - 2):
        if sequence[i:i+3] == "ATG":
            # Szukamy stopu w tej samej ramce odczytu (co 3 znaki)
            for j in range(i + 3, len(sequence) - 2, 3):
                if sequence[j:j+3] in stops:
                    length = j + 3 - i
                    if length >= min_len:
                        found.append((i, length))
                    break
    return found

def insert_name(sequence: str, name: str) -> str:
    """Wstawia imię (małymi literami) w losowe miejsce sekwencji."""
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name.lower() + sequence[pos:]

def format_fasta(seq_id: str, description: str, sequence: str, width: int = 80) -> str:
    """Formatuje tekst zgodnie z zasadami pliku FASTA (nagłówek i linie po 80 znaków)."""
    lines = [f">{seq_id} {description}".strip()]
    for i in range(0, len(sequence), width):
        lines.append(sequence[i : i + width])
    return "\n".join(lines)

def main():
    print("=== Generator FASTA s30972_2026 ===")

    # 1. Wejście danych
    length = validate_positive_int("Podaj długość sekwencji (1-100000): ")
    weights = get_nucleotide_weights()

    while True:
        seq_id = input("Podaj ID sekwencji (bez spacji): ")
        if seq_id and not any(c.isspace() for c in seq_id):
            break
        print("Błąd: ID musi być jednym słowem bez spacji.")

    description = input("Podaj opis sekwencji: ")
    user_name = input("Podaj swoje imię: ")

    # 2. Logika i obliczenia
    pure_seq = generate_sequence(length, weights)
    stats = calculate_stats(pure_seq)

    # Analiza okna i ORF
    win_size = validate_positive_int(f"Szerokość okna GC (1-{length}): ", 1, length)
    gc_data = sliding_window_gc(pure_seq, win_size)

    min_orf = validate_positive_int("Minimalna długość ORF (np. 30): ", 3, length)
    orfs = find_orfs(pure_seq, min_orf)

    # 3. Zapisywanie plików
    # FASTA
    modified_seq = insert_name(pure_seq, user_name)
    with open(f"{seq_id}.fasta", "w", encoding="utf-8") as f:
        f.write(format_fasta(seq_id, description, modified_seq))

    # CSV
    with open(f"{seq_id}_gc_analysis.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pozycja", "GC_Content"])
        writer.writerows(gc_data)

    # Wykres
    plt.figure(figsize=(10, 4))
    plt.plot([x[0] for x in gc_data], [x[1] for x in gc_data], color="blue")
    plt.title(f"Analiza GC dla {seq_id}")
    plt.xlabel("Pozycja w sekwencji")
    plt.ylabel("GC %")
    plt.savefig(f"{seq_id}_wykres.png")
    plt.close()

    # 4. Raport końcowy
    print(f"\nSukces! Wygenerowano pliki: {seq_id}.fasta, {seq_id}_gc_analysis.csv, {seq_id}_wykres.png")
    print(f"Statystyki: A:{stats['A']:.2f}% C:{stats['C']:.2f}% G:{stats['G']:.2f}% T:{stats['T']:.2f}%")
    print(f"Całkowity GC-content: {stats['GC']:.2f}%")
    print(f"Liczba znalezionych ORF: {len(orfs)}")

if __name__ == "__main__":
    main()
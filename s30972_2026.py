import random

def generate_sequence(length: int) -> str:
    """Zwraca losową sekwencję DNA o zadanej długości."""
    # Używamy random.choices, aby wylosować 'length' znaków z podanej listy
    return "".join(random.choices(['A', 'C', 'G', 'T'], k=length))

def calculate_stats(sequence: str) -> dict:
    """Zwraca słownik ze statystykami sekwencji.
    Klucze: "A", "C", "G", "T" (wartości float, %),
           "GC" (wartość float, %)."""
    length = len(sequence)
    if length == 0:
        return {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0, "GC": 0.0}

    # Zliczamy wystąpienia każdego nukleotydu
    counts = {
        "A": sequence.count("A"),
        "C": sequence.count("C"),
        "G": sequence.count("G"),
        "T": sequence.count("T")
    }

    # Obliczamy procenty
    stats = {k: (v / length) * 100 for k, v in counts.items()}
    # GC-content to suma procentowej zawartości G i C
    stats["GC"] = stats["G"] + stats["C"]

    return stats

def insert_name(sequence: str, name: str) -> str:
    """Wstawia imię w losową pozycję sekwencji.
    Imię zapisane małymi literami."""
    if not sequence:
        return name.lower()

    # Losujemy pozycję od 0 do długości sekwencji
    insert_pos = random.randint(0, len(sequence))

    # Rozcinamy sekwencję na dwie części i wstawiamy imię w środek
    return sequence[:insert_pos] + name.lower() + sequence[insert_pos:]

def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Zwraca sformatowany rekord FASTA jako string."""
    # Budujemy linię nagłówkową
    header = f">{seq_id}"
    if description:
        header += f" {description}"

    lines = [header]

    # Tniemy sekwencję na kawałki o maksymalnej szerokości 'line_width'
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])

    # Łączymy linie używając znaku nowej linii
    return "\n".join(lines)

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Pobiera od użytkownika liczbę całkowitą z zakresu.
    W przypadku błędu powtarza pytanie."""
    while True:
        user_input = input(prompt)
        try:
            value = int(user_input)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Błąd: wartość musi być liczbą całkowitą z zakresu [{min_val}, {max_val}].")
        except ValueError:
            print(f"Błąd: wartość musi być liczbą całkowitą z zakresu [{min_val}, {max_val}].")

def main():
    """Główna funkcja programu odpowiedzialna za interakcję i logikę."""
    # 1. Pobieranie zwalidowanej długości
    length = validate_positive_int("Podaj długość sekwencji: ")

    # 2. Pobieranie i walidacja ID sekwencji (brak białych znaków)
    while True:
        seq_id = input("Podaj ID sekwencji: ")
        # Sprawdzamy czy w stringu są spacje, tabulatory lub czy jest pusty
        if not seq_id or any(char.isspace() for char in seq_id):
            print("Błąd: ID sekwencji nie może być puste i nie może zawierać białych znaków.")
        else:
            break

    # 3. Pobieranie opisu i imienia
    description = input("Podaj opis sekwencji: ")
    name = input("Podaj imię: ")

    # 4. Generowanie logiki biologicznej
    pure_sequence = generate_sequence(length)
    stats = calculate_stats(pure_sequence)

    # 5. Modyfikacja wygenerowanej sekwencji (Wyzwanie)
    modified_sequence = insert_name(pure_sequence, name)

    # 6. Formatowanie do standardu FASTA
    fasta_content = format_fasta(seq_id, description, modified_sequence)

    # 7. Zapis do pliku
    filename = f"{seq_id}.fasta"
    with open(filename, "w", encoding="utf-8") as file:
        file.write(fasta_content)

    # 8. Wyświetlanie wyników w konsoli
    print(f"\nSekwencja zapisana do pliku: {filename}\n")
    print(f"Statystyki sekwencji (n={length}):")
    print(f"  A: {stats['A']:.2f}%")
    print(f"  C: {stats['C']:.2f}%")
    print(f"  G: {stats['G']:.2f}%")
    print(f"  T: {stats['T']:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")

if __name__ == "__main__":
    main()
# Projekt: Beam-Focusing for Physical Layer Security in Massive MIMO Systems

**Skład zespołu:** Antoni, Jakub, Michał, Dawid, Filip

## 📌 O projekcie
Niniejsze repozytorium zawiera symulacje badawcze z zakresu bezpieczeństwa warstwy fizycznej (Physical Layer Security - PLS) dla systemów Massive MIMO w sieciach 6G. Projekt skupia się na porównaniu wydajności technik precodingu w paśmie **Sub-6 GHz** oraz paśmie fal milimetrowych **mmWave (28 GHz)** w obecności pasywnych i aktywnych podsłuchiwaczy (Eavesdroppers).

---

## 📂 Struktura Plików i Architektura

### 1. Silnik Fizyczny (Baza)
📄 **`Ghz6_band_vs_mmWave_band.m`**
* **Opis:** Główny skrypt konfiguracyjny (Phased Array System Toolbox).
* **Funkcja:** Definiuje fizykę propagacji (FSPL) oraz parametry szyków ULA (32 anteny dla 6 GHz, 512 anten dla 28 GHz).

### 2. Scenariusze Badawcze (Folder `/simulations`)

📄 **`sim_spatial_correlation.m`**
* **Analiza:** Wpływ korelacji przestrzennej na precoding ZF. Udowadnia, że przy wysokim $\rho$ sumaryczna przepustowość drastycznie spada, mimo prób zachowania sprawiedliwości wektorowej.

📄 **`sim_colluding_eavesdroppers_v2.m`**
* **Analiza:** Odporność 6 GHz vs 28 GHz na grupę współpracujących podsłuchiwaczy (MRC). Pokazuje przewagę ultra-wąskich wiązek mmWave.

📄 **`sim_artificial_noise.m`**
* **Analiza:** Skuteczność Sztucznego Szumu (AN) w ochronie pasma 6 GHz. Udowadnia, że poświęcenie części mocy na zakłócenia w przestrzeni zerowej (Null Space) drastycznie zwiększa Secrecy Rate.

📄 **`sim_pilot_contamination.m`**
* **Analiza:** Aktywny atak typu "Pilot Contamination". Pokazuje fizyczne zjawisko "porwania wiązki" i spadek bezpieczeństwa do zera przy silnym ataku hakerskim.

### 3. Galeria Wyników (Folder `/results`)
Folder zawiera automatycznie generowane wykresy w formacie `.png`. Każdy skrypt symulacyjny po zakończeniu pracy zapisuje tam zaktualizowany wykres z białym tłem, gotowy do wykorzystania w raporcie lub prezentacji.

---

## 💡 Jak uruchamiać symulacje?
1. Sklonuj repozytorium i otwórz folder w MATLAB.
2. Uruchom wybrany plik z prefiksem `sim_`.
3. Wykresy zostaną wyświetlone w nowym oknie i automatycznie zapisane w folderze `/results`.

---

## 🚀 Planowane kolejne kroki
* Badanie efektu **CSI Aging** (przesunięcie Dopplera dla poruszających się użytkowników).
* Implementacja **Filtru Wienera** (predykcja kanału) jako metody przeciwdziałania starzeniu się informacji o kanale.
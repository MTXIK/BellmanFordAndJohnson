# Алгоритмы Беллмана-Форда и Джонсона для работы с графами

## Описание

Данная программа реализует алгоритмы Беллмана-Форда и Джонсона для нахождения кратчайших путей в графах. Программа читает граф из бинарного файла, представленного в виде матрицы смежности с 16-битными целыми числами (`int16`). После загрузки графа, программа выполняет следующие операции:

1. **Обнаружение циклов с отрицательным весом** с помощью алгоритма Беллмана-Форда.
2. **Вычисление кратчайших путей между всеми парами вершин** с использованием алгоритма Джонсона (если граф не содержит циклов с отрицательным весом).
3. **Определение диаметра и радиуса графа**, а также **находжение центральных и периферийных вершин** на основе матрицы расстояний.

Результаты работы программы записываются в текстовый файл, который содержит информацию о наличии отрицательных циклов, вектор расстояний от нулевой вершины, а также графические характеристики графа.

## Требования

- **Компилятор C++**: Поддерживающий стандарт C++11 или выше (например, `g++`, `clang++`).
- **Операционная система**: Windows, Linux, macOS или любая другая ОС с поддержкой C++.

## Сборка

1. **Скачайте исходный код** программы и сохраните его в файл, например, `main.cpp`.

2. **Откройте терминал** или командную строку и перейдите в директорию с исходным кодом.

3. **Скомпилируйте программу** с помощью компилятора C++. Например, используя `g++`:

   ```bash
   g++ -o BellmanFordAndJohnson main.cpp
   ```

   Где `main.cpp` — имя файла с исходным кодом, а `BellmanFordAndJohnson` — имя создаваемого исполняемого файла.

## Использование

Запустите программу, передав в качестве обязательного аргумента имя входного бинарного файла. Опционально можно указать имя выходного файла с помощью ключа `-o`.

### Синтаксис

```bash
./BellmanFordAndJohnson <input_file.bm> [-o <output_file.txt>]
```

- `<input_file.bm>`: Имя входного бинарного файла, содержащего матрицу смежности графа.
- `-o <output_file.txt>`: (Опционально) Имя выходного текстового файла. Если ключ `-o` не указан, результаты записываются в файл `output.txt` по умолчанию.

### Пример использования

**Без указания выходного файла (результаты будут записаны в `output.txt`):**

```bash
./BellmanFordAndJohnson graph.bm
```

**С указанием выходного файла:**

```bash
./BellmanFordAndJohnson graph.bm -o results.txt
```

## Формат входного файла

Входной файл должен быть бинарным и содержать матрицу смежности графа с 16-битными целыми числами (`int16`). Структура файла следующая:

1. **Первое число (`int16`)**: Размер матрицы `n` (количество вершин в графе).
2. **Далее следует `n x n` чисел (`int16`)**: Веса ребер между вершинами. Если вес равен `0`, то предполагается отсутствие ребра между соответствующими вершинами.

### Пример создания бинарного файла

Вот пример того, как можно создать бинарный файл с матрицей смежности на языке C++:

```cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

int main() {
    int16_t n = 4; // Число вершин
    // Матрица смежности (4x4)
    std::vector<int16_t> adjMatrix = {
        0, 1, 0, 4,
        1, 0, 2, 0,
        0, 2, 0, 3,
        4, 0, 3, 0
    };
    
    std::ofstream outFile("graph.bm", std::ios::binary);
    if (!outFile) {
        std::cerr << "Не удалось открыть файл для записи." << std::endl;
        return 1;
    }
    
    outFile.write(reinterpret_cast<char*>(&n), sizeof(int16_t));
    outFile.write(reinterpret_cast<char*>(adjMatrix.data()), adjMatrix.size() * sizeof(int16_t));
    
    outFile.close();
    return 0;
}
```

Этот код создаст бинарный файл `graph.bm`, содержащий матрицу смежности для графа с 4 вершинами.

## Структура выходного файла

Выходной файл представляет собой текстовый файл, в котором содержится следующая информация:

1. **Наличие цикла с отрицательным весом:**
   - Если такой цикл обнаружен, выводится сообщение: `Graph contains a negative-weight cycle.`
   
2. **Вектор расстояний от нулевой вершины до всех остальных:**
   - Формат: `0 - <vertex>: <distance>`
   - Если расстояние до вершины невозможно (достижение невозможно), выводится `INF`.

3. **Диаметр и радиус графа:**
   - **Диаметр**: Максимальное эксцентриситет среди всех вершин.
   - **Радиус**: Минимальный эксцентриситет среди всех вершин.

4. **Список центральных и периферийных вершин:**
   - **Центральные вершины**: Вершины с эксцентриситетом, равным радиусу.
   - **Периферийные вершины**: Вершины с эксцентриситетом, равным диаметру.

### Пример выходного файла

```
Shortest paths lengths:
0 - 0: 0
0 - 1: 1
0 - 2: 3
0 - 3: 4

Diameter: 4
Radius: 2
Central vertices: 1 2
Peripheral vertices: 3
```

## Функциональность

- **Чтение из файла:** Программа открывает бинарный файл, считывает размер матрицы и саму матрицу смежности графа.

- **Обнаружение отрицательных циклов:** С помощью алгоритма Беллмана-Форда определяется наличие циклов с отрицательным весом.

- **Вычисление кратчайших путей:** Если граф не содержит отрицательных циклов, применяется алгоритм Джонсона для нахождения кратчайших путей между всеми парами вершин.

- **Графические характеристики графа:**
  - **Диаметр:** Наибольшее расстояние между любой парой вершин.
  - **Радиус:** Наименьшее расстояние, такое что каждая вершина находится на расстоянии не более этого значения от некоторой центральной вершины.
  - **Центральные вершины:** Вершины, обладающие минимальным эксцентриситетом (равным радиусу).
  - **Периферийные вершины:** Вершины, обладающие максимальным эксцентриситетом (равным диаметру).

- **Запись результатов:** Все результаты записываются в указанный выходной файл в текстовом формате.

## Тесты

Тестовые файлы находятся в папке `task1tests`. В этой папке вы найдете примеры бинарных файлов с графами различных типов, включая графы с и без циклов с отрицательным весом. Эти тесты помогут проверить корректность работы программы.

### Пример использования тестов

1. **Перейдите в папку с тестами:**

   ```bash
   cd task1tests
   ```

2. **Запустите программу с одним из тестовых файлов:**

   ```bash
   ../BellmanFordAndJohnson test1.bm -o result1.txt
   ```

3. **Просмотрите результаты в выходном файле:**

   ```bash
   cat result1.txt
   ```

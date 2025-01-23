#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <climits>
#include <cstring>
#include <queue>
#include <set>
#include <utility>

// Структура для хранения ребра графа
struct Edge {
    // Исходная вершина (source)
    int src;
    // Конечная вершина (destination)
    int dest;
    // Вес ребра (weight)
    int weight;
};

// Функция для чтения графа из бинарного файла
// Параметры:
// - inputFileName: имя входного файла
// - adjMatrix: ссылка на матрицу смежности, которая будет заполнена
bool readGraph(const std::string& inputFileName, std::vector< std::vector<int> >& adjMatrix) {
    // Открываем бинарный файл для чтения
    std::ifstream inputFile(inputFileName.c_str(), std::ios::binary);
    if (!inputFile) {
        // Если файл не удалось открыть, выводим сообщение об ошибке
        std::cerr << "Не удалось открыть входной файл." << std::endl;
        return false;
    }

    // Считываем размер матрицы (количество вершин)
    int16_t n;
    inputFile.read(reinterpret_cast<char*>(&n), sizeof(int16_t));
    int numVertices = n;

    // Инициализируем матрицу смежности графа размером numVertices x numVertices
    adjMatrix.resize(numVertices, std::vector<int>(numVertices));

    // Считываем веса ребер из бинарного файла и заполняем матрицу смежности
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            int16_t weight;
            inputFile.read(reinterpret_cast<char*>(&weight), sizeof(int16_t));
            adjMatrix[i][j] = weight;
        }
    }

    // Закрываем файл после чтения
    inputFile.close();
    return true;
}

// Алгоритм Беллмана-Форда для нахождения кратчайших путей от источника до всех остальных вершин
// Параметры:
// - numVertices: количество вершин в графе
// - edges: список всех ребер графа
// - source: индекс исходной вершины
// - distance: вектор, в который будут записаны расстояния до каждой вершины
// - negativeCycle: индикатор наличия цикла с отрицательным весом
bool bellmanFord(int numVertices, const std::vector<Edge>& edges, int source, std::vector<int>& distance, bool& negativeCycle) {
    // Инициализируем расстояния до всех вершин как бесконечные (INT_MAX)
    distance.resize(numVertices, INT_MAX);
    distance[source] = 0; // Расстояние до исходной вершины равно 0
    negativeCycle = false;

    // Выполняем релаксацию всех ребер numVertices-1 раз
    for (int i = 1; i <= numVertices - 1; ++i) {
        for (size_t j = 0; j < edges.size(); ++j) {
            int u = edges[j].src;
            int v = edges[j].dest;
            int w = edges[j].weight;
            // Обновляем расстояние, если найден более короткий путь
            if (distance[u] != INT_MAX && distance[u] + w < distance[v]) {
                distance[v] = distance[u] + w;
            }
        }
    }

    // Проверка на наличие циклов с отрицательным весом
    for (size_t j = 0; j < edges.size(); ++j) {
        int u = edges[j].src;
        int v = edges[j].dest;
        int w = edges[j].weight;
        if (distance[u] != INT_MAX && distance[u] + w < distance[v]) {
            negativeCycle = true;
            return false; // Цикл с отрицательным весом найден
        }
    }

    return true; // Циклов с отрицательным весом нет
}

// Алгоритм Джонсона для нахождения кратчайших путей между всеми парами вершин
// Параметры:
// - adjMatrix: матрица смежности исходного графа
// - dist: матрица, в которую будут записаны кратчайшие расстояния между всеми парами вершин
// - diameter: диаметр графа (максимальное эксцентриситет среди вершин)
// - radius: радиус графа (минимальный эксцентриситет среди вершин)
// - centralVertices: список центральных вершин графа
// - peripheralVertices: список периферийных вершин графа
bool johnsonAlgorithm(const std::vector< std::vector<int> >& adjMatrix, std::vector< std::vector<int> >& dist, int& diameter, int& radius, std::vector<int>& centralVertices, std::vector<int>& peripheralVertices) {
    int numVertices = adjMatrix.size();
    int newVertex = numVertices; // Индекс новой вершины, добавляемой в граф
    std::vector<Edge> edges;
    std::vector<Edge> newEdges;

    // Создание списка ребер из матрицы смежности
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            int weight = adjMatrix[i][j];
            if (weight != 0) { // Предполагается, что вес 0 означает отсутствие ребра
                Edge edge;
                edge.src = i;
                edge.dest = j;
                edge.weight = weight;
                edges.push_back(edge);
            }
        }
    }

    // Добавляем фиктивную вершину и новые ребра с весом 0 от этой вершины ко всем другим вершинам
    newEdges = edges; // Копируем существующие ребра
    for (int i = 0; i < numVertices; ++i) {
        Edge edge;
        edge.src = newVertex;
        edge.dest = i;
        edge.weight = 0;
        newEdges.push_back(edge);
    }

    // Инициализация потенциалов вершин (h) как бесконечные
    std::vector<int> h(numVertices + 1, INT_MAX);
    h[newVertex] = 0; // Потенциал новой вершины равен 0

    // Применение алгоритма Беллмана-Форда для вычисления потенциалов вершин
    for (int i = 1; i <= numVertices - 1; ++i) {
        for (size_t j = 0; j < newEdges.size(); ++j) {
            int u = newEdges[j].src;
            int v = newEdges[j].dest;
            int w = newEdges[j].weight;
            if (h[u] != INT_MAX && h[u] + w < h[v]) {
                h[v] = h[u] + w;
            }
        }
    }

    // Проверка на наличие циклов с отрицательным весом после добавления фиктивной вершины
    for (size_t j = 0; j < newEdges.size(); ++j) {
        int u = newEdges[j].src;
        int v = newEdges[j].dest;
        int w = newEdges[j].weight;
        if (h[u] != INT_MAX && h[u] + w < h[v]) {
            return false; // Цикл с отрицательным весом найден
        }
    }

    // Корректировка весов ребер с использованием потенциалов
    // Создаем список смежности с корректированными весами
    std::vector< std::vector< std::pair<int, int> > > adjList(numVertices);
    for (size_t i = 0; i < edges.size(); ++i) {
        int u = edges[i].src;
        int v = edges[i].dest;
        int w = edges[i].weight + h[u] - h[v]; // Корректировка веса
        adjList[u].emplace_back(std::make_pair(v, w));
    }

    // Применение алгоритма Дейкстры для нахождения кратчайших путей между всеми парами вершин
    dist.resize(numVertices, std::vector<int>(numVertices, INT_MAX)); // Инициализация матрицы расстояний
    for (int src = 0; src < numVertices; ++src) {
        std::vector<int> d(numVertices, INT_MAX); // Временный вектор расстояний от src
        d[src] = 0; // Расстояние до самого себя равно 0

        // Используем множество для реализации приоритетной очереди (на основе расстояния)
        std::set< std::pair<int, int> > pq;
        pq.insert(std::make_pair(0, src));

        // Основной цикл алгоритма Дейкстры
        while (!pq.empty()) {
            std::pair<int, int> top = *pq.begin(); // Выбираем вершину с минимальным расстоянием
            pq.erase(pq.begin()); // Удаляем выбранную вершину из очереди

            int u_dist = top.first; // Текущее минимальное расстояние до вершины u
            int u = top.second;     // Индекс вершины u

            if (u_dist > d[u])
                continue; // Если найденное расстояние больше текущего, пропускаем

            // Проходимся по всем соседям вершины u
            for (size_t i = 0; i < adjList[u].size(); ++i) {
                int v = adjList[u][i].first;   // Соседняя вершина
                int weight = adjList[u][i].second; // Вес ребра после корректировки
                // Проверяем возможность улучшения расстояния до вершины v
                if (d[u] != INT_MAX && d[u] + weight < d[v]) {
                    if (d[v] != INT_MAX) {
                        // Если вершина уже в очереди, удаляем ее для обновления расстояния
                        pq.erase(std::make_pair(d[v], v));
                    }
                    d[v] = d[u] + weight; // Обновляем расстояние до вершины v
                    pq.insert(std::make_pair(d[v], v)); // Вставляем вершину с обновленным расстоянием
                }
            }
        }

        // Корректируем расстояния обратно, учитывая потенциалы
        for (int i = 0; i < numVertices; ++i) {
            if (d[i] != INT_MAX) {
                dist[src][i] = d[i] - h[src] + h[i];
            }
        }
    }

    // Вычисление эксцентриситетов вершин, диаметра и радиуса графа
    std::vector<int> eccentricity(numVertices, 0); // Вектор эксцентриситетов
    diameter = 0; // Инициализация диаметра графа
    radius = INT_MAX; // Инициализация радиуса графа

    for (int i = 0; i < numVertices; ++i) {
        int maxDist = 0; // Максимальное расстояние от вершины i до остальных
        for (int j = 0; j < numVertices; ++j) {
            if (dist[i][j] != INT_MAX && dist[i][j] > maxDist) {
                maxDist = dist[i][j];
            }
        }
        eccentricity[i] = maxDist; // Записываем эксцентриситет вершины i
        if (maxDist > diameter) {
            diameter = maxDist; // Обновляем диаметр графа, если текущий эксцентриситет больше
        }
        if (maxDist < radius) {
            radius = maxDist; // Обновляем радиус графа, если текущий эксцентриситет меньше
        }
    }

    // Определение центральных и периферийных вершин графа
    for (int i = 0; i < numVertices; ++i) {
        if (eccentricity[i] == radius) {
            centralVertices.push_back(i); // Вершина является центральной
        }
        if (eccentricity[i] == diameter) {
            peripheralVertices.push_back(i); // Вершина является периферийной
        }
    }

    return true; // Алгоритм успешно завершен без обнаружения отрицательных циклов
}

// Функция для записи результатов в выходной файл
// Параметры:
// - outputFileName: имя выходного файла
// - negativeCycle: наличие отрицательного цикла (булевый флаг)
// - distance: вектор расстояний от источника
// - negativeCycleInJohnson: наличие отрицательного цикла в алгоритме Джонсона
// - diameter: диаметр графа
// - radius: радиус графа
// - centralVertices: список центральных вершин
// - peripheralVertices: список периферийных вершин
void writeResults(const std::string& outputFileName, bool negativeCycle, const std::vector<int>& distance,
                  bool negativeCycleInJohnson, int diameter, int radius,
                  const std::vector<int>& centralVertices, const std::vector<int>& peripheralVertices) {
    // Открываем файл для записи
    std::ofstream outputFile(outputFileName.c_str());
    if (!outputFile) {
        std::cerr << "Не удалось открыть выходной файл." << std::endl;
        return;
    }

    if (negativeCycle) {
        // Если найден отрицательный цикл, выводим соответствующее сообщение
        outputFile << "Graph contains a negative-weight cycle." << std::endl;
    } else {
        // Иначе выводим длины кратчайших путей
        outputFile << "Shortest paths lengths:" << std::endl;
        for (size_t i = 0; i < distance.size(); ++i) {
            if (distance[i] == INT_MAX) {
                // Если расстояние до вершины i бесконечно, выводим INF
                outputFile << "0 - " << i << ": INF" << std::endl;
            } else {
                // Иначе выводим фактическое расстояние
                outputFile << "0 - " << i << ": " << distance[i] << std::endl;
            }
        }
        outputFile << std::endl;

        if (negativeCycleInJohnson) {
            // Если в алгоритме Джонсона найден отрицательный цикл, выводим сообщение
            outputFile << "Graph contains a negative-weight cycle." << std::endl;
        } else {
            // Иначе выводим диаметр, радиус, центральные и периферийные вершины
            outputFile << "Diameter: " << diameter << std::endl;
            outputFile << "Radius: " << radius << std::endl;

            // Выводим список центральных вершин
            outputFile << "Central vertices:";
            for (size_t i = 0; i < centralVertices.size(); ++i) {
                outputFile << " " << centralVertices[i];
            }
            outputFile << std::endl;

            // Выводим список периферийных вершин
            outputFile << "Peripheral vertices:";
            for (size_t i = 0; i < peripheralVertices.size(); ++i) {
                outputFile << " " << peripheralVertices[i];
            }
            outputFile << std::endl;
        }
    }

    // Закрываем файл после записи
    outputFile.close();
}

// Главная функция программы
int main(int argc, char* argv[]) {
    std::string inputFileName;
    std::string outputFileName = "output.txt"; // Имя выходного файла по умолчанию

    // Проверяем, указано ли имя входного файла в аргументах командной строки
    if (argc < 2) {
        std::cerr << "How to use: \n   " << argv[0] << " <input_file.bm> \n   "<< argv[0] << " -o <input_file.bm> <output_file.txt>" << std::endl;
        return 1; // Завершаем программу с кодом ошибки 1
    }

    inputFileName = argv[1]; // Первый аргумент после имени программы - имя входного файла

    // Обработка дополнительных параметров командной строки для указания выходного файла
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            // Если найден флаг -o, следующий аргумент считается именем выходного файла
            outputFileName = argv[i + 1];
            i++; // Пропускаем следующий аргумент, так как он уже обработан
        }
    }

    // Чтение графа из входного файла
    std::vector< std::vector<int> > adjMatrix; // Матрица смежности графа
    if (!readGraph(inputFileName, adjMatrix)) {
        // Если чтение графа не удалось, завершаем программу с кодом ошибки 1
        return 1;
    }

    int numVertices = adjMatrix.size(); // Определяем количество вершин в графе

    // Построение списка ребер из матрицы смежности
    std::vector<Edge> edges; // Список ребер графа
    bool hasNegativeEdge = false; // Флаг наличия ребра с отрицательным весом
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            int weight = adjMatrix[i][j];
            if (weight != 0) { // Предполагается, что вес 0 означает отсутствие ребра
                Edge edge;
                edge.src = i;
                edge.dest = j;
                edge.weight = weight;
                edges.push_back(edge); // Добавляем ребро в список
                if (weight < 0) {
                    hasNegativeEdge = true; // Обнаружено ребро с отрицательным весом
                }
            }
        }
    }

    // Применение алгоритма Беллмана-Форда для обнаружения отрицательных циклов
    std::vector<int> distance; // Вектор расстояний от источника
    bool negativeCycle = false; // Флаг наличия отрицательного цикла
    if (!bellmanFord(numVertices, edges, 0, distance, negativeCycle)) {
        // Если алгоритм вернул false, значит найден отрицательный цикл
        negativeCycle = true;
    }

    // Применение алгоритма Джонсона для нахождения кратчайших путей между всеми парами вершин
    std::vector< std::vector<int> > dist; // Матрица расстояний между всеми парами вершин
    int diameter = 0; // Диаметр графа
    int radius = 0; // Радиус графа
    std::vector<int> centralVertices; // Список центральных вершин
    std::vector<int> peripheralVertices; // Список периферийных вершин
    bool negativeCycleInJohnson = false; // Флаг наличия отрицательного цикла в алгоритме Джонсона
    if (!negativeCycle) {
        // Если отрицательный цикл не найден алгоритмом Беллмана-Форда, продолжаем с алгоритмом Джонсона
        if (!johnsonAlgorithm(adjMatrix, dist, diameter, radius, centralVertices, peripheralVertices)) {
            // Если алгоритм Джонсона вернул false, значит найден отрицательный цикл
            negativeCycleInJohnson = true;
        }
    }

    // Запись результатов в выходной файл
    writeResults(outputFileName, negativeCycle, distance, negativeCycleInJohnson, diameter, radius, centralVertices, peripheralVertices);

    return 0; // Успешное завершение программы
}

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

// funcion para mapear nucleotidos a indices
int nucleotidoAIndice(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:
            cerr << "error: caracter no valido en nucleotido: " << c << endl;
            exit(1);
    }
}

// funcion para leer una secuencia desde el archivo
string leerSecuencia(const string& archivo) {
    ifstream inFile(archivo);
    if (!inFile.is_open()) {
        cerr << "error al abrir el archivo: " << archivo << endl;
        exit(1);
    }
    string secuencia;
    inFile >> secuencia;
    inFile.close();
    return secuencia;
}

// funcion para leer la matriz de similitud desde el archivo
vector<vector<int>> leerMatriz(const string& archivo, int alfabetoSize) {
    ifstream inFile(archivo);
    if (!inFile.is_open()) {
        cerr << "error al abrir el archivo: " << archivo << endl;
        exit(1);
    }

    vector<vector<int>> matriz(alfabetoSize, vector<int>(alfabetoSize));
    for (int i = 0; i < alfabetoSize; ++i) {
        for (int j = 0; j < alfabetoSize; ++j) {
            inFile >> matriz[i][j];
        }
    }
    inFile.close();
    return matriz;
}

// funcion para inicializar la matriz de puntajes
vector<vector<int>> inicializarMatriz(int n, int m, int penalizacion) {
    vector<vector<int>> matriz(n + 1, vector<int>(m + 1, 0));
    for (int i = 1; i <= n; ++i) matriz[i][0] = i * penalizacion;
    for (int j = 1; j <= m; ++j) matriz[0][j] = j * penalizacion;
    return matriz;
}

// funcion para calcular la matriz de puntajes
void llenarMatriz(vector<vector<int>>& matriz, const string& S, const string& T,
                  const vector<vector<int>>& U, int penalizacion) {
    for (size_t i = 1; i <= S.size(); ++i) {
        for (size_t j = 1; j <= T.size(); ++j) {
            int indiceS = nucleotidoAIndice(S[i - 1]);
            int indiceT = nucleotidoAIndice(T[j - 1]);

            int match = matriz[i - 1][j - 1] + U[indiceS][indiceT];
            int deleteGap = matriz[i - 1][j] + penalizacion;
            int insertGap = matriz[i][j - 1] + penalizacion;

            matriz[i][j] = max({match, deleteGap, insertGap});
        }
    }
}

// funcion para calcular el porcentaje de coincidencia
double calcularPorcentajeCoincidencia(const string& alineamientoS, const string& alineamientoT) {
    int matches = 0;
    int total = alineamientoS.size();  // longitud del alineamiento (incluye gaps)

    for (size_t i = 0; i < alineamientoS.size(); ++i) {
        if (alineamientoS[i] == alineamientoT[i] && alineamientoS[i] != '-') {
            ++matches;  // contar coincidencias
        }
    }

    return (static_cast<double>(matches) / total) * 100.0;
}

// funcion para generar un archivo Graphviz (.dot) con el porcentaje de coincidencia
void generarGraphviz(const string& alineamientoS, const string& alineamientoT, const string& archivoSalida, double porcentajeCoincidencia) {
    ofstream outFile(archivoSalida);
    if (!outFile.is_open()) {
        cerr << "Error al crear el archivo Graphviz: " << archivoSalida << endl;
        exit(1);
    }

    outFile << "digraph Alineamiento {" << endl;
    outFile << "    rankdir=LR;" << endl;
    outFile << "    node [shape=box, style=filled, fontsize=16, fontname=\"Arial\"];" << endl;

    // Leyenda
    outFile << "    subgraph cluster_legend {" << endl;
    outFile << "        label = \"Leyenda\";" << endl;
    outFile << "        fontsize = 20;" << endl;
    outFile << "        fontname = \"Arial\";" << endl;
    outFile << "        LegendMatch [label=\"Match (A|A)\", fillcolor=\"lightgreen\", style=filled];" << endl;
    outFile << "        LegendMismatch [label=\"Mismatch (A|T)\", fillcolor=\"lightcoral\", style=filled];" << endl;
    outFile << "        LegendGap [label=\"Gap (A|-)\", fillcolor=\"lightblue\", style=filled];" << endl;
    outFile << "        Porcentaje [label=\"Porcentaje de coincidencia: " << fixed << setprecision(2) << porcentajeCoincidencia << "%\", shape=plaintext];" << endl;
    outFile << "    }" << endl;

    // nodos del alineamiento
    const int nodosPorFila = 10;  // cantidad de nodos por fila
    int filaActual = 0;

    for (size_t i = 0; i < alineamientoS.size(); ++i) {
        string color;
        if (alineamientoS[i] == alineamientoT[i]) {
            color = "lightgreen";  // coincidencia
        } else if (alineamientoS[i] == '-' || alineamientoT[i] == '-') {
            color = "lightblue";  // gap
        } else {
            color = "lightcoral";  // desajuste
        }

        outFile << "    Nodo" << i << " [label=\"" << alineamientoS[i] << "|" << alineamientoT[i] << "\", fillcolor=\"" << color << "\"];" << endl;

        // agregar conexion entre nodos
        if (i > 0) {
            outFile << "    Nodo" << i - 1 << " -> Nodo" << i << ";" << endl;
        }

        // dividir en filas
        if ((i + 1) % nodosPorFila == 0) {
            outFile << "    { rank = same; ";
            for (int j = filaActual * nodosPorFila; j <= i; ++j) {
                outFile << "Nodo" << j << " ";
            }
            outFile << "}" << endl;
            filaActual++;
        }
    }

    // finalizar ultima fila
    if (alineamientoS.size() % nodosPorFila != 0) {
        outFile << "    { rank = same; ";
        for (int j = filaActual * nodosPorFila; j < alineamientoS.size(); ++j) {
            outFile << "Nodo" << j << " ";
        }
        outFile << "}" << endl;
    }

    outFile << "}" << endl;
    outFile.close();
}

// funcion para reconstruir el alineamiento
pair<string, string> reconstruirAlineamiento(const vector<vector<int>>& matriz, const string& S,
                                             const string& T, const vector<vector<int>>& U,
                                             int penalizacion) {
    string alineamientoS, alineamientoT;
    int i = S.size(), j = T.size();

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int indiceS = nucleotidoAIndice(S[i - 1]);
            int indiceT = nucleotidoAIndice(T[j - 1]);

            if (matriz[i][j] == matriz[i - 1][j - 1] + U[indiceS][indiceT]) {
                alineamientoS = S[i - 1] + alineamientoS;
                alineamientoT = T[j - 1] + alineamientoT;
                --i, --j;
                continue;
            }
        }
        if (i > 0 && matriz[i][j] == matriz[i - 1][j] + penalizacion) {
            alineamientoS = S[i - 1] + alineamientoS;
            alineamientoT = "-" + alineamientoT;
            --i;
        } else if (j > 0 && matriz[i][j] == matriz[i][j - 1] + penalizacion) {
            alineamientoS = "-" + alineamientoS;
            alineamientoT = T[j - 1] + alineamientoT;
            --j;
        } else {
            cerr << "error durante la reconstruccion del alineamiento." << endl;
            exit(1);
        }
    }

    return {alineamientoS, alineamientoT};
}

// Funcion principal
int main(int argc, char* argv[]) {
    if (argc != 9) {
        cerr << "Uso: ./programa -C1 cad1.txt -C2 cad2.txt -U matriz.txt -V penalizacion" << endl;
        return 1;
    }

    string archivoS = argv[2];
    string archivoT = argv[4];
    string archivoU = argv[6];
    int penalizacion;

    try {
        penalizacion = stoi(argv[8]);
    } catch (exception& e) {
        cerr << "error: La penalizacion debe ser un numero entero." << endl;
        return 1;
    }

    cout << "leyendo secuencias y matriz..." << endl;
    string S = leerSecuencia(archivoS);
    string T = leerSecuencia(archivoT);
    vector<vector<int>> U = leerMatriz(archivoU, 4);  // alfabeto de ADN: A, C, G, T
    cout << "lectura completada." << endl;

    cout << "inicializando matriz de puntajes..." << endl;
    vector<vector<int>> matriz = inicializarMatriz(S.size(), T.size(), penalizacion);

    cout << "llenando matriz de puntajes..." << endl;
    llenarMatriz(matriz, S, T, U, penalizacion);

    cout << "reconstruyendo alineamiento..." << endl;
    auto alineamiento = reconstruirAlineamiento(matriz, S, T, U, penalizacion);

    // calcula el porcentaje de coincidencia
    double porcentajeCoincidencia = calcularPorcentajeCoincidencia(alineamiento.first, alineamiento.second);

    cout << "generando archivo Graphviz..." << endl;
    generarGraphviz(alineamiento.first, alineamiento.second, "alineamiento.dot", porcentajeCoincidencia);

    cout << "Secuencia 1 alineada: " << alineamiento.first << endl;
    cout << "Secuencia 2 alineada: " << alineamiento.second << endl;
    cout << "Puntaje máximo: " << matriz[S.size()][T.size()] << endl;
    cout << "Porcentaje de coincidencia: " << fixed << setprecision(2) << porcentajeCoincidencia << "%" << endl;
    cout << "Archivo alineamiento.dot generado correctamente." << endl;

    return 0;
}

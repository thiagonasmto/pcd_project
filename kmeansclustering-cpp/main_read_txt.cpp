#include "./lib/kMeansCalc.cpp"
#include <cstdlib> // Para usar atoi

int main(int argc, char **argv){
    // Verificar se o tamanho do problema foi fornecido
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <dimSize>\n";
        return 1;
    }

    int dimSize = std::atoi(argv[1]); // Receber o tamanho do problema como argumento
    int clusterCount = 20;
    int iterationCount = 100;
    double threshold = 1;

    kMeansCalc<double> k {"./data/mnist_train.txt", dimSize};
    std::cout << "Avg fitness: " << k.doubleFindOptimalClusters(clusterCount, iterationCount, threshold, 1) << "\n";

    k.voidWritePointsToFile("out.txt");
    k.voidPrintSummary();
    k.voidPrintClusterSummary();

    return 0;
}


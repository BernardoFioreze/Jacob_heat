// Mede a latência e a velocidade de comunicação entre dois processos MPI (ping-pong).
//
// O processo 0 manda uma mensagem ao 1, que a devolve; isso se repete várias vezes
// para cada tamanho de mensagem. O tempo de ida é metade do tempo de ida e volta.
// Rodando os dois processos em nós diferentes, mede a rede entre os nós; no mesmo nó,
// mede a comunicação pela memória compartilhada.
//
// A mensagem de 8192 bytes é do tamanho de uma linha de borda do solver com N = 1024
// (1024 doubles), que é o que cada processo troca com cada vizinho a cada iteração.
//
// Uso: mpirun -np 2 ./pingpong

#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#include <mpi.h>
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        if (rank == 0) std::fprintf(stderr, "Use exatamente 2 processos.\n");
        MPI_Finalize();
        return 1;
    }

    // Onde está cada processo.
    char name[MPI_MAX_PROCESSOR_NAME] = {0};
    int len = 0;
    MPI_Get_processor_name(name, &len);
    char other[MPI_MAX_PROCESSOR_NAME] = {0};
    MPI_Sendrecv(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 1 - rank, 0,
                 other, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 1 - rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (rank == 0) {
        std::printf("Processo 0 em %s, processo 1 em %s (%s)\n", name, other,
                    std::string(name) == other ? "mesmo nó" : "nós diferentes");
        std::printf("%10s %12s %14s %14s\n", "bytes", "ida (us)", "MB/s", "Mb/s");
    }

    const std::vector<long> sizes = {8, 64, 512, 1024, 4096, 8192, 32768,
                                     131072, 524288, 1048576, 4194304};
    std::vector<char> buf(sizes.back(), 1);

    for (long bytes : sizes) {
        // Cerca de 64 MB trafegados por tamanho, entre 50 e 2000 repetições.
        int reps = static_cast<int>(std::clamp(64L * 1024 * 1024 / bytes, 50L, 2000L));
        int warmup = std::max(5, reps / 10);

        double t0 = 0.0;
        for (int r = -warmup; r < reps; ++r) {
            if (r == 0) {
                MPI_Barrier(MPI_COMM_WORLD);
                t0 = MPI_Wtime();
            }
            if (rank == 0) {
                MPI_Send(buf.data(), static_cast<int>(bytes), MPI_CHAR, 1, 1, MPI_COMM_WORLD);
                MPI_Recv(buf.data(), static_cast<int>(bytes), MPI_CHAR, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(buf.data(), static_cast<int>(bytes), MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(buf.data(), static_cast<int>(bytes), MPI_CHAR, 0, 1, MPI_COMM_WORLD);
            }
        }
        double oneWay = (MPI_Wtime() - t0) / reps / 2.0;  // segundos, só ida

        if (rank == 0) {
            double mbps = bytes / oneWay / 1e6;  // megabytes por segundo
            std::printf("%10ld %12.1f %14.1f %14.1f%s\n", bytes, oneWay * 1e6, mbps, mbps * 8.0,
                        bytes == 8192 ? "   <- linha de borda do solver (N=1024)" : "");
            std::printf("PINGPONG bytes=%ld ida_us=%.2f MBps=%.1f Mbps=%.1f\n",
                        bytes, oneWay * 1e6, mbps, mbps * 8.0);
        }
    }

    MPI_Finalize();
    return 0;
}

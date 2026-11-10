// poisson.cpp
// Solver da equacao de Poisson  -Laplaciano(u) = f  em [0,1]^d (d = 2 ou 3)
// com condicao de contorno de Dirichlet, diferencas finitas (5 ou 7 pontos)
// e SOR vermelho-preto. Exporta os resultados em .vti / .vtk para o ParaView.
//
// Casos disponiveis:
//   seno  : solucao manufaturada u = prod sin(pi x_i)  (tem solucao exata -> erro)
//   gauss : fonte gaussiana no centro, u = 0 no contorno (sem solucao exata)
//
// Compilar: make        Executar: ./poisson --help
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <sys/stat.h>  // mkdir (POSIX: Linux e macOS)

#include "vtk_writer.hpp"

static const double PI = 3.14159265358979323846;

struct Opcoes {
    int n;                // pontos por direcao (inclui contorno)
    int dim;              // 2 ou 3
    double tol;           // tolerancia do residuo relativo
    int maxit;            // maximo de iteracoes
    double omega;         // fator de sobre-relaxacao (<=0 => otimo teorico)
    std::string caso;     // seno | gauss
    std::string saida;    // pasta de saida
    std::string formato;  // vti | vtk | ambos
    Opcoes()
        : n(33), dim(3), tol(1e-9), maxit(20000), omega(0.0), caso("seno"),
          saida("saida"), formato("ambos") {}
};

static void ajuda() {
    std::cout <<
        "Uso: ./poisson [opcoes]\n"
        "  --n N          pontos por direcao (padrao 33)\n"
        "  --dim D        2 ou 3 (padrao 3)\n"
        "  --caso C       seno | gauss (padrao seno)\n"
        "  --tol T        tolerancia do residuo relativo (padrao 1e-9)\n"
        "  --maxit M      maximo de iteracoes (padrao 20000)\n"
        "  --omega W      fator SOR; 0 = otimo teorico (padrao 0)\n"
        "  --saida DIR    pasta de saida (padrao saida)\n"
        "  --formato F    vti | vtk | ambos (padrao ambos)\n";
}

static bool lerArgs(int argc, char** argv, Opcoes& o) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--help" || a == "-h") { ajuda(); return false; }
        if (i + 1 >= argc) { std::cerr << "Falta valor para " << a << "\n"; return false; }
        std::string v = argv[++i];
        if (a == "--n") o.n = std::atoi(v.c_str());
        else if (a == "--dim") o.dim = std::atoi(v.c_str());
        else if (a == "--caso") o.caso = v;
        else if (a == "--tol") o.tol = std::atof(v.c_str());
        else if (a == "--maxit") o.maxit = std::atoi(v.c_str());
        else if (a == "--omega") o.omega = std::atof(v.c_str());
        else if (a == "--saida") o.saida = v;
        else if (a == "--formato") o.formato = v;
        else { std::cerr << "Opcao desconhecida: " << a << "\n"; return false; }
    }
    if (o.n < 3 || (o.dim != 2 && o.dim != 3) ||
        (o.caso != "seno" && o.caso != "gauss")) {
        std::cerr << "Parametros invalidos.\n";
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    Opcoes o;
    if (!lerArgs(argc, argv, o)) return 1;

    const int nx = o.n, ny = o.n, nz = (o.dim == 3) ? o.n : 1;
    const double h = 1.0 / (o.n - 1);
    const long N = static_cast<long>(nx) * ny * nz;
    const double h2 = h * h;
    const double hd = std::pow(h, o.dim);  // peso de integracao (h^d)
    mkdir(o.saida.c_str(), 0755);  // cria a pasta de saida (ignora se ja existe)
    if (o.omega <= 0.0) o.omega = 2.0 / (1.0 + std::sin(PI * h));  // omega otimo

#define IDX(i, j, k) ((i) + nx * ((j) + ny * (k)))

    std::vector<double> u(N, 0.0), f(N, 0.0), uex(N, 0.0);
    const bool temExata = (o.caso == "seno");

    // --- Termo fonte e solucao exata ---
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const double x = i * h, y = j * h, z = (o.dim == 3) ? k * h : 0.0;
                const long p = IDX(i, j, k);
                if (o.caso == "seno") {
                    double s = std::sin(PI * x) * std::sin(PI * y);
                    if (o.dim == 3) s *= std::sin(PI * z);
                    uex[p] = s;
                    f[p] = o.dim * PI * PI * s;
                } else {
                    const double sig = 0.08;
                    double r2 = (x - .5) * (x - .5) + (y - .5) * (y - .5);
                    if (o.dim == 3) r2 += (z - .5) * (z - .5);
                    f[p] = 100.0 * std::exp(-r2 / (2 * sig * sig));
                }
            }
    // Contorno: u = 0 (ja e a condicao inicial; no caso seno a exata tambem e 0 la).

    // --- Norma de f para o residuo relativo ---
    double nf = 0.0;
    for (long p = 0; p < N; ++p) nf += f[p] * f[p];
    nf = std::sqrt(nf * hd);

    const int kmin = (o.dim == 3) ? 1 : 0, kmax = (o.dim == 3) ? nz - 1 : 1;
    const double diag = 2.0 * o.dim;

    auto residuo = [&]() {
        double s = 0.0;
        for (int k = kmin; k < kmax; ++k)
            for (int j = 1; j < ny - 1; ++j)
                for (int i = 1; i < nx - 1; ++i) {
                    double viz = u[IDX(i - 1, j, k)] + u[IDX(i + 1, j, k)] +
                                 u[IDX(i, j - 1, k)] + u[IDX(i, j + 1, k)];
                    if (o.dim == 3) viz += u[IDX(i, j, k - 1)] + u[IDX(i, j, k + 1)];
                    const double r = f[IDX(i, j, k)] - (diag * u[IDX(i, j, k)] - viz) / h2;
                    s += r * r;
                }
        return std::sqrt(s * hd) / (nf > 0 ? nf : 1.0);
    };

    // --- SOR vermelho-preto ---
    std::ofstream hist((o.saida + "/convergencia.csv").c_str());
    if (!hist) {
        std::cerr << "Crie a pasta de saida primeiro: mkdir -p " << o.saida << "\n";
        return 1;
    }
    hist << "iteracao,residuo_relativo\n";
    int it = 0;
    double res = residuo();
    hist << 0 << "," << res << "\n";
    while (it < o.maxit && res > o.tol) {
        for (int cor = 0; cor < 2; ++cor)
            for (int k = kmin; k < kmax; ++k)
                for (int j = 1; j < ny - 1; ++j)
                    for (int i = 1; i < nx - 1; ++i) {
                        if (((i + j + k) & 1) != cor) continue;
                        double viz = u[IDX(i - 1, j, k)] + u[IDX(i + 1, j, k)] +
                                     u[IDX(i, j - 1, k)] + u[IDX(i, j + 1, k)];
                        if (o.dim == 3) viz += u[IDX(i, j, k - 1)] + u[IDX(i, j, k + 1)];
                        const double gs = (viz + h2 * f[IDX(i, j, k)]) / diag;
                        u[IDX(i, j, k)] += o.omega * (gs - u[IDX(i, j, k)]);
                    }
        ++it;
        res = residuo();
        hist << it << "," << res << "\n";
    }
    hist.close();

    // --- Pos-processamento: erro e gradiente ---
    std::vector<double> err(N, 0.0), grad(3 * N, 0.0);
    double l2 = 0.0, linf = 0.0;
    if (temExata) {
        for (long p = 0; p < N; ++p) {
            err[p] = u[p] - uex[p];
            l2 += err[p] * err[p];
            linf = std::max(linf, std::fabs(err[p]));
        }
        l2 = std::sqrt(l2 * hd);
    }
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const long p = IDX(i, j, k);
                const int im = std::max(i - 1, 0), ip = std::min(i + 1, nx - 1);
                const int jm = std::max(j - 1, 0), jp = std::min(j + 1, ny - 1);
                grad[3 * p + 0] = (u[IDX(ip, j, k)] - u[IDX(im, j, k)]) / ((ip - im) * h);
                grad[3 * p + 1] = (u[IDX(i, jp, k)] - u[IDX(i, jm, k)]) / ((jp - jm) * h);
                if (o.dim == 3) {
                    const int km = std::max(k - 1, 0), kp = std::min(k + 1, nz - 1);
                    grad[3 * p + 2] = (u[IDX(i, j, kp)] - u[IDX(i, j, km)]) / ((kp - km) * h);
                }
            }

    double umin = u[0], umax = u[0];
    for (long p = 0; p < N; ++p) {
        umin = std::min(umin, u[p]);
        umax = std::max(umax, u[p]);
    }

    // --- Exportacao VTK ---
    vtkio::Grade g = {nx, ny, nz, 0, 0, 0, h, h, h};
    std::vector<vtkio::Campo> campos;
    campos.push_back(vtkio::Campo("u_numerico", 1, u));
    campos.push_back(vtkio::Campo("f", 1, f));
    if (temExata) {
        campos.push_back(vtkio::Campo("u_exato", 1, uex));
        campos.push_back(vtkio::Campo("erro", 1, err));
    }
    campos.push_back(vtkio::Campo("grad_u", 3, grad));

    if (o.formato == "vti" || o.formato == "ambos")
        vtkio::escreverVTI(o.saida + "/solucao.vti", g, campos);
    if (o.formato == "vtk" || o.formato == "ambos")
        vtkio::escreverVTKLegacy(o.saida + "/solucao.vtk", g, campos,
                                 "Poisson " + o.caso + " " + std::to_string(o.dim) + "D");

    std::ofstream r((o.saida + "/resumo.txt").c_str());
    r << std::setprecision(10)
      << "caso=" << o.caso << "\ndim=" << o.dim << "\nn=" << o.n << "\nh=" << h
      << "\nomega=" << o.omega << "\niteracoes=" << it << "\nresiduo_relativo=" << res
      << "\nu_min=" << umin << "\nu_max=" << umax << "\n";
    if (temExata) r << "erro_L2=" << l2 << "\nerro_Linf=" << linf << "\n";
    r.close();

    std::cout << std::setprecision(6) << "Poisson " << o.dim << "D, caso '" << o.caso
              << "', n=" << o.n << ", omega=" << o.omega << "\n"
              << "Iteracoes: " << it << "  residuo relativo: " << res << "\n"
              << "u em [" << umin << ", " << umax << "]\n";
    if (temExata) std::cout << "Erro L2: " << l2 << "  Linf: " << linf << "\n";
    std::cout << "Arquivos gravados em '" << o.saida << "/'\n";
    return (res <= o.tol) ? 0 : 2;
}

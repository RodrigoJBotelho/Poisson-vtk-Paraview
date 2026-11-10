// vtk_writer.hpp
// Escritores de arquivos VTK para grades estruturadas uniformes (ImageData):
//   - .vti : formato XML moderno (ParaView abre direto)
//   - .vtk : formato "legacy" STRUCTURED_POINTS (ASCII)
// Compativel com C++11. Sem dependencias externas.
#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace vtkio {

// Um campo definido nos pontos da grade.
// ncomp = 1 (escalar) ou 3 (vetor). Ordem dos dados: x varia mais rapido,
// depois y, depois z  ->  indice = i + nx*(j + ny*k)  (convencao VTK).
struct Campo {
    std::string nome;
    int ncomp;
    std::vector<double> dados;
    Campo(const std::string& n, int c, const std::vector<double>& d)
        : nome(n), ncomp(c), dados(d) {}
};

struct Grade {
    int nx, ny, nz;
    double ox, oy, oz;  // origem
    double dx, dy, dz;  // espacamento
};

inline void escreverVTI(const std::string& caminho, const Grade& g,
                        const std::vector<Campo>& campos) {
    std::ofstream f(caminho.c_str());
    if (!f) throw std::runtime_error("nao foi possivel abrir " + caminho);
    f << std::setprecision(9) << std::scientific;
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <ImageData WholeExtent=\"0 " << g.nx - 1 << " 0 " << g.ny - 1 << " 0 "
      << g.nz - 1 << "\" Origin=\"" << g.ox << " " << g.oy << " " << g.oz
      << "\" Spacing=\"" << g.dx << " " << g.dy << " " << g.dz << "\">\n";
    f << "    <Piece Extent=\"0 " << g.nx - 1 << " 0 " << g.ny - 1 << " 0 "
      << g.nz - 1 << "\">\n";
    f << "      <PointData>\n";
    for (size_t c = 0; c < campos.size(); ++c) {
        const Campo& k = campos[c];
        f << "        <DataArray type=\"Float64\" Name=\"" << k.nome
          << "\" NumberOfComponents=\"" << k.ncomp << "\" format=\"ascii\">\n";
        for (size_t i = 0; i < k.dados.size(); ++i) {
            f << k.dados[i] << ((i + 1) % 6 == 0 ? '\n' : ' ');
        }
        f << "\n        </DataArray>\n";
    }
    f << "      </PointData>\n";
    f << "      <CellData>\n      </CellData>\n";
    f << "    </Piece>\n  </ImageData>\n</VTKFile>\n";
}

inline void escreverVTKLegacy(const std::string& caminho, const Grade& g,
                              const std::vector<Campo>& campos,
                              const std::string& titulo) {
    std::ofstream f(caminho.c_str());
    if (!f) throw std::runtime_error("nao foi possivel abrir " + caminho);
    f << std::setprecision(9) << std::scientific;
    const long n = static_cast<long>(g.nx) * g.ny * g.nz;
    f << "# vtk DataFile Version 3.0\n" << titulo << "\nASCII\n";
    f << "DATASET STRUCTURED_POINTS\n";
    f << "DIMENSIONS " << g.nx << " " << g.ny << " " << g.nz << "\n";
    f << "ORIGIN " << g.ox << " " << g.oy << " " << g.oz << "\n";
    f << "SPACING " << g.dx << " " << g.dy << " " << g.dz << "\n";
    f << "POINT_DATA " << n << "\n";
    for (size_t c = 0; c < campos.size(); ++c) {
        const Campo& k = campos[c];
        if (k.ncomp == 1) {
            f << "SCALARS " << k.nome << " double 1\nLOOKUP_TABLE default\n";
            for (size_t i = 0; i < k.dados.size(); ++i)
                f << k.dados[i] << ((i + 1) % 6 == 0 ? '\n' : ' ');
        } else {
            f << "VECTORS " << k.nome << " double\n";
            for (size_t i = 0; i < k.dados.size(); i += 3)
                f << k.dados[i] << " " << k.dados[i + 1] << " " << k.dados[i + 2] << "\n";
        }
        f << "\n";
    }
}

}  // namespace vtkio

#endif

#!/usr/bin/env python3
"""
posprocessar.py - Pos-processador do solver de Poisson (sem precisar do ParaView).

Le o arquivo .vti (ASCII) gerado pelo solver C++, extrai os campos e gera:
  - fatias (planos medios) de u, erro e gradiente
  - perfil central numerico x exato
  - historico de convergencia
  - relatorio de estatisticas (relatorio.txt)

Dependencias: numpy, matplotlib.
Uso:  python3 scripts/posprocessar.py saida/solucao.vti --destino docs/images
"""
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")  # sem janela: funciona em servidor, Colab etc.
import matplotlib.pyplot as plt
import numpy as np


def ler_vti(caminho):
    """Parser minimo de .vti ASCII (ImageData). Retorna (campos, dims, espacamento)."""
    with open(caminho, "r") as f:
        texto = f.read()
    ext = [int(v) for v in re.search(r'WholeExtent="([^"]+)"', texto).group(1).split()]
    esp = [float(v) for v in re.search(r'Spacing="([^"]+)"', texto).group(1).split()]
    nx, ny, nz = ext[1] - ext[0] + 1, ext[3] - ext[2] + 1, ext[5] - ext[4] + 1
    campos = {}
    padrao = re.compile(
        r'<DataArray[^>]*Name="([^"]+)"[^>]*NumberOfComponents="(\d+)"[^>]*>(.*?)</DataArray>',
        re.S)
    for nome, nc, corpo in padrao.findall(texto):
        v = np.array(corpo.split(), dtype=float)
        nc = int(nc)
        # ordem VTK: x mais rapido -> reshape (nz, ny, nx[, nc])
        campos[nome] = v.reshape((nz, ny, nx, nc)) if nc > 1 else v.reshape((nz, ny, nx))
    return campos, (nx, ny, nz), esp


def ler_resumo(pasta):
    r = {}
    caminho = os.path.join(pasta, "resumo.txt")
    if os.path.exists(caminho):
        for linha in open(caminho):
            if "=" in linha:
                k, v = linha.strip().split("=", 1)
                r[k] = v
    return r


def fatia_media(campo):
    return campo[campo.shape[0] // 2]  # plano z = meio (ou o unico plano, em 2D)


def painel(ax, dados, titulo, cmap="viridis", rotulo=""):
    im = ax.imshow(dados, origin="lower", extent=[0, 1, 0, 1], cmap=cmap)
    ax.set_title(titulo)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.colorbar(im, ax=ax, shrink=0.85, label=rotulo)


def tres_planos_3d(u, destino):
    """Vista 3D com tres planos ortogonais (aproxima a renderizacao do ParaView)."""
    nz, ny, nx = u.shape
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    z = np.linspace(0, 1, nz)
    norma = plt.Normalize(u.min(), u.max())
    cmap = plt.cm.viridis
    fig = plt.figure(figsize=(6.5, 5.5))
    ax = fig.add_subplot(111, projection="3d")
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y, np.full_like(X, 0.5), rstride=1, cstride=1, shade=False,
                    facecolors=cmap(norma(u[nz // 2, :, :])))
    X, Z = np.meshgrid(x, z)
    ax.plot_surface(X, np.full_like(X, 0.5), Z, rstride=1, cstride=1, shade=False,
                    facecolors=cmap(norma(u[:, ny // 2, :])))
    Y, Z = np.meshgrid(y, z)
    ax.plot_surface(np.full_like(Y, 0.5), Y, Z, rstride=1, cstride=1, shade=False,
                    facecolors=cmap(norma(u[:, :, nx // 2])))
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title("u em tres planos ortogonais")
    m = plt.cm.ScalarMappable(norm=norma, cmap=cmap)
    fig.colorbar(m, ax=ax, shrink=0.65, label="u")
    fig.tight_layout()
    fig.savefig(os.path.join(destino, "tres_planos_3d.png"), dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vti", help="arquivo .vti gerado pelo solver")
    ap.add_argument("--destino", default="docs/images", help="pasta das imagens")
    args = ap.parse_args()

    os.makedirs(args.destino, exist_ok=True)
    pasta = os.path.dirname(args.vti) or "."
    campos, (nx, ny, nz), (hx, hy, hz) = ler_vti(args.vti)
    resumo = ler_resumo(pasta)
    u = campos["u_numerico"]
    tem_exata = "u_exato" in campos
    dim = 3 if nz > 1 else 2

    # 1) Solucao, erro e |grad u| no plano medio
    n_pan = 3 if tem_exata else 2
    fig, eixos = plt.subplots(1, n_pan, figsize=(5 * n_pan, 4.4))
    painel(eixos[0], fatia_media(u), "u numerico (plano medio)", "viridis", "u")
    g = np.linalg.norm(campos["grad_u"], axis=-1)
    painel(eixos[1], fatia_media(g), "|grad u| (plano medio)", "magma", "|grad u|")
    if tem_exata:
        painel(eixos[2], fatia_media(campos["erro"]), "erro (numerico - exato)", "coolwarm", "erro")
    fig.suptitle("Poisson %dD - n=%d" % (dim, nx))
    fig.tight_layout()
    fig.savefig(os.path.join(args.destino, "fatias.png"), dpi=140)
    plt.close(fig)

    if dim == 3:
        tres_planos_3d(u, args.destino)

    # 2) Perfil central numerico x exato (reta y = z = 1/2, x variando)
    x = np.linspace(0, 1, nx)
    jc, kc = ny // 2, nz // 2
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x, u[kc, jc, :], "o", ms=3, label="numerico")
    if tem_exata:
        ax.plot(x, campos["u_exato"][kc, jc, :], "-", label="exato")
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.set_title("Perfil central (y = z = 0.5)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(args.destino, "perfil_central.png"), dpi=140)
    plt.close(fig)

    # 3) Convergencia do solver iterativo
    csv = os.path.join(pasta, "convergencia.csv")
    if os.path.exists(csv):
        d = np.loadtxt(csv, delimiter=",", skiprows=1)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.semilogy(d[:, 0], d[:, 1])
        ax.set_xlabel("iteracao (SOR vermelho-preto)")
        ax.set_ylabel("residuo relativo")
        ax.set_title("Convergencia do solver")
        ax.grid(alpha=0.3, which="both")
        fig.tight_layout()
        fig.savefig(os.path.join(args.destino, "convergencia.png"), dpi=140)
        plt.close(fig)

    # 4) Relatorio
    linhas = ["Relatorio de pos-processamento", "arquivo: %s" % args.vti,
              "grade: %d x %d x %d  (h = %.5g)" % (nx, ny, nz, hx),
              "u_min = %.6g   u_max = %.6g" % (u.min(), u.max()),
              "|grad u|_max = %.6g" % g.max()]
    if tem_exata:
        e = campos["erro"]
        peso = hx ** dim
        linhas += ["erro L2   = %.6e" % np.sqrt(np.sum(e ** 2) * peso),
                   "erro Linf = %.6e" % np.abs(e).max()]
    for k in ("iteracoes", "residuo_relativo", "omega"):
        if k in resumo:
            linhas.append("%s = %s" % (k, resumo[k]))
    texto = "\n".join(linhas)
    with open(os.path.join(pasta, "relatorio.txt"), "w") as f:
        f.write(texto + "\n")
    print(texto)
    print("Imagens salvas em", args.destino)


if __name__ == "__main__":
    main()

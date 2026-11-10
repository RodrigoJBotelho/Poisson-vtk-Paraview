#!/usr/bin/env python3
"""
estudo_convergencia.py - Verifica a ordem de precisao do solver.

Roda o executavel ./poisson para varias malhas (caso 'seno', que tem solucao exata),
le o erro L2/Linf do resumo.txt e plota erro x h em escala log-log.
Para diferencas finitas centradas esperamos ordem 2.

Uso:  python3 scripts/estudo_convergencia.py --dim 3 --malhas 9 17 33 49 65
"""
import argparse
import os
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="./poisson")
    ap.add_argument("--dim", type=int, default=3)
    ap.add_argument("--malhas", type=int, nargs="+", default=[9, 17, 33, 49])
    ap.add_argument("--destino", default="docs/images")
    args = ap.parse_args()

    hs, l2s, lis = [], [], []
    for n in args.malhas:
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run([args.exe, "--dim", str(args.dim), "--n", str(n), "--caso", "seno",
                            "--saida", tmp, "--formato", "vti", "--tol", "1e-11"],
                           check=True, stdout=subprocess.DEVNULL)
            r = dict(l.strip().split("=") for l in open(os.path.join(tmp, "resumo.txt")))
        hs.append(float(r["h"]))
        l2s.append(float(r["erro_L2"]))
        lis.append(float(r["erro_Linf"]))
        print("n=%3d  h=%.4f  L2=%.3e  Linf=%.3e" % (n, hs[-1], l2s[-1], lis[-1]))

    hs, l2s, lis = map(np.array, (hs, l2s, lis))
    ordem_l2 = np.polyfit(np.log(hs), np.log(l2s), 1)[0]
    ordem_li = np.polyfit(np.log(hs), np.log(lis), 1)[0]
    print("Ordem observada: L2 = %.2f   Linf = %.2f" % (ordem_l2, ordem_li))

    os.makedirs(args.destino, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.loglog(hs, l2s, "o-", label="erro L2 (ordem %.2f)" % ordem_l2)
    ax.loglog(hs, lis, "s-", label="erro Linf (ordem %.2f)" % ordem_li)
    ax.loglog(hs, l2s[0] * (hs / hs[0]) ** 2, "k--", lw=1, label="referencia h^2")
    ax.set_xlabel("h")
    ax.set_ylabel("erro")
    ax.set_title("Ordem de convergencia - Poisson %dD" % args.dim)
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.destino, "ordem_convergencia.png"), dpi=140)


if __name__ == "__main__":
    main()

#!/usr/bin/env pvpython
"""
renderizar_paraview.py - Gera renderizacoes 2D/3D profissionais no ParaView, sem abrir a GUI.

Uso (precisa do ParaView instalado; pvpython vem junto):
    pvpython --force-offscreen-rendering scripts/renderizar_paraview.py saida/solucao.vti docs/images

Imagens geradas na pasta de destino (prefixo pv_):
    pv_fatia.png          plano medio colorido por u
    pv_isosuperficies.png isosuperficies de u com contorno do dominio (so 3D)
    pv_tres_planos.png    tres planos ortogonais (so 3D)
    pv_gradiente.png      glifos (setas) do gradiente sobre o plano medio
    pv_erro.png           plano medio colorido pelo erro (se existir)

Nota: foi escrito para a API paraview.simple (ParaView 5.x). Ajustes de nomes de
propriedades podem ser necessarios entre versoes; veja a Python Shell do ParaView
(View > Python Shell) ou 'Tools > Start Trace' para descobrir o equivalente.
"""
import os
import sys

from paraview.simple import *  # noqa: F401,F403

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

arquivo = sys.argv[1]
destino = sys.argv[2] if len(sys.argv) > 2 else "docs/images"
os.makedirs(destino, exist_ok=True)

LARGURA, ALTURA = 1600, 1000

# --- Leitura (.vti usa o leitor XML; .vtk usa o leitor legacy) ---
if arquivo.lower().endswith(".vti"):
    leitor = XMLImageDataReader(FileName=[arquivo])
else:
    leitor = LegacyVTKReader(FileNames=[arquivo])
leitor.UpdatePipeline()

info = leitor.GetDataInformation()
ext = info.GetExtent()  # (xmin, xmax, ymin, ymax, zmin, zmax)
limites = info.GetBounds()
eh3d = (ext[5] - ext[4]) > 0
centro = [(limites[0] + limites[1]) / 2.0, (limites[2] + limites[3]) / 2.0,
          (limites[4] + limites[5]) / 2.0]
campos = [leitor.PointData.GetArray(i).GetName() for i in range(len(leitor.PointData))]
print("Campos encontrados:", campos)

# --- Vista ---
vista = CreateView("RenderView")
vista.ViewSize = [LARGURA, ALTURA]
vista.OrientationAxesVisibility = 1
try:  # ParaView >= 5.10: fundo branco para publicacao
    vista.UseColorPaletteForBackground = 0
except Exception:
    pass
vista.Background = [1.0, 1.0, 1.0]


def esconder_tudo():
    for fonte in GetSources().values():
        Hide(fonte, vista)


def colorir(exibicao, campo, titulo=None, paleta="Viridis (matplotlib)"):
    ColorBy(exibicao, ("POINTS", campo))
    lut = GetColorTransferFunction(campo)
    try:
        lut.ApplyPreset(paleta, True)
    except Exception:
        pass
    exibicao.RescaleTransferFunctionToDataRange(True, False)
    exibicao.SetScalarBarVisibility(vista, True)
    barra = GetScalarBar(lut, vista)
    barra.Title = titulo or campo
    barra.ComponentTitle = ""
    barra.TitleColor = [0, 0, 0]
    barra.LabelColor = [0, 0, 0]
    barra.TitleFontSize = 20
    barra.LabelFontSize = 16
    return lut


def camera_3d():
    vista.CameraPosition = [centro[0] + 2.2, centro[1] - 2.6, centro[2] + 1.8]
    vista.CameraFocalPoint = centro
    vista.CameraViewUp = [0, 0, 1]
    vista.ResetCamera()


def camera_2d():
    vista.CameraPosition = [centro[0], centro[1], centro[2] + 5.0]
    vista.CameraFocalPoint = centro
    vista.CameraViewUp = [0, 1, 0]
    vista.ResetCamera()


def salvar(nome):
    Render(vista)
    caminho = os.path.join(destino, nome)
    SaveScreenshot(caminho, vista, ImageResolution=[LARGURA, ALTURA])
    print("salvo:", caminho)


# 1) Plano medio colorido por u ---------------------------------------------
esconder_tudo()
if eh3d:
    plano = Slice(Input=leitor)
    plano.SliceType = "Plane"
    plano.SliceType.Origin = centro
    plano.SliceType.Normal = [0, 0, 1]
    fonte_plano = plano
else:
    fonte_plano = leitor
rep = Show(fonte_plano, vista)
colorir(rep, "u_numerico", "u")
camera_2d()
salvar("pv_fatia.png")

# 2) Isosuperficies + contorno do dominio (3D) ------------------------------
if eh3d:
    esconder_tudo()
    contorno = Outline(Input=leitor)
    rc = Show(contorno, vista)
    rc.AmbientColor = [0, 0, 0]
    rc.DiffuseColor = [0, 0, 0]
    rng = leitor.PointData["u_numerico"].GetRange()
    iso = Contour(Input=leitor)
    iso.ContourBy = ["POINTS", "u_numerico"]
    iso.Isosurfaces = [rng[0] + f * (rng[1] - rng[0]) for f in (0.2, 0.4, 0.6, 0.8)]
    ri = Show(iso, vista)
    ri.Opacity = 0.55
    colorir(ri, "u_numerico", "u")
    camera_3d()
    salvar("pv_isosuperficies.png")

    # 3) Tres planos ortogonais --------------------------------------------
    esconder_tudo()
    Show(contorno, vista)
    for normal in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
        p = Slice(Input=leitor)
        p.SliceType = "Plane"
        p.SliceType.Origin = centro
        p.SliceType.Normal = normal
        rp = Show(p, vista)
        colorir(rp, "u_numerico", "u")
    camera_3d()
    salvar("pv_tres_planos.png")

# 4) Gradiente como setas sobre o plano medio ---------------------------------
if "grad_u" in campos:
    esconder_tudo()
    base = fonte_plano
    rb = Show(base, vista)
    colorir(rb, "u_numerico", "u", "Cool to Warm")
    gl = Glyph(Input=base, GlyphType="Arrow")
    gl.OrientationArray = ["POINTS", "grad_u"]
    gl.ScaleArray = ["POINTS", "No scale array"]
    gl.ScaleFactor = 0.04
    gl.GlyphMode = "Every Nth Point"
    gl.Stride = 3
    rg = Show(gl, vista)
    rg.AmbientColor = [0.1, 0.1, 0.1]
    rg.DiffuseColor = [0.1, 0.1, 0.1]
    ColorBy(rg, None)
    camera_2d()
    salvar("pv_gradiente.png")

# 5) Erro numerico ---------------------------------------------------------
if "erro" in campos:
    esconder_tudo()
    re_ = Show(fonte_plano, vista)
    colorir(re_, "erro", "erro (numerico - exato)", "Cool to Warm")
    camera_2d()
    salvar("pv_erro.png")

print("Pronto. Imagens em", destino)

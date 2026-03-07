# =========================================================
# ANALISIS DE ARCHIVOS NIFTI EN R
# 1) T1 estructural:
#    - Tamaño del volumen
#    - Vistas ortogonales
#    - Histograma
#    - Intervalo aproximado de tonos de gris del cuerpo calloso
#
# 2) BOLD fMRI:
#    - Tamaño del volumen 3D
#    - Numero de volumenes
# =========================================================

# -----------------------------
# Instalar/cargar paquete
# -----------------------------
if (!requireNamespace("RNifti", quietly = TRUE)) install.packages("RNifti")
library(RNifti)

# =========================================================
# PARTE 1. T1 ESTRUCTURAL
# =========================================================

cat("=====================================================\n")
cat("PARTE 1. ANALISIS T1 ESTRUCTURAL\n")
cat("=====================================================\n\n")

# Ruta del archivo T1
nii_path_t1 <- "sub-01_T1w.nii.gz"

# Leer imagen
img_t1 <- readNifti(nii_path_t1)

# Dimensiones y tamaño de voxel
dim_img_t1 <- dim(img_t1)              # (X, Y, Z)
vox_t1 <- voxelSize(img_t1)            # tamaño voxel (mm)
phys_mm_t1 <- dim_img_t1 * vox_t1      # tamaño fisico aproximado (mm)

cat("Dimensiones (voxeles) [X,Y,Z]: ", paste(dim_img_t1, collapse = " x "), "\n", sep = "")
cat("Tamaño voxel (mm) [dx,dy,dz]: ", paste(round(vox_t1, 4), collapse = " x "), "\n", sep = "")
cat("Tamaño fisico aprox (mm): ", paste(round(phys_mm_t1, 2), collapse = " x "), "\n\n", sep = "")

# -----------------------------
# 1) Tamaño del volumen
# -----------------------------
n_vox_t1 <- prod(dim_img_t1)
vox_vol_mm3_t1 <- prod(vox_t1)               # mm^3 por voxel
total_vol_mm3_t1 <- n_vox_t1 * vox_vol_mm3_t1
total_vol_ml_t1 <- total_vol_mm3_t1 / 1000   # 1 mL = 1000 mm^3

cat("1) TAMAÑO DEL VOLUMEN\n")
cat("Numero total de voxeles: ", format(n_vox_t1, big.mark = ","), "\n", sep = "")
cat("Volumen por voxel (mm^3): ", round(vox_vol_mm3_t1, 4), "\n", sep = "")
cat("Volumen total (mm^3): ", format(round(total_vol_mm3_t1, 2), big.mark = ","), "\n", sep = "")
cat("Volumen total (mL): ", round(total_vol_ml_t1, 2), "\n\n", sep = "")

# -----------------------------
# 2) Vistas ortogonales en (i0, j0, k0)
# -----------------------------
i0 <- round(dim_img_t1[1] / 2)
j0 <- round(dim_img_t1[2] / 2)
k0 <- round(dim_img_t1[3] / 2)

cat("2) VISTAS ORTOGONALES\n")
cat("Mostrando en (i0, j0, k0) = ", i0, ", ", j0, ", ", k0, "\n\n", sep = "")

orthographic(img_t1, xyz = c(i0, j0, k0))

# -----------------------------
# 3) Histograma del volumen
# -----------------------------
vals_t1 <- as.numeric(img_t1)
vals_t1 <- vals_t1[is.finite(vals_t1)]
vals_t1 <- vals_t1[vals_t1 > 0]   # quitar fondo

# Si son demasiados voxeles, muestreo para acelerar
set.seed(123)
max_n <- 2e6
if (length(vals_t1) > max_n) {
  vals_plot_t1 <- sample(vals_t1, max_n)
} else {
  vals_plot_t1 <- vals_t1
}

hist(
  vals_plot_t1,
  breaks = 200,
  main = "Histograma de intensidades (voxeles > 0) - T1",
  xlab = "Intensidad (tono de gris)",
  ylab = "Frecuencia"
)

cat("3) HISTOGRAMA DEL VOLUMEN\n")
cat("Se genero el histograma de intensidades del volumen T1.\n\n")

# -----------------------------
# 4) Intervalo aproximado de tonos de gris del cuerpo calloso
# -----------------------------
# Rebanada sagital media
i_mid <- round(dim_img_t1[1] / 2)

# ROI aproximada en el plano sagital medio
# j = anterior-posterior
# k = inferior-superior
j1 <- round(0.35 * dim_img_t1[2])
j2 <- round(0.65 * dim_img_t1[2])
k1 <- round(0.45 * dim_img_t1[3])
k2 <- round(0.70 * dim_img_t1[3])

roi_t1 <- img_t1[i_mid, j1:j2, k1:k2]
roi_vals_t1 <- as.numeric(roi_t1)
roi_vals_t1 <- roi_vals_t1[is.finite(roi_vals_t1) & roi_vals_t1 > 0]

p_t1 <- quantile(roi_vals_t1, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)

cat("4) INTERVALO APROXIMADO DE TONOS DE GRIS DEL CUERPO CALLOSO\n")
cat("ROI (sagital medio) usada:\n")
cat("  i = ", i_mid, "\n", sep = "")
cat("  j = ", j1, " a ", j2, "\n", sep = "")
cat("  k = ", k1, " a ", k2, "\n\n", sep = "")
cat("Tonos de gris aproximados en esta ROI:\n")
cat("  p05 = ", round(p_t1[1], 2), "\n", sep = "")
cat("  p50 = ", round(p_t1[2], 2), "\n", sep = "")
cat("  p95 = ", round(p_t1[3], 2), "\n", sep = "")
cat("Intervalo sugerido ~ [p05, p95] = [", round(p_t1[1], 2), ", ", round(p_t1[3], 2), "]\n\n", sep = "")

hist(
  roi_vals_t1,
  breaks = 100,
  main = "Histograma ROI aprox. Cuerpo Calloso - sagital medio",
  xlab = "Intensidad (tono de gris)",
  ylab = "Frecuencia"
)

# =========================================================
# PARTE 2. BOLD fMRI
# =========================================================

cat("=====================================================\n")
cat("PARTE 2. ANALISIS BOLD fMRI\n")
cat("=====================================================\n\n")

# Ruta del archivo BOLD
nii_path_bold <- "sub-01_task-restaurative_bold.nii.gz"

# Leer imagen
img_bold <- readNifti(nii_path_bold)

# Dimensiones y voxel
dim_img_bold <- dim(img_bold)        # (X, Y, Z, T)
vox_bold <- voxelSize(img_bold)

cat("Dimensiones del archivo BOLD: ", paste(dim_img_bold, collapse = " x "), "\n", sep = "")
cat("Tamaño voxel (mm): ", paste(round(vox_bold, 4), collapse = " x "), "\n\n", sep = "")

# -----------------------------
# 1) Tamaño del volumen
# -----------------------------
dim_spatial_bold <- dim_img_bold[1:3]
n_vox_3D_bold <- prod(dim_spatial_bold)
vox_vol_mm3_bold <- prod(vox_bold[1:3])
vol_3D_mm3_bold <- n_vox_3D_bold * vox_vol_mm3_bold
vol_3D_ml_bold <- vol_3D_mm3_bold / 1000

cat("1) TAMAÑO DEL VOLUMEN BOLD (por cada volumen 3D)\n")
cat("Dimensiones espaciales (X,Y,Z): ", paste(dim_spatial_bold, collapse = " x "), "\n", sep = "")
cat("Numero de voxeles por volumen 3D: ", format(n_vox_3D_bold, big.mark = ","), "\n", sep = "")
cat("Volumen total por cada 3D (mm^3): ", format(round(vol_3D_mm3_bold, 2), big.mark = ","), "\n", sep = "")
cat("Volumen total por cada 3D (mL): ", round(vol_3D_ml_bold, 2), "\n\n", sep = "")

# -----------------------------
# 2) Numero de volumenes
# -----------------------------
if (length(dim_img_bold) == 4) {
  n_vol_bold <- dim_img_bold[4]
} else {
  n_vol_bold <- 1
}

cat("2) NUMERO DE VOLUMENES\n")
cat("Numero de volumenes (dimension temporal): ", n_vol_bold, "\n", sep = "")

cat("\n=====================================================\n")
cat("ANALISIS COMPLETADO\n")
cat("=====================================================\n")

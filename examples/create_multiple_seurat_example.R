# Exemple : Remplacement du code manuel par CreateMultipleSeurat
# 
# Cet exemple montre comment remplacer le workflow manuel par la nouvelle fonction

library(Seurat)

# ==========================================
# AVANT : Code manuel (ce que vous faisiez)
# ==========================================

# data_dir <- "/path/to/cellranger/outputs/"
# list_data <- c("sample1", "sample2", "sample3")
# 
# paths <- paste0(data_dir, list_data, "/filtered_feature_bc_matrix/")
# NeuroProx <- 
#   mclapply(X = paths, 
#          FUN =  function(paths){
#            CreateSeuratObject(counts = Read10X(data.dir = paths), 
#                               project = str_split(paths, pattern = "/")[[1]][7])
#          },
#          mc.cores = detectCores()
#          )
# names(NeuroProx) <- str_split(paths, pattern = "/", simplify = T)[, 7]
# qsave(x = NeuroProx, file = "20250528_NeuroProx_full.qs", nthreads = 10)
# NeuroProx <- qread(file = "20250528_NeuroProx_full.qs", nthreads = 10)

# ==========================================
# MAINTENANT : Version simplifiée
# ==========================================

# Définir les paramètres
data_dir <- "/path/to/cellranger/outputs/"
list_data <- c("sample1", "sample2", "sample3")

# Créer les objets Seurat en une seule ligne
NeuroProx <- CreateMultipleSeurat(
  data_dir = data_dir,
  list_data = list_data,
  save_file = "20250528_NeuroProx_full.qs",
  mc.cores = parallel::detectCores(),
  nthreads = 10
)

# Charger plus tard (remplace qread)
NeuroProx <- LoadMultipleSeurat("20250528_NeuroProx_full.qs", nthreads = 10)

# ==========================================
# Exemple avec des paramètres personnalisés
# ==========================================

# Avec filtrage des cellules et gènes
NeuroProx_filtered <- CreateMultipleSeurat(
  data_dir = data_dir,
  list_data = list_data,
  min.cells = 3,      # Arguments passés à CreateSeuratObject
  min.features = 200,
  save_file = "NeuroProx_filtered.qs"
)

# Avec un autre sous-répertoire de matrice
NeuroProx_raw <- CreateMultipleSeurat(
  data_dir = data_dir,
  list_data = list_data,
  matrix_subdir = "raw_feature_bc_matrix",  # Au lieu de filtered
  save_file = "NeuroProx_raw.qs"
)

# ==========================================
# Avantages de la nouvelle fonction
# ==========================================

# 1. Plus simple : 1 ligne au lieu de 6+
# 2. Gestion d'erreurs intégrée
# 3. Validation des chemins automatique
# 4. Compatible avec votre workflow existant
# 5. Même résultat exact que le code manuel

cat("Conversion terminée ! NeuroProx créé avec la nouvelle fonction.\n")
cat("Nombre d'échantillons :", length(NeuroProx), "\n")
cat("Noms des échantillons :", paste(names(NeuroProx), collapse = ", "), "\n")

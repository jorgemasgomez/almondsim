library(dplyr)
library(igraph)
library(RColorBrewer)

df <- read.table(
  "parent_gs_tracking.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

df_scheme <- df %>%
  filter(scheme == unique(scheme)[1], rep == 1) %>%
  distinct(id, mother, father, year, rep, scheme, .keep_all = TRUE) %>%
  arrange(year)

years <- sort(unique(df_scheme$year))

folder_name <- unique(df_scheme$scheme)
if(!dir.exists(folder_name)) dir.create(folder_name)

for (y in years) {
  
  df_y <- df_scheme %>% filter(year <= y)
  
  edges <- bind_rows(
    df_y %>% select(from = father, to = id),
    df_y %>% select(from = mother, to = id)
  ) %>%
    filter(!is.na(from), from != "", from != "0") %>%
    distinct(from, to)  # <- elimina duplicados, una línea por conexión
  
  nodes <- unique(c(edges[["from"]], edges[["to"]]))
  nodes <- nodes[!is.na(nodes) & nodes != "" & nodes != "0"]
  
  g <- graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = data.frame(name = nodes)
  )
  
  # Los 30 padres de este año
  parents_y <- df_scheme %>% filter(year == y) %>% slice(1:30)
  
  # Crear grupos de full-sibs y half-sibs
  sib_groups <- list()
  group_id <- 1
  used <- rep(FALSE, nrow(parents_y))
  
  for (i in seq_len(nrow(parents_y))) {
    if(used[i]) next
    
    full_sibs <- which(
      parents_y$father == parents_y$father[i] &
        parents_y$mother == parents_y$mother[i]
    )
    half_sibs <- which(
      (parents_y$father == parents_y$father[i] |
         parents_y$mother == parents_y$mother[i]) &
        !(parents_y$father == parents_y$father[i] &
            parents_y$mother == parents_y$mother[i])
    )
    
    sib_groups[[group_id]] <- list(
      full = full_sibs,
      half = half_sibs
    )
    used[c(full_sibs, half_sibs)] <- TRUE
    group_id <- group_id + 1
  }
  
  # Generar colores base
  n_groups <- length(sib_groups)
  base_colors <- brewer.pal(min(12, n_groups), "Set3")
  if(n_groups > 12) base_colors <- rep(base_colors, length.out = n_groups)
  
  vertex_colors <- rep("grey80", length(nodes))
  names(vertex_colors) <- nodes
  
  for (i in seq_along(sib_groups)) {
    full_ids <- parents_y$id[sib_groups[[i]]$full]
    half_ids <- parents_y$id[sib_groups[[i]]$half]
    
    vertex_colors[nodes %in% full_ids] <- base_colors[i]          # full sibs
    vertex_colors[nodes %in% half_ids] <- adjustcolor(base_colors[i], alpha.f = 0.5) # half sibs más claros
  }
  
  V(g)$color <- vertex_colors
  V(g)$size <- ifelse(V(g)$name %in% parents_y$id, 8, 5)
  
  roots <- V(g)$name[degree(g, mode = "in") == 0]
  
  png_filename <- paste0(folder_name, "/Pedigree_", y, ".png")
  png(png_filename, width = 1200, height = 800)
  
  plot(
    g,
    layout = layout_as_tree(g, root = roots),
    vertex.label.cex = 0.6,
    main = paste("Pedigree dinámico – Año", y)
  )
  
  dev.off()
}

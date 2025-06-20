# app.R

library(shiny)
library(dplyr)
library(tidyr)
library(pheatmap)
library(vegan)
library(ggplot2)
library(plotly)
library(tools)

# -------------------------------------------------
# 1) Fonctions pour lire un CSV / CSV2 sans se tromper
# -------------------------------------------------
guess_sep <- function(path) {
  ligne1 <- readLines(path, n = 1, encoding = "UTF-8")
  if (grepl(";", ligne1, fixed = TRUE)) {
    return(";")
  } else {
    return(",")
  }
}

read_flexible_csv <- function(path) {
  sep_utilise <- guess_sep(path)
  if (sep_utilise == ";") {
    read.csv(
      path,
      sep = ";",
      dec = ",",
      stringsAsFactors = FALSE,
      encoding = "UTF-8",
      check.names = FALSE
    )
  } else {
    read.csv(
      path,
      sep = ",",
      dec = ".",
      stringsAsFactors = FALSE,
      encoding = "UTF-8",
      check.names = FALSE
    )
  }
}

# -------------------------------------------------
# 2) Fonction de lecture + filtrage d’un rapport CZID
# -------------------------------------------------
#   - lvl : 1 pour espèce, 2 pour genre
read_czid_report <- function(
    path,
    metric = "nt_bpm",
    lvl = 1,
    category_list = c("viruses", "bacteria"),
    min_nt_bpm = 1,
    min_nr_bpm = 1,
    min_nt_contigs = 0,
    min_nr_contigs = 0
) {
  df <- read_flexible_csv(path)
  df[is.na(df)] <- 0
  
  # Convertir en numérique les colonnes de seuils
  cols_num <- c("nt_bpm", "nr_bpm", "nt_contig_b", "nr_contig_b")
  for (col in cols_num) {
    if (col %in% colnames(df)) {
      df[[col]] <- as.numeric(gsub(",", ".", df[[col]]))
    }
  }
  
  # 1) Filtre sur niveau taxonomique
  if ("tax_level" %in% colnames(df)) {
    df <- df %>% filter(.data$tax_level == lvl)
  }
  
  # 2) Filtre catégories
  if ("category" %in% colnames(df)) {
    df <- df %>%
      mutate(category = trimws(tolower(category))) %>%
      filter(category %in% tolower(category_list))
  }
  
  # 3) Seuils
  if ("nt_bpm" %in% colnames(df)) {
    df <- df %>% filter(nt_bpm >= min_nt_bpm)
  }
  if ("nr_bpm" %in% colnames(df)) {
    df <- df %>% filter(nr_bpm >= min_nr_bpm)
  }
  if ("nt_contig_b" %in% colnames(df)) {
    df <- df %>% filter(nt_contig_b >= min_nt_contigs)
  }
  if ("nr_contig_b" %in% colnames(df)) {
    df <- df %>% filter(nr_contig_b >= min_nr_contigs)
  }
  
  return(df)
}

# -------------------------------------------------
# 3) Fonction pour tracer la heatmap
# -------------------------------------------------
plot_czid_heatmap <- function(
    df_joined,
    metric = "nt_bpm",
    top_n = 10,
    log_transform = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE
) {
  df_agg <- df_joined %>%
    group_by(name, display_sample) %>%
    summarise(valeur = sum(.data[[metric]], na.rm = TRUE), .groups = "drop")
  
  if (nrow(df_agg) == 0) {
    plot.new()
    text(0.5, 0.5, "Aucune donnée à afficher\navec ces filtres.", col = "red", cex = 1.2)
    return()
  }
  
  df_wide <- df_agg %>%
    pivot_wider(names_from = display_sample, values_from = valeur, values_fill = 0)
  
  mat_all <- as.matrix(df_wide[ , -1, drop = FALSE ])
  rownames(mat_all) <- df_wide$name
  
  somme_taxa <- rowSums(mat_all, na.rm = TRUE)
  n_effectif <- min(top_n, length(somme_taxa))
  top_taxa <- names(sort(somme_taxa, decreasing = TRUE))[1:n_effectif]
  mat_top <- mat_all[top_taxa, , drop = FALSE]
  
  if (log_transform) {
    mat_top <- log2(mat_top + 1)
  }
  
  mat_top <- mat_top[order(rownames(mat_top)), , drop = FALSE]
  
  colors <- colorRampPalette(c("white", "orange", "red"))(50)
  
  pheatmap(
    mat_top,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    fontsize_row = ifelse(show_rownames, 8, 0),
    fontsize_col = ifelse(show_colnames, 8, 0),
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    color = colors
  )
}

# -------------------------------------------------
# 4) UI
# -------------------------------------------------
ui <- fluidPage(
  titlePanel("Heatmap CZID + Analyses Diversité / PCoA / Stacked Bar"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "metadonnees",
        "Charger métadonnées (CSV ou CSV2)\n(Colonnes attendues : 'sample' + 'sampleid')",
        accept = c(".csv", ".txt")
      ),
      fileInput(
        "reports",
        "Charger rapports CZID (CSV ou CSV2)",
        accept = c(".csv", ".txt"),
        multiple = TRUE
      ),
      
      tags$hr(),
      
      # Paramètres Heatmap
      conditionalPanel(
        condition = "input.mainTabs == 'Heatmap'",
        radioButtons(
          "tax_choice",
          "Niveau taxonomique",
          choices = list("Espèce (tax_level = 1)" = 1, "Genre (tax_level = 2)" = 2),
          selected = 1
        ),
        checkboxGroupInput(
          "category_list",
          "Catégories à inclure",
          choices = list("viruses" = "viruses", "bacteria" = "bacteria",
                         "archaea" = "archaea", "eukaryota" = "eukaryota"),
          selected = c("viruses", "bacteria")
        ),
        numericInput("min_nt_bpm", "Seuil min nt_bpm", value = 1, min = 0),
        numericInput("min_nr_bpm", "Seuil min nr_bpm", value = 0, min = 0),
        numericInput("min_nt_contigs", "Seuil min nt_contig_b", value = 0, min = 0),
        numericInput("min_nr_contigs", "Seuil min nr_contig_b", value = 0, min = 0),
        selectInput(
          "metric",
          "Metric pour heatmap",
          choices = c("nt_bpm", "nr_bpm", "nt_contig_b", "nr_contig_b"),
          selected = "nt_bpm"
        ),
        numericInput("top_n", "Top N taxa", value = 10, min = 1),
        checkboxInput("log_transform", "Appliquer log2(x+1)", value = FALSE),
        checkboxInput("show_rownames", "Afficher noms des lignes", value = TRUE),
        checkboxInput("show_colnames", "Afficher noms des colonnes", value = TRUE),
        checkboxInput("cluster_rows", "Cluster lignes", value = FALSE),
        checkboxInput("cluster_cols", "Cluster colonnes", value = FALSE),
        downloadButton("downloadHeatmap", "Télécharger Heatmap (PDF)")
      ),
      
      # Paramètres Alpha Diversity
      conditionalPanel(
        condition = "input.mainTabs == 'Alpha Diversity'",
        selectInput(
          "alpha_index",
          "Index de diversité",
          choices = c("Shannon" = "shannon", "Simpson" = "simpson", "Richness" = "richness"),
          selected = "shannon"
        ),
        uiOutput("groupSelect"),
        downloadButton("downloadAlpha", "Télécharger Alpha (PNG)")
      ),
      
      # Paramètres PCoA
      conditionalPanel(
        condition = "input.mainTabs == 'PCoA'",
        uiOutput("pcoaGroupSelect"),
        checkboxInput("show_centroid", "Afficher centroïdes", value = FALSE),
        checkboxInput("show_segments", "Afficher segments vers centroïdes", value = FALSE),
        checkboxInput("show_ellipse", "Afficher ellipses (95%)", value = FALSE),
        downloadButton("downloadPCoA", "Télécharger PCoA (PNG)")
      ),
      
      # Paramètres Stacked Bar
      conditionalPanel(
        condition = "input.mainTabs == 'Stacked Bar'",
        sliderInput("num_stacks", "Nombre de taxa", min = 2, max = 40, value = 10),
        uiOutput("stackGroupSelect"),
        checkboxGroupInput(
          "sb_category_list",
          "Filtrer par catégorie",
          choices = list("viruses" = "viruses", "bacteria" = "bacteria",
                         "archaea" = "archaea", "eukaryota" = "eukaryota"),
          selected = c("viruses", "bacteria", "archaea", "eukaryota")
        ),
        downloadButton("downloadStacked", "Télécharger Stacked Bar (PNG)")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "mainTabs",
        
        tabPanel("Aperçu métadonnées", tableOutput("previewMeta")),
        tabPanel("Aperçu rapports", uiOutput("previewReports")),
        
        tabPanel("Heatmap", plotOutput("heatmapPlot", height = "700px")),
        
        tabPanel(
          "Alpha Diversity",
          plotOutput("alphaPlot", height = "500px"),
          tags$br(),
          h4("Table Alpha Diversity"),
          tableOutput("alphaTable")
        ),
        
        tabPanel("PCoA", plotOutput("pcoaPlot", height = "500px")),
        
        tabPanel("Stacked Bar", plotlyOutput("stacked_bar", height = "600px"))
      )
    )
  )
)

# -------------------------------------------------
# 5) Server
# -------------------------------------------------
server <- function(input, output, session) {
  
  # 5.1) Métadonnées réactives
  meta_df <- reactive({
    req(input$metadonnees)
    df <- tryCatch({
      read_flexible_csv(input$metadonnees$datapath)
    }, error = function(e) {
      showNotification("Erreur lecture métadonnées.", type = "error")
      return(NULL)
    })
    if (!all(c("sample", "sampleid") %in% colnames(df))) {
      showNotification("Colonnes 'sample' et 'sampleid' requises.", type = "error")
      return(NULL)
    }
    return(df)
  })
  
  # 5.2) Aperçu métadonnées
  output$previewMeta <- renderTable({
    head(meta_df(), 10)
  })
  
  # 5.3) Lecture rapports
  reports_list <- reactive({
    req(input$reports)
    paths <- input$reports$datapath
    names(paths) <- input$reports$name
    return(paths)
  })
  
  # 5.4) Aperçu rapports
  output$previewReports <- renderUI({
    req(reports_list())
    lapply(names(reports_list()), function(nf) {
      chemin <- reports_list()[[nf]]
      df <- tryCatch({ read_flexible_csv(chemin) }, error = function(e) NULL)
      if (is.null(df)) return(NULL)
      tagList(tags$h4(nf), renderTable({ head(df, 5) }))
    }) %>% do.call(tagList, .)
  })
  
  # 5.5) base_data rassemble df_plot et df_full (avec category)
  base_data <- reactive({
    req(meta_df(), reports_list())
    meta <- meta_df()
    paths <- reports_list()
    
    dfs_filtres <- lapply(seq_along(paths), function(i) {
      chemin <- paths[i]
      nf <- names(paths)[i]
      df_raw <- read_czid_report(
        path = chemin,
        metric = input$metric,
        lvl = as.numeric(input$tax_choice),
        category_list = input$category_list,
        min_nt_bpm = input$min_nt_bpm,
        min_nr_bpm = input$min_nr_bpm,
        min_nt_contigs = input$min_nt_contigs,
        min_nr_contigs = input$min_nr_contigs
      )
      base <- file_path_sans_ext(nf)
      base_clean <- sub("_taxon_report$", "", base)
      df_raw$samplename <- base_clean
      return(df_raw)
    })
    
    big_df <- bind_rows(dfs_filtres)
    
    df_joined <- big_df %>%
      left_join(meta %>% select(sample, sampleid), by = c("samplename" = "sample")) %>%
      mutate(display_sample = ifelse(is.na(sampleid), samplename, sampleid))
    
    df_plot <- df_joined %>% select(name, display_sample, all_of(input$metric))
    return(list(df_plot = df_plot, df_full = df_joined, meta = meta))
  })
  
  # 5.6) matrice taxa × échantillons
  abundance_matrix <- reactive({
    df_plot <- base_data()$df_plot
    req(df_plot)
    df_agg <- df_plot %>%
      group_by(name, display_sample) %>%
      summarise(valeur = sum(.data[[input$metric]], na.rm = TRUE), .groups = "drop")
    if (nrow(df_agg) == 0) return(NULL)
    df_wide <- df_agg %>%
      pivot_wider(names_from = display_sample, values_from = valeur, values_fill = 0)
    mat_all <- as.matrix(df_wide[ , -1, drop = FALSE ])
    rownames(mat_all) <- df_wide$name
    return(mat_all)
  })
  
  # 5.7) Heatmap plot réactif
  output$heatmapPlot <- renderPlot({
    mat_all <- abundance_matrix()
    df_plot <- base_data()$df_plot
    req(mat_all, df_plot)
    plot_czid_heatmap(
      df_joined = base_data()$df_full,
      metric = input$metric,
      top_n = input$top_n,
      log_transform = input$log_transform,
      show_rownames = input$show_rownames,
      show_colnames = input$show_colnames,
      cluster_rows = input$cluster_rows,
      cluster_cols = input$cluster_cols
    )
  })
  
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste0("heatmap_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      mat_all <- abundance_matrix()
      df_plot <- base_data()$df_plot
      req(mat_all, df_plot)
      pdf(file, width = 8, height = 8)
      plot_czid_heatmap(
        df_joined = base_data()$df_full,
        metric = input$metric,
        top_n = input$top_n,
        log_transform = input$log_transform,
        show_rownames = input$show_rownames,
        show_colnames = input$show_colnames,
        cluster_rows = input$cluster_rows,
        cluster_cols = input$cluster_cols
      )
      dev.off()
    }
  )
  
  # 5.8) UI dynamique : variable de regroupement pour Alpha
  output$groupSelect <- renderUI({
    choix <- setdiff(colnames(base_data()$meta), c("sample", "sampleid"))
    selectInput("alpha_group", "Grouper par (metadata)",
                choices = c("Aucun" = "None", choix), selected = "None")
  })
  
  # 5.9a) Calcul Alpha Diversity
  alpha_df <- reactive({
    mat_all <- abundance_matrix()
    req(mat_all)
    index <- input$alpha_index
    if (index == "richness") {
      div_vals <- apply(mat_all, 2, function(x) vegan::specnumber(x))
    } else {
      div_vals <- apply(mat_all, 2, function(x) vegan::diversity(x, index = index))
    }
    df_div <- data.frame(display_sample = colnames(mat_all), diversity = div_vals,
                         stringsAsFactors = FALSE)
    meta2 <- base_data()$meta %>%
      mutate(display_sample = ifelse(is.na(sampleid), sample, sampleid))
    df_div <- df_div %>% left_join(meta2, by = "display_sample")
    return(df_div)
  })
  
  # 5.9b) Affichage plot Alpha
  output$alphaPlot <- renderPlot({
    df_div <- alpha_df()
    req(df_div)
    index <- input$alpha_index
    grp_var <- input$alpha_group
    if (grp_var == "None") {
      ggplot(df_div, aes(x = display_sample, y = diversity)) +
        geom_jitter(width = 0.2, size = 2, shape = 21, fill = "white", color = "black") +
        labs(x = "Échantillon", y = paste("Indice", index),
             title = paste0("Alpha Diversity (", index, ")")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              plot.title = element_text(hjust = 0.5))
    } else {
      df_div[[grp_var]] <- as.factor(df_div[[grp_var]])
      ggplot(df_div, aes_string(x = grp_var, y = "diversity", fill = grp_var)) +
        geom_boxplot(color = "black") +
        geom_jitter(width = 0.2, size = 2, shape = 21, fill = "white", color = "black") +
        labs(x = grp_var, y = paste("Indice", index),
             title = paste0("Alpha Diversity (", index, ") par ", grp_var),
             fill = grp_var) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 12),
              plot.title = element_text(hjust = 0.5))
    }
  })
  
  output$alphaTable <- renderTable({
    df_div <- alpha_df()
    req(df_div)
    df_div %>% select(display_sample, diversity) %>% arrange(display_sample)
  }, rownames = TRUE)
  
  output$downloadAlpha <- downloadHandler(
    filename = function() {
      paste0("alpha_diversity_", Sys.Date(), ".png")
    },
    content = function(file) {
      df_div <- alpha_df()
      req(df_div)
      png(file, width = 1000, height = 700, res = 150)
      index <- input$alpha_index
      grp_var <- input$alpha_group
      if (grp_var == "None") {
        print(
          ggplot(df_div, aes(x = display_sample, y = diversity)) +
            geom_jitter(width = 0.2, size = 2, shape = 21, fill = "white", color = "black") +
            labs(x = "Échantillon", y = paste("Indice", index),
                 title = paste0("Alpha Diversity (", index, ")")) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  plot.title = element_text(hjust = 0.5))
        )
      } else {
        df_div[[grp_var]] <- as.factor(df_div[[grp_var]])
        print(
          ggplot(df_div, aes_string(x = grp_var, y = "diversity", fill = grp_var)) +
            geom_boxplot(color = "black") +
            geom_jitter(width = 0.2, size = 2, shape = 21, fill = "white", color = "black") +
            labs(x = grp_var, y = paste("Indice", index),
                 title = paste0("Alpha Diversity (", index, ") par ", grp_var),
                 fill = grp_var) +
            theme_bw() +
            theme(axis.text.x = element_text(size = 12),
                  plot.title = element_text(hjust = 0.5))
        )
      }
      dev.off()
    }
  )
  
  # 5.10) UI dynamique : choix du groupe pour PCoA
  output$pcoaGroupSelect <- renderUI({
    choix <- setdiff(colnames(base_data()$meta), c("sample", "sampleid"))
    selectInput("pcoa_group", "Colorer/par (metadata)",
                choices = c("Aucun" = "None", choix), selected = choix[1])
  })
  
  # 5.11) Affichage PCoA
  output$pcoaPlot <- renderPlot({
    mat_all <- abundance_matrix()
    req(mat_all)
    dist_bc <- vegdist(t(mat_all), method = "bray")
    pcoa_res <- cmdscale(dist_bc, eig = TRUE, k = 2)
    coords <- as.data.frame(pcoa_res$points)
    colnames(coords) <- c("Axis1", "Axis2")
    coords$display_sample <- rownames(coords)
    meta2 <- base_data()$meta %>%
      mutate(display_sample = ifelse(is.na(sampleid), sample, sampleid))
    df_plot <- coords %>% left_join(meta2, by = "display_sample")
    
    grp_var <- input$pcoa_group
    show_centroid <- input$show_centroid
    show_segments <- input$show_segments
    show_ellipse <- input$show_ellipse
    
    if (grp_var == "None") {
      ggplot(df_plot, aes(x = Axis1, y = Axis2)) +
        geom_point(size = 3) +
        labs(x = paste0("PCoA1 (", round(100 * pcoa_res$eig[1] / sum(pcoa_res$eig), 1), "%)"),
             y = paste0("PCoA2 (", round(100 * pcoa_res$eig[2] / sum(pcoa_res$eig), 1), "%)"),
             title = "PCoA (Bray-Curtis)") +
        theme_minimal()
    } else {
      df_plot[[grp_var]] <- as.factor(df_plot[[grp_var]])
      centroids <- df_plot %>%
        group_by(.data[[grp_var]]) %>%
        summarise(cx = mean(Axis1, na.rm = TRUE), cy = mean(Axis2, na.rm = TRUE), .groups = "drop")
      colnames(centroids)[1] <- grp_var
      
      p <- ggplot(df_plot, aes_string(x = "Axis1", y = "Axis2", color = grp_var)) +
        geom_point(size = 3)
      if (show_segments) {
        p <- p + geom_segment(
          data = df_plot %>% left_join(centroids, by = grp_var),
          aes(xend = cx, yend = cy), alpha = 0.5
        )
      }
      if (show_centroid) {
        p <- p + geom_point(
          data = centroids,
          aes_string(x = "cx", y = "cy"),
          shape = 21, fill = "white", color = "black", size = 4
        )
      }
      if (show_ellipse) {
        p <- p + stat_ellipse(aes_string(color = grp_var), type = "norm", level = 0.95)
      }
      p <- p +
        labs(x = paste0("PCoA1 (", round(100 * pcoa_res$eig[1] / sum(pcoa_res$eig), 1), "%)"),
             y = paste0("PCoA2 (", round(100 * pcoa_res$eig[2] / sum(pcoa_res$eig), 1), "%)"),
             color = grp_var, title = "PCoA (Bray-Curtis)") +
        theme_minimal()
      print(p)
    }
  })
  
  output$downloadPCoA <- downloadHandler(
    filename = function() {
      paste0("pcoa_", Sys.Date(), ".png")
    },
    content = function(file) {
      mat_all <- abundance_matrix()
      req(mat_all)
      dist_bc <- vegdist(t(mat_all), method = "bray")
      pcoa_res <- cmdscale(dist_bc, eig = TRUE, k = 2)
      coords <- as.data.frame(pcoa_res$points)
      colnames(coords) <- c("Axis1", "Axis2")
      coords$display_sample <- rownames(coords)
      meta2 <- base_data()$meta %>%
        mutate(display_sample = ifelse(is.na(sampleid), sample, sampleid))
      df_plot <- coords %>% left_join(meta2, by = "display_sample")
      
      grp_var <- input$pcoa_group
      show_centroid <- input$show_centroid
      show_segments <- input$show_segments
      show_ellipse <- input$show_ellipse
      
      png(file, width = 1000, height = 700, res = 150)
      if (grp_var == "None") {
        print(
          ggplot(df_plot, aes(x = Axis1, y = Axis2)) +
            geom_point(size = 3) +
            labs(x = paste0("PCoA1 (", round(100 * pcoa_res$eig[1] / sum(pcoa_res$eig), 1), "%)"),
                 y = paste0("PCoA2 (", round(100 * pcoa_res$eig[2] / sum(pcoa_res$eig), 1), "%)"),
                 title = "PCoA (Bray-Curtis)") +
            theme_minimal()
        )
      } else {
        df_plot[[grp_var]] <- as.factor(df_plot[[grp_var]])
        centroids <- df_plot %>%
          group_by(.data[[grp_var]]) %>%
          summarise(cx = mean(Axis1, na.rm = TRUE), cy = mean(Axis2, na.rm = TRUE), .groups = "drop")
        colnames(centroids)[1] <- grp_var
        
        p <- ggplot(df_plot, aes_string(x = "Axis1", y = "Axis2", color = grp_var)) +
          geom_point(size = 3)
        if (show_segments) {
          p <- p + geom_segment(
            data = df_plot %>% left_join(centroids, by = grp_var),
            aes(xend = cx, yend = cy), alpha = 0.5
          )
        }
        if (show_centroid) {
          p <- p + geom_point(
            data = centroids,
            aes_string(x = "cx", y = "cy"),
            shape = 21, fill = "white", color = "black", size = 4
          )
        }
        if (show_ellipse) {
          p <- p + stat_ellipse(aes_string(color = grp_var), type = "norm", level = 0.95)
        }
        p <- p +
          labs(x = paste0("PCoA1 (", round(100 * pcoa_res$eig[1] / sum(pcoa_res$eig), 1), "%)"),
               y = paste0("PCoA2 (", round(100 * pcoa_res$eig[2] / sum(pcoa_res$eig), 1), "%)"),
               color = grp_var, title = "PCoA (Bray-Curtis)") +
          theme_minimal()
        print(p)
      }
      dev.off()
    }
  )
  
  
  # 5.13) UI dynamique : choix du groupe pour Stacked Bar
  output$stackGroupSelect <- renderUI({
    choix <- setdiff(colnames(base_data()$meta), c("sample", "sampleid"))
    selectInput("stack_group", "Facer par (metadata)",
                choices = c("Aucun" = "None", choix), selected = "None")
  })
  
  # 5.14) Stacked Bar avec filtre catégories
  output$stacked_bar <- renderPlotly({
    df_full <- base_data()$df_full
    df_md <- base_data()$meta
    grp <- input$stack_group
    n_tax <- input$num_stacks
    cats <- input$sb_category_list
    
    # Filtrer catégories
    df_filt <- df_full %>% filter(category %in% cats)
    
    df_agg0 <- df_filt %>%
      group_by(name, display_sample, category) %>%
      summarise(count = sum(.data[[input$metric]], na.rm = TRUE), .groups = "drop")
    
    abund_mat <- df_agg0 %>%
      group_by(name, display_sample) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      pivot_wider(names_from = display_sample, values_from = count, values_fill = 0)
    
    df_long0 <- abund_mat %>%
      rename(Taxon = name) %>%
      pivot_longer(cols = -Taxon, names_to = "SampleID", values_to = "Count") %>%
      left_join(df_md %>% rename(SampleID = sampleid), by = "SampleID")
    
    total_per_sample <- df_long0 %>%
      group_by(SampleID) %>%
      summarise(TotalSample = sum(Count, na.rm = TRUE), .groups = "drop")
    
    long1 <- df_long0 %>%
      left_join(total_per_sample, by = "SampleID") %>%
      mutate(RelAbundance = ifelse(TotalSample > 0, Count / TotalSample, 0))
    
    total_ab_catalog <- long1 %>%
      group_by(Taxon) %>%
      summarise(TotalAcross = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(TotalAcross))
    
    top_taxa <- head(total_ab_catalog$Taxon, n_tax)
    
    if (grp == "None") {
      long2 <- long1 %>%
        mutate(Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
        group_by(SampleID, Taxon2) %>%
        summarise(RelAbund2 = sum(RelAbundance, na.rm = TRUE), .groups = "drop")
    } else {
      long2 <- long1 %>%
        mutate(Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
        group_by(SampleID, Taxon2, .data[[grp]]) %>%
        summarise(RelAbund2 = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
        mutate(!!grp := as.factor(.data[[grp]]))
    }
    
    taxon_levels <- c(top_taxa, "Other")
    long2$Taxon2 <- factor(long2$Taxon2, levels = taxon_levels)
    
    pal <- viridis::viridis(length(taxon_levels))
    colors_named <- setNames(pal, taxon_levels)
    
    if (grp == "None") {
      p <- ggplot(long2, aes(x = SampleID, y = RelAbund2, fill = Taxon2)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_fill_manual(values = colors_named) +
        labs(
          x = "Échantillon", y = "Abondance relative",
          title = "Stacked Bar (pas de faceting)", fill = "Taxon"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          panel.grid.major.x = element_blank()
        )
    } else {
      p <- ggplot(long2, aes(x = SampleID, y = RelAbund2, fill = Taxon2)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_fill_manual(values = colors_named) +
        facet_wrap(as.formula(paste("~", grp)), scales = "free_x", nrow = 1) +
        labs(
          x = grp, y = "Abondance relative",
          title = paste0("Stacked Bar (échantillons par ", grp, ")"),
          fill = "Taxon"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          strip.text.x = element_text(size = 12),
          panel.grid.major.x = element_blank()
        )
    }
    
    ggplotly(p, tooltip = c("y", "fill"))
  })
  
  output$downloadStacked <- downloadHandler(
    filename = function() {
      paste0("stacked_bar_", Sys.Date(), ".png")
    },
    content = function(file) {
      df_full <- base_data()$df_full
      df_md <- base_data()$meta
      grp <- input$stack_group
      n_tax <- input$num_stacks
      cats <- input$sb_category_list
      
      df_filt <- df_full %>% filter(category %in% cats)
      
      df_agg0 <- df_filt %>%
        group_by(name, display_sample, category) %>%
        summarise(count = sum(.data[[input$metric]], na.rm = TRUE), .groups = "drop")
      
      abund_mat <- df_agg0 %>%
        group_by(name, display_sample) %>%
        summarise(count = sum(count), .groups = "drop") %>%
        pivot_wider(names_from = display_sample, values_from = count, values_fill = 0)
      
      df_long0 <- abund_mat %>%
        rename(Taxon = name) %>%
        pivot_longer(cols = -Taxon, names_to = "SampleID", values_to = "Count") %>%
        left_join(df_md %>% rename(SampleID = sampleid), by = "SampleID")
      
      total_per_sample <- df_long0 %>%
        group_by(SampleID) %>%
        summarise(TotalSample = sum(Count, na.rm = TRUE), .groups = "drop")
      
      long1 <- df_long0 %>%
        left_join(total_per_sample, by = "SampleID") %>%
        mutate(RelAbundance = ifelse(TotalSample > 0, Count / TotalSample, 0))
      
      total_ab_catalog <- long1 %>%
        group_by(Taxon) %>%
        summarise(TotalAcross = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(TotalAcross))
      
      top_taxa <- head(total_ab_catalog$Taxon, n_tax)
      
      if (grp == "None") {
        long2 <- long1 %>%
          mutate(Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
          group_by(SampleID, Taxon2) %>%
          summarise(RelAbund2 = sum(RelAbundance, na.rm = TRUE), .groups = "drop")
      } else {
        long2 <- long1 %>%
          mutate(Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
          group_by(SampleID, Taxon2, .data[[grp]]) %>%
          summarise(RelAbund2 = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
          mutate(!!grp := as.factor(.data[[grp]]))
      }
      
      taxon_levels <- c(top_taxa, "Other")
      long2$Taxon2 <- factor(long2$Taxon2, levels = taxon_levels)
      
      pal <- viridis::viridis(length(taxon_levels))
      colors_named <- setNames(pal, taxon_levels)
      
      png(file, width = 1000, height = 700, res = 150)
      if (grp == "None") {
        p <- ggplot(long2, aes(x = SampleID, y = RelAbund2, fill = Taxon2)) +
          geom_bar(stat = "identity", width = 0.8) +
          scale_fill_manual(values = colors_named) +
          labs(
            x = "Échantillon", y = "Abondance relative",
            title = "Stacked Bar (pas de faceting)", fill = "Taxon"
          ) +
          theme_minimal(base_size = 14) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            panel.grid.major.x = element_blank()
          )
      } else {
        p <- ggplot(long2, aes(x = SampleID, y = RelAbund2, fill = Taxon2)) +
          geom_bar(stat = "identity", width = 0.8) +
          scale_fill_manual(values = colors_named) +
          facet_wrap(as.formula(paste("~", grp)), scales = "free_x", nrow = 1) +
          labs(
            x = grp, y = "Abondance relative",
            title = paste0("Stacked Bar (échantillons par ", grp, ")"),
            fill = "Taxon"
          ) +
          theme_minimal(base_size = 14) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            strip.text.x = element_text(size = 12),
            panel.grid.major.x = element_blank()
          )
      }
      print(p)
      dev.off()
    }
  )
  
}

shinyApp(ui = ui, server = server)

# app.R (完整替换版本)
#Sys.setenv(http_proxy="http://172.16.196.1:7890")
#Sys.setenv(https_proxy="http://172.16.196.1:7890")
#Sys.setenv(all_proxy="socks5://172.16.196.1:7890")

library(DNAcopy)
library(GenomicAlignments)
library(patchwork)
library(shiny)
library(DT)
library(jsonlite)
library(shinyWidgets)
library(CNVision)
library(ggplot2)
library(ggpubr)
library(gghalves)
library(ggsci)
library(tidyverse)
library(scales)
library(future)
library(promises)
library(uuid)   # 用于生成 UUID

plan(multisession)  # 并发后台任务，按服务器资源调整 workers 参数

options(shiny.maxRequestSize = 1024 * 1024 * 1024) # 1GB 上限示例，按需调整

# 工作目录与临时目录
setwd('~/msCNVS_web_app')
app_dir <- getwd()
TEMP_DIR <- file.path(app_dir, "TMP")
if (!dir.exists(TEMP_DIR)) dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)

# 允许通过 URL 访问 TMP 下的静态文件（Task Viewer 使用）
addResourcePath("TMP", TEMP_DIR)

log2_range <- c(-0.45, 0.3)

# ---------------- UI ----------------
ui <- fluidPage(
  navbarPage(
    title = div(
      tags$img(src="sequmed.png", height="30px"),
      span("Sequmed CNV Web App")
    ),
    
    # WIKI tab (保留你的样式，但去掉固定高度)
    tabPanel(
      "WIKI",
      tags$head(
        tags$style(HTML("
          /* Sidebar styling */
          #wiki_sidebar {
            background: #f5f7fa;
            border-right: 1px solid #d9d9d9;
            padding-top: 25px;
            padding-left: 18px;
            padding-right: 18px;
            min-height: 880px;
          }
          #wiki_sidebar h4 { font-size: 20px; font-weight: 600; color: #1f4e79; margin-bottom: 18px; }
          #wiki_sidebar ul { list-style-type: none; padding-left: 0; margin-top: 0; }
          #wiki_sidebar li { margin-bottom: 12px; }
          #wiki_sidebar a { font-size: 15px; color: #1f4e79; text-decoration: none; font-weight: 500; }
          #wiki_sidebar a:hover { color: #163a5c; text-decoration: underline; }
          #wiki_main { padding: 30px 35px; background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-top: 20px; }
        "))
      ),
      sidebarLayout(
        sidebarPanel(
          id = "wiki_sidebar",
          h4("Documentation"),
          tags$ul(
            tags$li(actionLink("wiki_intro", "What is this tool?")),
            tags$li(actionLink("wiki_principles", "Principles")),
            tags$li(actionLink("wiki_analysis", "Analysis")),
            tags$li(actionLink("wiki_cite", "How to Cite"))
          ),
          width = 3
        ),
        mainPanel(
          div(id = "wiki_main", uiOutput("wiki_content")),
          width = 9
        )
      )
    ),
    
    # Analysis tab
    tabPanel(
      title = "Analysis",
      tags$head(
        tags$style(HTML("
          /* Sidebar styling */
          #analysis_sidebar { background: #f5f7fa; border-right: 1px solid #d9d9d9; padding: 25px 20px; }
          #analysis_sidebar h3 { font-size: 20px; font-weight: 600; color: #1f4e79; margin-top: 25px; margin-bottom: 12px; }
          #analysis_sidebar label { font-weight: 500; color: #333; }
          #analysis_sidebar .form-control, #analysis_sidebar .selectize-input { border-radius: 6px; }
          #analysis_sidebar hr { border-top: 1px solid #c9c9c9; margin-top: 20px; margin-bottom: 20px; }
          #analysis_main { padding: 30px 35px; background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-top: 20px; }
          #analysis_main h3 { font-size: 22px; color: #1f4e79; margin-top: 10px; margin-bottom: 15px; }
          .btn-primary { background-color: #1f4e79 !important; border-color: #1f4e79 !important; border-radius: 6px; font-weight: 600; }
          .btn-primary:hover { background-color: #163a5c !important; border-color: #163a5c !important; }
          .download-btn { margin-top: 10px; }
        "))
      ),
      sidebarLayout(
        sidebarPanel(
          id = "analysis_sidebar",
          width = 4,
          h3("Upload & Directory"),
          fileInput("file", "Upload ZIP File", accept = ".zip"),
          hr(),
          h3("LoadBins Parameters"),
          numericInput("resolution", "Resolution", value = 600),
          numericInput("length", "Length", value = 150),
          selectInput("genome", "Genome", choices = c("hg19", "hg38"), selected = "hg19"),
          hr(),
          h3("Maskbins"),
          checkboxGroupInput("ctypes", "Mask Types",
                             choices = c("centromere", "clone", "contig", "heterochromatin", "scaffold", "short_arm"),
                             selected = c("centromere", "clone")),
          hr(),
          h3("Segmentation"),
          numericInput("alpha", "Alpha", value = 0.01),
          numericInput("nperm", "nperm", value = 1000),
          numericInput("SD", "undo.SD", value = 1),
          numericInput("undo_min_width", "min.width", value = 5),
          hr(),
          h3("InferPloidy"),
          numericInput("m_start", "m_start", value = 1.5),
          numericInput("m_end", "m_end", value = 2.5),
          numericInput("m_step", "m_step", value = 0.1),
          numericInput("b_start", "b_start", value = -0.2),
          numericInput("b_end", "b_end", value = 0.2),
          numericInput("b_step", "b_step", value = 0.01),
          hr(),
          h3("Peak Detection"),
          numericInput("minpeakdistance", "minpeakdistance", value = 0.3),
          hr(),
          h3("PlotPeaks Output"),
          numericInput("plot_width", "Width", value = 3.7),
          numericInput("plot_height", "Height", value = 0.55),
          selectInput("plot_unit", "Unit", choices = c("in", "cm", "mm"), selected = "in"),
          numericInput("plot_dpi", "DPI", value = 300),
          hr(),
          h3("PlotCNV Output"),
          numericInput("plot_width_2", "Width", value = 3.7),
          numericInput("plot_height_2", "Height", value = 0.55),
          selectInput("plot_unit_2", "Unit", choices = c("in", "cm", "mm"), selected = "in"),
          numericInput("plot_dpi_2", "DPI", value = 300),
          hr(),
          actionButton("run_analysis", "Run Analysis", class = "btn-primary")
        ),
        mainPanel(
          width = 8,
          div(
            id = "analysis_main",
            h3("Peaks Plot"),
            imageOutput("peaks_plot"),  # 自适应高度
            br(),
            downloadButton("download_peaks", "Download Peaks", class = "download-btn"),
            hr(),
            h3("CNV Plot"),
            imageOutput("cnv_plot"),
            br(),
            downloadButton("download_cnv", "Download CNV", class = "download-btn")
          )
        )
      )
    ),
    
    # Tasks tab
    tabPanel(
      title = "Tasks",
      tags$head(
        tags$style(HTML("
    /* Tasks 面板响应式高度与滚动 */
    .tasks-column {
      min-height: 60vh;           /* 左右两栏最小高度为视窗的 60% */
      max-height: 80vh;           /* 最大高度限制，避免过高 */
      overflow: auto;             /* 超出时显示滚动条 */
      padding-right: 12px;
      padding-left: 12px;
    }

    /* Task Viewer 内图片自适应并限制高度 */
    .task-viewer img {
      max-width: 100%;
      height: auto;
      display: block;
      margin-bottom: 12px;
      max-height: 60vh;           /* 单张图片最大高度，防止撑开页面 */
      object-fit: contain;
    }

    /* 日志区域样式：固定高度并可滚动 */
    .task-log {
      background: #fafafa;
      border: 1px solid #e6e6e6;
      padding: 10px;
      border-radius: 6px;
      max-height: 40vh;
      overflow: auto;
      white-space: pre-wrap;
      font-family: monospace;
      font-size: 12px;
      color: #333;
    }

    /* 在小屏幕上调整高度 */
    @media (max-width: 768px) {
      .tasks-column { min-height: 50vh; max-height: 70vh; }
      .task-log { max-height: 30vh; }
    }
  "))
      ),
        fluidRow(
          column(
            width = 6,
            class = "tasks-column",
            h3("Task List"),
            DT::DTOutput("tasks_table")
          ),
          column(
            width = 6,
            class = "tasks-column",
            h3("Task Viewer"),
            div(class = "task-viewer-container",
                uiOutput("task_viewer")
            ),
            hr(),
            h4("Task Log"),
            div(class = "task-log",
                verbatimTextOutput("task_log")
            )
          )
        )
      )
      ,
  
  # Footer
  tags$footer(
    tags$div(
      style = "text-align:center; padding: 14px 10px; background: #f5f7fa; border-top: 1px solid #d9d9d9; margin-top: 30px;",
      tags$div("Sequmed CNV Web App is a cutting-edge tool for analyzing copy number variations in genomic sequencing data, providing accurate and efficient results for researchers in the field of genomics.", style = "font-size:13px; margin-bottom:6px; color:#333;"),
      tags$div("© 2026 Sequmed Medical Diagnosis Center | Developed by Zexin Zheng, Yusong Shi | Contact: yusong_shi@hotmail.com", style = "font-size:12px; color:#666;")
    )
  )
)
)
# ---------------- Server ----------------
server <- function(input, output, session) {
  
  # WIKI content switching
  current_section <- reactiveVal("wiki_intro")
  observeEvent(input$wiki_intro, { current_section("wiki_intro") })
  observeEvent(input$wiki_principles, { current_section("wiki_principles") })
  observeEvent(input$wiki_analysis, { current_section("wiki_analysis") })
  observeEvent(input$wiki_cite, { current_section("wiki_cite") })
  
  output$wiki_content <- renderUI({
    sec <- current_section()
    if (sec == "wiki_intro") {
      tags$iframe(src = "wiki_intro.html", style = "width:100%; min-height:800px; border:0;")
    } else if (sec == "wiki_principles") {
      tags$iframe(src = "wiki_principles.html", style = "width:100%; min-height:800px; border:0;")
    } else if (sec == "wiki_analysis") {
      tags$iframe(src = "wiki_analysis.html", style = "width:100%; min-height:800px; border:0;")
    } else if (sec == "wiki_cite") {
      tags$iframe(src = "wiki_cite.html", style = "width:100%; min-height:800px; border:0;")
    } else {
      HTML("<p>Click the left menu to view content</p>")
    }
  })
  
  # Reactive values for current session
  rv <- reactiveValues(
    io_dir = NULL,
    zip_file = NULL,
    peaks_file = NULL,
    cnv_file = NULL,
    last_target_dir = NULL,
    task_id = NULL
  )
  
  # 上传时创建目录并复制 ZIP（上传即建 task，使用 uuid）
  observeEvent(input$file, {
    req(input$file)
    
    task_id <- paste0("task_", as.integer(Sys.time()), "_", UUIDgenerate())
    rv$task_id <- task_id
    
    io_dir <- normalizePath(file.path(TEMP_DIR, task_id), mustWork = FALSE)
    dir.create(io_dir, recursive = TRUE, showWarnings = FALSE)
    
    safe_name <- gsub("[^A-Za-z0-9._-]", "_", input$file$name)
    zip_target <- file.path(io_dir, safe_name)
    file.copy(input$file$datapath, zip_target, overwrite = TRUE)
    
    task_info <- list(
      id = task_id,
      file = safe_name,
      time_uploaded = as.character(Sys.time()),
      status = "uploaded",
      message = ""
    )
    write_json(task_info, file.path(io_dir, "task.json"), auto_unbox = TRUE, pretty = TRUE)
    writeLines(character(0), con = file.path(io_dir, "log.txt"))
    
    rv$io_dir <- io_dir
    rv$zip_file <- zip_target
    
    showNotification(paste("Upload complete. Task created:", task_id), type = "message")
  })
  
  # run_segment: 纯函数，写日志并返回结构化结果
  run_segment <- function(zip_path, io_dir,
                          resolution = 600, length = 150, genome = "hg19",
                          ctypes = c("centromere","clone"), nperm = 1000, SD = 1, undo_min_width = 5, alpha = 0.01,
                          m_start = 1.5, m_end = 2.5, m_step = 0.1,
                          b_start = -0.2, b_end = 0.2, b_step = 0.01,
                          minpeakdistance = 0.3,
                          plot_width = 3.7, plot_height = 0.55, plot_unit = "in", plot_dpi = 300,
                          plot_width_2 = 3.7, plot_height_2 = 0.55, plot_unit_2 = "in", plot_dpi_2 = 300) {
    log_file <- file.path(io_dir, "log.txt")
    tryCatch({
      
      cat(paste0("[", Sys.time(), "] Starting analysis\n"), file = log_file, append = TRUE)
      
      file_name <- tools::file_path_sans_ext(basename(zip_path))
      target_dir <- file.path(io_dir, file_name)
      if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
      
      
      cat(paste0("[", Sys.time(), "] Unzipping to\n"), file = log_file, append = TRUE)
      
      unzip(zip_path, exdir = target_dir)
      cat(paste0("[", Sys.time(), "] Running CNVision pipeline\n"), file = log_file, append = TRUE)

      a <- CNVision(dir = target_dir)
      b <- LoadBins(a, resolution = resolution, length = length, genome = genome)
      c <- CountRead(b)
      d <- Maskbins(c, mask_types = ctypes)
      e <- NormalizeData(d, method = "Normalize")
      f <- Segment(e, alpha = alpha, undo.splits = "sdundo", nperm = nperm, undo.SD = SD, min.width = undo_min_width)
      g <- InferPloidy(f, m_start = m_start, m_end = m_end, m_step = m_step,
                       b_start = b_start, b_end = b_end, b_step = b_step)
      cell <- g@config$cells[1]
      

      cat(paste0("[", Sys.time(), "] Detecting peaks\n"), file = log_file, append = TRUE)
      h <- detect_peaks(g, minpeakdistance = minpeakdistance, plot = TRUE)
      set.seed(123)
      m <- laplaceMM(h, core_prob = 0.5, dist_type = "gaussian")
      
      # 保存图像
      p1 <- plotPeaks(m)
      peaks_file <- file.path(target_dir, "Peaks.png")
      ggsave(peaks_file, plot = p1, width = plot_width, height = plot_height, units = plot_unit, dpi = plot_dpi)
      
      p2 <- plotCNV(m, cell = cell, without_x = TRUE)
      cnv_file <- file.path(target_dir, paste0(cell, "_CNV.png"))
      ggsave(cnv_file, plot = p2, width = plot_width_2, height = plot_height_2, units = plot_unit_2, dpi = plot_dpi_2)
      cat(paste0("[", Sys.time(), "] Analysis finished successfully\n"), file = log_file, append = TRUE)
      
      list(status = "finished", target_dir = target_dir, peaks = peaks_file, cnv = cnv_file, error = NULL)
    }, error = function(e) {
      cat(paste0("[", Sys.time(), "] ERROR:\n ", e$message), file = log_file, append = TRUE)

      list(status = "failed", target_dir = NULL, peaks = NULL, cnv = NULL, error = e$message)
    })
  }
  
  # Run analysis in background using future/promises
  observeEvent(input$run_analysis, {
    req(rv$io_dir, rv$zip_file)
    io_dir_local <- rv$io_dir 
    zip_file_local <- rv$zip_file
    task_json <- file.path(io_dir_local, "task.json")
    task_info <- read_json(task_json)
    task_info$status <- "running"
    task_info$time_started <- as.character(Sys.time())
    write_json(task_info, task_json, auto_unbox = TRUE, pretty = TRUE)
    genome_local <- input$genome
    resolution_local <- input$resolution
    length_local <- input$length
    ctypes_local <- input$ctypes
    nperm_local <- input$nperm
    SD_local <- input$SD
    undo_min_width_local <- input$undo_min_width
    alpha_local <- input$alpha
    m_start_local <- input$m_start
    m_end_local <- input$m_end
    m_step_local <- input$m_step
    b_start_local <- input$b_start
    b_end_local <- input$b_end
    b_step_local <- input$b_step
    minpeakdistance_local <- input$minpeakdistance
    plot_width_local <- input$plot_width
    plot_height_local <- input$plot_height
    plot_unit_local <- input$plot_unit
    plot_dpi_local <- input$plot_dpi
    plot_width_2_local <- input$plot_width_2
    plot_height_2_local <- input$plot_height_2
    plot_unit_2_local <- input$plot_unit_2
    plot_dpi_2_local <- input$plot_dpi_2
    future(seed = TRUE, {
      run_segment(
        zip_path = zip_file_local,
        io_dir = io_dir_local,
        resolution = resolution_local,
        length = length_local,
        genome = genome_local,
        ctypes = ctypes_local,
        nperm = nperm_local,
        SD = SD_local,
        undo_min_width = undo_min_width_local,
        alpha = alpha_local,
        m_start = m_start_local,
        m_end = m_end_local,
        m_step = m_step_local,
        b_start = b_start_local,
        b_end = b_end_local,
        b_step = b_step_local,
        minpeakdistance = minpeakdistance_local,
        plot_width = plot_width_local,
        plot_height = plot_height_local,
        plot_unit = plot_unit_local,
        plot_dpi = plot_dpi_local,
        plot_width_2 = plot_width_2_local,
        plot_height_2 = plot_height_2_local,
        plot_unit_2 = plot_unit_2_local,
        plot_dpi_2 = plot_dpi_2_local
      )
    })    %...>% (function(res) {
      task_info <- read_json(task_json)
      task_info$status <- res$status
      task_info$time_finished <- as.character(Sys.time())
      if (!is.null(res$error)) task_info$message <- res$error
      write_json(task_info, task_json, auto_unbox = TRUE, pretty = TRUE)
      
      if (res$status == "finished") {
        rv$peaks_file <- res$peaks
        rv$cnv_file <- res$cnv
        rv$last_target_dir <- res$target_dir
        showNotification("Analysis finished", type = "message")
      } else {
        showNotification(paste("Analysis failed:", res$error), type = "error")
      }
    }) %...!% (function(e) {
      task_info <- read_json(task_json)
      task_info$status <- "failed"
      task_info$message <- e$message
      write_json(task_info, task_json, auto_unbox = TRUE, pretty = TRUE)
      cat(paste0("[", Sys.time(), "] Starting analysis\n"), file = log_file, append = TRUE)
      showNotification("Background job failed", type = "error")
    })
  })
  
  # 列出任务函数（Task Manager 使用）
  list_tasks <- function() {
    dirs <- list.dirs(TEMP_DIR, recursive = FALSE, full.names = TRUE)
    tasks <- lapply(dirs, function(d) {
      jf <- file.path(d, "task.json")
      if (file.exists(jf)) {
        info <- tryCatch(jsonlite::read_json(jf), error = function(e) NULL)
        if (!is.null(info)) {
          data.frame(
            id = ifelse(!is.null(info$id), info$id, basename(d)),
            file = ifelse(!is.null(info$file), info$file, NA),
            status = ifelse(!is.null(info$status), info$status, NA),
            time_uploaded = ifelse(!is.null(info$time_uploaded), info$time_uploaded, NA),
            dir = basename(d),
            stringsAsFactors = FALSE
          )
        } else NULL
      } else NULL
    })
    df <- do.call(rbind, Filter(Negate(is.null), tasks))
    if (is.null(df) || nrow(df) == 0) df <- data.frame(id=character(0), file=character(0), status=character(0), time_uploaded=character(0), dir=character(0), stringsAsFactors = FALSE)
    df
  }
  
  # 渲染任务表（初次）
  output$tasks_table <- DT::renderDT({
    df <- list_tasks()
    DT::datatable(df, selection = 'single', rownames = FALSE, options = list(pageLength = 10))
  })
  
  # 选择任务时显示 viewer 与日志
  observeEvent(input$tasks_table_rows_selected, {
    sel <- input$tasks_table_rows_selected
    if (length(sel) == 0) {
      output$task_viewer <- renderUI({ HTML("<p>No task selected</p>") })
      output$task_log <- renderText({ "" })
      return()
    }
    
    df <- list_tasks()
    row <- df[sel, , drop = FALSE]
    task_dir_name <- row$dir
    task_dir <- file.path(TEMP_DIR, task_dir_name)
    
    peaks_path <- file.path(task_dir, "Peaks.png")
    cnv_files <- list.files(task_dir, pattern = "_CNV\\.png$", full.names = TRUE)
    
    output$task_viewer <- renderUI({
      items <- list()
      if (file.exists(peaks_path)) {
        items <- c(items, tags$div(tags$strong("Peaks Plot:"), tags$br(), tags$img(src = paste0("TMP/", task_dir_name, "/Peaks.png"), style = "max-width:100%; height:auto;")))
      } else {
        items <- c(items, tags$p("Peaks plot not available yet."))
      }
      if (length(cnv_files) > 0) {
        items <- c(items, tags$div(tags$strong("CNV Plot:"), tags$br(), tags$img(src = paste0("TMP/", task_dir_name, "/", basename(cnv_files[1])), style = "max-width:100%; height:auto;")))
      } else {
        items <- c(items, tags$p("CNV plot not available yet."))
      }
      items <- c(items,
                 tags$div(style="margin-top:8px;",
                          if (file.exists(peaks_path)) tags$a(href = paste0("TMP/", task_dir_name, "/Peaks.png"), "Download Peaks", target = "_blank") else NULL,
                          HTML("&nbsp;&nbsp;"),
                          if (length(cnv_files) > 0) tags$a(href = paste0("TMP/", task_dir_name, "/", basename(cnv_files[1])), "Download CNV", target = "_blank") else NULL,
                          HTML("&nbsp;&nbsp;"),
                          if (file.exists(file.path(task_dir, "log.txt"))) tags$a(href = paste0("TMP/", task_dir_name, "/log.txt"), "Download Log", target = "_blank") else NULL
                 )
      )
      do.call(tagList, items)
    })
    
    output$task_log <- renderText({
      logf <- file.path(task_dir, "log.txt")
      if (!file.exists(logf)) return("No log available.")
      txt <- readLines(logf, warn = FALSE)
      n <- length(txt)
      if (n > 500) txt <- txt[(n-499):n]
      paste(txt, collapse = "\n")
    })
  })
  
  # 自动刷新任务表（每 10 秒）
  observe({
    invalidateLater(10000, session)
    df <- list_tasks()
    output$tasks_table <- DT::renderDT({
      DT::datatable(df, selection = 'single', rownames = FALSE, options = list(pageLength = 10))
    })
  })
  
  # 渲染 Peaks / CNV 图（当前会话最近运行的任务）
  output$peaks_plot <- renderImage({
    req(rv$peaks_file)
    list(src = rv$peaks_file, contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)
  
  output$cnv_plot <- renderImage({
    req(rv$cnv_file)
    list(src = rv$cnv_file, contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)
  
  # 下载 handlers
  output$download_peaks <- downloadHandler(
    filename = function() { if (!is.null(rv$peaks_file)) basename(rv$peaks_file) else "Peaks.png" },
    content = function(file) { req(rv$peaks_file); file.copy(rv$peaks_file, file, overwrite = TRUE) },
    contentType = "image/png"
  )
  
  output$download_cnv <- downloadHandler(
    filename = function() { if (!is.null(rv$cnv_file)) basename(rv$cnv_file) else "CNV.png" },
    content = function(file) { req(rv$cnv_file); file.copy(rv$cnv_file, file, overwrite = TRUE) },
    contentType = "image/png"
  )
  
  # 自动清理：删除超过 24 小时的任务目录（放在 server 内）
  clean_old_tasks <- function() {
    dirs <- list.dirs(TEMP_DIR, recursive = FALSE, full.names = TRUE)
    now <- Sys.time()
    for (d in dirs) {
      info <- file.info(d)
      if (is.na(info$mtime)) next
      if (difftime(now, info$mtime, units = "hours") > 24) {
        unlink(d, recursive = TRUE, force = TRUE)
      }
    }
  }
  
  observe({
    invalidateLater(3600 * 1000, session)  # 每小时执行一次
    clean_old_tasks()
  })
}

# Run the application
shinyApp(ui = ui, server = server)

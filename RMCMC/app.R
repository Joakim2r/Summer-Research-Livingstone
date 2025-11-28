library(shiny)

# Helper: parse comma-separated inputs into vectors
parse_num <- function(x) {
  x <- trimws(unlist(strsplit(x, ",")))
  x[x == ""] <- NA
  y <- suppressWarnings(as.numeric(ifelse(tolower(x) == "na", NA, x)))
  y
}
parse_chr <- function(x) {
  x <- trimws(unlist(strsplit(x, ",")))
  x[x != ""]
}
parse_lgl <- function(x) {
  x <- trimws(unlist(strsplit(x, ",")))
  if (!length(x)) return(logical())
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

ui <- fluidPage(
  titlePanel("Sampler diagnostics (quick & dirty)"),
  sidebarLayout(
    sidebarPanel(
      textInput("posterior", "posterior_name",
                "eight_schools-eight_schools_noncentered"),
      textInput("samplers", "samplers (comma-separated)",
                "barker, amh, hmc, hmc_nuts"),
      textInput("adapters", "adapters (comma-separated; per run or grid for non-NUTS)",
                "adapter1, adapter2"),
      numericInput("chains", "chains", 2, min = 1),
      numericInput("warmup", "warmup_iter", 200, min = 0),
      numericInput("main", "main_iter", 2000, min = 1),
      numericInput("thin", "thin", 1, min = 1),
      textInput("nstep", "HMC n_step (comma-separated; use NA for no-fixed-step)",
                "5, 10, NA"),
      textInput("aux", "HMC sample_auxiliary (comma-separated TRUE/FALSE)",
                "TRUE, FALSE"),
      actionButton("run", "Run"),
      hr(),
      uiOutput("method_picker"),
      uiOutput("var_picker")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary table", tableOutput("tbl")),
        tabPanel("Curves",
                 h4("Running MSE (mean)"), plotOutput("plot_mse", height = 250),
                 h4("Running DT (log-variance distance)"), plotOutput("plot_dt", height = 250),
                 h4("Running ESJD"), plotOutput("plot_esjd", height = 250)
        )
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(eval = NULL, methods = NULL, vars = NULL)

  observeEvent(input$run, {
    withProgress(message = "Running samplers…", value = 0, {
      samplers <- parse_chr(input$samplers)
      adapters <- parse_chr(input$adapters)
      n_step   <- parse_num(input$nstep)
      aux      <- parse_lgl(input$aux)
      if (!length(samplers)) samplers <- "hmc"
      if (!length(adapters)) adapters <- "adapter1"
      if (!length(aux)) aux <- TRUE

      # Call your evaluate_samplers() (the “grid” logic for HMC/adapters lives inside)
      incProgress(0.1, detail = "Evaluating…")
      ev <- evaluate_samplers(
        posterior_name = input$posterior,
        samplers       = samplers,
        adapter        = adapters,
        chains         = input$chains,
        main_iter      = input$main,
        warmup_iter    = input$warmup,
        thin           = input$thin,
        n_step         = n_step,
        sample_auxiliary = aux,
        plot           = FALSE
      )
      rv$eval <- ev

      # discover method labels and variables
      methods <- names(ev$basic_all)
      rv$methods <- methods
      updateCheckboxGroupInput(
        session, "methods",
        label = "Which runs to show?",
        choices = methods, selected = methods
      )

      # Pick variable names from the first available summary
      any_tbl <- ev$basic_all[[methods[1]]]
      var_col <- if ("variable" %in% names(any_tbl)) "variable" else names(any_tbl)[1]
      vars <- sort(unique(any_tbl[[var_col]]))
      rv$vars <- vars
      updateSelectizeInput(
        session, "vars",
        choices = vars, selected = vars,
        server = TRUE
      )
    })
  })

  output$method_picker <- renderUI({
    req(rv$methods)
    checkboxGroupInput("methods", "Which runs to show?",
                       choices = rv$methods, selected = rv$methods)
  })

  output$var_picker <- renderUI({
    req(rv$vars)
    selectizeInput("vars", "Filter parameters (multiselect)", choices = rv$vars,
                   selected = rv$vars, multiple = TRUE)
  })

  # Summary table filtered by selected methods and parameters
  output$tbl <- renderTable({
    req(rv$eval, input$methods, input$vars)
    ev <- rv$eval
    out <- lapply(input$methods, function(m) {
      df <- ev$basic_all[[m]]
      var_col <- if ("variable" %in% names(df)) "variable" else names(df)[1]
      df <- df[df[[var_col]] %in% input$vars, , drop = FALSE]
      cbind(Sampler = m, df)
    })
    do.call(rbind, out)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", width = "100%")

  # Generic plotter for curves matrix (columns = methods)
  plot_curves <- function(mat, selected) {
    if (is.null(mat)) { plot.new(); title("No data"); return() }
    keep <- intersect(selected, colnames(mat))
    if (!length(keep)) { plot.new(); title("No selected runs"); return() }
    y <- mat[, keep, drop = FALSE]
    matplot(y, type = "l", lty = 1, xlab = "Iteration", ylab = "", main = "",
            col = seq_len(ncol(y)))
    legend("topright", legend = colnames(y), lty = 1, col = seq_len(ncol(y)), cex = 0.7)
  }

  output$plot_mse  <- renderPlot({ req(rv$eval, input$methods); plot_curves(rv$eval$mse_curves,  input$methods) })
  output$plot_dt   <- renderPlot({ req(rv$eval, input$methods); plot_curves(rv$eval$dt_curves,   input$methods) })
  output$plot_esjd <- renderPlot({ req(rv$eval, input$methods); plot_curves(rv$eval$esjd_curves, input$methods) })
}

shinyApp(ui, server)

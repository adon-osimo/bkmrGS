library(shiny)

ui <- fluidPage(
  
  tags$head(
    tags$script(HTML("
    MathJax = {
      tex: { inlineMath: [['$','$'], ['\\\\(','\\\\)']] }
    };

    Shiny.addCustomMessageHandler('typeset', function(message) {
      if (window.MathJax) {
        MathJax.typesetPromise();
      }
    });
  "))
  ),
  
  withMathJax(),
  
  titlePanel("BKMR Plot Generator"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      textInput(
        "fit_name",
        "Fit Object Name",
        placeholder = "e.g. fitkm"
      ),
      
      hr(),
      
      selectInput(
        "exposure",
        "(What is changing?) Exposure Variation",
        choices = c("Overall", "Single", "Bivariate")
      ),
      
      #selectInput(
      #  "analysis",
      #  "(What is being compared?) Difference Type:",
      #  choices = c("Group Specific", "Between Groups", "Response Function")
      #),
      
      uiOutput("analysis_ui"),
      uiOutput("exposure_Overall"),
      uiOutput("analysis_GS"),
      uiOutput("analysis_BG"),
      uiOutput("exposure_Single_1"),
      uiOutput("exposure_Single_2"),
      uiOutput("analysis_er"),
      uiOutput("interpretation"),
      
      actionButton(
        "prev",
        "Create Preview Plot"
      ),
      
      actionButton(
        "go",
        "Create Full Plot"
      )
      
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel(
          "Preview Plot",
          br(),
          plotOutput("preview")
        ),
        
        tabPanel(
          "Full Plot",
          br(),
          plotOutput("full")
        ),
        
        tabPanel(
          "Generated Code",
          br(),
          tags$pre(
            verbatimTextOutput("generated_code")
          )
        ),
        
        tabPanel(
          "Results",
          br(),
          tableOutput("results_table")
        ),
        
        tabPanel(
          "Interpretation",
          br(),
          uiOutput("interpretation")
        )
        
      )
    )
  )
)

server <- function(input, output, session){
  
  #---------------------------------------
  # Retrieve fit object
  #---------------------------------------
  
  fit <- reactive({
    
    req(input$fit_name)
    
    validate(
      need(
        exists(input$fit_name, envir = .GlobalEnv),
        "Fit object not found."
      )
    )
    
    get(
      input$fit_name,
      envir = .GlobalEnv
    )
    
  })
  
  output$analysis_ui <- renderUI({
    
    req(input$exposure)
    
    if(input$exposure == "Overall"){
      
      selectInput(
        "analysis",
        "(What is being compared?) Analysis Type",
        choices = c(
          "Group Specific",
          "Between Groups"
        )
      )
      
    } else if(input$exposure == "Single"){
      
      selectInput(
        "analysis",
        "(What is being compared?) Analysis Type",
        choices = c(
          "Group Specific",
          "Between Groups",
          "Response Function"
        )
      )
      
    } else if(input$exposure == "Bivariate"){
      
      selectInput(
        "analysis",
        "(What is being compared?) Analysis Type",
        choices = c(
          "Response Function"
        )
      )
      
    }
    
  })
  
  output$exposure_Overall <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(input$exposure != "Overall"){
      return(NULL)
    }
    
    tagList(
      
      textInput(
        "qs",
        "(qs) exposure Quantiles (e.g. c(0.25, 0.5, 0.75) or seq(0.25, 0.75, by = 0.05))",
        value = "seq(0.25, 0.75, by = 0.05)"
      ),
      
      numericInput(
          "q.fixed",
          "(q.fixed) Quantile Fixed At",
          value = 0.5,
          min = 0,
          max = 1,
          step = 0.05
        )

    )
    
  })
  
  output$analysis_GS <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(input$analysis != "Group Specific"){
      return(NULL)
    }
    
    tagList(
      
      selectInput(
        "m.fixed",
        "(m.fixed) Fixed Modifier",
        choices = c(unique(as.character(fit()$modifier)), "All")
      )
      
    )
    
  })
  
  output$analysis_BG <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(input$analysis != "Between Groups"){
      return(NULL)
    }
    
    tagList(
      
      textInput(
        "mod.diff",
        "(mod.diff) Modifier Values to Compare (e.g. c('Group_1', 'Group_2'))",
        value = paste0("c('", unique(as.character(fit()$modifier))[1], "' ,'", unique(as.character(fit()$modifier))[2], "' )")
      )
      
    )
    
  })  
  
  output$exposure_Single_1 <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(!(input$analysis == "Group Specific" && input$exposure == "Single")){
      return(NULL)
    }
    
    tagList(
      
      textInput(
        "qs.diff",
        "(qs.diff) Quantiles to Compare",
        value = "c(0.25, 0.75)"
      ), 
      
      numericInput(
        "q.fixed",
        "(q.fixed) Quantile Fixed At",
        value = 0.5,
        min = 0,
        max = 1,
        step = 0.05
      )
      
    )
    
  }) 
  
  output$exposure_Single_2 <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(!(input$analysis == "Between Groups" && input$exposure == "Single")){
      return(NULL)
    }
    
    tagList(
      
      textInput(
        "qs.diff",
        "(qs.diff) Quantiles to Compare",
        value = "c(0.25, 0.75)"
      ), 
      
      numericInput(
        "q.fixed",
        "(q.fixed) Quantile Fixed At",
        value = 0.5,
        min = 0,
        max = 1,
        step = 0.05
      )
      
    )
    
  }) 
  
  output$analysis_er <- renderUI({
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(!(input$analysis == "Response Function" && input$exposure == "Single")){
      return(NULL)
    }
    
    tagList(
    
      numericInput(
        "ngrid",
        "(ngrid) Ngrid",
        value = 50,
        min = 1,
        max = 100,
        step = 1
      ),
      
      numericInput(
        "q.fixed",
        "Quantile Fixed At",
        value = 0.5,
        min = 0,
        max = 1,
        step = 0.05
      )
    )
    
  })
  
  output$interpretation <- renderUI({
    
    interpretations <- list(
      list("Overall", "Group Specific"),
      list("Single", "Group Specific"),
      list("Overall", "Between Groups"),
      list("Single", "Between Groups"),
      list("Single", "Response Function")
    )
    
    tagList(
      lapply(interpretations, function(x) {
        
        active <- input$exposure == x[[1]] && input$analysis == x[[2]]
        
        div(
          style = paste(
            "color:", if (active) "black" else "lightgray",
            "; margin-bottom: 15px;",
            if (active) "font-weight: bold;" else ""
          ),
          
          h4(paste(x[[1]], " Exposure with ", x[[2]], " Movement type")),
          
          div(
            style = "white-space: pre-wrap;",
            withMathJax(HTML(generate_interpretation(x[[1]], x[[2]])))
          )
        )
      })
    )
  })
  
  observe({
    input$exposure
    input$analysis
    session$sendCustomMessage("typeset", list())
  })

  generated_code_preview <- eventReactive(input$prev, {
    generate_code(
      fit_name = input$fit_name,
      exposure = input$exposure,
      analysis = input$analysis,
      centered = input$centered,
      sel = '1:10',
      m.fixed = input$m.fixed,
      qs = input$qs,
      q.fixed = input$q.fixed,
      qs.diff = input$qs.diff, 
      mod.diff = input$mod.diff,
      ngrid = input$ngrid
    )
  })
  
  generated_code <- reactive({
    
    req(input$fit_name)
    
    generate_code(
      fit_name = input$fit_name,
      exposure = input$exposure,
      analysis = input$analysis,
      centered = input$centered,
      m.fixed = input$m.fixed,
      qs = input$qs,
      q.fixed = input$q.fixed,      
      qs.diff = input$qs.diff, 
      mod.diff = input$mod.diff,
      ngrid = input$ngrid
    )
  })
  
  plot_preview <- eventReactive(input$prev, {
    code_text <- generated_code_preview()
    prev <- new.env(parent = globalenv())
    eval(parse(text = code_text), envir = prev)
    prev
  })
  
  plot_env <- eventReactive(input$go,{
    code_text <- generated_code()
    env <- new.env(parent = globalenv())
    eval(parse(text = code_text), envir = env)
    env
  })
  
  output$preview <- renderPlot({
    env <- plot_preview()
    env$plot_obj
    
  })
  
  output$full <- renderPlot({
    env <- plot_env()
    env$plot_obj
    
  })
  
  output$generated_code <- renderText({
    
    generated_code()
  })
  
  output$results_table <- renderTable({
    
    if(input$go > input$prev){
      env <- plot_env()
    } else if(input$prev > 0){
      env <- plot_preview()
    } else{
      return(NULL)
    }
    
    env$results
  
  })
  
}

shinyApp(ui, server)

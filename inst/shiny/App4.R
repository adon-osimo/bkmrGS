library(shiny)

ui <- fluidPage(
  
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
        "comparison",
        "Exposure Structure",
        choices = c("Overall", "Single", "Bivariate")
      ),
      
      selectInput(
        "movement",
        "Movement Type",
        choices = c("Overall", "Group Specific")
      ),
      
      uiOutput("comparison_Overall"),
      uiOutput("movement_Overall"),
      uiOutput("movement_GS"),
      uiOutput("comparison_Single_1"),
      uiOutput("comparison_Single_2"),
      
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
  
  output$comparison_Overall <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(input$comparison == "Single"){
      return(NULL)
    }
    
    tagList(
      
      numericInput(
        "q.fixed",
        "(q.fixed) Quantile Fixed At",
        value = 0.5,
        min = 0,
        max = 1,
        step = 0.05
      ),
    
      textInput(
        "qs",
        "(qs) Comparison Quantiles (e.g. c(0.25, 0.5, 0.75) or seq(0.25, 0.75, by = 0.05))",
        value = "seq(0.25, 0.75, by = 0.05)"
      )
    )
    
  })
  
  output$movement_Overall <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(input$movement == "Group Specific"){
      return(NULL)
    }
    
    tagList(
      
      selectInput(
        "m.fixed",
        "(m.fixed) Fixed Modifier",
        choices = unique(as.character(fit()$modifier))
      )
      
    )
    
  })
  
  output$movement_GS <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(input$movement == "Overall"){
      return(NULL)
    }
    
    tagList(
      
      textInput(
        "mod.diff",
        "(mod.diff) Modifier Values to Compare (e.g. c('Group_1', 'Group_2'))",
        value = "c('', '')"
      )
      
    )
    
  })  
  
  output$comparison_Single_1 <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(!(input$movement == "Overall" && input$comparison == "Single")){
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
  
  output$comparison_Single_2 <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier)){
      return(NULL)
    }
    
    if(!(input$movement == "Group Specific" && input$comparison == "Single")){
      return(NULL)
    }
    
    tagList(
      
      textInput(
        "qs.diff",
        "(qs.diff) Quantiles to Compare",
        value = "c(0.25, 0.75)"
      ), 
      
      textInput(
        "qs.fixed",
        "(qs.fixed) Quantiles to Be Fixed",
        value = "c(0.5, 0.5)"
      )
      
    )
    
  }) 
  
  #---------------------------------------
  # Generate code
  #---------------------------------------
  
  generated_code_preview <- eventReactive(input$prev, {
    generate_code(
      fit_name = input$fit_name,
      comparison = input$comparison,
      movement = input$movement,
      centered = input$centered,
      sel = '5:10',
      m.fixed = input$m.fixed,
      qs = input$qs,
      q.fixed = input$q.fixed,
      qs.diff = input$qs.diff, 
      mod.diff = input$mod.diff, 
      qs.fixed = input$qs.fixed
    )
  })
  
  generated_code <- reactive({
    
    req(input$fit_name)
    
    generate_code(
      fit_name = input$fit_name,
      comparison = input$comparison,
      movement = input$movement,
      centered = input$centered,
      m.fixed = input$m.fixed,
      qs = input$qs,
      q.fixed = input$q.fixed,      
      qs.diff = input$qs.diff, 
      mod.diff = input$mod.diff, 
      qs.fixed = input$qs.fixed
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
    
    env <- plot_env()
    env$results
  })
  
}

shinyApp(ui, server)

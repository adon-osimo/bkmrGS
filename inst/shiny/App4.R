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
      
      uiOutput("modifier_controls"),
      
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
  
  output$modifier_controls <- renderUI({
    
    req(fit())
    
    if(is.null(fit()$modifier))
      return(helpText("No modifier found in fit object."))
    
    tagList(
      
      selectInput(
        "m.fixed",
        "Fixed Modifier",
        choices = unique(as.character(fit()$modifier))
      ),
      
      numericInput(
        "q.fixed",
        "Quantile Fixed At",
        value = 0.5,
        min = 0,
        max = 1,
        step = 0.05
      ),
      
      textInput(
        "qs",
        "Comparison Quantiles (e.g. c(0.25, 0.5, 0.75) or seq(0.25, 0.75, by = 0.05))",
        value = "seq(0.25, 0.75, by = 0.05)"
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
      q.fixed = input$q.fixed
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
      q.fixed = input$q.fixed
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

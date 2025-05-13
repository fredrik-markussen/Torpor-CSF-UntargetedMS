library(shiny)
library(plotly)
library(dplyr)
library(reshape2)
library(ggplot2)


rsconnect::setAccountInfo(name='60afyc-fredrik0markussen', token='EC4C048045D615675794BE7BC5A6B8FF', secret='HlfchM40J2wWcwks8zWV7hI7iXJQpXvqVUKHwThc')



# Load your data
followtb <- read.csv("./metabolitesfollowstb.csv")
pretb <- read.csv("./metabolitesPREtb.csv")
synchtb <- read.csv("./metabolitesSYNCHtb.csv")



# UI
ui <- fluidPage(
  titlePanel("Metabolite Dynamics vs Body Temperature"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("metGroup", "Choose Metabolite Timing Group:",
                  choices = c("Follows Tb" = "followtb",
                              "Precedes Tb" = "pretb",
                              "Synchronous with Tb" = "synchtb"),
                  selected = "followtb")
    ),
    
    mainPanel(
      plotlyOutput("corPlot", height = "800px")
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive data selection
  selectedData <- reactive({
    switch(input$metGroup,
           followtb = followtb,
           pretb = pretb,
           synchtb = synchtb)
  })
  
  output$corPlot <- renderPlotly({
    df <- selectedData()
    
    # Create the plot
    p <- ggplot(df, aes(x = sample_nr*2, y = value, color = variable, group = variable)) +
      geom_line(data = df, aes(x = sample_nr*2, y = tb_mean/10-2), size=2, color = "gray20")+
      geom_point(alpha = 0.8) +
      geom_line() +
      geom_smooth(
        aes(group = variable),
        method = "loess", span = 0.3, se = FALSE,
        linetype = "solid", size = 0.9
      ) +
      scale_y_continuous(
        name = "Metabolite Abundance (z-score)",
        limits = c(-1.3, 1.6),
        sec.axis = sec_axis(~ . * 10 + 20, name = "Body Temperature (°C)")
      ) +
      scale_x_continuous(
        name = "Time (hours)",
        breaks = seq(0, 220, 4)
      ) +
      theme_classic()  +
      labs(
        title = paste("Metabolite Group:", input$metGroup),
        color = "Metabolite"
      )
    
    ggplotly(p)
  })
}

# Run the app
shinyApp(ui, server)






#-----------------------------------------
# Shiny web application to view the 
# enteric antibody response distributions
# for the different datasets
#
#-----------------------------------------

#-----------------------------------------
# preamble
#-----------------------------------------
library(here)
here()

library(shiny)
library(tidyverse)

# bright color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"

#-----------------------------------------
# load the underlying datasets for viewing
#-----------------------------------------
d_leogane <- readRDS(here("data","haiti_analysis2.rds"))
d_kongwa <- readRDS(here("data","kongwa_analysis2.rds"))
d_asembo <- readRDS(here("data","asembo_analysis2.rds"))


# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Enteric Antibody Distribution Viewer"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Selector for choosing dataset ----
      selectInput(inputId = "dataset",
                  label = "Choose a dataset:",
                  choices = c("Kongwa, Tanzania", "Leogane, Haiti","Asembo, Kenya")),
      
      # Input: Select an antibody ----
      uiOutput("d_antigens"),
      
      # Input: Select ages ----
      uiOutput("d_ages")
      
    ),
  
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Formatted text for caption ----
      h3(textOutput("selected_ages", container = span)),
      
      h5("If available, the figure includes seropositivity cutoffs based on an ROC analysis (solid line) or mixture model (dashed line)."),
      
      
      # Output: HTML table with requested number of observations ----
      plotOutput("fig",width="100%",height="600px")
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  # Return the requested dataset ----
  # By declaring datasetInput as a reactive expression we ensure
  # that:
  #
  # 1. It is only called when the inputs it depends on changes
  # 2. The computation and result are shared by all the callers,
  #    i.e. it only executes a single time
  datasetInput <- reactive({
    switch(input$dataset,
           "Kongwa, Tanzania" = d_kongwa,
           "Leogane, Haiti" = d_leogane,
           "Asembo, Kenya" = d_asembo
           )
  })
  
  # Output: List of antigens included in the dataset
  output$d_antigens <- renderUI({
    dataset <- datasetInput()
    antigens <- unique(dataset$antigenf)
    selectInput("antigens", "Choose an antibody:", antigens)
  })
  
  # Output: Age range in the dataset
  output$d_ages <- renderUI({
    dataset <- datasetInput()
    agemin <- floor(min(dataset$age))
    agemax <- ceiling(max(dataset$age))
    sliderInput("ages","Age range",
                min=agemin,max=agemax,
                value=c(agemin,agemax)
                )
  })
  
  # Output: Age range description
  output$selected_ages <- renderText({
    if(input$ages[1]!=input$ages[2]) {
      paste(input$dataset,", Ages ",input$ages[1]," to ",input$ages[2],sep="")
    } else{
      paste(input$dataset,", Ages ",input$ages[1],sep="")
    }

  })
  
  
  # Show the distribution of the data ----
  # The output$view depends on the databaseInput reactive
  # expression, the antigen selected, and the age range
  # so it will be re-estimated when any one of them changes
  output$fig <- renderPlot({
    dl <- datasetInput() %>%
      filter(age >= input$ages[1] & age <= input$ages[2]) %>%
      filter(antigenf %in% input$antigens)
    
    # custom color blind color palette is in the preamble chunck
    # pcols <- c(cred,corange,cchartr,cgreen,cteal,cblue,cmagent)
    plotantigen <- dl$antigen[1]
    pcols <- cblue
    if(plotantigen %in% c("vsp3","vsp5")) pcols <- cblue
    if(plotantigen %in% c("cp17","cp23")) pcols <- cteal
    if(plotantigen %in% c("leca")) pcols <- cgreen
    if(plotantigen %in% c("salb","sald")) pcols <- corange
    if(plotantigen %in% c("norogi","norogii")) pcols <- cred
    if(plotantigen %in% c("etec","cholera")) pcols <- cmagent
    
    # custom log labels
    log10labs <- c( 
      expression(10^0),
      expression(10^1),
      expression(10^2),
      expression(10^3),
      expression(10^4)
    )
    
    p <- ggplot(data=dl,aes(x=logmfi,group=pathogen,color=pathogen,fill=pathogen)) +
      # facet_wrap(~antigenf,nrow=6,ncol=2,scales="free_y") +
      geom_histogram(aes(y=..density..),bins=50,alpha=0.7) +
      geom_density(aes(fill=NULL),color=cgrey) +
      geom_vline(aes(xintercept=min(roccut)),size=1.5)+
      geom_vline(aes(xintercept=min(mixcut)),lty="dashed",size=1.5)+
      # geom_vline(aes(xintercept=mixcut),linetype="dashed",size=1.5)+
      scale_x_continuous(limits = c(0,4.5),breaks = 0:4,labels=log10labs) +
      # coord_cartesian(ylim=c(0,1))+
      scale_fill_manual(values=pcols) +
      scale_color_manual(values=pcols) +
      labs(x="Luminex Response (MFI-bg)",y="Density",title=paste(input$antigens)) +
      theme_minimal(base_size=20) +
      theme(legend.position = "none")
    
    p
    
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)


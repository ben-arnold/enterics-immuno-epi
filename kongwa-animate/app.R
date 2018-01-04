
#-----------------------------------------
# Shiny web application to make an
# interactive animation of the enteric
# antibody distributions in the Kongwa
# Tanzania study
#
#-----------------------------------------

#-----------------------------------------
# preamble
#-----------------------------------------
rm(list=ls())
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
# load the underlying dataset for viewing
#-----------------------------------------
d_kongwa <- readRDS("~/enterics-immuno-epi/data/kongwa_analysis2.rds")

# drop antigens that were not measured in >1 year
d_kongwa <- d_kongwa %>%
  filter(!antigen %in% c("salb","sald","p18","p39"))


# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Kongwa, Tanzania Study Enteric Antibody Distributions"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select an antibody ----
      uiOutput("d_antigens"),
      
      # Input: Select ages ----
      uiOutput("d_ages")
      
    ),
  
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Formatted text for caption ----
      h3(textOutput("selected_ages", container = span)),
      
      h5("If available, the figure includes seropositivity cutoffs based on an ROC analysis (solid line)."),
      
      # Output: HTML table with requested number of observations ----
      plotOutput("fig",width="100%",height="600px")
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  # Output: List of antigens included in the dataset
  output$d_antigens <- renderUI({
    antigens <- unique(d_kongwa$antigenf)
    selectInput("antigens", "Choose an antibody:", antigens)
  })
  
  # Output: Age range in the dataset
  output$d_ages <- renderUI({
    agemin <- floor(min(d_kongwa$age))
    agemax <- ceiling(max(d_kongwa$age))
    sliderInput("ages","Age in years (press play to animate)",
                min=agemin,max=agemax,value=agemin,
                step=1,animate=TRUE
                )
  })
  
  # Output: Age selected
  output$selected_ages <- renderText({ 
    paste("Kongwa, Age ",input$ages, " years",sep="")
  })
  
  
  # Show the distribution of the data ----
  # The output$view depends on the databaseInput reactive
  # expression, the antigen selected, and the age range
  # so it will be re-estimated when any one of them changes
  output$fig <- renderPlot({
    dl <- d_kongwa %>%
      filter(floor(age) >= input$ages & floor(age) <= input$ages) %>%
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
      geom_vline(aes(xintercept=log10(min(roccut))),size=1.3)+
      # geom_vline(aes(xintercept=mixcut),linetype="dashed",size=1.3)+
      geom_text(aes(x=0.5,y=0.9,label=paste(input$ages)),size=14,color=cgrey) +
      scale_x_continuous(limits = c(0,4.5),breaks = 0:4,labels=log10labs) +
      coord_cartesian(ylim=c(0,1))+
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


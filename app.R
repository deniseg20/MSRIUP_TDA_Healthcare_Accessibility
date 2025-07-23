#installing packages

options(repos = c(CRAN = "https://cran.rstudio.com/"))
#install.packages("htmltools")
#install.packages('shinydashboard')

library(htmltools)
library(flexdashboard)
library(shiny)        # For interactive elements
library(ggplot2)      # For visualizations
library(dplyr)        # For data manipulation
library(plotly)       # For interactive plots
library(readr)        # For data loading
library(tidyr)        # For data reshaping
library(scales)
library(viridis)
library(leaflet)
library(shinydashboard)


#Creating Dashboard aesthetic
ui<-dashboardPage(skin='purple',
                  dashboardHeader(title='Defunding Sexual healthcare'),
                  dashboardSidebar(sidebarMenu(
                    menuItem("Intro", tabName = "Intro", icon = icon("default")),
                    menuItem("Healthcare Gaps", tabName = "gaps", icon=icon("default")),
                    menuItem("What Now?",tabName="end",icon=icon('default'))
                  )
                  ),
                  dashboardBody(
                    tabItems(
                      
                      #first tab
                      
                      tabItem(tabName = "Intro",
                              fluidPage(
                                #title of the dashboard and explanation of dashboard
                                h1("Defunding Sexual Healthcare: Visualizing the Problem"),
                                HTML("<center><img src='PPHC_Timline.png', height='600px',width='1100px'></center>"),
                                br(),
                                h3('In the last 7 years, threats against the federal funding of Planned Parenthood have only escalated. As an organization that heavily depends on federal funding, specifically Medicaid reimbursements, this can disrupt sexual healthcare access for Medicaid beneficiaries. While Federally Qualified Health Centers (FQHCs) have been considered possible replacements because they also provide covered services to Medicaid beneficiaries, this plan ignores the locational accessibility of sexual healthcare when Planned Parenthood Health Centers (PPHCs) are no longer in option. Here, we visual this dilemma.'),
                                #caption describing the heat maps
                                br(),
                                align="center",
                                h3("The heatmaps below describe the concentration of healthcare facilities in the State of California. Switching between maps, one can observe that the removal of PPPHCs, representing a Medicaid user's inability to receive care in the facilities if they are federally defunded, decreases sexual health care opportunities. We examine California because it is one of the states with the highest percentage of Medicaid beneficiaries."),
                                br(),
                                #inputs for the heatmap
                                selectInput("selected_map", "Choose Facility Type:",
                                            choices = c('Only PPHCs','Only FQHCs','Both PPHCs and FQHCs'),
                                            selected = "Only PPHCs"),
                                
                                #laying out the heatmaps
                                leafletOutput("map")
                              )
                      ),
                      #laying out the holes from the persistence diagrams onto maps of california
                      
                      
                      #second tab
                      
                      tabItem(tabName = "gaps",
                            align="center",
                            h1("Persistant Homology to Represent Gaps in Sexual Healthcare Coverage"),
                              h4("Persistant Homology is a sect of of topological data analysis that assesses the significance holes in data. Our data is FQHCs and PPHCs in California, a state with one of the highest concentrations of Medicaid beneficiaries, and the distance between them, minimum time is takes to get from one location to another. The holes in this data represent low sexual health coverage. We present diagrams where FQHCs and PPHCs are present and diagrams where PPHCs are absent to model a Medicaid beneficiaries lack of access to PPHCs."),
                              br(),
                              br(),
                           
                            fluidRow(
                              align='center',
                                column(width = 6,
                                       box(
                                         align='center',
                                         title = "0D Homology for FQHC and PPHC Locations", width = NULL, status = "primary", solidheader=TRUE,
                                         img(src="0D_All_Locations.png", align = "right", height='600px',width='600px')
                                       ),
                                       box(
                                         align='center',
                                         title = "1D Homology for FQHC and PPHC Locations", width = NULL, status = "primary",solidheader = TRUE,
                                         img(src="1D_All_Locations.png", align = "right", height='600px',width='600px')
                                       )),
                            fluidRow(
                              align='center',
                                column(width = 6,
                                       box(
                                         align='center',
                                         title = "0D Homology for FQHC Locations", width = NULL, status = "primary", solidheader=TRUE,
                                         img(src="OD_FQHC_Locations.png", align = "right", height='600px',width='600px') 
                                       ),
                                       box(
                                         align='center',
                                         title = "1D Homology for FQHC Locations", width = NULL, status = "primary", solidheader = TRUE,
                                         img(src="1D_FQHC_Locations.png", align = "right", height='600px',width='600px')
                                       )
                                )
                              ),
                            h3("These 0D homology diagrams, homology referring to the holes, shows where there is low sexual healthcare resources. Each line is connected to two healthcare facilities. The colors of the lines between the locations shows the amount of time it takes for one to get to a healthcare center in that region. When PPHCs are removed, sexual healthcare coverage worsens as it takes more time to reach a facilitiy."),
                            br(),
                            h3("These 1D homoloy diagrams, similarly to the 0D homology diagrams, show regions of poor sexual healthcare coverage. The color of each hole, a triangulation between three healthcare locations, shows the amount of time in minutes it takes for an individual in the region to access a health center. When PPHCs are removed, new gaps in coverage are created, adding to pre-existent coverage gaps from when PPHCs were present."),  
                            
                              )
                      
                              ),
                      
                      #third tab
                      tabItem(tabName='end',
                              align='center',
                              fluidRow(
                              align='center',
                              h1('We Found the Problem... Now What?'),
                              h3("While we found out that there are new holes in coverage with the 1D homology diagrams, we have to see if the greater distance between facilities shown in the 0D homology diagram, after removing the PPHCs, is significant. We do this using the box plot below and significant testing. After completing a t-test, we realize the difference is signiciant."),
                              br(),
                              
                              align='center',
                              img(src="Box_Plot.png",height='600px',width='600px'),

                               h3("Removing PPHCs both creates new gaps in sexual healthcare coverage and makes it harder to reach sexual healthcare facilities.Therefore, the government attacks against Planned Parenthood, under the guise of protecting life, is endangering the lives of Medicaid beneficiaries. Planned Parenthood is a key source of sexual healthcare for Medicaid beneficiaries, potentially reflecting the experiences of Medicaid beneficiaries across the county. Either it is necessary that this fact is recognized, or the government puts equal effort into expand affordable sexual healthcare."),
                              ))
                              
                        
                  )
)
)



#running the server

server<-function(input,output){
  #Heatmap
  output$map <- renderLeaflet({
    # 1. Filter your data based on input
    if (input$selected_map == "Both PPHCs and FQHCs") {
      map<-leaflet_map3
    } else if (input$selected_map == "Only PPHCs") {
      map<-leaflet_map1
    } else {
      map<-leaflet_map2
    }
    map
    
  })
  #images
  #output$ZD_All_Locations <- renderUI({
   # tags$img(src='0D_All_Locations.png',
    #         style = "width: 100%; height: auto;")
  #})
}

shinyApp(ui, server)


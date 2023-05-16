#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(Seurat)
library(ape)
library(limma)
library(tidyverse)
#library(dplyr)
library(shiny)
library(ggplot2)


# Define UI for Hippocampus and SUB-PROS seurat
ui <- fluidPage(

    # Application title
    titlePanel("HIP and SUB data"),

    # Sidebar with several tab separated side pannel
    #   1. maps
    #       umap/tsne:          choose from umap or tsne plot
    #   2. differential expression
    #       ident reference:    differential expression for main ident 
    #       Decide:             activated button for ident variable 
    #   3. plots
    #       gene name:          user input gene name for feature/violin plot
    #       decide:             activated button for gene name
    #   4. compare two cluster
    #       ident1:             first ident
    #       ident2:             second ident
    #       dif:                different between ident1 and ident2
    #       go2:                activated button for all other 3 variables
    sidebarLayout(
        sidebarPanel(
            tabsetPanel(
                tabPanel(
                    "map",
                    selectInput("maps","Umap/Tsne",
                                choices = c("umap","tsne"),
                                selected = "umap")
                ),
                tabPanel(
                    "differential expression",
                    selectInput("ident","indent reference",
                                choices = levels(HipSub.s4obj.cluster),
                                selected = "CA1d"),
                    actionButton("go","Decide"),
                ),
                tabPanel(
                  "plots",
                  textInput("gene","gene name",value = "Dcn"),
                  actionButton("go1","decide")
                ),
                tabPanel(
                    "compare two  cluster",
                    selectInput("ident1","First Ident",
                                choices = levels(HipSub.s4obj.cluster),
                                selected = "CA1d"),
                    selectInput("ident2","Second Ident",
                                choices = levels(HipSub.s4obj.cluster),
                                selected = "CA1v"),
                    sliderInput("dif","difference",
                                min = 0,max = 0.9,value = 0,step = 0.1),
                    actionButton("go2","decide")
                ),
            ),
        ),

        # Show several tab separated plot 
        #   1. maps: plot umap  or tsne plot 
        #   2. ident: ident for differential expression
        #   3. feature: feature plot 
        #      violin: violin plot
        #   4. compare: compare two specific idents
        #       
        mainPanel(
            tabsetPanel(
                tabPanel("map",
                         plotOutput("maps"),
                         ),
                tabPanel("differential expression",
                         dataTableOutput("ident"),
                         ),
                tabPanel("plots",
                         plotOutput("feature"),
                         plotOutput("violin")
                         ),
                tabPanel("compare two cluster",
                         dataTableOutput("compare")
                         ),
            ),
        )
    )
)

# Define server logic required to display several tab panel for mainPanel
server <- function(input, output) {

    output$maps <- renderPlot({
        DimPlot(HipSub.s4obj.cluster,reduction = input$maps)
    })
    
    active<- eventReactive(input$go, {
        
        input$ident
    })
    
    output$ident <- renderDataTable({
        HipSub<-FindMarkers(HipSub.s4obj.cluster,
                            ident.1 = active(),
                            only.pos = T)
        HipSub<-mutate(HipSub,gene= rownames(HipSub))
        HipSub
    })
    
    gene<-eventReactive(input$go1,{
        input$gene
    })
    
    output$feature <- renderPlot({
        FeaturePlot(HipSub.s4obj.cluster,features = gene())
    })
    
    output$violin <- renderPlot({
        VlnPlot(HipSub.s4obj.cluster,features = gene())
    })
    
    ident1<-eventReactive(input$go2,{
        input$ident1
    })
    
    ident2<-eventReactive(input$go2,{
        input$ident2
    })
    
    dif<-eventReactive(input$go2,{
        input$dif
    })
    output$compare <- renderDataTable({
        HipSub<-FindMarkers(HipSub.s4obj.cluster,
                            ident.1 = ident1(),
                            ident.2 = ident2(),
                            only.pos = T,
                            min.diff.pct = dif())
        HipSub<-mutate(HipSub,gene=rownames(HipSub))
        HipSub
    })
}

# Run the application 
shinyApp(ui = ui, server = server)


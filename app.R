## for sig, bar and violin plots ##
#20180920

library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(googlesheets) # for loading the GeoPT data from Google drive
library(readr)
library(MASS) # for hubers estimate
library(modeest) # for mode calculations
library(stringr) #for combining letters and numbers for lab selection
library(plotflow) # needed for reorder in plots


loadData <- function() {
  # Grab the Google Sheet
  sheet <- gs_title("GeoPTall") 
  # Read the data
  gs_read_csv(sheet, col_types = cols(
    lab = col_character(),
    analyte = col_character(),
    method = col_character(),
    prep = col_character(),
    mass = col_double() ,
    measurand = col_double(),
    round = col_character()
  ))
}



GeoPT <- loadData() # file in Google Drive needs to be named GeoPTall and not GeoPTall.csv
                    # need to check if headers are available
                    # need to upload latest version of "~/GitHub/GeoPT02app/GeoPTall.csv" then open with Google Sheets
                    # and rename GeoPTall.csv to GeoPTall
#GeoPT <- read.csv("~/GitHub/GeoPT02app/GeoPTall.csv")

GeoPT <- as.data.table(GeoPT)
mort <- c("SiO2", "TiO2", "Al2O3", "Fe2O3T", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "LOI", "Fe(II)O", "CO2", "H2O+")  # needed for assigning T for trace or M for major element as is used for unit and z-score calcuations

GeoPT$MorT <- "T" # setting all to trace element "T
GeoPT[analyte %in% mort, MorT := "M"] # exchanging the major element analytes with "M" for major

GeoPT$unit <- "mg/kg" # assigning unit used for ggplot2
GeoPT$unit[GeoPT$MorT == "M"] <- "g/100g"

## introduced to simplify the digestion and measurement methods ##

setkey(GeoPT, prep)


GeoPT[prep == "FM+AD", prep := "FM_AD"] #need to change for scale_fill_manual colours
GeoPT[prep == "AD+FM", prep := "AD_FM"]  #need to change for scale_fill_manual colours
GeoPT[method == "ICP-OES/AES",method := "ICP"]  #need to change for scale_fill_manual colours
GeoPT[method == "ICP-MS", method := "ICP_MS"]  #need to change for scale_fill_manual colours
GeoPT[method == "GRAV+VOL", method := "GRAV_VOL"] #need to change for scale_fill_manual colours


anal <- unique(GeoPT$analyte) # needed for ui selector
rnd <- unique(GeoPT$round) # needed for ui selector

findrangeDT <- function(measurand) {
  ##test.f <- test.f[order(measurand)] # sorting by increasing value
  df <- as.data.table(measurand)
  df <- df[order(measurand)]
  #measurand <- order(measurand)
  nx <- length(measurand) # getting length of values
  i <- 1:(nx/2) # i to half the values needed to find minimum halfrange
  #df.f1 <- test.f[,measurand] # getting list
  
  test.range <- lapply(i, function(x){ # testing all ranges for minimum (max-min)
    max(df[x:(x+nx/2)]) - min(df[x:(x+nx/2)])
  })
  # df.range <- as.data.table(test.range)
  df.min.position <- which.min(test.range) # getting the position of the mininum range value
  xa <- df[df.min.position:(df.min.position+nx/2)] # getting the range
  xa <- xa[, measurand]
  return(xa)
}

lientz.mode <- function(x, bt, type, bw1){ # bt is  n for resampling bw1 is bandwidth
  #M <- mlv(findrangeDT(x), method = "lientz", bw = 0.5, boot = TRUE, R=300)
  TorF <- ifelse(type %in% c("lientz", "hsm", "kernel"), TRUE, FALSE)
  
  M <- mlv(x, method = type, bw = as.numeric(bw1), boot = TorF, R=bt)
  data <- M[["boot.M"]]
  mean <- mean(data)
  median <- median(data)
  sd <- sd(data)
  se <- sd/sqrt(length(x))
  semode <- median(abs(data-median),na.rm=TRUE)*1.5
  estimates <- list(mean, median, se, semode, data)
  #return(M[1])
  return(estimates)
}

newmode <- function(x, h){ # h is a smoothing parameter # x is column of data, no missing values
  x <- na.omit(x)
  nx <- length(x) #nx is the count of data in x
  if(nx > 3){
    top <- max(x)
    bottom <- min(x)
    increment <- (top-bottom)/20
    xpoint <- seq(bottom, top, increment) #Define points for calculating kernel density in column xpoint
    df <- as.data.frame(xpoint)
    # df$den <- 0 # initialise corresponding densities as zero in column den
    srob <- h #median(abs(x-median(x)))*1.5 
    df$d <- dnorm(df$xpoint, mean = median(x), sd = srob/2)
    df$den <- pnorm(df$xpoint, mean = median(x), sd = srob/2)
    dmax <- max(df$d)	#find highest density
    center <- which.max(df$d) # row index with highest density
    # find points at three higher and three lower than highest density
    low <- center -3
    high <- center + 3
    if (low<1) {low <- 1}	# this sequence prevents selection of undefined values of xpoint
    # ATTENTION I change this values from 5 to 1
    if (high > 21) {high <- 21}
    x2df <- df[low:high,] #Select seven xpoints and corresponding densities 
    x2df$subxsq <- x2df$xpoint^2
    model <- lm(x2df$d ~ poly(x2df$xpoint, 2, raw=TRUE))
    #pol2 <- function(x) model$coefficient[3]*x^2 + model$coefficient[2]*x + model$coefficient[1]
    #plot(x2df$xpoint, x2df$d, type="l", lwd=3)  
    #points(x2df$xpoint, predict(model), type="l", col="tomato1", lwd=2 )
    #  curve(pol2, col="tomato1", lwd=2)
    mode <- -as.numeric(model$coefficient[2]/(2*model$coefficient[3]))
    #plot(density(x))
    #M <- mlv(x, method = "asselin", bw = 2/3)
    #plot(M)
    #plot(x)
    return(mode)
  } else
    mode <- median(x)
  
}

estimator.all <- function(round.select, anal, bt, type, bw){ # xcol = analyte, radio = bt == # boot resampling
  #selectedData <- GeoPT[round == "40A" & analyte == "TiO2", measurand, by = .(lab, analyte, quality, round, method, prep)][order(round, measurand)]
  selectedData <- reactive({GeoPT[round == round.select & analyte == anal, measurand, by = .(lab, analyte, quality, round, method, prep)][order(round, measurand)]})
  
  xx <- selectedData()$measurand # extracting the data only
  xx <- na.omit(xx) # removing NA
  nx <- length(xx) # getting n
  
  ## Mike's mode
  xa <- findrangeDT(xx) # using function to define the data set window for mode calculation
  srob <- mad(xx)
  Mode <-  signif(newmode(xa, srob),5)  # Mode based on newmode window
  
  modeNM <- rep(0,bt)
  if(nx > 7) {
    for (i in 1:bt){ # bootstrapping by replacing one value of the entire data range xx
      sam <- sample(xx, nx, replace =TRUE)
      srob <- mad(sam)
      samp <- findrangeDT(sam)
      if(srob != 0) {modeNM[i] <- newmode(samp, srob)}
      
    }
    Mode_NM <- signif(median(modeNM, na.rm=TRUE),5)
    Semode <- signif(median(abs(modeNM-Mode_NM),na.rm=TRUE)*1.5,3)
    Mike.mode <- c(Mode_NM, Semode)
  } else {
    Mode_NM <- NA
    Semode <- NA
    Mike.mode <- c(Mode_NM, Semode)
  }
  
  ## median calculations for lines
  med.sig <- signif(median(xx),5)
  median_se <- signif(mad(xx)/sqrt(length(xx)),3)
  
  ## robust calculation for lines
  z <- hubers(xx)
  robhub <- signif(z$mu,5)# the robust estimate
  robhubs <- signif(z$s,3) # the robust standard deviation
  robhubse <- signif(z$s/sqrt(nx),3) # the robust standard error
  
  ## Lientz mode    
  estimators <- lientz.mode(xx, bt, type, bw)
  #mode_lientz_all <- signif(as.numeric(lientz.mode(xx)[2]),5) # median
  mode_lientz_mean <- signif(estimators[[1]],5) # mean
  mode_lientz_med <- signif(estimators[[2]],5) # median
  mode_lientz_sd <- signif(estimators[[3]],3) # standard deviation
  mode_lientz_se <- signif(estimators[[4]],3) # mad
  mode_lientz_data <- estimators[[5]]
  
  all.data <- list(anal, nx, med.sig, median_se, 
                   Mode, mode_lientz_data, mode_lientz_mean, 
                   mode_lientz_med, mode_lientz_sd, mode_lientz_se, 
                   robhub, robhubs, robhubse, Mike.mode[1], Mike.mode[2],  mode_lientz_data)
  return(all.data)
}



# Define UI for application that draws a histogram
ui <- fluidPage(
  navbarPage("Chart Selection",
             tabPanel("Sigmoidal",
                      pageWithSidebar(
                        headerPanel('GeoPT charts'),
                        sidebarPanel(
                          selectInput('round.select', 'GeoPT Round', rnd, selected=rnd[[49]]),
                          selectInput('xcol', 'analyte', anal,
                                      selected=anal[[1]]),
                          actionButton(inputId = "go", label = "Update calculations",style = "background-color:#E69F00"),
                          #h4("restore excluded values"),
                          actionButton("exclude_reset", "Reset"),
                          helpText("Note: Press Reset button before changing the analyte. In case an error msg occurs, press Reset button"),
                          radioButtons("sig_colour", label = "select point grouping",
                                       choices = list("quality", "prep", "method"), 
                                       selected = "quality"),
                          numericInput("num", label = "type lab# for highlighting", value = 0),
                          radioButtons("assigned", label = "select estimator for Horwitz function",
                                       choices = list("median", "robust", "Rmode median", "Rmode mean","Mike's mode", "Mike's mode boot"),
                                       selected = "robust"),
                          checkboxGroupInput("checkGroup", label = "select estimators", 
                                             choices = list("median" = 1, "robust" = 2, "Rmode median" = 3, "Rmode mean" = 4, "Mike's mode" = 5, "Mike's mode boot" = 6)
                                             #,selected = c(1,2,3,4)
                          ),
                          radioButtons("radio", label = "bootstrap resampling values",
                                       choices = list("100" = 100, "300" = 300, "500" = 500, "1000" = 1000), 
                                       selected = 100)
                        ),
                        mainPanel(
                          h4("sigmoidal chart plot - click to remove single values"),
                          plotOutput('plot2', height = "600px",
                                     click = "plot1_click"),
                          
                          h4("estimators"),
                          tableOutput("view1"),      
                          h4("estimators uncertainties"),
                          tableOutput("view2")
                        )
                      )
             ),
             tabPanel("Bar",
                      h4("bar chart plot - click to remove single values"),
                      plotOutput('plot3', height = "600px",
                                 click = "plot1_click")
             ),
             tabPanel("Violin",
                      h4("Violin plot - outlier to be removed in bar or sigmoidal plots"),
                      plotOutput('plot4', height = "600px"
                                 ),
                      plotOutput('plot5', height = "600px"
                                 )
             )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  selectedData <- reactive({GeoPT[round == input$round.select & analyte == input$xcol, measurand, by = .(lab, analyte, quality, round, method, prep, MorT, unit)][order(round, measurand)]})
  for.rows <- reactive({GeoPT[round == input$round.select & analyte == input$xcol, .(lab)]})
  
  # doing the calcuations with function estimator.all
  
  estimators <- eventReactive(input$go, {estimator.all(input$round.select, input$xcol, input$radio, "lientz", 0.5)}) # starting function estimator.all and passing on information
  
  # needed for point removal
  
  vals <- reactiveValues(keeprows = rep(TRUE, length(for.rows))) #needed for point removal
  ranges <- reactiveValues(x = NULL, y = NULL) #needed for point removal
  
  # assigning estimators for plot
  
  estimator.plot <- reactive({
    switch(input$assigned,
           "median" = estimators()[[3]],
           "robust" = estimators()[[11]],
           "Rmode median" = estimators()[[8]],
           "Rmode mean" = estimators()[[7]],
           "Mike's mode" = estimators()[[5]],
           "Mike's mode boot" = estimators()[[14]]
    )
  })
  
  
  #all.data <- list(anal, nx, med.sig, median_se, # 1,2,3,4
  #                 Mode, mode_lientz_data, mode_lientz_mean,  #5,6,7
  #                 mode_lientz_med, mode_lientz_sd, mode_lientz_se, #8,9,10
  #                 robhub, robhubs, robhubse, Mike.mode[Mike.mode1], [2]) #11,12,13,14,15
  
  output$view1 <- renderTable({
    
    tab.est <- cbind(estimators()[[1]], estimators()[[2]], estimators()[[3]], 
                     estimators()[[11]], estimators()[[7]], estimators()[[8]], 
                     estimators()[[5]], estimators()[[14]])
    
    
    tab.est <- as.data.frame(tab.est)
    names(tab.est) <- c("analyte", "n", "median", "robust", "Rmode mean", "Rmode median", "Mike's mode", "Mike's mode boot")
    tab.est
    
  })
  
  output$view2 <- renderTable({
    
    tab.u <- cbind(estimators()[[1]], estimators()[[2]], estimators()[[4]], 
                   estimators()[[13]], estimators()[[9]], estimators()[[10]], 
                   estimators()[[15]])
    
    tab.u <- as.data.frame(tab.u)
    names(tab.u) <- c("analyte", "n", "median.se", "robust.se", "Rmode.sd", "Rmode.se", "Mike.se")
    tab.u})
  
  
  output$view2 <- renderTable({
    
    tab.u <- cbind(estimators()[[1]], estimators()[[2]], estimators()[[4]], 
                   estimators()[[13]], estimators()[[9]], estimators()[[10]], 
                   estimators()[[15]])
    
    tab.u <- as.data.frame(tab.u)
    names(tab.u) <- c("analyte", "n", "median.se", "robust.se", "Rmode.sd", "Rmode.se", "Mike.se")
    tab.u})
  
  output$plot2 <- renderPlot({                                       ### sigmoidal plot
    #all.data <- list(anal, nx, med.sig, median_se, # 1,2,3,4
    #                 Mode, mode_lientz_data, mode_lientz_mean,  #5,6,7
    #                 mode_lientz_med, mode_lientz_sd, mode_lientz_se, #8,9,10
    #                 robhub, robhubs, robhubse, Mike.mode[1], Mike.mode[2]) #11,12,13,14,15
    
    #keep    <<- selectedData()[ vals$keeprows, , drop = FALSE]
    #exclude <<- selectedData()[!vals$keeprows, , drop = FALSE]
    
    
    #xx <- keep$measurand # extracting the data only
    
    
    keep <- selectedData()[ vals$keeprows, , drop = FALSE]
    exclude <- selectedData()[!vals$keeprows, , drop = FALSE]
    
    df <- as.data.table(keep)
    df <- as.data.table(df)
    
    df$rank <- 1:length(df$measurand)

    #all.data <- list(anal, nx, med.sig, median_se, # 1,2,3,4
    #                 Mode, mode_lientz_data, mode_lientz_mean,  #5,6,7
    #                 mode_lientz_med, mode_lientz_sd, mode_lientz_se, #8,9,10
    #                 robhub, robhubs, robhubse, Mike.mode[1], Mike.mode[2]) #11,12,13,14,15
    
    x <- estimator.plot() # estimator for assigned value as basis for horwitz tramline
    
    ## horwitz value calculation
    MorT <- unique(df$MorT)
    Horwitz <- ifelse(MorT == "T", 1*0.01*(x/100/10000)^0.8495*10000, 0.01*(x/100)^0.8495)
    horwitz <- Horwitz*100 #getting the horwitz s for plotting range
    
    # lines showing up in plot
    #choices = list("median" = 1, "robust" = 2, "Rmode median" = 3, "Rmode mean" = 4, "Mike's mode" = 5, "Mike's mode boot" = 6)
    median <- ifelse("1" %in% input$checkGroup, estimators()[[3]], x)
    robhub <- ifelse("2" %in% input$checkGroup, estimators()[[11]], x)
    mode_lientz_med <- ifelse("3" %in% input$checkGroup, estimators()[[8]], x)
    mode_lientz_mean <- ifelse("4" %in% input$checkGroup, estimators()[[7]], x)
    Mode <- ifelse("5" %in% input$checkGroup, estimators()[[5]], x) # Mike without bootstrapping
    Mode_boot <- ifelse("6" %in% input$checkGroup, estimators()[[14]], x) # Mike's mode with bootstrapping
    
    
    # getting the y-axis label
    unit <- unique(selectedData()$unit)
    
    method <- df$method
    
    # needs to be selected before ggplot is executed. "select" is added to ggplot. Point removal does not work with print("plot)
    
    select <- switch(input$sig_colour,
                     "method" = geom_point(aes(fill=method), size=6, alpha=1, shape=21) ,
                     "prep" = geom_point(aes(fill=prep), size=6, alpha=1, shape=21) , 
                     "quality" = geom_point(size=6, alpha=1, aes(fill = quality), shape=21)
    )
    
    select.col <- switch(input$sig_colour,
                     "method" = scale_fill_manual(values=c(CSAN="#fdae61",
                                                           CVAAS="gray80",
                                                           EDXRF="#d73027",
                                                           ETAAS="#ffffbf",
                                                           FLAAS="#66bd63",
                                                           GFAAS="#56B4E9",
                                                           GRAV="#f46d43",
                                                           GRAV_VOL="black",
                                                           HYAAS ="#d9ef8b",
                                                           ICP="#1a9850", 
                                                           ICP_MS="#006837",
                                                           INAA="#dfc27d",
                                                           IRS="#a6d96a",
                                                           ISE="darkorchid2",
                                                           other="#80cdc1",
                                                           SPHOT="#8073ac",
                                                           UVS="gray50",
                                                           VOL="#66bd63",
                                                           WDXRF="#a50026")) ,
                     "prep" = scale_fill_manual(values=c(AD="#a50026",
                                                         AD_FM="#abd9e9",
                                                         CB="#fee090",
                                                         FA="black",
                                                         FD="#313695", 
                                                         FM="#4575b4",
                                                         FM_AD="#74add1",
                                                         NO="#dfc27d",
                                                         other="#80cdc1",
                                                         PF="gray50",
                                                         PP="#d73027",
                                                         PY="gray50",
                                                         SD="#f46d43",
                                                         SI="#fdae61")),
                     "quality" = scale_fill_manual(values=c(Applied ="mediumpurple",
                                                            Pure = "darkorange"))
    )
    
    # finding the position for lab highlighting
    
    round.ini.let <- substr(df$lab[1], 1,1) # getting the inital letter assigned to this round
    lab <- paste(round.ini.let, input$num, sep="") # pasting the initial letter with the lab number
    
    position <- df$rank[df$lab == lab] # finding the position for the plot
    
    if(lab %in% df$lab){
    rect <- geom_rect(xmin = position - 0.5, xmax = position +0.5, alpha = 0.9, ymin = -Inf, ymax = Inf, fill = "gray90")
    } else {
    rect <-   NULL
    }
    
    
    
    ggplot(df, aes(reorder(x=lab, rank), measurand)) +
      rect + # adding rectangle for selected lab (input)
      geom_point(size=1, alpha=.1) + # need to have this line for "rect" otherwise error msg
      annotate("rect", xmin=-Inf, xmax=+Inf, ymin= x - 4*horwitz, ymax=x + 4*horwitz, fill="gray", alpha=0.2) + # gray box for horwitz range
      geom_hline(yintercept = median, colour="darkgreen") + # median
      geom_hline(yintercept = robhub, linetype="dotdash", colour="khaki3") + #robhub
      geom_hline(yintercept = mode_lientz_mean, colour="skyblue3") + #Lientz mean
      geom_hline(yintercept = mode_lientz_med, colour="skyblue3", linetype = "dashed") +  # Lientz median      
      geom_hline(yintercept = Mode_boot, linetype="dashed", colour="tomato1") + # Mike mode boot
      geom_hline(yintercept = Mode, linetype="solid", colour="tomato1") + # Mike mode
      geom_hline(yintercept =  x + 2*horwitz, colour = "gray20", linetype="solid") + # horwitz
      geom_hline(yintercept =  x - 2*horwitz, colour="gray20", linetype="solid") + # horwitz
      geom_hline(yintercept =  x + 4*horwitz, colour = "gray50", linetype="dashed") + # horwitz*2
      geom_hline(yintercept =  x - 4*horwitz, colour="gray50", linetype="dashed") + # horwitz*2
      annotate("text", x=-Inf, y=Inf , hjust= -0.2, vjust= 2, label = 'solid = Lientz mean', size=4, colour="skyblue3") +
      annotate("text", x=-Inf, y=Inf , hjust= -0.16, vjust= 3.5, label = 'dashed = Lientz median', size=4, colour="skyblue3") +
      annotate("text", x=-Inf, y=Inf, hjust= -0.32, vjust= 5, label = 'solid = Mike', size = 4, colour="tomato1" ) +
      annotate("text", x=-Inf, y=Inf, hjust= -0.2, vjust= 6.5, label = 'dashed = Mike boot', size = 4, colour="tomato1" ) +
      #annotate("text", x=-Inf, y=Inf, hjust= -0.2, vjust= 3.5, label = 'dotted = se Lientz ', size = 4, colour="skyblue3") +
      annotate("text", x=-Inf, y=Inf, hjust= -0.28, vjust= 8, label = 'solid = median', size = 4, colour="darkgreen" ) +
      annotate("text", x=-Inf, y=Inf, hjust= -0.24, vjust= 9.5, label = 'dotdash = huber', size = 4, colour="khaki3" ) +
      #annotate("text", x=Inf, y=-Inf , hjust= +2, vjust= -3, label = 'Pure', size=4, colour="darkorange") +
      #annotate("text", x=Inf, y=-Inf , hjust= +1.5, vjust= -2, label = 'Applied', size=4, colour="mediumpurple") +
      annotate("text", x=-Inf, y=x+4*horwitz, hjust= -0.2, vjust= -0.3, label = "z'= 2", size=3) +
      annotate("text", x=-Inf, y=x+2*horwitz, hjust= -0.2, vjust= -0.35, label = "z = 2", size=3) +
      annotate("text", x=Inf, y=x-2*horwitz, hjust= +1.25, vjust= -0.3, label = "z'= -2", size=3) +
      annotate("text", x=Inf, y=x-4*horwitz, hjust= +1.2, vjust= -0.4, label = "z = -2", size=3) +
      xlab("laboratory code") +
      ylab(unit) +
      theme_bw() +
      select + # the selection of the groupings (quality, prep, method)
      select.col + 
      theme(legend.position = c(1,0), legend.justification = c(1,-0.05)) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.4, size = rel(1.4))) 


    
  })
  
  output$plot3 <- renderPlot({                                      # barchart plot
    
    
    keep <- selectedData()[ vals$keeprows, , drop = FALSE]
    exclude <- selectedData()[!vals$keeprows, , drop = FALSE]
    
    df <- as.data.table(keep)
    df <- as.data.table(df)
    
    df$rank <- 1:length(df$measurand)
    
    x <- estimator.plot() # estimator for assigned value as basis for horwitz tramline
    
    ## horwitz value calculation
    MorT <- unique(df$MorT)
    Horwitz <- ifelse(MorT == "T", 1*0.01*(x/100/10000)^0.8495*10000, 0.01*(x/100)^0.8495)
    horwitz <- Horwitz*100 #getting the horwitz s for plotting range
    
    ## horwitz value calculation
    df$hor[df$MorT == "T"] <- 1*0.01*(x/100/10000)^0.8495*10000
    df$hor[df$MorT == "M"] <- 0.01*(x/100)^0.8495
    df$z <- (df$measurand-x)/horwitz
    
    ## range testing
    df$colour <- "gray40"
    df$colour[df$z > 8] <- "gray80" # changing colour in plot if out of range
    df$colour[df$z < - 8] <- "gray80" # changing colour in plot if out of range
    #df$colour[df$lab == lab] <- "red"
    #df$z[df$z > 8] <- 8 # setting maximum to 8 if out of range
    #df$z[df$z < - 8] <- -8 # setting minimum to 8 if out of range
    #####
    
    #all.data <- list(anal, nx, med.sig, median_se, # 1,2,3,4
    #                 Mode, mode_lientz_data, mode_lientz_mean,  #5,6,7
    #                 mode_lientz_med, mode_lientz_sd, mode_lientz_se, #8,9,10
    #                 robhub, robhubs, robhubse, Mike.mode[1], Mike.mode[2]) #11,12,13,14,15
    
    ### assigning estimator to variable
    
    robhub <- estimators()[[11]] # robust huber
    med <-  estimators()[[3]] # median
    mode_lientz_boot <- estimators()[[8]] # median Lientz
    Mode <- estimators()[[14]] # Mike's mode boot
    
    ## calculations for bar plot
    
    robhub.z <- (robhub-x)/horwitz
    #robhub.z <- ifelse(robhub.z < -8 | robhub.z > 8,8.5, robhub.z) # not showing line when outside the -8 and +8 z limits
    med.z <- (med-x)/horwitz
    #med.z <- ifelse(med.z < -8 | med.z > 8,8.5, med.z) # not showing line when outside the -8 and +8 z limits
    mod.z <- (mode_lientz_boot-x)/horwitz
    #mod.z <- ifelse(mod.z < -8 | mod.z > 8,8.5, mod.z) # not showing line when outside the -8 and +8 z limits
    #mod.se.z <- (Semode_lientz_boot)/horwitz
    #mod.z <- ifelse(mod.z < -8 | mod.z > 8,8.5, mod.z) # not showing line when outside the -8 and +8 z limits
    #mod.se.z <- ifelse(mod.se.z < -8 | mod.se.z > 8,8.5, mod.se.z) # not showing line when outside the -8 and +8 z limits
    mike.mod.z <- (Mode-x)/horwitz
    #mike.mod.z <- ifelse(mike.mod.z < -8 | mike.mod.z > 8,8.5, mike.mod.z) # not showing line when outside the -8 and +8 z limits
    #mike.mod.z <- ifelse(is.na(mike.mod.z), mike.mod.z <- 0, mike.mod.z) # not showing line when outside the -8 and +8 z limits # need additional test as unplottable mike.mod.ze == "NA_real_" occurs otherwise
    
    colour <- df$colour # getting the colours for barcharts
    
    # getting the y-axis label
    unit <- unique(selectedData()$unit)
    
    # lines showing up in plot
    #choices = list("median" = 1, "robust" = 2, "Rmode median" = 3, "Rmode mean" = 4, "Mike's mode" = 5, "Mike's mode boot" = 6)
    med.z <- ifelse("1" %in% input$checkGroup, med.z, 0)
    robhub.z <- ifelse("2" %in% input$checkGroup, robhub.z, 0)
    mod.z <- ifelse("3" %in% input$checkGroup,  mod.z, 0)
    #mode_lientz_mean <- ifelse("4" %in% input$checkGroup, estimators()[[7]], 0)
    #Mode <- ifelse("5" %in% input$checkGroup, estimators()[[5]], 0) # Mike without bootstrapping
    mike.mod.z <- ifelse("6" %in% input$checkGroup, mike.mod.z, 0) # Mike's mode with bootstrapping
    
    
    ggplot(df, aes(reorder(x=lab, measurand), z)) +

      geom_hline(yintercept = med.z, colour="aquamarine3") + # median
      
      geom_hline(yintercept =  + 2, colour = "gray20", linetype="solid") + # horwitz "tramline"
      geom_hline(yintercept =  - 2, colour="gray20", linetype="solid") + # horwitz
      geom_hline(yintercept =  + 4, colour = "gray70", linetype="solid") + # horwitz*2 "tramline for applied labs"
      geom_hline(yintercept =  - 4, colour="gray70", linetype="solid") + # horwitz*2
      geom_bar(stat="identity", position="identity", width=0.6, colour=colour, fill=colour) +
      annotate("rect", xmin=-Inf, xmax=+Inf, ymin=  - 4, ymax= + 4, fill="gray", alpha=0.2) + # gray box for horwitz range
      
      ylab("implied z-score") +
      xlab("laboratory code") +
      geom_hline(yintercept = robhub.z, linetype="dotdash", colour="khaki3") + # robust estimator
      geom_hline(yintercept = mike.mod.z, linetype="dashed", colour="tomato1") + # mike's mode
      geom_hline(yintercept = mod.z,linetype="dotdash", colour="skyblue3") + # lientz_mode
      geom_hline(yintercept = med.z, linetype="dotted", colour="aquamarine3") +
      #ylim(-8, 8) + # restickting the range
      annotate("text", x=0, y=4, hjust= -0.2, vjust= -0.3, label = "z'= 2", size=3) +
      annotate("text", x=0, y=2, hjust= -0.2, vjust= -0.35, label = "z = 2", size=3) +
      annotate("text", x=Inf, y=-4, hjust= +1.2, vjust= -0.3, label = "z'= -2", size=3) +
      annotate("text", x=Inf, y=-2, hjust= +1.2, vjust= -0.4, label = "z = -2", size=3) +
      geom_hline(yintercept = + 8, linetype="dotted", colour="gray70") +
      geom_hline(yintercept = - 8, linetype="dotted", colour="gray70") +
      scale_y_continuous(breaks = c(-8, -4,-2,2,4,8)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.4, size = rel(.9)))
    
  })
  
  output$plot4 <- renderPlot({                                                  #### violin method plot

    keep <- selectedData()[ vals$keeprows, , drop = FALSE]
    exclude <- selectedData()[!vals$keeprows, , drop = FALSE]
    
    df <- as.data.table(keep)
    df <- as.data.table(df)
    #df <- na.omit(df, cols = "method")
    
    df$rank <- 1:length(df$measurand)

    x <- estimator.plot() # estimator for assigned value as basis for horwitz tramline
    
    ## horwitz value calculation
    MorT <- unique(df$MorT)
    Horwitz <- ifelse(MorT == "T", 1*0.01*(x/100/10000)^0.8495*10000, 0.01*(x/100)^0.8495)
    horwitz <- Horwitz*100 #getting the horwitz s for plotting range
    
    # lines showing up in plot
    #choices = list("median" = 1, "robust" = 2, "Rmode median" = 3, "Rmode mean" = 4, "Mike's mode" = 5, "Mike's mode boot" = 6)
    
    # getting the y-axis label
    unit <- unique(selectedData()$unit)

    median <- ifelse("1" %in% input$checkGroup, estimators()[[3]], x)
    
    
    ggplot(reorder_by(method, ~measurand, df, median), aes(method, measurand)) +
      geom_hline(yintercept =  x + 2*horwitz, colour = "gray20", linetype="solid") + # horwitz
      geom_hline(yintercept =  x -  2*horwitz, colour="gray20", linetype="solid") + # horwitz
      geom_hline(yintercept =  x + 4*horwitz, colour = "gray70", linetype="solid") + # horwitz*2
      geom_hline(yintercept =  x - 4*horwitz, colour="gray70", linetype="solid") + # horwitz*2
      geom_hline(yintercept = median, colour="aquamarine3") + # median
      geom_violin(scale = "count") +
      geom_jitter(aes(fill=prep), size = 5, alpha = .8, position = position_jitter(width = .15), shape=21) +
      stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=3.5) +

      ylab(unit) + 
      labs(colour = "prep") +
      scale_fill_manual(values=c(AD="#a50026",
                                 AD_FM="#abd9e9",
                                 CB="#fee090",
                                 FA="black",
                                 FD="#313695", 
                                 FM="#4575b4",
                                 FM_AD="#74add1",
                                 NO="#dfc27d",
                                 other="#80cdc1",
                                 PF="gray50",
                                 PP="#d73027",
                                 PY="gray50",
                                 SD="#f46d43",
                                 SI="#fdae61")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = rel(1.2))) +
      theme(legend.text=element_text(size=rel(1.2))) +
      theme(legend.position="top")
    
  })
  
  output$plot5 <- renderPlot({                                                  #### violin prep plot
    #all.data <- list(anal, nx, med.sig, median_se, # 1,2,3,4
    #                 Mode, mode_lientz_data, mode_lientz_mean,  #5,6,7
    #                 mode_lientz_med, mode_lientz_sd, mode_lientz_se, #8,9,10
    #                 robhub, robhubs, robhubse, Mike.mode[1], Mike.mode[2]) #11,12,13,14,15
    
    #keep    <<- selectedData()[ vals$keeprows, , drop = FALSE]
    #exclude <<- selectedData()[!vals$keeprows, , drop = FALSE]
    
    
    #xx <- keep$measurand # extracting the data only
    
    
    keep <- selectedData()[ vals$keeprows, , drop = FALSE]
    exclude <- selectedData()[!vals$keeprows, , drop = FALSE]
    
    df <- as.data.table(keep)
    df <- as.data.table(df)
    
    df$rank <- 1:length(df$measurand)
    
    
    #all.data <- list(anal, nx, med.sig, median_se, # 1,2,3,4
    #                 Mode, mode_lientz_data, mode_lientz_mean,  #5,6,7
    #                 mode_lientz_med, mode_lientz_sd, mode_lientz_se, #8,9,10
    #                 robhub, robhubs, robhubse, Mike.mode[1], Mike.mode[2]) #11,12,13,14,15
    
    x <- estimator.plot() # estimator for assigned value as basis for horwitz tramline
    
    ## horwitz value calculation
    MorT <- unique(df$MorT)
    Horwitz <- ifelse(MorT == "T", 1*0.01*(x/100/10000)^0.8495*10000, 0.01*(x/100)^0.8495)
    horwitz <- Horwitz*100 #getting the horwitz s for plotting range
    
    # lines showing up in plot
    #choices = list("median" = 1, "robust" = 2, "Rmode median" = 3, "Rmode mean" = 4, "Mike's mode" = 5, "Mike's mode boot" = 6)
    
    median <- ifelse("1" %in% input$checkGroup, estimators()[[3]], x)
    robhub <- ifelse("2" %in% input$checkGroup, estimators()[[11]], x)
    mode_lientz_med <- ifelse("3" %in% input$checkGroup, estimators()[[8]], x)
    mode_lientz_mean <- ifelse("4" %in% input$checkGroup, estimators()[[7]], x)
    Mode <- ifelse("5" %in% input$checkGroup, estimators()[[5]], x) # Mike without bootstrapping
    Mode_boot <- ifelse("6" %in% input$checkGroup, estimators()[[14]], x) # Mike's mode with bootstrapping
    
    
    # getting the y-axis label
    unit <- unique(selectedData()$unit)
    
    ggplot(reorder_by(method, ~measurand, df, median), aes(prep, measurand)) +
      geom_hline(yintercept =  x + 2*horwitz, colour = "gray20", linetype="solid") + # horwitz
      geom_hline(yintercept =  x -  2*horwitz, colour="gray20", linetype="solid") + # horwitz
      geom_hline(yintercept =  x + 4*horwitz, colour = "gray70", linetype="solid") + # horwitz*2
      geom_hline(yintercept =  x - 4*horwitz, colour="gray70", linetype="solid") + # horwitz*2
      geom_hline(yintercept = median, colour="aquamarine3") + # median
      geom_violin(scale = "count") +
      geom_jitter(aes(fill=method), size = 5, alpha = .8, position = position_jitter(width = .15), shape=21) +
      stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=3.5) +
      
      ylab(unit) + 
      labs(colour = "prep") +
      scale_fill_manual(values=c(CSAN="#fdae61",
                                 CVAAS="gray80",
                                 EDXRF="#d73027",
                                 ETAAS="#ffffbf",
                                 FLAAS="#66bd63",
                                 GFAAS="#56B4E9",
                                 GRAV="#f46d43",
                                 GRAV_VOL="black",
                                 HYAAS ="#d9ef8b",
                                 ICP="#1a9850", 
                                 ICP_MS="#006837",
                                 INAA="#dfc27d",
                                 IRS="#a6d96a",
                                 ISE="darkorchid2",
                                 other="#80cdc1",
                                 SPHOT="#8073ac",
                                 UVS="gray50",
                                 VOL="#66bd63",
                                 WDXRF="#a50026")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = rel(1.2))) +
      theme(legend.text=element_text(size=rel(1.2))) +
      theme(legend.position="top")
    
  })
  
  observeEvent(input$plot1_click, {
    x <- round(input$plot1_click$x)
    df <- selectedData()
    df$rank <- rep(1:length(df$measurand))
    res <- round(input$plot1_click$x) == df$rank
    vals$keeprows <- xor(vals$keeprows, res)
  })
  
  
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  # Reset all points
  observeEvent(input$exclude_reset, {
    vals$keeprows <-  rep(TRUE, length(for.rows))
  })
  
  output$out_table <- DT::renderDataTable({
    DT::datatable(selectedData(), options = list(lengthMenu = c(5, 30, 50), pageLength = 15))
    #DT::datatable(selectedData(), options = list(orderClasses = TRUE))
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)


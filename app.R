###############
# Autor: Valentin Koob
# Description: Shows the influence of assumption violations 
# on the t-test
###############





######### Ab hier nur noch Berechnungen und Plots #########
plotPopulation = function(sampleInput, 
                          sdInput1,
                          sdInput2){
  
  
  # find the min and max values of the population
  if (sampleInput == "norm") {
    xmin = floor(-max(sdInput1, sdInput2) * 4)
    xmax = ceiling(-xmin)
  } else if (sampleInput == "beta") { # varianz bei alpha = 2 = beta -> 0.05
    scale_1 = sqrt(sdInput1^2/0.05)
    scale_2 = sqrt(sdInput2^2/0.05)
    xmax = ceiling(max(scale_1, scale_2)/2) # da spaeter auf 0 zentriert 
    xmin = -xmax
  } else if (sampleInput == "unif") { # Var(unif) = 1/12 [0,1]
    scale_1 = sqrt(sdInput1^2/(1/12))
    scale_2 = sqrt(sdInput2^2/(1/12))
    xmax = ceiling(max(scale_1, scale_2)/2) # da spaeter auf 0 zentriert 
    xmin = -xmax
  }
  
  # calculate density distributions
  xs = seq(xmin, xmax, length.out = 10000)
  if (sampleInput == "norm"){
    dichte1 = dnorm(xs, mean = 0, sd = sdInput1)
    dichte2 = dnorm(xs, mean = 0, sd = sdInput2)
  } else if (sampleInput == "beta"){
    dichte1 = dbeta(seq(0, 1, length.out = 10000), 2, 2)/(scale_1)
    dichte2 = dbeta(seq(0, 1, length.out = 10000), 2, 2)/(scale_2)
  } else if(sampleInput == "unif"){
    dichte1 = dunif(xs, min = -scale_1/2, max = scale_1/2)
    dichte2 = dunif(xs, min = -scale_2/2, max = scale_2/2)
  }
  
  # plot the distributions
  if(sampleInput == "norm" | sampleInput == "unif"){
    plot(c(1,2) ~ c(1,1), xlim = c(xmin, xmax), ty = "l", col = "white",
         xlab = "moegliche Werte x", 
         ylab = "f(x)", ylim = c(0, max(c(dichte1, dichte2)) + max(c(dichte1, dichte2))*0.35),
         cex.lab = 1.25, cex.axis = 1.25)
    polygon(x = xs, y = dichte1, col = rgb(0,1,0, 0.2))
    polygon(x = xs, y = dichte2, col = rgb(1,0,0, 0.2))
  }else if(sampleInput == "beta"){
    plot(c(1,2) ~ c(1,1), xlim = c(xmin, xmax), ty = "l", col = "white",
         xlab = "moegliche Werte x", 
         ylab = "f(x)", ylim = c(0, max(c(dichte1, dichte2)) + max(c(dichte1, dichte2))*0.35))
    polygon(x = seq(0-scale_1/2, scale_1/2, length.out = 10000), 
            y = dichte1, col = rgb(0,1,0, 0.2))
    polygon(x = seq(0-scale_2/2, scale_2/2, length.out = 10000), 
            y = dichte2, col = rgb(1,0,0, 0.2))
    
  }
  
  legend("topright", legend = c("Population 1", "Population 2"), col = c(rgb(0,1,0, 0.8), rgb(1,0,0, 0.8)), pch = 15, 
         cex = 1.25)
  
  
}

plotSimulation = function(sampleInput, 
                          sdInput1,
                          sdInput2, 
                          n,
                          nSim, 
                          buttonState){
  
  if(!sampleInput %in% c("norm", "beta", "unif"))
    stop("NOT DEFINED DISTRIBUTION")
  

  buttonState = as.numeric(buttonState)
  
  # if buttonState > 0  run the simulation
  if (buttonState > 0) {
    
    # place holder while the user is waiting
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("Bitte warten :)"), 
         cex = 1.6, col = "black")
    
    # run the simulation
    res = 
    lapply(1:nSim, function(oneRep){
      # get the samples in one repetition
      if (sampleInput == "norm") {
        sample1 = rnorm(n, mean = 0, sd = sdInput1)
        sample2 = rnorm(n, mean = 0, sd = sdInput2)
      } else if (sampleInput == "beta") {
        scale_1 = sqrt(sdInput1 ^ 2 / 0.05) # Streckung der Daten
        scale_2 = sqrt(sdInput2 ^ 2 / 0.05) # Streckung der Daten
        sample1 = (rbeta(n, 2, 2) - 0.5) * scale_1
        sample2 = (rbeta(n, 2, 2) - 0.5) * scale_2
      } else if (sampleInput == "unif") {
        scale_1 = sqrt(sdInput1 ^ 2 / (1 / 12)) # obere Grenze finden
        scale_2 = sqrt(sdInput2 ^ 2 / (1 / 12)) # obere Grenze finden
        sample1 = runif(n, min = -scale_1 / 2, max = scale_1 / 2)
        sample2 = runif(n, min = -scale_2 / 2, max = scale_2 / 2)
      }    
      t_obj <- t.test(x = sample1,
                      y = sample2,
                      alternative =  "two.sided",
                      var.equal = TRUE) # t test
      t <- t_obj$statistic # t wert
      p <- t_obj$p.value # p wert
      return(matrix(c(t, p), nrow = 1, ncol = 2)) # Rueckgabe
    })
    res = do.call("rbind", res)
    
    #### plotten der theoretischen t-Verteilung 
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    dfs = (2*n) - 2
    xmin = -5
    xmax = 5
    xs = seq(xmin, xmax, length.out = 10000)
    dichteT = dt(xs, df = dfs)
    # hist(res[,1], xlim = c(xmin, xmax), prob = TRUE, breaks = 100,
    #      xlab = "moegliche t-Werte", 
    #      ylab = "f(t)", ylim = c(0, max(dichteT) + max(dichteT)*0.3), 
    #      main = "") 
    plot(c(1,2) ~ c(1,1), col = "white", xlim = c(xmin, xmax),
         xlab = "moegliche t-Werte", 
         ylab = "f(t)", ylim = c(0, max(dichteT) + max(dichteT)*0.35), 
         main = "", cex.lab = 1.25, cex.axis = 1.25) 
    dens = density(res[,1])
    polygon(x = dens$x, y = dens$y, col = "skyblue")
    points(dichteT ~ xs, col = "red", ty = "l", lwd = 2)
    legend("topright", legend = c(paste("theoretische t-Verteilung\ndf =", dfs)
                                        , "simulierte t-Werte"), col = c("red", "skyblue"), 
           lty = 1, cex = 1.25)
    
    nomAlpha = mean(res[,2] <= 0.05)
    mtext(side = 3, paste("% Fehlentscheidungen\n", round(nomAlpha, 4)*100), cex = 1.5, line = 0.5)
    t_krit = qt(0.975, dfs)
    segments(x0 = -t_krit, y0 = 0, x1 = -t_krit, y1 = 0.1)
    segments(x0 = t_krit, y0 = 0, x1 = t_krit, y1 = 0.1)
    text(x = -t_krit, y = 0.15, expression(t[krit]), cex = 1.25)
    text(x = t_krit, y = 0.15, expression(t[krit]), cex = 1.25)
    #text(x = 0, y = 0.45, paste("t(", dfs, ")", sep = "" ), cex = 1.25, col = "red")
    
  }else{
  
  
    # place holder while the user is waiting
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("Bitte Simulation starten :)"), 
         cex = 1.6, col = "black")
 
  } 
}



################## Shiny Implementierung


ui <- fluidPage(
  titlePanel("Auswirkung der Annahmensverletzungen auf das nominelle Alpha-Niveau (beim t-Test fuer 
             unabhaengige Stichproben"),
  sidebarLayout(
    
    sidebarPanel(
      
      sliderInput(
        "SdInput1",
        h4("Standardabweichung der 1. Population"),
        min = 5, max = 50, value = 5, step = 5),
      sliderInput(
        "SdInput2",
        h4("Standardabweichung der 2. Population"),
        min = 5, max = 50, value = 5, step = 5),
      sliderInput(
        "nsInput",
        h4("Stichprobengroeße"),
        min = 3, max = 15, value = 7, step = 1),
      selectInput(
        "SampleInput",
        h4("Woraus soll gezogen werden?"),
        list(Normalverteilung = "norm", 
             `skalierte Beta-Verteilung` = "beta", 
             Gleichverteilung = "unif"),
        "Normalverteilung"),
      sliderInput(
        "nSimInput",
        h4("Wie oft soll gezogen werden? (Dauer der Simulation)"),
        min = 5000, max = 15000, value = 5000, step = 5000),
      h4("Simulation starten"),
      actionButton(
        "buttonInput", 
        "go!")
      

    ),
    
    mainPanel(fluidPage(
      
      h4("Veranschaulichung der Populationen aus denen gezogen wird (unter der H0)"),
      fluidRow(plotOutput(outputId = "population")),
      
      h4("Theoretische vs. simulierte Verteilung der t-Werte"),
      fluidRow(plotOutput(outputId = "sim")),
      
    ))
  )
)




server <- function(input, output) {
  output$population <- renderPlot({
    plotPopulation(input$SampleInput, 
                   input$SdInput1,
                   input$SdInput2)
  })
  
  output$sim <- renderPlot({
    plotSimulation(input$SampleInput, 
                   input$SdInput1,
                   input$SdInput2, 
                   input$nsInput,
                   input$nSimInput, 
                   input$buttonInput)
  })
  
}

shinyApp(ui = ui, server = server)

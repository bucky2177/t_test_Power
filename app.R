###############
# Autor: Valentin Koob
###############

list.of.packages <- c("colorspace")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(colorspace)

######### Ab hier nur noch Berechnungen und Plots #########

plotMeans = function(means, sigma.2) {
  means = unlist(strsplit(means, "[ ]"))
  means = as.numeric(means)
  groups = 1:length(means)
  sigma.2 = as.numeric(sigma.2)
  
  plot(
    means ~ groups,
    ylim = c(min(means) - 2 - sqrt(sigma.2), max(means) + 2 + sqrt(sigma.2)),
    type = "b",
    axes = FALSE,
    ylab = "Mittelwerte",
    xlab = "Gruppen",
    pch = 19
  )
  
  axis(side = 1,
       at = groups)
  axis(
    side = 2,
    las = 2,
    at = seq(
      from = min(means) - 4 - ceil(sqrt(sigma.2)),
      to = max(means) + 4 + ceil(sqrt(sigma.2)),
      by = 2.5
    )
  )
  errbar(
    x = groups,
    y = means,
    yplus = means + sqrt(sigma.2),
    yminus = means - sqrt(sigma.2),
    add = TRUE
  )
}

plotNormalPops = function(means, sigma.2){
  
  means = unlist(strsplit(means, "[ ]"))
  means = stringr::str_replace(means, ",", ".")
  means = as.numeric(means)
  groups = 1:length(means)
  sigma.2 = unlist(strsplit(sigma.2, "[ ]"))
  sigma.2 = stringr::str_replace(sigma.2, ",", ".")
  sigma.2 = as.numeric(sigma.2)
  if(length(sigma.2) >= 2){
    warning("Mehr als eine Varianz spezifiziert")
    return(NULL)
  }
  if(length(means) > 2){
    warning("Mehr als zwei Mittelwerte spezifiziert")
    return(NULL)
  }
  
  par(mfrow = c(1, 2))
  
  plot(c(1,2) ~ 1, 
       col = "white", 
       xlab = "Werte in den Populationen", 
       ylab = "Dichte", 
       axes = FALSE, 
       cex.lab = 1.2,
       ylim = c(0, dnorm(means[1], mean = means[1], sd = sqrt(sigma.2)))*1.8,
       xlim = c(floor(min(means) - 6 * sqrt(sigma.2)), ceiling(max(means) + 6*sqrt(sigma.2))))
  
  axis(side = 1)
  axis(side = 2)
  x = seq(floor(min(means) - 3 * sqrt(sigma.2)),  ceiling(max(means) + 3*sqrt(sigma.2)), 0.01)
  for(i in 1:length(means)){
    #segments(means[i], 0, means[i], dnorm(means[i], mean = means[i], sd = sqrt(sigma.2)), col = "gray")
    dfs = dnorm(x, mean = means[i], sd = sqrt(sigma.2))
    points(dfs~ x, ty = "l", 
           col = rainbow_hcl(10)[i], lwd = 2)
    polygon(c(min(x),x, max(x)), c(0, dfs, 0), col = adjustcolor(rainbow_hcl(10)[i], alpha.f=0.2), border = NA)
    
  }
  
}



plotDensities = function(ns, means, sigma.2, alpha) {
  alpha = unlist(strsplit(alpha, "[ ]"))
  alpha = as.numeric(alpha)[1]
  if(length(alpha) >= 2){
    warning("Mehr als ein alpha spezifiziert")
    return(NULL)
  }
  means = unlist(strsplit(means, "[ ]"))
  means = stringr::str_replace(means, ",", ".")
  means = as.numeric(means)
 
  ns = unlist(strsplit(ns, "[ ]"))
  if(length(ns) >= 2){
    warning("Mehr als eine Stichprobengrößen spezifiziert")
    return(NULL)
  }
  ns = as.numeric(ns)[1]
  ns = rep(ns, 2)
  
  sigma.2 = unlist(strsplit(sigma.2, "[ ]"))
  sigma.2 = stringr::str_replace(sigma.2, ",", ".")
  sigma.2 = as.numeric(sigma.2)
  if(length(sigma.2) >= 2){
    warning("Mehr als eine Varianz spezifiziert")
    return(NULL)
  }
  
  par(mfrow = c(1, 2))
  x = seq(-5, 10, 0.01)
  dfs_null = dt(x, sum(ns) - 2)
  plot(
    x,
    rep(1, length(x)),
    col = "white",
    ylim = c(0, 1),
    ylab = "Dichte",
    xlab = "t",
    cex.lab = 1.2,
    axes = FALSE,
    main = "t-Verteilung unter der H0"
  )
  axis(side = 1)
  axis(side = 2)
  
  points(x, dfs_null, ty = "l", pch = 5)
  x_quantil = qt(1 - alpha, sum(ns) - 2)
  polygon(c(0, x[x < x_quantil], x_quantil),
          c(0, dfs_null[x < x_quantil], 0),
          col = "green3",
          border = NA)
  polygon(c(x_quantil, x[x > x_quantil], max(x)),
          c(0, dfs_null[x > x_quantil], 0),
          col = "red",
          border = NA)
  
  segments(x_quantil, 0, x_quantil, .6)
  t_krit = round(x_quantil, 2)
  text(x = x_quantil, y = 0.65, bquote(t[krit] ~ "=" ~ .(t_krit)), cex = 1.5)
  text(x_quantil + 0.5, .1, expression(alpha), cex = 1.5)
  
  
  
  plot(
    x,
    rep(1, length(x)),
    col = "white",
    ylim = c(0, 1),
    ylab = "Dichte",
    xlab = "t",
    cex.lab = 1.2,
    axes = FALSE,
    main = "t-Verteilung gegeben der wahren Situation \n (= unter der H1) \n mit gestrichelter H0"
  )
  axis(side = 1)
  axis(side = 2)
  
  ncp = ((means[2] - means[1]) / sqrt(sigma.2)) * sqrt(sum(ns)/4) #ncp (vgl. power.t.test)
  dfs_h1 = dt(x, sum(ns) - 2, ncp)
  points(x,
         dfs_h1,
         ty = "l",
         pch = 5,
         col = "black")
  polygon(
    c(0, x[x < x_quantil], x_quantil),
    c(0, dfs_h1[x < x_quantil], 0),
    col = "orange",
    border = NA,
    lty = 0
  )
  polygon(
    c(x_quantil, x[x > x_quantil], max(x)),
    c(0, dfs_h1[x > x_quantil], 0),
    col = "dodgerblue",
    border = NA,
    lty = 0
  )
  segments(x_quantil, 0, x_quantil, .6)
  text(x = x_quantil, y = 0.65, bquote(t[krit] ~ "=" ~ .(t_krit)), cex = 1.5)
  points(x,
         dfs_null,
         ty = "l",
         pch = 5,
         lty = 3)
  
  
  ##Griechische Buchstaben
  text(x_quantil - .9, .1, expression(beta), cex = 1.5)
  text(x_quantil + 1, .1, expression(1 - beta), cex = 1.5)
  
  
  ##Power
  power = 1 - pt(x_quantil, sum(ns) - 2, ncp)
  mtext(paste("Power = ", round(power, 2)), side=1, line=-15, at=8, cex = 1.5)
}




################## Shiny Implementierung


ui <- fluidPage(
  titlePanel("Power beim t-Test (unabh. Stichproben)"),
  sidebarLayout(
    sidebarPanel(
      textInput(
        "MittelwerteInput",
        h3("Erwartungswerte der Populationen"),
        value = "1 1.5"
      ),
      textInput(
        "Sigma.2Input",
        h3("Varianz innerhalb der Populationen"),
        value = "1"
      ),
      textInput("nsInput", h3("n pro Gruppe"),
                value = "10"),
      textInput("alphaInput", h3("alpha (Fehler 1. Art)"),
                value = "0.05")
    ),
    
    mainPanel(fluidPage(
      h3(
        "Veranschaulichung der wahren Situation auf Populationsebene"
      ),
      fluidRow(plotOutput(outputId = "means")),
      h3("Veranschaulichung der t-Verteilungen"),
      fluidRow(plotOutput(outputId = "densities")),
      
    ))
  )
)




server <- function(input, output) {
  output$densities <- renderPlot({
    plotDensities(input$nsInput,
                  input$MittelwerteInput,
                  input$Sigma.2Input,
                  input$alphaInput)
  })
  
  output$means <- renderPlot({
    plotNormalPops(input$MittelwerteInput, input$Sigma.2Input)
  })
  
}

shinyApp(ui = ui, server = server)

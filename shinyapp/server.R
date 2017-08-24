library(mvtnorm) #2D normal density
library(fields) #image.plot for heatmap
library(RColorBrewer) # Color palette

residual.heat <- function(x, y, res=NULL, n_1=50, n_2=n_1, border=0.01, cut = 5, size=1000, sigma=diag(c((max(x)-min(x))^2/size, (max(y)-min(y))^2/size)), output = "plot", ...){
    # if res is NULL show a scatter plot of the input in x
    if(is.null(res)){
        my.cols <- rev(brewer.pal(11, "RdYlBu"))
        #return (smoothScatter(x, y, nrpoints=0, colramp=colorRampPalette(my.cols), pch=19, cex=.8, add=F, ...))
        return (smoothScatter(x, y, ...))
    }
    
    # initialize grid 
    ax1 <- seq(min(x), max(x), length.out = n_1)
    ax2 <- seq(min(y), max(y), length.out = n_2)
    grids <- expand.grid(x=ax1, y=ax2)
    
    # calculate the weights
    # first column is the first grid point with weights for all residuals on this grid point
    # first row contains all weights for the first residual for all grid points (begin top right)
    matr <- matrix(nrow = length(x), ncol = nrow(grids))
    for(i in 1:length(x)){
        mvnorm_swap <- dmvnorm(grids, mean=c(x[i],y[i]), sigma = sigma)
        mvnorm_swap[mvnorm_swap <= border] <- 0
        matr[i,] <- mvnorm_swap 
    }
    
    # calculate values for each grid point
    w_error <- rep(NaN, n_1*n_2)
    for(i in 1:ncol(matr)){
        if(sum(matr[,i]>0) >= cut) w_error[i] <- sum(res*(matr[,i]/sum(matr[,i])))
    }
    
    # arrange grid points in heat map matrix
    z <- matrix(w_error, length(ax1), length(ax2), byrow=F)
    
    # plot heat map
    b <- list(x=ax1,y=ax2,z=z)
    if(output == "image") image(b, ...)
    else if(output == "plot") image.plot(b, ...)
    else if(output == "matrix") b
    else NULL
} 

server <- function(input, output){
    
    set.seed(42)
    resid_test <- data.frame(a=rnorm(1000), b=rnorm(1000), c=rnorm(1000))
    heatmat <- residual.heat(x=resid_test$a, y=resid_test$b, res=resid_test[,3], n_1=100, cut = 5, output = "matrix")
    
    # Plot time series chart 
    output$base <- renderPlotly({
        p <- plot_ly(source = "source") %>%
            add_heatmap(x=heatmat$x, y=heatmat$y , z = heatmat$z, colorscale = "Greys", name="Heat") %>%
            add_trace(data = resid_test, x = ~a, y = ~b, type = "scatter", name="scatter")%>%
            layout(
                title = "Double Y Axis"
            )
        p
    })
    
    # Coupled hover event
    output$select <- renderPlotly({
        
        # Read in hover data
        eventdata <- event_data("plotly_selected", source = "source")
        validate(need(!is.null(eventdata), "Hover over the time series chart to populate this heatmap"))
        # Get point number
        datapoints <- as.numeric(eventdata$pointNumber+1)
        selectedforheat <- resid_test[datapoints,]
        heatmatselected <- residual.heat(x=selectedforheat$a, y=selectedforheat$b, res=selectedforheat[,3], n_1=100, cut = 5, output = "matrix")
        plot_ly(source = "source") %>%
            add_heatmap(x=heatmatselected$x, y=heatmatselected$y , z = heatmatselected$z, colorscale = "Greys", name="Heat") %>%
            add_trace(data = selectedforheat, x = ~a, y = ~b, type = "scatter", name="scatter")%>%
            layout(
                title = "Double Y Axis"
            )
    })
    
}
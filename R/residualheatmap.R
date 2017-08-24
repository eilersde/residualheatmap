#' Calculate residual heat maps

library(mvtnorm) 
library(fields) 
library(RColorBrewer) 

residual_heatmap <- function(x, y, res, n_1=50, n_2=n_1, border=0.01, cut = 5, size=1000, sigma=diag(c((max(x)-min(x))^2/size, (max(y)-min(y))^2/size)), output = "plot", ...){
    # initialize grid 
    ax1 <- seq(min(x), max(x), length.out = n_1)
    ax2 <- seq(min(y), max(y), length.out = n_2)
    grids <- expand.grid(x=ax1, y=ax2)
    # calculate the weights
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
    # return heat map matrix
    list(x=ax1,y=ax2,z=z)
} 

plot_residual_heatmap <- function(heatmap_matrix, type = "plot", pal = brewer.pal(11, "RdYlBu"), divisions = 64, balanced = F, ...){
    breaks <- seq(min(heatmap_matrix[['z']], na.rm =T),max(heatmap_matrix[['z']], na.rm =T),length.out=divisions+1)
    max_absolute_value <- max(abs(heatmap_matrix[['z']]), na.rm =T)
    if(balanced) breaks <- seq(-max_absolute_value,max_absolute_value,length.out=divisions+1)
    if(type == "image") image(heatmap_matrix, col=colorRampPalette(pal)(divisions), breaks=breaks,...)
    else if(type == "plot") image.plot(heatmap_matrix, col=colorRampPalette(pal)(divisions), breaks=breaks,...)
    else image.plot(heatmap_matrix, col=colorRampPalette(pal)(divisions), breaks=breaks,...)
}

## fig 1: Scatterplot
set.seed(42)
resid_test <- data.frame(a=rnorm(500), b=rnorm(500), c=rnorm(500))
smoothScatter(resid_test$a, 
              resid_test$b, 
              nbin=300,
              xlab = "x",
              ylab = "y",
              add=F,
              colramp = colorRampPalette(c("white", brewer.pal(9, "Greens"))))

## fig 2: 
plot_residual_heatmap(
    residual_heatmap(x=resid_test$a,
                  y=resid_test$b,
                  res=resid_test$c, 
                  n_1=100,
                  cut = 9),
    xlab = "x",
    ylab = "y",
    pal = rev(brewer.pal(11, "RdYlBu")),
    balanced = T)
    

# fig 3: Bias
resid_test[resid_test[,1] < 0.2 & resid_test[,1] > -0.2, 3] <- resid_test[resid_test[,1] < 0.2 & resid_test[,1] > -0.2, 3] + 2
plot_residual_heatmap(
    residual_heatmap(x=resid_test$a,
                     y=resid_test$b,
                     res=resid_test$c, 
                     n_1=100,
                     cut = 9),
    xlab = "x",
    ylab = "y",
    pal = rev(brewer.pal(11, "RdYlBu")),
    balanced = T)

# fig 4: linear regression
set.seed(1233)
nonli <- function(a) a^2+a + rnorm(length(a), 0, 12)
a <- runif(500,-5,10) 
y <- nonli(a)
linear_model <- lm(y ~ a)
pred_linear <- predict(linear_model)
plot_residual_heatmap(residual_heatmap(x=a,
                                       y=y,
                                       res = (y-as.vector(pred_linear)),
                                       n_1=100),
                      xlab = "x",
                      ylab = "f(x)",
                      pal = rev(brewer.pal(11, "RdYlBu")),
                      balanced = T)
points(a, y)
abline(linear_model, lwd=3)

## fig 5: heaton(2016)
ratio <- function(a,b,c,d) (a-b)/(c-d) + rnorm(length(a), 0, 1)
a <- runif(2000,1,10) 
b <- runif(2000,1,10)
c <- runif(2000,1,10)
d <- runif(2000,1,10)
y <- ratio(a,b,c,d)
training <- data.frame(y=y, a=a, b=b, c=c, d=d)
library(h2o)
h2o.init(nthreads=-1, max_mem_size="2G")
h2o.removeAll()
train_h2o <- as.h2o(training)
y_dep <- 1
x_indep <- c(2:length(train_h2o))
ann_model <- h2o.deeplearning(x = x_indep,
                              y = y_dep,
                              training_frame = train_h2o,
                              epochs = 1000,
                              hidden = c(50,25),
                              activation = "Tanh")
pred_ann <- as.vector(h2o.predict(ann_model, train_h2o))
rows_select_ann <- abs(training$y-as.vector(pred_ann)) < quantile(abs(training$y-as.vector(pred_ann)), 0.95)
#a
plot_residual_heatmap(residual_heatmap(x=training[rows_select_ann,4],
                                       y=training[rows_select_ann,5],
                                       res = abs(training[,1]-as.vector(pred_ann))[rows_select_ann],
                                       n_1 = 100),
                      xlab = "x3",
                      ylab = "x4",
                      pal = brewer.pal(9, "YlOrRd"))
#b
set.seed(1234)
plot_residual_heatmap(residual_heatmap(x=training[rows_select_ann,4]-training[rows_select_ann,5],
                                       y=runif(length(training[rows_select_ann,4])),
                                       res = abs(training[,1]-as.vector(pred_ann))[rows_select_ann],
                                       n_1 = 100),
                      xlab = "x3-x4",
                      ylab = "independent random variable",
                      pal = brewer.pal(9, "YlOrRd"))
#c
plot_residual_heatmap(residual_heatmap(x=training[rows_select_ann,2],
                                       y=training[rows_select_ann,3],
                                       res = abs(training[,1]-as.vector(pred_ann))[rows_select_ann],
                                       n_1 = 100),
                      xlab = "x1",
                      ylab = "x2",
                      pal = brewer.pal(9, "YlOrRd"))
#d
plot_residual_heatmap(residual_heatmap(x=training[rows_select_ann,2],
                                       y=training[rows_select_ann,4],
                                       res = abs(training[,1]-as.vector(pred_ann))[rows_select_ann],
                                       n_1 = 100),
                      xlab = "x1",
                      ylab = "x3",
                      pal = brewer.pal(9, "YlOrRd"))


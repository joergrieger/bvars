pltBvarIrf <- function(bvarObj){
  # Initialize list to store impulse-response functions
  nLength <- length(bvarObj$varnames)
  pltList <- list()

  for(ii in 1:nLength){
    for(jj in 1:nLength){
      irf1     <- bvarObj$irfdraws[ii,jj,,1]
      irfUpper <- bvarObj$irfdraws[ii,jj,,2]
      irfLower <- bvarObj$irfdraws[ii,jj,,3]
      irfLength <- length(irf1)

      # Put all the information into a data frame
      irfDf <- data.frame(x = seq(1:irfLength), irf = irf1, Upper = irfUpper, Lower = irfLower)

      # Draw the irf plot
      p1 <- ggplot(data = irfDf) + geom_line(mapping = aes(x = x,y = irf)) +
        geom_line(mapping = aes(x = x,y = Upper,color = "red")) +
        geom_line(mapping = aes(x = x,y = Lower,color = "red")) +
        theme(legend.position = "none")
      if(ii == 1){
        p1 <- p1 + ylab(bvarObj$varnames[jj])
      }
      else{
        p1 <- p1 + theme(axis.title.y = element_blank())
      }
      if(jj == 1){
        p1 <- p1 + ggtitle(bvarObj$varnames[ii])
      }
      if(jj == nLength){
        p1 <- p1 + xlab("Horizon")
      }
      else{
        p1 <- p1 + theme(axis.title.x = element_blank())
      }
      # Store all plots in a list
      pltList[[(jj-1)*nLength+ii]] <- p1
    }
  }

  # Create the final plot
  do.call("grid.arrange",c(pltList,nrow=nLength))
}

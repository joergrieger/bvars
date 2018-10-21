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

pltTvarIrf <- function(tvarObj){

  # Initialize list to store impulse-response functions
  nLength <- length(tvarObj$varnames)
  pltListReg1 <- list()
  pltListReg2 <- list()

  for(ii in 1:nLength){
    for(jj in 1:nLength){
      irf1 <- tvarObj$irf[ii,jj,,1,1]
      irf2 <- tvarObj$irf[ii,jj,,2,1]

      irfLower1 <- tvarObj$irf[ii,jj,,1,2]
      irfLower2 <- tvarObj$irf[ii,jj,,2,2]

      irfUpper1 <- tvarObj$irf[ii,jj,,1,3]
      irfUpper2 <- tvarObj$irf[ii,jj,,2,3]

      irfLength <- length(irf1)

      irfDf1 <- data.frame(x = seq(1:irfLength), irf = irf1, Upper = irfUpper1, Lower = irfLower1)
      irfDf2 <- data.frame(x = seq(1:irfLength), irf = irf2, Upper = irfUpper2, Lower = irfLower2)

      # Draw the irf plot

      # First Regime
      p1 <- ggplot(data = irfDf1) + geom_line(mapping = aes(x = x,y = irf)) +
        geom_line(mapping = aes(x = x,y = Upper,color = "red")) +
        geom_line(mapping = aes(x = x,y = Lower,color = "red")) +
        theme(legend.position = "none")

      # Second Regime
      p2 <- ggplot(data = irfDf2) + geom_line(mapping = aes(x = x,y = irf)) +
        geom_line(mapping = aes(x = x,y = Upper,color = "red")) +
        geom_line(mapping = aes(x = x,y = Lower,color = "red")) +
        theme(legend.position = "none")

      if(ii == 1){
        p1 <- p1 + ylab(tvarObj$varnames[jj])
        p2 <- p2 + ylab(tvarObj$varnames[jj])
      }
      else{
        p1 <- p1 + theme(axis.title.y = element_blank())
        p2 <- p2 + theme(axis.title.y = element_blank())
      }
      if(jj == 1){
        p1 <- p1 + ggtitle(tvarObj$varnames[ii])
        p2 <- p2 + ggtitle(tvarObj$varnames[ii])
      }
      if(jj == nLength){
        p1 <- p1 + xlab("Horizon")
        p2 <- p2 + xlab("Horizon")
      }
      else{
        p1 <- p1 + theme(axis.title.x = element_blank())
        p2 <- p2 + theme(axis.title.x = element_blank())
      }

      # Store all plots in a list
      pltListReg1[[(jj-1)*nLength+ii]] <- p1
      pltListReg2[[(jj-1)*nLength+ii]] <- p2

    }
  }

  do.call("grid.arrange",c(pltListReg1,nrow=nLength))
  readline("Press [Enter] to continue")
  do.call("grid.arrange",c(pltListReg2,nrow=nLength))

}

pltmsvarirf <- function(msvarObj){

  # Initialize list to store impulse-response functions
  nLength <- length(msvarObj$varnames)

  for(kk in 1:msvarObj$NoRegimes){
    pltList <- list()
    # Plot over all regimes
    for(ii in 1:nLength){
      for(jj in 1:nLength){
        irf1     <- msvarObj$irf[ii,jj,,kk,1]
        irfUpper <- msvarObj$irf[ii,jj,,kk,2]
        irfLower <- msvarObj$irf[ii,jj,,kk,3]
        irfLength <- length(irf1)

        # Put all the information into a data frame
        irfDf <- data.frame(x = seq(1:irfLength), irf = irf1, Upper = irfUpper, Lower = irfLower)
        # Draw the irf plot
        p1 <- ggplot(data = irfDf) + geom_line(mapping = aes(x = x,y = irf)) +
          geom_line(mapping = aes(x = x,y = Upper,color = "red")) +
          geom_line(mapping = aes(x = x,y = Lower,color = "red")) +
          theme(legend.position = "none")
        if(ii == 1){
          p1 <- p1 + ylab(msvarObj$varnames[jj])
        }
        else{
          p1 <- p1 + theme(axis.title.y = element_blank())
        }
        if(jj == 1){
          p1 <- p1 + ggtitle(msvarObj$varnames[ii])
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
    readline("Press [Enter] to continue")

  }

}

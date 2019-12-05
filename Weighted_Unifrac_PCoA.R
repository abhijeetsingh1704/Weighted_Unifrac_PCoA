      ### Code by Abhijeet Singh
      #   email - abhijeetsingh.aau@gmail.com
      
      ### R Script fot the weighted unifrac with environmental parameter fitting
      
      ### code is intended to use together with the phyloseq visualization pipeline
      
      ### phyloseq data must contain tree element
      ##############################################
      #load R libraries 
      
      library(phyloseq)
      library(vegan)
      library(ggplot2)
      
      
      # If you are using real data, escape/avoid this section
      
      
      # loading sample data from phyloseq
      data("GlobalPatterns")
      ps <- GlobalPatterns
      ps
      # extract sample data 
      sample_df <- data.frame(sample_data(ps),
                              check.rows = T)
      # visualize the sample data
      head(sample_df)
      
      ############################################## X X X X CAUTION X X X X
      # do not do it if you have your own data
      #make new column /FAKE DATA
      TEMP <- sample(0:50, size=nsamples(ps), replace=T)
      DEPTH <- sample(100:500, size=nsamples(ps), replace=T)
      AGE <- sample(1:100, size=nsamples(ps), replace=T)
      # add fake data to sample data
      sample_df$Temperature=TEMP
      sample_df$Depth=DEPTH
      sample_df$Age=AGE
      # visualize the sample data
      head(sample_df)
      ############################################## X X X X CAUTION X X X X
      
      # If you have your own real data/phyloseq object, continue from here
      # Computing weighted unifrac from phyloseq object
      
      # visualize the details of phyloseq object
      ps
      
      # Write the phyloseq details to text file
      sink(file = "PCOA_weighted.txt")
      ps
      rep("#BRAKE#", 8)
      sink()
      
      # set random seed
      set.seed(10)
      
      ### weighted unifrac
      ps_UFw <- UniFrac(ps, weighted = T) 
      
      # visualize weighted unifrac
      ps_UFw
      
      # Write the weighted unifrac details to text file
      sink(file = "PCOA_weighted.txt", append = T)
      ps_UFw
      rep("#BRAKE#", 8)
      sink()
      
      ### computing PCOA from weighted unifrac 
      
      # set random seed
      set.seed(100)
      
      ### PCoA or Classical metric/multidimendional scaling
      ps_UFw_pcoa <- cmdscale(ps_UFw, eig=TRUE)
      
      # visualize weighted unifrac PCoA / CMD
      ps_UFw_pcoa
      
      # Write the weighted unifrac PCoA / CMD details to text file
      sink(file = "PCOA_weighted.txt", append = T)
      ps_UFw_pcoa
      rep("#BRAKE#", 8)
      sink()
      
      # getting the x and y aaxis values to label the graph
      X_Axis <- paste("PCoA 1 - ", round(ps_UFw_pcoa$eig[1]/sum(ps_UFw_pcoa$eig)*100), "%")
      Y_Axis <- paste("PCoA 2 - ", round(ps_UFw_pcoa$eig[2]/sum(ps_UFw_pcoa$eig)*100), "%")
      
      # visualize
      X_Axis
      Y_Axis
      
      # MAKING DF FROM THE PCOA POINTS
      PCoA <- data.frame(PCoA1 = ps_UFw_pcoa$points[,1], 
                         PCoA2 = ps_UFw_pcoa$points[,2],
                         check.rows = T,
                         check.names = T)
      # visualize PCOA values/coordinates
      PCoA
      
      # Write the PCOA values/coordinates details to text file
      sink(file = "PCOA_weighted.txt", append = T)
      PCoA
      rep("#BRAKE#", 8)
      sink()
      
      # plot the priliminary plot, just to see the distribution of points
      plot(PCoA)
      
      # merging PCOA and sample data
      sample_df_PCOA <- merge(sample_df, PCoA, by = 0,
                        sort = TRUE)
      sample_df_PCOA
      
      
      ## fit environmental parameters on PCOA plot
      ps_UFw_pcoa_EFit <- envfit(ps_UFw_pcoa, sample_df)
      
      # visualize environment fit attributes
      ps_UFw_pcoa_EFit
      
      # Write environment fit attribute details to text file
      sink(file = "PCOA_weighted.txt", append = T)
      ps_UFw_pcoa_EFit
      rep("#BRAKE#", 8)
      sink()
      
      # calculating the scores 
      ps_UFw_pcoa_EFit_scores <- as.data.frame(scores(ps_UFw_pcoa_EFit, 
                                                      display = "vectors"))
      
      # giving name of arrow
      ps_UFw_pcoa_EFit_scores <- cbind(ps_UFw_pcoa_EFit_scores, 
                                       Species = rownames(ps_UFw_pcoa_EFit_scores))
      
      # Write environment fit scores details to text file
      sink(file = "PCOA_weighted.txt", append = T)
      paste("Singnificant parameters are listed below:")
      ps_UFw_pcoa_EFit_scores
      rep("#BRAKE#", 8)
      sink()
      
      
      #  plot the PCOA graph with environment fit parameters
      my_plot <- 
        
        # plotting
        ggplot()+
        
        # making intercept at zero position
        geom_hline(yintercept = 0, linetype="dashed",color="grey20")+
        geom_vline(xintercept = 0, linetype="dashed",color="grey20")+
        
        # lines on axis
        theme(axis.line.y.left = element_line(colour = "black"))+
        theme(axis.line.x.bottom = element_line(colour = "black"))+
        
        # axis text
        theme(axis.text.x = element_text(colour = "black", face = "bold", size = 10))+
        theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10))+
        
        # axis label style
        theme(axis.title.x = element_text(colour = "black",face = "bold",size = 12))+
        theme(axis.title.y = element_text(colour = "black",face = "bold",size = 12))+
        
        # making points
        geom_point(mapping=aes(x=sample_df_PCOA$PCoA1, y=sample_df_PCOA$PCoA2, 
                               color=sample_df_PCOA$SampleType), size=4)+
        
        #point lables
        geom_text(data= sample_df_PCOA,
                  x=sample_df_PCOA$PCoA1, 
                  y=sample_df_PCOA$PCoA2,
                  mapping= aes(label = sample_df_PCOA$Temperature),
                  size=2)+
        
        # making arrows
        geom_segment(data = ps_UFw_pcoa_EFit_scores, 
                     mapping = aes(x=0, y=0, xend=Dim1, yend=Dim2),
                     size=0.5, 
                     linetype="solid",
                     lineend = "round", 
                     color="red",
                     arrow=arrow(type="open",
                                 angle = 20,
                                 length = unit(0.3, "cm")))+
        
        # arrow lables
        geom_text(data=ps_UFw_pcoa_EFit_scores,
                  aes(x=Dim1,y=Dim2,label=Species),
                  vjust=-1, 
                  #hjust=1, lwd=5, 
                  col="black", 
                  size=4) + 
        
        
        # Legend on right
        theme(legend.position = "right")+
        # removing legend title
        theme(legend.title=element_blank())+
        
        # x and y axis label
        xlab(X_Axis)+
        ylab(Y_Axis)+
         
        # x and y axis limits           # edit if required
        #xlim(-0.4,0.4)+
        #ylim(-0.4,0.4)+
        
        #
        ggtitle("Weighted Unifrac Principal Coordinate Analysis")
        
      # rendering plot
      my_plot  
      
      
      # saving plot as image
      tiff("pcoa_ggplot.tif", width = 8, height = 6, units = "in", res = 250)
      my_plot
      dev.off()
      
    
      
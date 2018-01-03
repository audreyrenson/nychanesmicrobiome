get_edgeR_results <- function(formla, pseq=NYC_HANES, method=c("BH","IHW"),
                              coef=2, alph=0.01, filtering=method[1] == "BH", 
                              countMinimum = 8, 
                              percentMinimumHaveCount = NULL, 
                              nMinimumHaveCount = 3, NA.RM=TRUE) {
  ###################################################################
  # Get a GLM-like object filtered by FDR p values less than 'alph' #
  # for one comparison, specified by 'coef'                         #
  ###################################################################
  
  #coef=2 sets default to test first coefficient
  # If filtering = TRUE, user must provide either percentMinimumHaveCount or 
  # nMinimumHaveCount, indicating how many samples must meet the countMinimum criteria
  # If method="IHW", no filtering is performed and and filtering arguments are ignored.
  
  suppressPackageStartupMessages({
    require(edgeR)
    require(magrittr)
    require(phyloseq)
  })
  
  
  #create model matrix
  dsg.mtrx <- model.matrix(formla, data=data.frame(sample_data(pseq)))
  #drop na's
  if(NA.RM) pseq <- prune_samples(samples = rownames(dsg.mtrx), x = pseq) 
  
  #extract otu table to use as a count matrix
  otus <- as(otu_table(pseq), "matrix")
  
  #initialize DGEList object
  dge <- DGEList(counts = otus)
  
  #filtering: keep OTU if minimum number of samples has minimum count
  if(filtering & !method[1] == "IHW") {
    if(!is.null(percentMinimumHaveCount)) {
      nMinimumHaveCount = ncol(otus) * percentMinimumHaveCount
    } else if(is.null(nMinimumHaveCount)) {
      stop("User must enter either nMinimumHaveCount or percentMinimumHaveCount.")
    }
    keep <- rowSums(dge$counts >= countMinimum ) >= nMinimumHaveCount
    dge <- dge[keep, , FALSE]
  }
  

  
  #set up the glm fit
  disp <- estimateDisp(dge, design = dsg.mtrx)
  fit <- glmFit(disp, dsg.mtrx)
  lrt <- glmLRT(fit, coef = coef)
  results <- topTags(lrt, n=nrow(lrt$table))  
  
  if(method[1]=="IHW") {
    require("IHW")
    lrt$table$mean_abundance <- rowSums(otus) / ncol(otus)
    ihwResult <- ihw(PValue ~ mean_abundance, data=lrt$table, alpha=alph)
    results$table$ihw_pvalue <- ihwResult@df$weighted_pvalue
    #filter out results that didn't pass alpha
    results <- results[results$table$ihw_pvalue < alph, ]
  } else if(method[1]=="BH") {
    #filter out results that didn't pass alpha
    results <- results[results$table$FDR < alph, ]
  } else stop("Argument 'method' must be either 'BH' or 'IHW'.")

  # attach taxonomy names
  if(nrow(results$table) > 0) { #if any results pass the filter
    results$table %<>% 
      cbind(as.matrix(tax_table(pseq)[rownames(results$table),]))
  }
  

  
  results
  
}




edger_get_p <- function(formla, pseq=NYC_HANES, alph=1, ...)
{
  
  #######################################################
  # Get a vector of raw p-values for regression of all  #
  # OTUs against all levels of one metadata variable    #
  #######################################################
  vect <- sample_data(pseq)[,as.character(formla)[2]][[1]] # exctract the vector to get information about it
  
  n_coefs <- ifelse(class(vect) %in% c("character","factor"),  
                    length(levels(factor(vect))) - 1, #the number of coeficients to test will be n levels - 1
                    1) # or 1 if it's a numeric etc.
  
  list_fits <- lapply(seq_len(n_coefs) + 1, # adding one because the coefficients start at 2 (skip the intercept)
                      function(coef)  
                        get_edgeR_results(formla, pseq, coef=coef, alph=alph)  )
  
  p_vect <- unlist(sapply(seq_along(list_fits), function(i) list_fits[[i]]$table[,"PValue"]))
  p_vect
}


plot_edgeR <- function(formla, ttle=NULL, varname=NULL,
                       pseq=NYC_HANES, coef=2, printIfEmpty=FALSE, 
                       invisbl=FALSE, color="FDR", sortby=NULL, ...) {
  #########################################################################
  # Get a dot plot of log-fold-change of all significant OTUs against one #
  # factor level specified by 'coef'. Also invisibly returns a data frame #
  # of edgeR results equivalent to get_edgeR_results()$table.             #
  #########################################################################
  
  # formla and ... are passed to get_edgeR_results
  suppressPackageStartupMessages({
    require(ggplot2)
    require(dplyr)
  })
  
  DGELRT = get_edgeR_results(formla, pseq, coef = coef, ...)
  
  if(!printIfEmpty & nrow(DGELRT$table) == 0) {
    return("No results")
  } 
  
  #extract the main results table
  dat <- data.frame(DGELRT$table)
  
  #drop NA Genus factor levels
  dat <- dat[!is.na(dat$Genus), ]
  
  
  #drop NA
  
  #sort genus factor levels by log fold change
  if(missing(sortby)) {
    dat$Genus %<>% as.character %>% factor(levels = unique(dat$Genus[order(dat$logFC, decreasing = TRUE)]))
  } else {
    dat$Genus %<>% as.character %>% factor(levels = unique(dat$Genus[order(dat[[sortby]], dat$logFC, decreasing=TRUE)]))
  }
  
  
  
  #create title
  if(is.null(ttle)) {
    metadata <- data.frame(sample_data(pseq))
    #extract the variable name from the formula
    colname <- strsplit(as.character(as.formula(formla)), "\\ ")[[2]][1]
    if(is.null(varname)) varname <- colname
    if(class(metadata[,colname]) %in% c("character","factor")) {
      #if it's a categorical variable
      ttle=paste0(varname, ": ", levels(metadata[,colname])[coef], #comparing the level specified by 'coef'
                  " vs. ", levels(metadata[,colname])[1])        #to the reference level
    } else {
      #if it's a quantitative variable
      ttle=paste0(varname, ": dy/dx")
    }
  }
  
  phyla_colors <- c('#00B2FFFF','#FF9900FF','#FF0000FF','#0DFF00FF','#00FF9FFF','#A600FFFF','#B9FF00FF','#FF00ACFF','#0006FFFF')
  
  #plot
  p <- dat %>%
        ggplot(aes(x=logFC, y=Genus, color=eval(as.name(color)))) +
        geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
        geom_point(size=3) +
        guides(color=guide_legend(title=color)) + 
        theme_light() +
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
              axis.text.y = element_text(face = "italic"),
              axis.text=element_text(size=12), 
              legend.text = element_text(size=12),
              legend.title = element_text(size=14),
              title = element_text(size=15)) +
        scale_y_discrete(position="right") +
        labs(x="Log10-fold change") +
        ggtitle(ttle)
  #add the colors set permanently to phyla (for multiple plots)
  if(color=="Phylum") p <- p + scale_color_manual(values=phyla_colors)
  if(!invisbl) print(p)
  invisible(dat)
}

edger_summary_plot <- function(list_models, ttle="Significant genera by sociodemographic characteristics",
         top_n_genera=NULL) {
  
  
  require(magrittr)
  require(dplyr)
  list_models %<>% nullify_nondataframes %>% drop_nulls_from_list
  
  #create a data frame with a column of all genera appearing in any of the models
  all_sig_OTUs <- data.frame(genus =
                               na.omit(unique(unlist(lapply(list_models, function(model)
                                 unlist(lapply(model, function(dataframe)
                                   as.data.frame(dataframe)[, "Genus"])))))))
  
  
  #function to count how many of the models (including all factor levels) contain a particular OTU
  count_of_variables_with_OTU <- function(otu_vect)  {
    sapply(otu_vect, function(otu)
      sum(unlist(lapply(list_models, function(model)
        any( unlist(lapply(model, function(df)
          otu %in% df[,"Genus"]))))  ))
    )
  }
  
  #create a column in the data frame with the count for each genus
  all_sig_OTUs$count <- count_of_variables_with_OTU(all_sig_OTUs$genus)
  
  #return only top genera if top_n_genera is a number
  if(!is.null(top_n_genera)) {
    stopifnot(is.numeric(top_n_genera))
    all_sig_OTUs <- all_sig_OTUs[order(all_sig_OTUs$count, decreasing = TRUE), ][1:top_n_genera,]
  }
  
  #sort models by number of significant genera
  n_genera <- sapply(names(list_models),
                     function(i) length(unique(unlist( lapply(list_models[[i]], function(j) j$Genus)))))
  n_genera <- n_genera[order(n_genera, decreasing = TRUE)]
  list_models <- list_models[ names(n_genera) ]
  
  #update list names to include number of significant genera (not doing is if
  #top_n_genera because numbers don't make sense)
  if(is.null(top_n_genera)) names(list_models) <- paste0(names(n_genera), " (", n_genera, ")")
  
  #generate a matrix indicating which OTU appeared in which model
  matrix_OTU_appears_in_model <-
    t(
      sapply( all_sig_OTUs$genus, function(otu)
        unlist(
          lapply(list_models, function(model)
            as.numeric(
              any(
                unlist(
                  lapply(model, function(df)
                    otu %in% df[,"Genus"]
                  )
                )
              )
            )
          )
        )
      )
    )
  
  #if these two procedures produced the same result,
  if( assertthat::are_equal(all_sig_OTUs$count, rowSums(matrix_OTU_appears_in_model))) {
    #add the rownames
    rownames(matrix_OTU_appears_in_model) <- all_sig_OTUs$genus
    
    #bind it to the data frame
    all_sig_OTUs %<>%
      cbind(as.data.frame(matrix_OTU_appears_in_model)) %>%
      arrange(desc(count)) %>% # sort by count of appearances
      select(-count) # remove the count variable, not necessary
    
  }
  
  melted_appearances <- reshape2::melt(matrix_OTU_appears_in_model)
  melted_appearances$X1 <- factor(melted_appearances$X1,
                                  levels = all_sig_OTUs$genus)
  melted_appearances %>%
    ggplot(aes(x=X1, y=X2, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(name="Appeared") +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust = 1),
          legend.position = "none") +
    labs(y="Sociodemographic Variable", x="Genus", title=ttle)
  
}

#functions to drop nulls from list of models
nullify_nondataframes <- function(lst)
  lapply(lst, function(l) lapply(l, function(k) if(is.data.frame(k)) {k} else {NULL}) )
drop_nulls_from_list <- function(lst){
  lst <- lst[ !sapply( lst, is.null ) ]
  if( !is.data.frame(lst) ){
    lst <- lapply( lst, drop_nulls_from_list)
  }
  return(lst[sapply(lst, function(i) length(i) > 0)  ])
}




edger_logfc_tileplot <- function(list_list_models) {
  
  dfs <- lapply(seq_along(list_list_models), function(i) edger_logfc_tileplot_df(list_list_models[[i]], oligo_name = names(list_list_models)[i]))
  df_tile <- do.call(rbind, dfs)
  
  
  df_tile %>%
    ggplot(aes(x=OTU, y=coef, fill=logFC)) +
    geom_tile() + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust = 1)) +
    scale_fill_gradient2(limits=c(-3,3.5),low = "seagreen4", mid = "white", high = "purple3") +
    facet_wrap(~oligo_name,ncol=1, dir="h", scales="free")
}

edger_logfc_tileplot_df <- function(list_models, oligo_name) {
  
  var_df <- function(edger_list) {
    list_edger_results <- lapply(seq_along(edger_list), function(i)
      data.frame(logFC = edger_list[[i]]@.Data[[1]]$logFC,
                 OTU = rownames(edger_list[[i]]@.Data[[1]]),
                 coef = edger_list[[i]]@.Data[[3]][1])
    )
    df <- do.call(rbind, list_edger_results)
    # df$OTU <- factor(df$OTU, labels = paste0("Oligotype ", levels(as.factor(df$OTU))))
    df <- aggregate(logFC~OTU+coef,df,function(x) x[which.max(abs(x))]) #when more than one OTU for a given genus was significant, select the max logFC to display
    df
  }
  
  
  for(i in seq_along(list_models)) list_models[[i]][sapply(list_models[[i]], nrow)==0] <- NULL 
  list_models[sapply(list_models, length)==0] <- NULL
  
  df_tile <- lapply(list_models, var_df)
  for(i in seq_along(df_tile)) df_tile[[i]]$variable <- names(df_tile)[i]
  df_tile <- do.call(rbind, df_tile)
  df_tile$OTU <- factor(df_tile$OTU, labels = paste0("O", 1:nlevels(df_tile$OTU)))
  df_tile$oligo_name <- oligo_name

  
  df_tile
}


edger_list_to_data.frame <- function(list_models) {
  #@param edger_list A list of objects each returned by 
  #                          lapply(2:length(levels(x)), function(i) get_edgeR_results(~x,coef=i))
  var_df <- function(edger_list) {
    list_edger_results <- lapply(seq_along(edger_list), function(i)
      data.frame(logFC = edger_list[[i]]@.Data[[1]]$logFC,
                 OTU = rownames(edger_list[[i]]@.Data[[1]]),
                 coef = edger_list[[i]]@.Data[[3]][1],
                 FDR = edger_list[[i]]@.Data[[1]]$FDR)
    )
    df <- do.call(rbind, list_edger_results)
    # df$OTU <- factor(df$OTU, labels = paste0("Oligotype ", levels(as.factor(df$OTU))))
    df <- aggregate(logFC~OTU+coef+FDR,df,function(x) x[which.max(abs(x))]) #when more than one OTU for a given genus was significant, select the max logFC to display
    df
  }
  for(i in seq_along(list_models)) list_models[[i]][sapply(list_models[[i]], nrow)==0] <- NULL 
  list_models[sapply(list_models, length)==0] <- NULL
  
  df_models <-  do.call(rbind, lapply(seq_along(list_models), function(i) data.frame(var_df(list_models[[i]]), variable=names(list_models)[i])) )
  df_models
}

get_all_edgeR_models <-  function(vars, varlabels, adjusted_for, to.data.frame=TRUE, ...) {
  as.form <- function(vars) as.formula(paste0("~", paste(vars, collapse="+")))
  
  if(missing(adjusted_for)) {
    models = lapply(seq_along(vars), 
                  function(i) get_edger_results_all_levels(as.form(vars[i]), ...))
  } else {
    models = lapply(seq_along(vars), 
                    function(i) get_edger_results_all_levels(as.form(c(vars[i], adjusted_for)), ...))
  }
  names(models) <- varlabels
  if(to.data.frame) edger_list_to_data.frame(models) else models
}

get_edger_results_all_levels <- function(formla, ...) {
  #the first variable in formula should be the one you want all levels of
  var1 <- strsplit(as.character(formla)[2], " \\+ ")[[1]][1]
  lapply(2:length(levels(metadata[[var1]])), 
         function(i) get_edgeR_results(formla = formla,coef=i, ...))
}



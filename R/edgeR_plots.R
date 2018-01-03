#' Dot plots based on edgeR results for microbiome data
#'
#' Plots results of \code{get_edgeR_results_all_levels()} using the \code{ggplot2} engine.
#'
#' @param formla formula. specifies the design matrix used by \code{edgeR::glmFit}.
#' @param ttle character. Title of plot.
#' @param varname character. Label for independent variable, for plotting.
#' @param pseq object of class \code{phyloseq}
#' @param coef integer. Specifies which linear model coefficient to test (default 2).
#' @param printIfEmpty logical. If no results pass alpha, whether or not to plot an empty plot.
#' @param invisbl logical. Whether or not to invisibly return the \code{topTags} object returned by \code{get_edgeR_results}.
#' @param color character. Variable  to color by, from either the \code{tax_table} or the \code{topTags$table}.
#' @param sortby character. Variable to sort the plot by, from either the \code{tax_table} or the \code{topTags$table}.
#' @param ... further arguments passed to \code{get_edgeR_results}.
#' @export
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

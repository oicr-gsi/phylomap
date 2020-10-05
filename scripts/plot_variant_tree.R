#
# Plot a phylogenetic tree with associated mutations
#

#
# this method of aligning the tree and tiles is from:
# https://thackl.github.io/ggtree-composite-plots
#



require(optparse)

option_list = list(
  make_option(c("-t", "--tree"), type="character", default="tree.nwk", 
              help="newick tree file", metavar="character"),
  make_option(c("-v", "--variants"), type="character", default="alleles.tsv", 
              help="variants file", metavar="character"),
  make_option(c("-c", "--covariates"), type="character", default=NULL, 
              help="covariate file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="plots", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

require(ggtree)
require(ggplot2)
require(patchwork)

# fix the alignment of the tree panel
scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
  scale_y_continuous(expand=expand, ...)
}

### theme for the covariate plots
themes<-list()
themes$covariates<-theme(axis.line = element_blank(), 
                       axis.title.y = element_blank(), 
                       #axis.title.x = element_blank(),
                       axis.title.x = element_text(angle = 90,hjust=0.5),
                       axis.ticks.y = element_blank(), 
                       axis.ticks.x = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       # panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       legend.position = "top",
                       plot.margin=unit(c(0,0,0,0),"cm"))


themes$variants<-theme(axis.line = element_blank(), 
                       axis.title.y = element_blank(), 
                       axis.title.x = element_blank(),
                       axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                       axis.text.x = element_text(angle = 90, hjust = 1),
                       legend.position = "top")

plot_tree<-function(t)
{
  g <- ggtree(t) + 
    #geom_tiplab(align=TRUE) + 
    #xlim(0, 0.00050) +
    scale_y_tree() +
    theme_tree2()
 return(g)
}

plot_variants<-function(a)
{
  a$pos = factor(alleles$pos)
  bases = c("A", "C", "G", "T", "N", "X")
  a$alt_allele = factor(a$alt_allele, bases)
  ambiguous_bases = !(a$alt_allele %in% bases)
  a$alt_allele[ambiguous_bases] = "X"

  cols <- c("blue", "red", "green3", "purple3", "lightgrey", "black")
  g <- ggplot(a, aes(x=pos, y=name)) + 
    geom_tile(aes(fill=alt_allele), color="white") +
    ylim(tip.order) +
    theme_bw() + 
    themes$variants + 
    scale_fill_manual(name="Variant", values=cols, drop=FALSE)
  
  return(g) 
}

plot_variants.2<-function(a)
{
  a$pos = factor(alleles$pos)
  bases = c("A", "C", "G", "T", "N")
  a$alt_allele = factor(a$alt_allele, bases)
  ambiguous_bases = !(a$alt_allele %in% bases)
  a$alt_allele[ambiguous_bases] = "X"
  
  cols <- c("blue", "red", "green3", "purple3", "lightgrey")
  g <- ggplot(a, aes(x=pos, y=name)) + 
    geom_tile(aes(fill=alt_allele), color="white") +
    ylim(tip.order) +
    theme_bw() + 
    themes$variants + 
    scale_fill_manual(name="Variant", values=cols, drop=FALSE)
  
  return(g) 
}



plot_covariate<-function(c,i)
{
  g<- ggplot(c,aes(x=1,y=name)) + 
    geom_tile(aes_string(fill=i)) + 
    xlab(i) +
    themes$covariate +
    scale_fill_discrete(na.translate=FALSE)
  return(g)
}

plot_covariate_text<-function(c,i)
{
  g<- ggplot(c,aes(x=1,y=name)) + 
    #geom_tile(aes(fill="white",colour="white")) + 
    geom_tile(colour="white",fill="white") + 
    geom_text(aes_string(label=i),size=1.6) +
    xlab(i) +
    themes$covariate +
    theme(legend.position = "none") +
    scale_fill_discrete(na.translate=FALSE)

  return(g)
}

#### MAIN
tree_path<-opt$tree
covariates_path<-opt$covariates
alleles_path<-opt$variants
prefix<-opt$prefix


### a list of plot panels, ordered from left to right
panels<-list()

### get the tree plot		
tree <- read.tree(tree_path)
# get order of tips in the tree
# from: https://groups.google.com/forum/#!topic/bioc-ggtree/LqRDK78m3U4
d = fortify(tree)
d = subset(d, isTip)
tip.order = with(d, label[order(y, decreasing=F)])
panels$tree<-plot_tree(tree)


### get the covariate plots, each a heatmap with no margins
if(! is.null(covariates_path))
{
  #refdf<-data.frame(name="MN908947.3",location=NA,site=NA,batch=NA)
  covariates<-read.table(covariates_path, header=T,as.is=T)
  covariate.ids<-names(covariates)[-1] 

  ### only use covariates that are in tip.order
  covariates<-covariates[covariates$name %in% tip.order,]
  ### add covariates that are missing from tip.order
  missing<-tip.order[! tip.order %in% covariates$name]
  #missing.df<-covariates[0,]
  missing.df<-data.frame(matrix(NA,ncol=ncol(covariates),nrow=length(missing)))
  colnames(missing.df)<-colnames(covariates)
  missing.df$name<-missing
  covariates<-rbind(covariates,missing.df)
  
  covariates$name = factor(covariates$name, levels=tip.order)
  
  for (id in covariate.ids){
    if(id=="CollectionDate"){
      panels[[id]]<-plot_covariate_text(covariates,id)
    }else{
      panels[[id]]<-plot_covariate(covariates,id)
    }
  }
}else{
  covariate.ids<-NULL
}


#### get the allele plot, a heatmap
alleles <- read.table(alleles_path, header=T)
alleles$name = factor(alleles$name, levels=tip.order)
#panels$variants<-plot_variants(alleles)   
panels$variants<-plot_variants.2(alleles) 





### set the panel width, 1 for the tree, 0.1 for each covariate and 4 for the variants
### the proportion for the variant panel should adjust based on the number of variants
panel.widths<-c(1,rep(0.1,length(covariate.ids)),4)
p<-wrap_plots(panels) + plot_layout(widths = panel.widths, guides="collect")

## save the plot    
#plot_path = sprintf("plots/%s_variant_tree.pdf", prefix)
plot_path = opt$out

# count number of samples, for scaling the plot
num_samples = nrow(d)
cat("number of samples for this heatmap is ",num_samples,"\n")
pdf_height = 0.125 * num_samples
pdf_height = max(8, pdf_height)
## wdith shoud be a factor of the number of variants
ggsave(plot_path, p, height=pdf_height, width=40, limitsize=F)

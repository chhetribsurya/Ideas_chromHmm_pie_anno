#library(plotrix)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)

### Dnase Whole genome piechart:
args <-  commandArgs(TRUE)
input_file <- args[1]
tf_name <- args[2]
out_dir_name <- args[3]

read_file <- fread("~/Dropbox/encode_3/tf_hotspots_dnase/files/merged_d_100/final_dnase_wholegenome_piechart_distribution_data.txt", sep="\t")
output_dir <-  "~/Dropbox/encode_3/tf_hotspots_dnase/plots"
tf_name <- "Dnase whole genome"
names(read_file) <- c("states", "intersect_count")
read_file

ideas_ordered_file <- fread("~/Dropbox/ideas_table.txt", sep="\t")
names(ideas_ordered_file) <-  c("State_index", "Mnemonics", "Rationale", "COLORCODE" )
target_vector <- ideas_ordered_file$Mnemonics

read_file <- merge(read_file, ideas_ordered_file, by.x="states", by.y="Mnemonics", all.x=FALSE, all.y=TRUE)
read_file <- read_file %>% select(states, intersect_count, COLORCODE)
read_file[is.na(read_file$intersect_count)]$intersect_count <- 0

#Order data frame rows according to a target vector that specifies the desired order:
read_file <- read_file[match(target_vector, read_file$states),]

#read_file_df  <- read_file %>% mutate(intersect_count=ifelse(is.na(intersect_count), 0, intersect_count))
# Changing the default alphabetical order for the discrete x-axis values and making them show up in the same order as in the dataframe...
# Maintaining the order of the dataframe...like before...
read_file$states <- factor(read_file$states, levels=read_file$states)
factor(read_file$states)
state_level <- levels(factor(read_file$states))
state_level
read_file$pct <- percent(read_file$intersect_count/sum(read_file$intersect_count)) #pct <- with(read_file, round(intersect_count/sum(intersect_count)*100))


# For extracting the color codes, copy and paste the data table from the browser in sublime, and name it as chromHMM_table.txt
ideas_ordered_file$COLORCODE 
read_file$rgb_code <- ideas_ordered_file$COLORCODE
cols <- sapply(strsplit(as.character(read_file$rgb_code), ","), function(x) {
  rgb(x[1], x[2], x[3], m=255) #rgb(rval, bval, gval, max)
}) 

#read_file$hex_code <- cols

prom_col <- colorRampPalette(c("red4","red"))
prom_flank_col <- colorRampPalette(c("firebrick2","firebrick1"))
prom_weak_col <- colorRampPalette(c("indianred2","indianred1"))
genebody_col <- colorRampPalette(c("darkgreen", "lightgreen"))
elonw_col <- colorRampPalette(c("seagreen"))
strong_enh_col <- colorRampPalette(c("orange3", "orange2"))
weak_enh_col <- colorRampPalette(c("yellow3", "yellow"))
dnase_col <- colorRampPalette(c("gold2","gold1"))
faire_col <-  colorRampPalette(c("khaki2","khaki1"))
ctcf_col <- colorRampPalette(c("purple4", "purple1"))
repr_col  <- colorRampPalette(c("gray28", "gray50"))
low_col <- colorRampPalette(c("gray67", "gray72" ))
repeats_region <- colorRampPalette(c("slateblue4", "slateblue1"))
zero_region <-  colorRampPalette("gainsboro")

custom_col <- c(prom_col(4), prom_flank_col(2), prom_weak_col(2), genebody_col(5),
    elonw_col(1),strong_enh_col(2), weak_enh_col(4), dnase_col(2), faire_col(2), ctcf_col(2),
    repr_col(4), low_col(2), repeats_region(3), zero_region(1) )

read_file$hex_code <- custom_col
read_file$rgb_code <- NULL


#changing from white to grey for the quiescent state:
#cols[25] <- rgb(127, 255, 212, m=255) ## find hex color code if rgb code given
#a <-colorRampPalette(c("grey")) ## find hex color code if color name given
#cols[25] <- a(1)

# Removing axes ticks & titles and modifying text size & color.....
blank_theme <-  theme_minimal()+
				theme(
					axis.ticks=element_blank(), # the axis ticks
					axis.title=element_blank(), # the axis labels
					axis.text.y=element_blank(),  # the 0.75, 1.00, 1.25 labels
					axis.text.x=element_blank(),   #numbers alongside the circle that is x axis
					#axis.text.x=element_text(color='black',size=15))
					legend.title=element_text(face="bold", hjust = 0.5),   #legend title
					legend.text=element_text(colour="black", size = 5), #color of the legend text
					legend.position="right",
					legend.background = element_rect(),
					panel.border = element_blank(),
  					panel.grid=element_blank(),
  					#plot.background = element_rect(),
  					plot.title=element_text(size=14, face="bold", hjust = 0.6)
				) 

#output_file_name <- paste0("~/Dropbox/local_miscellaneous_data/misc", "/", tf_name, "_ideas_piechart.pdf")				
output_file_name <- file.path(output_dir, "Dnase_whole_genome_ideas_piechart.pdf")				
pdf(output_file_name)
# The stacked bar plot..
bar <- ggplot(read_file, aes(x="", intersect_count, fill=states)) + 
			geom_bar(width = 1, stat = "identity", colour="black") +
			ggtitle(paste(tf_name,"distribution with IDEAS segmentation")) +
			##override black diagonal line from legend
			guides(fill=guide_legend(override.aes=list(colour=NA)))
			#guides(fill = guide_legend(keywidth = 3, keyheight = 1)) # legend/key size
			#guides(fill = guide_legend(title = "LEFT", title.position = "left")) # legend title position

# Converting the stacked bar plot to a pie chart..			
pie_plot <- bar + coord_polar("y", start=0) + 
			scale_fill_manual(name="IDEAS States", values = read_file$hex_code, na.value="grey50") + 
			blank_theme #guides(fill = guide_legend(title = "chromHMM State"))
print(pie_plot)
dev.off()






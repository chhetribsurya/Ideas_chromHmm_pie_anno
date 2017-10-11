#library(plotrix)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)

#read_file <- fread("~/Dropbox/local_miscellaneous_data/test_data/ENCODE_poster/chromHMM_overlap_encode/intersectbed_counts_with_chromHMM.txt", sep="\t")
#tf_name <- "RNA POLII"
args <-  commandArgs(TRUE)
input_file <- args[1]
tf_name <- args[2]
out_dir_name <- args[3]

#print(paste(input_file, "filename...."))
#print(paste(tf_name, "tf name...."))
#print(paste(out_dir_name, "output dir name...."))
#input_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chromHMM_overlap_encode/FOXO1_intersectbed_counts_with_chromHMM.txt"

read_file <- fread(input_file, sep="\t")
names(read_file) <- c("states", "intersect_count")
read_file <- read_file[mixedorder(read_file$states), ]
# Changing the default alphabetical order for the discrete x-axis values and making them show up in the same order as in the dataframe...
# Maintaining the order of the dataframe...like before...
read_file$states <- factor(read_file$states, levels=read_file$states)
factor(read_file$states)
state_level <- levels(factor(read_file$states))
state_level
read_file$pct <- percent(read_file$intersect_count/sum(read_file$intersect_count)) #pct <- with(read_file, round(intersect_count/sum(intersect_count)*100))


# For extracting the color codes, copy and paste the data table from the browser in sublime, and name it as chromHMM_table.txt
df <- fread("~/Dropbox/chromHMM_table.txt")
names(df) <-  str_replace(names(df), " ", "")
df$COLORCODE 
read_file$rgb_code <- df$COLORCODE

cols <- sapply(strsplit(as.character(read_file$rgb_code), ","), function(x) {
  rgb(x[1], x[2], x[3], m=255) #rgb(rval, bval, gval, max)
}) 

#changing from white to grey for the quiescent state:
#cols[25] <- rgb(127, 255, 212, m=255) ## find hex color code if rgb code given
a <-colorRampPalette(c("grey")) ## find hex color code if color name given
cols[25] <- a(1)
read_file$hex_code <- cols

# Removing axes ticks & titles and modifying text size & color.....
blank_theme <-  theme_minimal()+
				theme(
					axis.ticks=element_blank(), # the axis ticks
					axis.title=element_blank(), # the axis labels
					axis.text.y=element_blank(),  # the 0.75, 1.00, 1.25 labels
					axis.text.x=element_blank(),   #numbers alongside the circle that is x axis
					#axis.text.x=element_text(color='black',size=15))
					legend.title=element_text(face="bold"),   #legend title
					legend.text=element_text(colour="black", size = 8), #color of the legend text
					legend.position="right",
					legend.background = element_rect(),
					panel.border = element_blank(),
  					panel.grid=element_blank(),
  					#plot.background = element_rect(),
  					plot.title=element_text(size=14, face="bold", hjust = 0.6)
				) 

output_file_name <- paste0(out_dir_name, "/", tf_name, "_chromHMM_piechart.pdf")				
pdf(output_file_name)
# The stacked bar plot..
bar <- ggplot(read_file, aes(x="", intersect_count, fill=states)) + 
			geom_bar(width = 1, stat = "identity", colour="black") +
			ggtitle(paste(tf_name,"distribution with ChromHMM model")) +
			##override black diagonal line from legend
			guides(fill=guide_legend(override.aes=list(colour=NA)))
			#guides(fill = guide_legend(keywidth = 3, keyheight = 1)) # legend/key size
			#guides(fill = guide_legend(title = "LEFT", title.position = "left")) # legend title position

# Converting the stacked bar plot to a pie chart..			
pie_plot <- bar + coord_polar("y", start=0) + 
			scale_fill_manual(name="ChromHMM States", values = read_file$hex_code, na.value="grey50") + 
			blank_theme #guides(fill = guide_legend(title = "chromHMM State"))
print(pie_plot)
dev.off()






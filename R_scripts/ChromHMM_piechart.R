#library(plotrix)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)


read_file <- fread("~/Dropbox/local_miscellaneous_data/test_data/ENCODE_poster/chromHMM_overlap_encode/intersectbed_counts_with_chromHMM.txt", sep="\t")
tf_name <- "RNA POLII"
#print(paste(tf_name, "tf name...."))

read_file <- fread("~/Dropbox/local_miscellaneous_data/test_data/ENCODE_poster/chromHMM_overlap_encode/intersectbed_counts_with_chromHMM.txt", sep="\t")
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

output_file_name <- paste0("~/Dropbox/local_miscellaneous_data/misc", "/", tf_name, "_chromHMM_piechart.pdf")				
pdf(output_file_name)

# The stacked bar plot..
bar <- ggplot(read_file, aes(x="", intersect_count, fill=states)) + 
			geom_bar(width = 1, stat = "identity", colour="black") +
			ggtitle(paste(tf_name,"distribution with ChromHMM model"))+
			##override black diagonal line from legend
			guides(fill=guide_legend(override.aes=list(colour=NA)))
			
# Converting the stacked bar plot to a pie chart..			
pie_plot <- bar + coord_polar("y", start=0) + 
			scale_fill_manual(name="ChromHMM States", values = read_file$hex_code, na.value="grey50") + 
			blank_theme #guides(fill = guide_legend(title = "chromHMM State"))
print(pie_plot)
dev.off()


# In order for the text to appear in the middle of each arc representing a particular piece of the pie, a new variable called midpoint is created that identifies the middle of each arc. Using that and changing the numbers into % using the "scales"" library discussed before...Plot 12...
#midpoint <- cumsum(read_file$intersect_count) - read_file$intersect_count/2
#pie_plot <- pie_plot + scale_y_continuous(breaks=midpoint, labels=percent(read_file$intersect_count/sum(read_file$intersect_count)))


##########################################
##########################################
##########################################


# read_file <- fread("~/Dropbox/local_miscellaneous_data/test_data/ENCODE_poster/chromHMM_overlap_encode/intersectbed_counts_with_chromHMM.txt", sep="\t")
# names(read_file) <- c("states", "intersect_count")

# read_file$states <- str_replace(read_file$states, "\\d+_", "") 
# factor(read_file$states)
# levels(read_file$states)

# slices <- read_file$intersect_count
# lbls <- read_file$states
# pct <- with(read_file, round(intersect_count/sum(intersect_count)*100))
# pct_lbls <- paste(lbls, pct)
# lbls <-  paste(pct_lbls, "%", sep="")

# pdf("~/Dropbox/local_miscellaneous_data/test_data/ENCODE_poster/chromHMM_overlap_encode/Pol2_dist_chromHMM.pdf")
# pie3D(slices, explode=0.1, labels=lbls, labelcex=0.6, radius=0.9, labelrad=1, minsep=0.5, color=rainbow(length(slices)), main="Pie chart distribution for POLII with chromHMM genomic segmentation")
# legend(-0.5,1, read_file$states, cex=0.6, fill=rainbow(length(slices)))
# dev.off()


# read_file <- fread("~/Dropbox/local_miscellaneous_data/test_data/ENCODE_poster/ideas_overlap_encode/intersectbed_counts_with_ideas.txt", sep="\t")
# names(read_file) <- c("states", "intersect_count")

# read_file$states <- str_replace(read_file$states, "\\d+_", "") 
# read_file
# slices <- read_file$intersect_count
# lbls <- read_file$states
# pct <- with(read_file, round(intersect_count/sum(intersect_count)*100))
# pct_lbls <- paste(lbls, pct)
# lbls <-  paste(pct_lbls, "%", sep="")


# pdf("~/Dropbox/local_miscellaneous_data/test_data/ENCODE_poster/ideas_overlap_encode/Pol2_dist_ideas.pdf")
# pie3D(slices, explode=0.2, labels=lbls, labelcex=0.6, radius=0.9, labelrad=1, minsep=0.5, color=rainbow(length(slices)), main="Pie chart distribution for POLII with ideas genomic segmentation")
# legend(-0.5,1, read_file$states, cex=0.6, fill=rainbow(length(slices)))
# dev.off()

#line = DataFrame({"onset": 30.0, "length": 1.3}, index=[3])
#df2 = concat([df.ix[:2], line, df.ix[3:]]).reset_index(drop=True)






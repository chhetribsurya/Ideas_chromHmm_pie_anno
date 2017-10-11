library(ggplot2)
library(data.table)
library(stringr)
library(scales)
library(dplyr)
library(gtools)

chromHMM_hepg2_file <- "~/Dropbox/local_miscellaneous_data/test_data/chrom_impute_chromHMM_data/corrected/E118_25_imputed12marks_dense.bed"
chromHMM_file <- fread(chromHMM_hepg2_file, sep = "\t")
colnames(chromHMM_file) <- c("chrom", "start", "end", "states", "V5", "V6", "V7", "V8", "colorcode")

head(chromHMM_file)

#chromHMM_file$states <- str_replace(chromHMM_file$V4, "\\d+_", "") 

chromHMM_file$elm_lgth <- chromHMM_file$end - chromHMM_file$start
chromHMM_file

#state_levels <- chromHMM_file$states %>% as.factor %>% levels
#state_levels
elem_sum_table <- chromHMM_file %>% group_by(states) %>% summarize(elem_sum = sum(elm_lgth))
Total_bp = sum(elem_sum_table$elem_sum)
elem_sum_table$pct <- with(elem_sum_table, elem_sum/sum(elem_sum) *100)

ordered_chromHMM_file <- elem_sum_table[order(elem_sum_table$elem_sum),]
ordered_chromHMM_file

ordered_chromHMM_file <- ordered_chromHMM_file[mixedorder(ordered_chromHMM_file$states), ]
# Changing the default alphabetical order for the discrete x-axis values and making them show up in the same order as in the dataframe...
# Maintaining the order of the dataframe...like before...
ordered_chromHMM_file$states <- factor(ordered_chromHMM_file$states, levels=ordered_chromHMM_file$states)
factor(ordered_chromHMM_file$states)
state_level <- levels(factor(ordered_chromHMM_file$states))
state_level

df <- fread("~/Dropbox/chromHMM_table.txt")
names(df) <-  str_replace(names(df), " ", "")
df$COLORCODE 
ordered_chromHMM_file$rgb_code <- df$COLORCODE
ordered_chromHMM_file
cols <- sapply(strsplit(as.character(ordered_chromHMM_file$rgb_code), ","), function(x) {
  rgb(x[1], x[2], x[3], m=255) #rgb(rval, bval, gval, max)
}) 

#changing from white to grey for the quiescent state:
#cols[25] <- rgb(127, 255, 212, m=255) ## find hex color code if rgb code given
a <-colorRampPalette(c("grey")) ## find hex color code if color name given
cols[25] <- a(1)
ordered_chromHMM_file$hex_code <- cols

# Removing axes ticks & titles and modifying text size & color.....
blank_theme <-  theme_minimal() +
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
output_file_name <- paste0("/home/surya/Dropbox/local_miscellaneous_data/chip_hepg2_tf/pie_charts", "/", "whole_genome_chromHMM_piechart.pdf")
model_name = "chromHMM model"				
pdf(output_file_name)

# The stacked bar plot..
bar <- ggplot(ordered_chromHMM_file, aes(x="", pct, fill=states)) + 
			geom_bar(width = 1, stat = "identity", colour="black") +
			ggtitle(paste("Whole Genome distribution of", model_name))+
			##override black diagonal line from legend
			guides(fill=guide_legend(override.aes=list(colour=NA)))
			
# Converting the stacked bar plot to a pie chart..			
pie_plot <- bar + coord_polar("y", start=0) + 
			scale_fill_manual(name="ChromHMM States", values = ordered_chromHMM_file$hex_code, na.value="grey50") + 
			blank_theme #guides(fill = guide_legend(title = "chromHMM State"))
print(pie_plot)
dev.off()



#############################################
#############################################
#############################################



input_file1 = "/home/surya/Dropbox/local_miscellaneous_data/genomic_models/hepg2_ideas"

ideas_file <- fread(input_file1, sep = "\t")
colnames(ideas_file) <- c("chrom", "start", "end", "states", "V5", "V6", "V7", "V8")
head(ideas_file)

#ideas_file$states <- str_replace(ideas_file$V4, "\\d+_", "") 

ideas_file$elm_lgth <- ideas_file$end - ideas_file$start
ideas_file

#state_levels <- ideas_file$states %>% as.factor %>% levels
#state_levels
elem_sum_table <- ideas_file %>% group_by(states) %>% summarize(elem_sum = sum(elm_lgth))
Total_bp = sum(as.numeric(elem_sum_table$elem_sum))
elem_sum_table$pct <- with(elem_sum_table, elem_sum/sum(as.numeric(elem_sum)) *100)

elem_sum_table
ideas_file <- elem_sum_table[order(elem_sum_table$elem_sum),]
ideas_file
#sum(ordered_ideas_file$pct)


ideas_ordered_file <- fread("~/Dropbox/ideas_table.txt", sep="\t")
names(ideas_ordered_file) <-  c("State_index", "Mnemonics", "Rationale", "COLORCODE" )
target_vector <- ideas_ordered_file$Mnemonics

#df <- merge(read_file_ideas, read_file, by.x="Mnemonics", by.y="states")
#Order data frame rows according to a target vector that specifies the desired order:
ideas_file <- ideas_file[match(target_vector, ideas_file$states),]
ideas_file

# Changing the default alphabetical order for the discrete x-axis values and making them show up in the same order as in the dataframe...
# Maintaining the order of the dataframe...like before...
ideas_file$states <- factor(ideas_file$states, levels=ideas_file$states)
factor(ideas_file$states)
state_level <- levels(factor(ideas_file$states))
state_level
ideas_file$percent <- percent(ideas_file$elem_sum/sum(as.numeric(ideas_file$elem_sum))) #pct <- with(ideas_file, round(intersect_count/sum(intersect_count)*100))
ideas_file
# For extracting the color codes, copy and paste the data table from the browser in sublime, and name it as chromHMM_table.txt
ideas_ordered_file$COLORCODE 
ideas_file$rgb_code <- ideas_ordered_file$COLORCODE

cols <- sapply(strsplit(as.character(ideas_file$rgb_code), ","), function(x) {
  rgb(x[1], x[2], x[3], m=255) #rgb(rval, bval, gval, max)
}) 

ideas_file$hex_code <- cols
ideas_file
#changing from white to grey for the quiescent state:
#cols[25] <- rgb(127, 255, 212, m=255) ## find hex color code if rgb code given
#a <-colorRampPalette(c("grey")) ## find hex color code iiiiif color name given
#cols[25] <- a(1)

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

output_file_name <- paste0("/home/surya/Dropbox/local_miscellaneous_data/chip_hepg2_tf/pie_charts", "/", "whole_genome_ideas_piechart.pdf")
model_name = "IDEAS segmentation"				

pdf(output_file_name)
# The stacked bar plot..
bar <- ggplot(ideas_file, aes(x="", pct, fill=states)) + 
			geom_bar(width = 1, stat = "identity", colour="black") +
			ggtitle(paste("Whole Genome distribution of", model_name)) +
			##override black diagonal line from legend
			guides(fill=guide_legend(override.aes=list(colour=NA)))
			#guides(fill = guide_legend(keywidth = 3, keyheight = 1)) # legend/key size
			#guides(fill = guide_legend(title = "LEFT", title.position = "left")) # legend title position

# Converting the stacked bar plot to a pie chart..			
pie_plot <- bar + coord_polar("y", start=0) + 
			scale_fill_manual(name="IDEAS States", values = ideas_file$hex_code, na.value="grey50") + 
			blank_theme #guides(fill = guide_legend(title = "chromHMM State"))
print(pie_plot)
dev.off()








# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size = 4
# this assumes the BioGeoBears script has been run to the point of importing the data and 
# generating the "tipranges" object
areas = getareas_from_tipranges_object(tipranges)
#areas = c("A", "B", "C", "D")

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, 
                                                    include_null_range=TRUE)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based))
{    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list

#====================
areanames = names(tipranges@df)
areanames

include_null_range = TRUE

states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, 
                                                          include_null_range=include_null_range)

statenames = areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, 
                                           include_null_range=include_null_range, split_ABC=FALSE)
statenames

#create default color matrix based on the number of areanames (e.g., 4 areas)
colors_matrix = get_colors_for_numareas(length(areanames))

#mix colors for states
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, 
                                               plot_null_range=include_null_range)
colors_list_for_states
#lists colors for all 16 states as 6-digit codes

#info on colors from http://bxhorn.com/r-color-basics/

#convert RGB colors to hex triplicates
rgb(255, 165, 0, maxColorValue = 255)
#[1] "#FFA500"
#r, g, b, colors values are from 0 to 256.

# Find 3-digit RGB codes using color name

col2rgb("yellow")
#[,1]
#red    255
#green  255
#blue     0

col2rgb(c(R = "red", O = "orange", Y = "yellow"))
#R   O   Y
#red   255 255 255
#green 0   165 255
#blue  0   0   0

#function to display both hex and decimal values of a given color name.
GetColorHexAndDecimal <- function(color)
{
  c <- col2rgb(color)
  sprintf("#%02X%02X%02X %3d %3d %3d", c[1],c[2],c[3], c[1], c[2], c[3])
}
GetColorHexAndDecimal("yellow")
#[1] "#FFFF00 255 255 0"

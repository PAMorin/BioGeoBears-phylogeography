plot_BioGeoBEARS_results <-
function(results_object, analysis_titletxt=NULL, addl_params=list(), plotwhat="text", label.offset=NULL, tipcex=0.8, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, plotlegend=FALSE, legend_ncol=NULL, legend_cex=1, cornercoords_loc="manual", tr=NULL, tipranges=NULL, if_ties="takefirst", pie_tip_statecex=0.7, juststats=FALSE, xlab="Millions of years ago", manual_ranges_txt=NULL, root.edge=TRUE, colors_list_for_states=NULL, skiptree=FALSE, show.tip.label=TRUE, tipcol="black", dej_params_row=NULL, plot_max_age=NULL, skiplabels=FALSE, plot_stratum_lines=TRUE, include_null_range=NULL, plot_null_range=FALSE, simplify_piecharts=FALSE, tipboxes_TF=TRUE, tiplabel_adj=c(0.5), no.margin=FALSE, xlims=NULL, ylims=NULL)
	{
	
	junk='
	scriptdir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/a_scripts/"
	plot_BioGeoBEARS_results(results_object, analysis_titletxt=NULL, addl_params=list(), plotwhat="text", label.offset=NULL, tipcex=0.8, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=NULL, tipranges=NULL)
	
	# Defaults
	addl_params=list("j"); plotwhat="text"; label.offset=0.45; tipcex=0.7; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges; juststats = FALSE; plotlegend=FALSE; 	xlab="Millions of years ago"; if_ties="takefirst"
	
	
	# Setup
results_object = resDEC
analysis_titletxt ="BioGeoBEARS DEC on Mariana M1v4_unconstrained"
addl_params=list("j"); plotwhat="text"; label.offset=0.45; tipcex=0.7; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges
juststats=FALSE; plotlegend=FALSE; 	xlab="Millions of years ago"; if_ties="takefirst"
show.tip.label=TRUE
tipcol="black"; dej_params_row=NULL; plot_max_age=NULL; skiplabels=FALSE; 
colors_list_for_states=NULL
skiptree=FALSE
include_null_range=NULL
plot_stratum_lines=TRUE
	plot_null_range = FALSE
	' # endjunk
	
	# Default; no longer used
	if (is.null(include_null_range) == TRUE)
		{
		include_null_range = results_object$inputs$include_null_range
		}
	
	#######################################################
	# User modifications to border color (externally, using
	# 'par(fg=NA)' or whatever
	#######################################################
	# border color (for modifying this)
	tmp_fg = par("fg")
	par(fg="black")	# set to default for most things

	
	#######################################################
	# Plot ancestral states - DEC
	#######################################################


	# Setup
	#results_object = resDEC_strat
	BioGeoBEARS_run_object = results_object$inputs

	# Read the tree from file, if needed
	if (is.null(tr))
		{
		tr = read.tree(BioGeoBEARS_run_object$trfn)
		}
	tr_pruningwise = reorder(tr, "pruningwise")
	
	# Basic tree info
	tips = 1:length(tr_pruningwise$tip.label)
	nodes = (length(tr_pruningwise$tip.label)+1):(length(tr_pruningwise$tip.label)+tr_pruningwise$Nnode)


	
	# Read the tipranges from file, if needed.
	if (is.null(tipranges))
		{
		tipranges = getranges_from_LagrangePHYLIP(BioGeoBEARS_run_object$geogfn)
		}
	
	# Basic areas info
	areas = getareas_from_tipranges_object(tipranges)
	areas

	numareas = length(areas)
	numareas
	
	if (!is.na(results_object$inputs$max_range_size))
		{
		max_range_size = results_object$inputs$max_range_size
		} else {
		max_range_size = length(areas)
		}
	max_range_size


	if (is.null(results_object$inputs$states_list))
		{
		numstates = numstates_from_numareas(numareas=length(areas), maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range)
		numstates
		states_list_areaLetters = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range)
		#states_list
		states_list_0based_index = rcpp_areas_list_to_states_list(areas, maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range)
		#states_list_0based_index
		} else {
		states_list_0based_index = results_object$inputs$states_list
		#states_list = states_list_indexes_to_areastxt(states_list=states_list_0based_index, areanames=areas, counting_base=0, concat=FALSE, sep="")
		}


	# calculate the ML marginal probs of states at the base of each branch
	# above each split (local, non-joint probs, global model)
	# 2014-05-15_NJM: Used to add:
	# results_object$ML_marginal_prob_each_split_at_branch_bottom_BELOW_node = ML_marginal_prob_each_split_at_branch_bottom_BELOW_node / rowSums(ML_marginal_prob_each_split_at_branch_bottom_BELOW_node)
	# ... but this is now totally pointless since this is done automatically
	#results_object = get_MLsplitprobs_from_results(results_object)
	#names(results_object)

	# Extract ML parameter values, and LnL
	# This will work with optim, optimx2012, or optimx2013
	
	# Handy summary outputs
	param_ests = extract_params_from_BioGeoBEARS_results_object(results_object, returnwhat="table", addl_params=addl_params, paramsstr_digits=4)

	
	# If you want to skip the plotting and just extract
	# the parameter values
	if (juststats == TRUE)
		{
		return(param_ests)		
		} else {
		paramstr = extract_params_from_BioGeoBEARS_results_object(results_object, returnwhat="string", addl_params=addl_params, paramsstr_digits=4)
		} # if (juststats == TRUE)

	# Get the parameter names
	param_names = extract_params_from_BioGeoBEARS_results_object(results_object, returnwhat="param_names", addl_params=addl_params, paramsstr_digits=4)

		
	# Major title (short description)
	if (is.null(analysis_titletxt))
		{
		tmptxt = results_object$inputs$description
		if (any(is.null(tmptxt), tmptxt=="", tmptxt=="defaults", tmptxt=="default"))
			{
			analysis_titletxt = ""
			} else {
			analysis_titletxt = results_object$inputs$description
			}
		}
	
	
	if (is.null(dej_params_row))
		{
		# Default text for an inference or stochastic mapping
		analysis_titletxt = paste(analysis_titletxt, "\n", "ancstates: global optim, ", max_range_size, " areas max. ", paramstr, sep="")
		analysis_titletxt
		} else {
		# Text for an SSE simulation
		dej_params_row

		brate_col_TF = names(dej_params_row) == "brate"
		brate_col = (1:length(dej_params_row))[brate_col_TF]
		biogeog_params = dej_params_row[1:(brate_col-1)]
		biogeog_param_names = names(dej_params_row)[1:(brate_col-1)]
		equals_col = "="
		
		tmpcols = cbind(biogeog_param_names, equals_col, unlist(biogeog_params))
		tmpcols
		txtrows = apply(X=tmpcols, MARGIN=1, FUN=paste, sep="", collapse="")
		txtrows
		biogeog_params_txt = paste(txtrows, sep="", collapse="; ")
		biogeog_params_txt
			
		titletxt2 = bquote(paste(.(max_range_size), " areas max., ", .(biogeog_params_txt), "; ", lambda, "=", .(dej_params_row$brate), "; ",  mu, "=", .(dej_params_row$drate), "; ", alpha, "=", .(dej_params_row$rangesize_b_exponent), "; ", omega, "=", .(dej_params_row$rangesize_d_exponent), "", sep=""))
		
		#print(titletxt2)
		
		} # END if (is.null(dej_params_row))




	#######################################################
	# Get the marginal probs of the splits (global ML probs, not local)
	# (These are marginal, rather than joint probs; but not local optima)
	#######################################################
	leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr_pruningwise)

	# This gets you the prob. of each state at the left base above each node, and
	# at the right base above each node
	marprobs = results_object$ML_marginal_prob_each_state_at_branch_bottom_below_node
	left_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 2], ]
	right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 1], ]
	right_ML_marginals_by_node




	#######################################################
	# Extract the outputs ancestral states at nodes, and plot!
	#######################################################
	relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
	relprobs_matrix_for_internal_states = relprobs_matrix[nodes,]	# subset to just internal nodes
	relprobs_matrix
	
	if (is.null(states_list_0based_index))
		{
		statenames = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range, split_ABC=FALSE)
		ranges_list = as.list(statenames)
		statenames
		} else {
		ranges_list = states_list_0based_to_ranges_txt_list(state_indices_0based=states_list_0based_index, areanames=areas)
		ranges_list
		statenames = unlist(ranges_list)
		statenames
		}


	MLprobs = get_ML_probs(relprobs_matrix)
	MLprobs
	MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties=if_ties)


	# Set up colors for each state
	if (is.null(colors_list_for_states))
		{
		# Fix plot_null_range to FALSE (don't want to plot that color)
		colors_matrix = color_matrix_4areas
		colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, plot_null_range=results_object$inputs$include_null_range)
		colors_list_for_states
		} # END if (is.null(colors_list_for_states))

	# Set up colors by possible ranges
	if (is.null(ranges_list))
		{
		possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=results_object$inputs$include_null_range)
		} else {
		possible_ranges_list_txt = ranges_list
		} # if (is.null(ranges_list))
	#possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)

# 	if (plot_null_range == FALSE)
# 		{
# 		possible_ranges_list_txt[possible_ranges_list_txt == "_"] = NULL
# 		possible_ranges_list_txt[possible_ranges_list_txt == ""] = NULL
# 		possible_ranges_list_txt[possible_ranges_list_txt == "NA"] = NULL
# 		possible_ranges_list_txt[is.na(possible_ranges_list_txt)] = NULL
# 		possible_ranges_list_txt[is.null(possible_ranges_list_txt)] = NULL
# 		}

	cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

	# Legend, if desired
	if (plotlegend == TRUE)
		{
		colors_legend(possible_ranges_list_txt, colors_list_for_states, legend_ncol=legend_ncol, legend_cex=legend_cex)
		}
	
	
	# Put in a 0 for root.edge
	if (root.edge == FALSE)
		{
		tr$root.edge = 0
		} # END if (root.edge == FALSE)
	if (root.edge == TRUE)
		{
		if (is.null(tr$root.edge) == TRUE)
			{
			tr$root.edge = 0
			} # END if (is.null(tr$root.edge) == TRUE)
		} # END if (root.edge == TRUE)

	# Default label offset
	if (is.null(label.offset))
		{
		label.offset = 0.007 * (get_max_height_tree(tr) + tr$root.edge)
		}
		
	# Plot to screen
	if (show.tip.label == TRUE)
		{
		# If max_x is not specified
		if (is.null(plot_max_age))
			{
			max_x = 1.25*(get_max_height_tree(tr) + tr$root.edge)
			min_x = 0
			} else {
			nontree_part_of_x = plot_max_age - (get_max_height_tree(tr) + tr$root.edge)
			max_x = 1.25*(get_max_height_tree(tr) + tr$root.edge)
			min_x = -1 * nontree_part_of_x
			}
		} else {
		if (is.null(plot_max_age))
			{
			max_x = 1.05*(get_max_height_tree(tr) + tr$root.edge)
			min_x = 0
			} else {
			nontree_part_of_x = plot_max_age - (get_max_height_tree(tr) + tr$root.edge)
			max_x = 1.05*(get_max_height_tree(tr) + tr$root.edge)
			min_x = -1 * nontree_part_of_x
			}
		}


	###################################################
	# Calculate x-axis ticks (axisPhylo alternative)
	###################################################	
	max_tree_x = 1.0 * (get_max_height_tree(tr) + tr$root.edge)
	
	# Plot min/max
	if (is.null(xlims))
		{
		xlims = c(min_x, max_x)
		} else {
		xlims = xlims
		}
	#print(xlims)

	# Tree min/max
	nodecoords = node_coords(tr, tmplocation=cornercoords_loc, root.edge=root.edge)
	max_tree_x = max(nodecoords$x)
	
	# AxisPhylo() min/max
	if (is.null(plot_max_age))
		{
		xticks_desired_lims = c(0, max_tree_x)
		} else {
		xticks_desired_lims = c(0, plot_max_age)
		}

	xticks_desired = pretty(xticks_desired_lims)
	
	# Translate into plot coordinates
	xaxis_ticks_locs = max_tree_x - xticks_desired
	
	#print(xticks_desired)
	#print(xaxis_ticks_locs)
	###################################################	
		
	# Skip tree plotting if it has already been done:
	if (skiptree != TRUE)
		{
		#######################################################
		# Otherwise, plot the tree!!
		#######################################################
		plot(tr_pruningwise, x.lim=xlims, y.lim=ylims, show.tip.label=FALSE, label.offset=label.offset, cex=tipcex, no.margin=no.margin, root.edge=root.edge)
		
		if (show.tip.label == TRUE)
			{
			tiplabels_to_plot = sapply(X=tr_pruningwise$tip.label, FUN=substr, start=1, stop=30)
			if (skiplabels == FALSE)
				{
				tiplabels(text=tiplabels_to_plot, tip=tips, cex=tipcex, adj=0, bg="white", frame="n", pos=4, offset=label.offset, col=tipcol)	# pos=4 means labels go to the right of coords
				} # END if (skiplabels == FALSE)
			} # END if (show.tip.label == TRUE)
		
		#axisPhylo(cex.axis=titlecex)
		axis(side=1,  at=xaxis_ticks_locs, label=xticks_desired)
		# Plot the title
		mtext(text=xlab, side=1, line=2, cex=titlecex)
		}
	
	# Add states / piecharts
	if (plotwhat == "text")
		{
		# Use statecex for pie chart size at nodes AND states at tips
		par(fg=tmp_fg)	# so that user can change border externally
		
		if (skiplabels == FALSE)
			{
			nodelabels(text=MLstates[nodes], node=nodes, bg=cols_byNode[nodes], cex=statecex)		
			tiplabels(text=MLstates[tips], tip=tips, bg=cols_byNode[tips], cex=statecex)
			} # END if (skiplabels == FALSE)
		
		par(fg="black")	# set to default for most things
		}
	if (plotwhat == "pie")
		{
		# Use statecex for pie chart size at nodes BUT for states at tips,
		# use "pie_tip_statecex"
		par(fg=tmp_fg)	# so that user can change border externally

		if (skiplabels == FALSE)
			{
			# DOSIMPLIFY PIE CHARTS
			if (simplify_piecharts == TRUE)
				{
				# columns to keep in the final
				colnums_to_keep_in_probs = NULL
				
				# Probs table
				probs = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
				probs2 = probs
				maxprob = rep(0, nrow(probs))
				other = rep(0, nrow(probs))
				num_to_keep = 1
				cat("\nSince simplify_piecharts==TRUE, reducing prob pie charts to (most probable, other)...\n")
				for (i in 1:nrow(probs))
					{
					cat(i, " ", sep="")
					tmprow = probs[i,]
					positions_highest_prob_to_lowest = rev(order(tmprow))
					# If there are ties, we take the first ones
					positions_to_keep = positions_highest_prob_to_lowest[1:num_to_keep]
					colnums_to_keep_in_probs = c(colnums_to_keep_in_probs, positions_to_keep)
					keepTF = rep(FALSE, length(tmprow))
					keepTF[positions_to_keep] = TRUE
	
					# Sum the others
					otherTF = keepTF == FALSE
					other[i] = sum(tmprow[otherTF])
					tmprow[otherTF] = 0
					probs2[i,] = tmprow
					} # END for (i in 1:nrow(probs))
				cat("\n")
				
				colnums_to_keep_in_probs_in_order = sort(unique(colnums_to_keep_in_probs))
				probs3 = cbind(probs2[,colnums_to_keep_in_probs_in_order], other)
				# Subset to just internal nodes
				probs3 = probs3[nodes,]
				
				newcols = c(colors_list_for_states[colnums_to_keep_in_probs_in_order], "white")

				# DO SIMPLIFY PIE CHARTS
				nodelabels(pie=probs3, node=nodes, piecol=newcols, cex=statecex)
				} else {
				# DON'T SIMPLIFY PIE CHARTS
				nodelabels(pie=relprobs_matrix_for_internal_states, node=nodes, piecol=colors_list_for_states, cex=statecex)
				} # END if (simplify_piecharts == TRUE)

			# Plot the tiplabels, if desired
			if (tipboxes_TF == TRUE)
				{
				tiplabels(text=MLstates[tips], tip=tips, bg=cols_byNode[tips], cex=pie_tip_statecex, adj=tiplabel_adj)
				} # END if (tipboxes_TF = TRUE)
			} # END if (skiplabels == FALSE)

		par(fg="black")	# set to default for most things
		}
	
	if (skiptree != TRUE)
		{
		if (titlecex > 0)
			{
			#par(ps = 12, cex = titlecex, cex.main = titlecex)
			par(cex.main = titlecex)
			title(analysis_titletxt)
			if (!is.null(dej_params_row))
				{
				# Subtitle for SSE simulation plots
				title(titletxt2, line=1)
				#print(titletxt2)
				}
			#par("font.main") = 2
			}
		}

	if (plotsplits == TRUE)
		{
		# Error check; users must specify the location of the function "plot_phylo3_nodecoords"
		if (cornercoords_loc == "manual")
			{
			stoptxt = cat("\nNOTE: To plot splits, this function needs to access the function 'plot_phylo3_nodecoords'.\n",
							"The function is modified from an APE function, and cannot be directly included in the package,\n",
							"due to some C code that does not meet CRAN standards. To solve this, give plot_BioGeoBEARS_results\n",
							"a 'cornercoords_loc' string that gives the directory of plot_phylo3_nodecoords.R.  Typically this\n",
							"can be found via: ", 'tmp=np(system.file("extdata/a_scripts", package="BioGeoBEARS"))\n',
							"then: list.files(tmp); print(tmp)\n", sep="")
			plotsplits = FALSE
			}
		}

	if (plotsplits == TRUE)
		{
		#######################################################
		# Also add the splits to the plot
		#######################################################
		# First, get the corner coordinates
		coords = corner_coords(tr, tmplocation=cornercoords_loc, root.edge=root.edge)
		coords

		# LEFT SPLITS
		relprobs_matrix = left_ML_marginals_by_node
		
		if (plotwhat == "text")
			{
			MLprobs = get_ML_probs(relprobs_matrix)
			MLprobs
			MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties=if_ties)
			MLstates
			length(MLstates)
		
			# Set up colors
			#possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
			cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)
			
			par(fg=tmp_fg)	# so that user can change border externally
			if (skiplabels == FALSE)
				{			
				cornerlabels(text=MLstates, coords=coords$leftcorns, bg=cols_byNode, cex=splitcex)
				} # END if (skiplabels == FALSE)

			par(fg="black")	# set to default for most things
			}
		
		if (plotwhat == "pie")
			{
			par(fg=tmp_fg)	# so that user can change border externally
			cornerpies(pievals=relprobs_matrix, coords$leftcorns, piecol=colors_list_for_states, cex=splitcex)
			par(fg="black")	# set to default for most things
			}



		# RIGHT SPLITS
		relprobs_matrix = right_ML_marginals_by_node

		if (plotwhat == "text")
			{
			MLprobs = get_ML_probs(relprobs_matrix)
			MLprobs
			MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties=if_ties)
			MLstates
			length(MLstates)

			# Set up colors
			#possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
			cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

			par(fg=tmp_fg)	# so that user can change border externally
			if (skiplabels == FALSE)
				{
				cornerlabels(text=MLstates, coords=coords$rightcorns, bg=cols_byNode, cex=splitcex)
				} # END if (skiplabels == FALSE)
			par(fg="black")	# set to default for most things
			}

		if (plotwhat == "pie")
			{
			par(fg=tmp_fg)	# so that user can change border externally
			cornerpies(pievals=relprobs_matrix, coords$rightcorns, piecol=colors_list_for_states, cex=splitcex)			
			par(fg="black")	# set to default for most things
			}
		}

	# Plot vertical dashed lines for timeperiods
	# Is it time-stratified? Plot lines if desired
	if ( ((is.null(BioGeoBEARS_run_object$timeperiods) == FALSE)) && (plot_stratum_lines==TRUE) )
		{
		timeperiods = BioGeoBEARS_run_object$timeperiods
		line_positions_on_plot = add_statum_boundaries_to_phylo_plot(tr, timeperiods=timeperiods, lty="dashed", col="gray50", plotlines=TRUE)
		} # END plot vertical dashed lines for timeperiods



	# Handy summary outputs
	param_ests = matrix(data=param_ests, nrow=1)
	param_ests = adf2(param_ests)
	
	param_ests = dfnums_to_numeric(param_ests)
	names(param_ests) = c("LnL", "nparams", param_names)
	
	return(param_ests)
	}

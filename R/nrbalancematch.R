nrbalancematch <-
function(cardata.fil, trt_labs, ctrl_labs, stratify, extmatch = NA, distmat, balanceValues, nmatch = 1){

		if (!requireNamespace("optmatch", quietly = TRUE)) {
			stop("Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.")
		}

		treated = seq(from=1, length = sum(trt_labs))
		treated.names = rownames(cardata.fil)[trt_labs]

		control = seq(from=length(treated)+1, length = sum(ctrl_labs))
		control.names = rownames(cardata.fil)[ctrl_labs]

		controlLabs = seq(from=length(treated)+length(control)+1, length = length(table(stratify)))
		controlLabs.names = names(table(stratify))

		nnodes = c(treated.names, control.names, controlLabs.names)
		
		controlLab_extended = unlist(sapply(control.names, function(x) which(nnodes == stratify[x])))
		names(controlLab_extended) = c()



		startn = c(rep(treated, length(control)), control, control)
		endn = c(rep(control, rep(length(treated), length(control))), controlLab_extended, rep(length(nnodes)+1, length(control)))
		ucap = rep(1, length(startn))
		b = c(rep(nmatch, length(treated)), rep(0, length(control)), -nmatch*balanceValues, -nmatch*(length(treated)-sum(balanceValues)))

		arc.names.start = c(rep(treated.names, length(control)), control.names, control.names)
		arc.names.end = c(rep(control.names, rep(length(treated), length(control))), rep("balance", length(controlLab_extended)), rep('outflow', length(control)))

		### cost
		cost = c()
				##initialize cost on the edges
		cost = c(rep(0, (length(startn)-length(control))), rep(2*max(distmat, na.rm=T), length(control)))
				##rearrange distance matrix 
		distmat = distmat[treated.names, control.names]
		cost[1:(length(treated)*length(control))] = as.numeric(distmat)
		
		
		# Now the exact match
		if(!is.na(extmatch)){
			extmatch <- outer(treated.names, control.names, FUN = function(x, y)
						cardata.fil[x, extmatch] == cardata.fil[y, extmatch])
		
			# vapply(control.names, function(x) cardata.fil[x, 
							# extmatch] == cardata.fil[treated.names, extmatch], logical(length(treated.names)))
		
			ucap[1:(length(treated)*length(control))] = as.numeric(extmatch)
		}
	
		my.expr <- parse(text = ".Fortran(\"relaxalg\", length(b), as.integer(length(startn)), \n    \t    as.integer(startn), as.integer(endn), as.integer(cost), \n    \t    as.integer(ucap), as.integer(b), x1 = integer(length(startn)), \n    \t    crash1 = as.integer(0), large1 = as.integer(.Machine$integer.max/4), \n    \t    feasible1 = integer(1), NAOK = FALSE, DUP = TRUE, PACKAGE = \"optmatch\")")
		res <- eval(my.expr)
		
		#res = callrelax(net)
		if(res$feasible1==0){
			#return(-99)
			print('NOT FEASIBLE')
		}
		
		res = res$x1
		arc.names.start.res = arc.names.start[res==1]
		arc.names.end.res = arc.names.end[res==1]

		#	return(list(arc.names.start.res, arc.names.end.res))
			
		match.control = arc.names.start.res[arc.names.end.res == "balance"]
		match.treated = sapply(match.control, function(x){
							trts = arc.names.start.res[arc.names.end.res==x]
							trts = trts[which.min(distmat[trts, x])]
							ifelse(length(trts) > 0, trts, NA)
						}
				)

		
		
		triplets = cbind(match.treated, match.control)
		#triplets = triplets[complete.cases(triplets),]
		triplets  #colnames(triplets) = c('WithSAB', 'WOSAB')
}

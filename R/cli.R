# Defines the CLI for R/bakeoff.R
# -d defines what prewarpings to use; the three discussed in the paper are bcd
# -p defines what test problems to use; de corresponds to Fig 4; abc Fig 5; g Fig 6.

process_cli <- function(args) {
    ### Custom dimension reduction methods
    dr_ind <- which(args == '-d')
    if (length(dr_ind) > 0) {
        dr_types <- strsplit(args[dr_ind+1], '*')[[1]]
        DR_METHODS <- c()
        for (dr_type in dr_types) {
            if (dr_type == 'a') {
                DR_METHODS <- c(DR_METHODS, 'identity')
            } else if (dr_type == 'b'){ 
                DR_METHODS <- c(DR_METHODS, 'lengthscale')
            } else if (dr_type == 'c'){ 
                DR_METHODS <- c(DR_METHODS, 'lebesgue')
            } else if (dr_type == 'd'){ 
                DR_METHODS <- c(DR_METHODS, 'empirical')
            } else if (dr_type == 'e') {
                DR_METHODS <- c(DR_METHODS, 'iso_lebesgue')
            } else if (dr_type == 'f') {
                DR_METHODS <- c(DR_METHODS, 'iso_empirical')
            } else {
                stop("Unrecognized dr_type (-d argument) ])")
            }
        }
    } else {
        ## Default DR
        DR_METHODS <- c('identity', 'lebesgue')
    }

    # Handle input problems
    prob_ind <- which(args == '-p')
    ### Specify which problems. 
    if (length(prob_ind) > 0) {
        problems <- strsplit(args[prob_ind+1], '*')[[1]]
        PROBLEMS <- c()
        for (prob in problems) {
            if (prob == 'a') {
                PROBLEMS <- c(PROBLEMS, 'borehole')
            } else if (prob == 'b'){ 
                PROBLEMS <- c(PROBLEMS, 'robot arm')
            } else if (prob == 'c'){ 
                PROBLEMS <- c(PROBLEMS, 'piston')
            } else if (prob == 'd'){ 
                PROBLEMS <- c(PROBLEMS, 'CNC')
            } else if (prob == 'e'){ 
                PROBLEMS <- c(PROBLEMS, 'Temperature')
            } else if (prob == 'f'){ 
                PROBLEMS <- c(PROBLEMS, 'CT Scans')
            } else if (prob == 'g'){ 
                PROBLEMS <- c(PROBLEMS, 'MOPTA')
            } else if (prob == 'h'){ 
                stop("Not implemented yet. ")
            } else {
                stop("Unrecognized prob (-d argument) ])")
            }
        }
    } else {
        # Default problems
        PROBLEMS <- c('borehole', 'robot arm', 'piston')
    }

    return(list(DR_METHODS = DR_METHODS, PROBLEMS = PROBLEMS))
}

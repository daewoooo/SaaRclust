#' Read SaaRclust configuration file
#'
#' Read SaaRclust configuration file into a list structure. 
#' The configuration file has to be specified in INI format. 
#' R expressions can be used and will be evaluated.
#'
#' @param configfile Path to the configuration file
#' @return A \code{list} with one entry for each element in \code{configfile}.
#' @author Aaron Taudt, David Porubsky
#' @importFrom utils read.table

readConfig <- function(configfile) {

    connection <- file(configfile) 
    Lines  <- readLines(connection) 
    close(connection) 

    Lines <- chartr("[]", "==", Lines) # change section headers 
    Lines <- gsub(" ", "", Lines) # no spaces

    connection <- textConnection(Lines) 
    data <- utils::read.table(connection, as.is = TRUE, sep = "=", fill = TRUE, quote="") 
    close(connection) 
    names(data) <- c('argument','value','section')

    L <- data$argument == "" # location of section breaks 
    data$section <- data$value[which(L)[cumsum(L)]]
    data <- data[data$argument!="",]

    configlist <- list() 
    ToParse <- paste0("configlist$", data$argument, " <- ", data$value)
#   ToParse  <- paste0("configlist$", data$section, "$",  data$argument, " <- ", data$value) # with sections

    eval(parse(text=ToParse)) 

    return(configlist) 
} 

#' Write SaaRclust configuration file
#'
#' Write a SaaRclust configuration file from a list structure.
#'
#' @param config A list structure with parameter values. Each entry will be written in one line.
#' @param configfile Filename of the outputfile.
#' @return \code{NULL}
#' @author Aaron Taudt, David Porubsky
#' @importFrom utils write.table

writeConfig <- function(config, configfile) {

    ## Printing function
    formatstring <- function(string) {
    if (is.character(string) & length(string)>1) {
        string <- paste0("c('",paste0(string,collapse="','"),"')")
    } else if (is.character(string) & length(string)==1) {
        string <- paste0("'",string,"'")
    } else if (is.numeric(string) & length(string)>1) {
        string <- paste0("c(",paste0(string,collapse=','),")")
    } else if (is.numeric(string) & length(string)==1) {
        string <- string
    } else if (is.null(string)) {
        string <- "NULL"
    }
    return(string)
    }
    
    f <- file(configfile, open='w')
    cat("#============== SaaRclust configuration file ===============#\n", file=f)
    # cat("\n[General]\n", file=f)
    # for (i1 in c('numCPU','reuse.existing.files')) {
    #     cat(i1," = ",formatstring(config[[i1]]),"\n", file=f)
    # }
    cat("\n[SaaRclust]\n", file=f)
    for (i1 in c('min.mapq','min.contig.size','min.region.to.order','chromosomes','pairedEndReads','bin.size','step.size','bin.method','store.data.obj','reuse.data.obj','num.clusters','desired.num.clusters','max.cluster.length.mbp','alpha','prob.th','ord.method','assembly.fasta','concat.fasta','z.limit','remove.always.WC','mask.regions','eval.ploidy','allow.contig.cuts')) {
        cat(i1," = ",formatstring(config[[i1]]),"\n", file=f)
    }
    close(f, type='w')
}

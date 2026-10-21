# Explicit fresh burn-in-only launcher; all data remain outside the repository.
files<-Filter(Negate(is.null),lapply(sys.frames(),function(f) f$ofile))
flag<-grep('^--file=',commandArgs(),value=TRUE)
script<-if(length(files)) tail(files,1)[[1]] else sub('^--file=','',flag[1])
root<-dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE)))
source(file.path(root,'redesigned_schemes.R'))
config<-breeding_config();config$future<-0L
coordinate_breeding_experiment(root,config,breeding_scenarios()[0,,drop=FALSE])

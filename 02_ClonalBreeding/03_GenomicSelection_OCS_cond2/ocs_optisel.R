# Functions for optisel
library(optiSel)
library(ECOSolveR)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
percentage_skin <- 0.05

create_phen_df <- function(obj, breed_name = "Almond", sex = NA, herd = NA, year = NA) {
  # Slot validation
  if(!all(c("id", "ebv") %in% slotNames(obj))) {
    stop("El objeto no tiene slots 'id' y 'ebv'")
  }
  
  # Born year
  if(length(year) == 1 && is.na(year)) {
    year <- rep(NA, length(obj@id))
  }
  
  
  df_list <- list(
    Indiv = obj@id,
    Born  = 1,
    Breed = breed_name,
    BV    = as.numeric(obj@ebv),  
    Sex   = sex,
    herd  = herd
  )
  
  
  df <- as.data.frame(df_list, stringsAsFactors = FALSE)
  
  return(df)
}

export_snp_by_chromosome <- function(pop_obj, out_dir = ".") {
  # pop_obj: object like ETC1 or similar
  # out_dir: 
  
  
  obj_name <- deparse(substitute(pop_obj))
  
  # pull haplo
  geno_mat <- pullSnpHaplo(pop_obj, chr = NULL, asRaw = FALSE, simParam = NULL)
  
  # Names
  snp_names <- colnames(geno_mat)
  row_ids <- rownames(geno_mat)
  
  # IDs base 
  ind_ids <- unique(sub("_[12]$", "", row_ids))
  
  # T (rows = SNPs)
  snp_df <- data.frame(
    I = "M",
    id = snp_names,
    t(geno_mat)
  )
  
  # Rename columns and duplicate per individual
  colnames(snp_df)[3:ncol(snp_df)] <- rep(ind_ids, each = 2)
  
  # Detect chromosomes (ej: "1_1")
  chroms <- as.integer(sub("_.*", "", snp_df$id))
  unique_chroms <- sort(unique(chroms))
  
  # Vector files
  file_paths <- character(length(unique_chroms))
  
  # Export
  for(i in seq_along(unique_chroms)) {
    chr <- unique_chroms[i]
    chr_rows <- which(chroms == chr)
    df_chr <- snp_df[chr_rows, , drop = FALSE]
    
    # filename: nameObject.ChrN.txt
    file_name <- file.path(out_dir, paste0(obj_name, ".Chr", chr, ".txt"))
    
    # save
    write.table(df_chr, file_name, sep = " ", row.names = FALSE, quote = FALSE)
    
    # saver
    file_paths[i] <- normalizePath(file_name)
  }
  
  message("Files export to ", normalizePath(out_dir))
  
  return(file_paths)
}

build_map_df <- function(genMap, recombination_rate, SP = NULL) {
  
  stopifnot(
    is.list(genMap),
    length(genMap) == length(recombination_rate)
  )
  
  # ---- Build map ----
  map <- do.call(
    rbind,
    lapply(seq_along(genMap), function(chr) {
      
      M <- as.numeric(genMap[[chr]])  # Morgans
      Position <- M / recombination_rate[chr]  # bp
      
      marker_names <- names(genMap[[chr]])
      if (is.null(marker_names)) {
        marker_names <- paste0(chr, "_", seq_along(M))
      }
      
      data.frame(
        Name     = marker_names,
        Chr      = chr,
        Position = Position,
        cM       = M * 100,
        Mb       = Position / 1e6,
        stringsAsFactors = FALSE
      )
    })
  )
  
  # ---- Filter by SNP chip  ----
  if (!is.null(SP)) {
    
    loci_map <- SP$snpChips[[1]]
    
    loci_per_chr <- loci_map@lociPerChr
    loci_loc <- loci_map@lociLoc
    
    # vector for lociLoc
    chr_vec <- rep(seq_along(loci_per_chr), times = loci_per_chr)
    
    # build "Chr_Pos"
    snp_names <- paste0(chr_vec, "_", loci_loc)
    
    # Filter
    map <- map[map$Name %in% snp_names, ]
    
    # keep order
    map <- map[match(snp_names, map$Name), ]
  }
  
  rownames(map) <- NULL
  return(map)
}

#columns hermaf
matings_hermaf <- function(phen, Kin, alpha=1, ub.n=NA, max=FALSE, solver="default", pop, SIPos, ...){
  if(!("Indiv" %in% colnames(phen))){stop("Column 'Indiv' is missing in phen.")}
  if(!("Sex"   %in% colnames(phen))){stop("Column 'Sex' is missing in phen.")}
  if(!("n"     %in% colnames(phen))){stop("Column 'n' is missing in phen.")}
  
  #CHANGE
  # phen <- optiSel:::checkphen(phen, columns=c("Indiv", "Sex"), quiet=TRUE, na.Sex=FALSE) 
  # if(!is.numeric(phen$n) | any(is.na(phen$n)) | any(phen$n<0) | any(floor(phen$n)!= phen$n)){
  #   stop("column 'n' must be numeric with non-negative integer values.")
  # }
  
  #CHANGE
  # if(sum(phen$n[phen$Sex=="male"])!=sum(phen$n[phen$Sex=="female"])){
  #   stop("Total numbers of matings for males amnd females must be equal.")
  # }
  
  # if("herd" %in% colnames(phen)){
  #   phen$herd <- as.character(phen$herd)
  #   phen$herd[phen$Sex=="male"] <- NA
  #   if(any(is.na(phen$herd[phen$Sex=="female"]))){
  #     stop("Column 'herd' contains NA for some females.")
  #   }
  # }
  
  phen <- phen[phen$n>0,]
  
  #CHANGE
  # Sire <- phen$Indiv[phen$Sex=="male"]
  # Dam  <- phen$Indiv[phen$Sex=="female"]
  Sire <- phen$Indiv
  Dam  <- phen$Indiv

  
  
  image(as.matrix(Kin), axes=FALSE, col=rev(heat.colors(100)))
  axis(1, at=seq(0,1,length.out=ncol(Kin)), labels=colnames(Kin), las=2)
  axis(2, at=seq(0,1,length.out=nrow(Kin)), labels=rownames(Kin), las=2)
  title("Heatmap  Kin")
  
  Kin <- addHaploPenalty(Kin=Kin,pop=pop, SIPos= SIPos, penalty = 1e6) # For no compatible crosses
  
  # diag(Kin) <- 1e6   # STRONG PENAL for autopollination
  
  
  image(as.matrix(Kin), axes=FALSE, col=rev(heat.colors(100)))
  axis(1, at=seq(0,1,length.out=ncol(Kin)), labels=colnames(Kin), las=2)
  axis(2, at=seq(0,1,length.out=nrow(Kin)), labels=rownames(Kin), las=2)
  title("Heatmap Kin")
  
  
  if(is.null(rownames(Kin))||is.null(colnames(Kin))){
    stop("Row names and column names of Matrix 'Kin' must be the individual IDs.")
  }
  
  if(!all(Sire %in% rownames(Kin))){
    stop("All male candidates must appear in the row names of matrix Kin.")
  }
  
  if(!all(Dam %in% colnames(Kin))){
    stop("All female candidates must appear in the column names of matrix Kin.")
  }
  
  if(alpha<1){
    if("herd" %in% colnames(phen)){
      phen$herd <- as.character(phen$herd)
      ub.herd   <- alpha * tapply(phen$n, phen$herd,sum)
    }else{
      stop("Column 'herd' is needed if alpha<1.")
    }
  }
  
  Kin <- Kin[Sire, Dam]
  Zeros <- 0*Kin
  
  rhsM <- phen$n[phen$Indiv %in% Sire]
  ConM <- NULL
  for(k in 1:length(Sire)){
    Con <- Zeros
    Con[k,] <- 1
    ConM <- rbind(ConM, c(Con))
  }
  
  rhsF <- phen$n[phen$Indiv %in% Dam]
  ConF <- NULL
  for(k in 1:length(Dam)){
    Con <- Zeros
    Con[,k] <- 1
    ConF <- rbind(ConF, c(Con))
  }
  
  if(alpha<1){
    herds <- setdiff(unique(phen$herd), NA)
    rhsH  <- rep(ub.herd[as.character(herds)], each=length(Sire))
    ConH  <- NULL
    for(l in herds){
      for(k in 1:length(Sire)){
        Con <- Zeros
        Con[k,] <- phen[Dam,"herd"] %in% l  
        ConH <- rbind(ConH, c(Con))
      }
    }
  }
  
  nVar <- nrow(Kin)*ncol(Kin)
  
  G <- NULL
  h <- NULL
  
  if(identical(solver, "default")){
    G <- rbind(G, -diag(nVar))
    h <- c(h, rep(0, nVar))
  }
  
  if(identical(solver, "default") & !is.na(ub.n)){
    G <- rbind(G, diag(nVar))
    h <- c(h, rep(ub.n, nVar))
  }
  
  if(alpha<1){
    G <- rbind(G,  ConH)
    h <- c(h, rhsH)
  }
  
  A <- rbind(ConF, ConM)
  b <- c(rhsF, rhsM)
  
  if(identical(solver, "default")){
    opt <- list(...)
    if("maxit"        %in% names(opt)){opt$maxit        <- as.integer(opt$maxit)}
    if("verbose"      %in% names(opt)){opt$verbose      <- as.integer(opt$verbose)}
    if("mi_max_iters" %in% names(opt)){opt$mi_max_iters <- as.integer(opt$mi_max_iters)}
    opt <- do.call(ecos.control, opt)
    
    dims <- list(l=length(h), q=NULL, e=0L)
    #A    <- Matrix(A, sparse=TRUE)
    #G    <- Matrix(G, sparse=TRUE)
    sig  <- ifelse(max,-1,1)
    
    res    <- ECOS_csolve(c=sig*c(Kin), G=G, h=h, dims=dims, A=A, b=b, int_vars=as.integer(1:nVar), control=opt)
    Mating <- round(res$x, 0)
    info   <- res$infostring
    objval <- sum(c(Kin)*res$x)/sum(res$x)
  }else{ #use Rsymphony_solve_LP
    if(is.na(ub.n)){
      bounds <- NULL
    }else{
      bounds <- list(upper=list(ind=1:(length(Sire)*length(Dam)), val=rep(ub.n, length(Sire)*length(Dam))))
    }
    
    Dir    <- c(rep("==", length(b)), rep("<=", length(h)))
    res    <- solver(obj=c(Kin), mat=rbind(A, G), dir=Dir, rhs=c(b, h), types="I", bounds=bounds, max=max, ...)
    Mating <- res$solution
    if(res$status==0L){
      info  <- "Optimum solution found"
    }else{
      info <- "No solution found"
    }
    objval <- sum(c(Kin)*res$solution)/sum(res$solution)
  }
  
  
  Mating <- matrix(Mating, nrow=nrow(Zeros), ncol=ncol(Zeros))
  Mating <- as.data.table(Mating)
  colnames(Mating) <- Dam
  Mating$Sire <- Sire
  Matings <- melt(Mating, id.vars="Sire", variable.name="Dam", value.name="n")
  Matings$Dam <- as.character(Matings$Dam)
  Matings <- Matings[Matings$n>0,]
  setDF(Matings)
  if(alpha<1){
    herds <- phen$herd
    names(herds) <- phen$Indiv
    Matings$herd <- herds[Matings$Dam]
  }
  attributes(Matings)$objval <- objval
  attributes(Matings)$info <- info
  message(paste(info,"\n"))
  Matings
}


addHaploPenalty <- function(Kin, pop, SIPos, penalty = 1e6) {
  # convert
  Kin <- as.matrix(sKin)
  mode(Kin) <- "numeric"
  
  
  # pull haplotypes
  haplo1 <- pullMarkerHaplo(pop, markers = prms$SIPos, haplo = 1)
  haplo2 <- pullMarkerHaplo(pop, markers = prms$SIPos, haplo = 2)
  
  haplo1_id <- apply(haplo1, 1, paste, collapse = "_")
  haplo2_id <- apply(haplo2, 1, paste, collapse = "_")
  
  all_haplos <- c(haplo1_id, haplo2_id)
  
  # unique haplotypes
  unique_haplos <- unique(all_haplos)
  
  
  haplo_numbers <- setNames(seq_along(unique_haplos), unique_haplos)
  
  
  length(unique_haplos)
  
  haplo1_num <- haplo_numbers[haplo1_id]
  haplo2_num <- haplo_numbers[haplo2_id]
  
  # Extract IDs base
  indivs <- unique(sub("_.$", "", rownames(haplo1)))
  
  haplo_df <- data.frame(
    Individual = indivs,
    Haplo_1 = haplo1_num,
    Haplo_2 = haplo2_num,
    row.names = indivs
  )
  
  
  # Reorder haplo_df for Kin
  haplo_df <- haplo_df[rownames(Kin), ]
  
  haplo_list <- lapply(seq_len(nrow(haplo_df)), function(i) {
    c(haplo_df$Haplo_1[i], haplo_df$Haplo_2[i])
  })
  names(haplo_list) <- rownames(haplo_df)
  
  n <- nrow(haplo_df)
  shared_haplo <- matrix(FALSE, n, n,
                         dimnames = list(rownames(haplo_df),
                                         rownames(haplo_df)))
  
  for (i in seq_len(n)) {
    h_i <- sort(c(haplo_df$Haplo_1[i], haplo_df$Haplo_2[i]))
    for (j in i:n) {
      h_j <- sort(c(haplo_df$Haplo_1[j], haplo_df$Haplo_2[j]))
      if (all(h_i == h_j)) {   # both haplotypes coincide
        shared_haplo[i, j] <- TRUE
        shared_haplo[j, i] <- TRUE
      }
    }
  }
  
  Kin[shared_haplo] <- Kin[shared_haplo] + penalty
  
  return(Kin)
}

getImpossibleCrosses <- function(pop, SIPos) {
  # pull
  haplo1 <- pullMarkerHaplo(pop, markers = SIPos, haplo = 1)
  haplo2 <- pullMarkerHaplo(pop, markers = SIPos, haplo = 2)
  
  # convert
  haplo1_id <- apply(haplo1, 1, paste, collapse = "_")
  haplo2_id <- apply(haplo2, 1, paste, collapse = "_")
  
  # number
  all_haplos <- c(haplo1_id, haplo2_id)
  unique_haplos <- unique(all_haplos)
  haplo_numbers <- setNames(seq_along(unique_haplos), unique_haplos)
  
  haplo1_num <- haplo_numbers[haplo1_id]
  haplo2_num <- haplo_numbers[haplo2_id]
  
  # df
  indivs <- unique(sub("_.$", "", rownames(haplo1)))
  haplo_df <- data.frame(
    Individual = indivs,
    Haplo_1 = haplo1_num,
    Haplo_2 = haplo2_num,
    row.names = indivs
  )
  
  n <- nrow(haplo_df)
  impossible_crosses <- data.frame(Parent1 = character(), Parent2 = character(), stringsAsFactors = FALSE)
  
  # review pairs
  for (i in seq_len(n)) {
    h_i <- sort(c(haplo_df$Haplo_1[i], haplo_df$Haplo_2[i]))
    for (j in i:n) {
      h_j <- sort(c(haplo_df$Haplo_1[j], haplo_df$Haplo_2[j]))
      if (all(h_i == h_j)) {  # both coincides → cross impossible
        pair <- sort(c(haplo_df$Individual[i], haplo_df$Individual[j]))
        impossible_crosses <- rbind(impossible_crosses,
                                    data.frame(Parent1 = pair[1],
                                               Parent2 = pair[2]))
      }
    }
  }
  
  # remove duplicates
  impossible_crosses <- unique(impossible_crosses)
  return(impossible_crosses)
}


#############################################################




map <- build_map_df(
  genMap = genMap,
  recombination_rate = recombination_rate_test,
  SP=SP
)
cont  <- data.frame(
  age   = c(1), 
  male  = c(0.5),
  female= c(0.5))

Ne <- 50
L  <- 8   # generation interval

# 
# # 
# # 
# #########################################################################
# 
# 
# 
# phen_df <- create_phen_df(Parents, breed_name = "Almond", sex = NA, herd = NA)
# head(phen_df)
# 
# # Exportar y guardar rutas de los archivos
# files <- export_snp_by_chromosome(Parents, out_dir = ".")
# 
# # Ver rutas de los archivos generados
# files
# 
# 
# map <- build_map_df(
#   genMap = genMap,
#   recombination_rate = recombination_rate_test,
#   SP=SP
# )
# 
# cont  <- data.frame(
#   age   = c(1),
#   male  = c(0.5),
#   female= c(0.5))
# 
# 
# sKin <- segIBD(files, map)
# files <- unlist(files)
# file.remove(files)
# 
# 
# cand  <- candes(phen=phen_df, sKin=sKin, cont=cont)
# 
# 
# Ne <- 50
# L  <- 8   # intervalo generacional en años
# 
# con <- list(
#   ub.sKin = 1 - (1 - cand$mean$sKin) * (1 - 1/(2*Ne))^(1/L)
# )
# 
# # calcular contribuciones óptimas
# Offspring <- opticont("max.BV", cand, con, trace=FALSE)
# 
# class(Offspring$parent)
# 
# Offspring$info
# 
# Offspring$obj.fun
# Offspring$mean
# 
# Offspring$parent[order(Offspring$parent$oc, decreasing = TRUE), ][1:3, ]
# 
# candidates_list <- Offspring
# 
# top3_indiv <- candidates_list$parent$Indiv[
#   order(candidates_list$parent$oc, decreasing = TRUE)
# ][1:3]
# 
# Candidate <- Offspring$parent[,  c("Indiv", "Sex", "oc", "herd")]
# Candidate[Candidate$oc>0.001,]
# 
# Candidate$Sex <- "male"
# Candidate$n <- noffspring(Candidate, N=200)$nOff
# Candidate$Sex <- NA
# 
# Mating <- matings_hermaf(Candidate, Kin=sKin,pop = Parents, SIPos = prms$SIPos)
# 
# 
# # Crear pares ordenados
# Parent1 <- pmin(Mating$Sire, Mating$Dam)
# Parent2 <- pmax(Mating$Sire, Mating$Dam)
# 
# # Crear dataframe con pares ordenados
# Mating_ord <- data.frame(Parent1, Parent2, n = Mating$n)
# 
# # Sumar los n por cada par
# Mating_acum <- aggregate(n ~ Parent1 + Parent2, data = Mating_ord, sum)
# 
# Mating_acum
# 
# impossible_crosses <- getImpossibleCrosses(Parents, SIPos= prms$SIPos)
# 
# impossible_crosses$Cross <- "Impossible"
# 
# 
# # Merge directo: quedarnos con todos los cruces de Mating_acum y añadir Cross si coincide
# Mating_acum <- merge(
#   Mating_acum,
#   impossible_crosses[, c("Parent1","Parent2","Cross")],
#   by = c("Parent1","Parent2"),
#   all.x = TRUE
# )
# 
# Mating_acum
# 
# 
# 
# 
# # Sacar haplotipos
# haplo1 <- pullMarkerHaplo(Parents, markers = prms$SIPos, haplo = 1)
# haplo2 <- pullMarkerHaplo(Parents, markers = prms$SIPos, haplo = 2)
# 
# # Convertir a IDs únicos
# haplo1_id <- apply(haplo1, 1, paste, collapse = "_")
# haplo2_id <- apply(haplo2, 1, paste, collapse = "_")
# 
# # Numerar haplotipos
# all_haplos <- c(haplo1_id, haplo2_id)
# unique_haplos <- unique(all_haplos)
# haplo_numbers <- setNames(seq_along(unique_haplos), unique_haplos)
# 
# haplo1_num <- haplo_numbers[haplo1_id]
# haplo2_num <- haplo_numbers[haplo2_id]
# 
# # Crear dataframe de haplotipos
# indivs <- unique(sub("_.$", "", rownames(haplo1)))
# haplo_df <- data.frame(
#   Individual = indivs,
#   Haplo_1 = haplo1_num,
#   Haplo_2 = haplo2_num,
#   row.names = indivs
# )
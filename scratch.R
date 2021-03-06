source("R_profile_helper.R", keep.source = T)

# 
# run_profile <- function(src_data,nuc.model,gamma.model,seed,nCores){
#   set.seed(seed)
#   cat(sprintf("Start: %s_%s_%s_%s_%i_%i\n",
#               src_data$input.key,
#               nuc.model,
#               gamma.model,
#               selac_release,
#               nCores,
#               seed))
#   profile_prefix=sprintf("%s_%s_%s_%s_%i_%i",
#                          src_data$input.key,
#                          nuc.model,
#                          gamma.model,
#                          selac_release,
#                          nCores,
#                          seed)
#   model.LL=NA
#   try({
#     prof_obj <- profvis({
#       model.LL=test_selac_std(src_data$phy,
#                               src_data$codon.data,
#                               nuc.model = nuc.model,
#                               gamma.type = gamma.model,
#                               nCores = nCores)
#     }, prof_output = file.path("profile_results",paste0(profile_prefix,".Rprof")))
#     htmlwidgets::saveWidget(prof_obj, file.path("profile_results",paste0(profile_prefix,".Rprofvis.html")))
#   })
#   cat(sprintf("End: %s_%s_%s_%s_%i_%i\tLL: %0.3f\n",
#               src_data$input.key,
#               nuc.model,
#               gamma.model,
#               selac_release,
#               nCores,
#               seed,
#               model.LL))
#   if(!file.exists(paste0(src_data$input.key,"_LL_log.csv")))
#     cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n",
#         file=paste0(src_data$input.key,"_LL_log.csv"),
#         append = T)
#   cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
#               src_data$input.key,
#               nuc.model,
#               gamma.model,
#               selac_release,
#               nCores,
#               seed,
#               model.LL),
#       file=paste0(src_data$input.key,"_LL_log.csv"),
#       append = T)
#   model.LL
# }

run_standard_profiles<-function(seed=4,nCores=1,ref="v1.6.1-rc1"){
  # ref="v1.6.1-rc1"
  nuc.model=c("GTR","UNREST","GTR","GTR")
  # seed=4
  gamma.model=c("none","none","median","quadrature")
  # nCores=1
  temp=as.data.frame(expand.grid(key=1:4,s=seed))
  nuc.model=nuc.model[temp$key]
  gamma.model=gamma.model[temp$key]
  seed=temp$s
  names(nuc.model) <- paste(nuc.model,gamma.model,seed,sep="_")
  setup_selac_for_profiling(ref=ref)
  src_data=load_rokasYeast()
  # run_profile(src_data = src_data, nuc.model = nuc.model, gamma.model = gamma.model,seed=seed,nCores=nCores)
  setwd("profile_results")
  mapply(FUN=run_profile, nuc.model = nuc.model, gamma.model = gamma.model,seed=seed,
         MoreArgs = list(
           nCores=nCores,src_data = src_data),USE.NAMES = T) -> res
  print(res)
}

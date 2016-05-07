#This script looks for a profile.txt file, if it has been erased, it produces one
#so that my BatchRunPhrapl.R script doesn't throw me an error
i=0
while(i==0){
	Sys.sleep(5)
	if(!file.exists("profile.txt")){
		outputfile="profile.txt"
		cat("sample.interval=990000
\"system\" \"systemPerl\" \"PipeMS\" \"FUN\" \"<Anonymous>\" \"apply.default\" \"apply\" \"lapply\" \"unlist\" \"unique\" \"simplify2array\" \"SearchContinuousModelSpaceNLoptr\" \"doTryCatch\" \"tryCatchOne\" \"tryCatchList\" \"tryCatch\" \"try\" \"ExhaustiveSearchNLoptr\" \"DoRun\" 
\"system\" \"systemPerl\" \"PipeMS\" \"FUN\" \"<Anonymous>\" \"apply.default\" \"apply\" \"lapply\" \"unlist\" \"unique\" \"simplify2array\" \"SearchContinuousModelSpaceNLoptr\" \"doTryCatch\" \"tryCatchOne\" \"tryCatchList\" \"tryCatch\" \"try\" \"ExhaustiveSearchNLoptr\" \"DoRun\"\n",
		file=outputfile)
}
}
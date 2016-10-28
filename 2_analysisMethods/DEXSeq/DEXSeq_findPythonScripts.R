pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
if(any(list.files(pythonScriptsDir)== "dexseq_count.py")){pythonScriptsDir
  }else{stop("Cannot find DEXSeq python script", call.=FALSE)}

# this function will read in pdb files, 
AnalyzePDB <- function(pdb) { # pdb code/unique ID 
  
  library(bio3d) # first load bio3d package
  
  s <- read.pdb(pdb)
  
  # remove other chains and elety in pdb file
  s.chainA <- trim.pdb(s,chain="A",elety="CA")
  
  # save b factors for each atom as a vector
  s.b <- s.chainA$atom$b
  
  # plot the b factors vs. residue number and print the plot
  plotb3(s.b, sse=s.chainA,type="l",ylab="Bfactor")
  
  print(paste("B-factor analysis of",pdb, "PDB file is completed.")) 
  } # end function
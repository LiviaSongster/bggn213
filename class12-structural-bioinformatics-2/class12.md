Class 12 - structural bioinformatics part 2
================
Livia Songster
11/8/2019

# Obtaining and inspecting our input structure

``` r
library(bio3d)
name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
hiv <- read.pdb(name)
hiv
```

    ## 
    ##  Call:  read.pdb(file = name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
# trim to make protein-only and ligand-only objects
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

# save the objects to my directory
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")

# check them in R
library(bio3d.view)
library(rgl)
view(prot)
```

    ## Computing connectivity from coordinates...

``` r
#rglwidget()

view(lig)
```

    ## Computing connectivity from coordinates...

``` r
#rglwidget()
```

# Inspecting my docking results

I added hydrogens and charges to the protein and ligand using
AutoDocTools from MGL, and then ran the protein and ligand through
Autodock Vina.

``` r
# convert .pdbqt to .pdb
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

Next I opened the results of the docking simulation and the original
protein and ligand crystal structure in VMD and rendered .bmp files for
each docking alignment (0-17 frames). I imported the image sequence into
ImageJ and save the result as a gif (which is pushed to github).

Most of the simulated docks were not very close to the crystal
structure, but the last one (17) looked pretty close to the crystal
structure.

# Assessing docking results quantitatively

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("ligand.pdbqt")
# calculate the RMSD for the results of the docking vs. ligand alone
rmsd(ori, res)
```

    ##  [1]  0.649  4.206 11.110 10.529  4.840 10.932 10.993  3.655 10.996 11.222
    ## [11] 10.567 10.372 11.019 11.338  8.390  9.063  8.254  8.978

``` r
which.min(rmsd(ori, res))
```

    ## [1] 1

``` r
mean(rmsd(ori, res))
```

    ## [1] 8.7285

The first docking simulation has the smallest RMSD, which is 0.649
angstroms from the crystal structure. The average RMSD is 8.73
angstroms.

Next let’s determine the RMSD for heavy atoms only.

``` r
# select atoms that are not hydrogen
res.noh <- atom.select(res,"noh",value=TRUE)
ori.noh <- atom.select(ori,"noh",value=TRUE)
# calculate RMSD
rmsd(ori.noh, res.noh)
```

    ##  [1]  0.506  4.310 11.022 10.359  4.781 10.956 10.918  3.704 10.905 10.994
    ## [11] 10.432 10.328 10.846 11.208  8.324  8.935  8.272  8.870

``` r
which.min(rmsd(ori.noh, res.noh))
```

    ## [1] 1

``` r
mean(rmsd(ori.noh, res.noh))
```

    ## [1] 8.648333

The first simulation still has the lowest RMSD, and the mean RMSD
without heavy atoms (no hydrogen) is 8.65 angstroms.

# Normal Mode Analysis (NMA)

Normal mode analysis (NMA) is one of the major simulation techniques
used to probe largescale motions in biomolecules. Typical application is
for the prediction of functional motions in proteins.

Normal mode analysis (NMA) of a single protein structure can be carried
out by providing a PDB object to the function nma(). In the code below
we first load the Bio3D package and then download an example structure
of hen egg white lysozyme (PDB id 1hel) with the function read.pdb().
Finally the function nma() is used perform the normal mode calculation:

``` r
library(bio3d)
# download hen egg white lysozyme
pdb <- read.pdb("1HEL")
```

    ##   Note: Accessing on-line PDB file

``` r
# Calculate normal mode analysis
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.02 seconds.
    ##  Diagonalizing Hessian...    Done in 0.07 seconds.

``` r
# plot the NMA results, where the pdb file provides sum of squared estimate of errors (SSE)
par(mar=c(5.1,4.1,4.1,2.1))
plot(modes, sse=pdb)
```

![](class12_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# visualize NMA results in VMD
# make trajectories of atomic displacements from NMA:
mktrj(modes, mode=7, file="nma_7.pdb")
```

I looked at the results in VMD and saved a gif in this folder.

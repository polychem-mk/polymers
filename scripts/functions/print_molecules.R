#              Print molecules     

# required packages:
# library(rcdk)

# arguments:
#             smiles    character, SMILES (representation of the structure of the molecule);
# prefix, row_labels    to make row labels, use one 'prefix' and 'row_labels' (may be a 
#                       character or numeric) ;
#             layout    numeric, how many molecules are printed in each row.

print_molecules <- function(smiles = rep("C(Cl)CCC(c1ccccc1)CCC(Cl)C", 3),
                            row_labels = NA, 
                            prefix = "",
                            layout = 3)
  {
  # set image dimensions
  depictor = get.depictor(width = 300,  height = 300,  zoom = 3)
  
  # set layout
  if(layout != 1){
    par(mfrow = c( ceiling(length(smiles)/layout),layout), mar = rep(0.3, 4) )
  }
  
  for(i in 1:length(smiles)){
    # Parse SMILES into molecules
    molecule <- parse.smiles(smiles[i]) 
    
    # get 2d coordinates and convert it to the raster object
    rstr = view.image.2d(molecule[[1]], depictor =depictor) 
    rstr = as.raster(rstr)
    
    plot(rstr)  # plot molecule
    
    # plot row label
    if(!sum(is.na(row_labels)) & i%%3 == 1){
      text(x = 5, y = 280, cex = 0.8,
           labels = paste( prefix, row_labels[ceiling(i/3)])) 
    }
    
    # add horizontal line between rows
    if(!sum(is.na(row_labels))){
      abline(0,0, col = "#E0E5E5")
    }
  }
  
  if(layout != 1){
    # set layout back to one plot
    par(mfrow = c(1,1), mar = rep(2, 4)) 
  }
}


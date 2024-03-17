module Defs

export NDIMS, CYCLIC, MAX_BOND_J, MAX_P, BINDING_RANGE, WORM_RANGE

const NDIMS = 3
#a CYCLIC list used for the curl function
const CYCLIC = [1 2 3 1 2 3] #max bond value
const MAX_BOND_J = 100 #current variable
const MAX_P = 50
const BINDING_RANGE = 1:9
const WORM_RANGE = 1:(2*NDIMS)

end

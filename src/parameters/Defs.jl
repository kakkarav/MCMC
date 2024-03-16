module Defs

export NDIMS, cyclic, MAX_BOND_J, MAX_P, odd, binding_range, worm_range

const NDIMS = 3
#a cyclic list used for the curl function
const cyclic = [1 2 3 1 2 3] #max bond value
const MAX_BOND_J = 100 #current variable
const MAX_P = 50
#Levi Civita symbol: return is a missing index from 1:3, e.g. odd[1,2] =3, odd[2,3] = 1
const odd = [0 3 2; 3 0 1; 2 1 0]
const binding_range = 1:9
const worm_range = 1:(2*NDIMS)

end

# Non dimentional Parameters
Re=10
lam_star=10
h_star=10

# Mesh parameters
Nelx=10
Nely=$(bc <<< "$h_star*$Nelx/$lam_star")

# Outdated, I think n will always be 1
n=1

# Locked in parameters
lam=10

# Calculated geometry/flow parameters
h=$(bc <<< "scale=10; $lam*($h_star - 1)/(2*$lam_star)")
wave_amp=$(bc <<< "scale=10; $lam/$lam_star")
U_b=$(bc <<< "scale=10; $Re/$h")




# This is the initial meshing stuff
L=$(bc <<< "scale=10; $n*$lam")
BigH=$(bc <<< "scale=10; 2*$h + $wave_amp")
echo $BigH

echo $L  
sed -i "9s/.*/-$Nelx  -$Nely                                    Nelx Nely/" Test.box
sed -i "10s/.*/0.0  $L  1.0                              x0 x1 ratio/" Test.box
sed -i "11s/.*/0.0  $BigH  1.0                              y0 y1 ratio/" Test.box

# This is the rebuild the mesh
genbox << EOF
Test.box
EOF
#
mv box.re2 channel.re2

genmap << EOF
channel

EOF

# Now we need to edit the .usr file
sed -i "/HarryIsTheBomb/c\	lam = $lam ! HarryIsTheBomb" channel.usr 
sed -i "/HarryIsTheBigDog/c\	wave_amp = $wave_amp ! HarryIsTheBigDog" channel.usr
sed -i "/BigDiggityHDizzle/c\	h = $h ! BigDiggityHDizzle" channel.usr
sed -i "/HarryDancesLikeAGod/c\	param(55) = $U_b ! HarryDancesLikeAGod" channel.usr
sed -i "/#define NUMBER_ELEMENTS_X/c\#define NUMBER_ELEMENTS_X $Nelx" channel.usr
sed -i "/#define NUMBER_ELEMENTS_Y/c\#define NUMBER_ELEMENTS_Y $Nely" channel.usr
# And the SIZE file


./compile.sh --all

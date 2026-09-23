#!/bin/bash

iter=6      # reference model m4
iter2=88    # max / m4
iter3=77    # m4 / min


# ============================================================
# Plot settings
# ============================================================

Hrange=-4/4
Vrange=-6/6
AreaSize=8/12

depths=(0 0.2 0.5 0.8 1.1 1.5 2.5)

ncols=${#depths[@]}
nrows=3


# ============================================================
# Colour palettes
# ============================================================

gmt makecpt -Cjet -Z -I -T0/4/0.2 > Res.cpt

# Same colour scale for spread below and spread above
gmt makecpt -Chot -Z -T0/2/0.1 > spread.cpt


# ============================================================
# GMT settings
# ============================================================

gmt set FONT_ANNOT_PRIMARY 14p,Helvetica,black
gmt set FONT_ANNOT_SECONDARY 12p,Helvetica,black
gmt set FONT_HEADING 14p,Helvetica,black
gmt set FONT_LABEL 14p,Helvetica,black

gmt set MAP_FRAME_WIDTH 0.1c
gmt set MAP_TICK_LENGTH_PRIMARY 0.05c
gmt set MAP_ANNOT_OFFSET_PRIMARY 0.2c
gmt set COLOR_NAN 255/255/255
gmt set MAP_FRAME_TYPE PLAIN
gmt set MAP_FRAME_PEN thin
gmt set MAP_LABEL_OFFSET 5p


# ============================================================
# Start GMT session
# ============================================================

gmt begin Model_spread_7maps_0-2.5kmdepth pdf

gmt subplot begin ${nrows}x${ncols} \
    -Fs$AreaSize \
    -M0.7c/0.4c \
    -R$Hrange/$Vrange \
    -JX$AreaSize


i=0


for depth in "${depths[@]}"; do

    # ========================================================
    # Define subplot numbers
    #
    # Row 1 : 0 ... ncols-1
    # Row 2 : ncols ... 2*ncols-1
    # Row 3 : 2*ncols ... 3*ncols-1
    # ========================================================

    ref_panel=$i
    below_panel=$((i + ncols))
    above_panel=$((i + 2 * ncols))

    ref_label=$((ref_panel + 1))
    below_label=$((below_panel + 1))
    above_label=$((above_panel + 1))


    echo
    echo "======================================================"
    echo "Depth = ${depth} km"
    echo "Reference panel     = ${ref_panel}"
    echo "Spread-below panel  = ${below_panel}"
    echo "Spread-above panel  = ${above_panel}"
    echo "======================================================"


    # ========================================================
    # REFERENCE RESISTIVITY SLICE
    # ========================================================

    cat << EOF > param_H_${depth}km.dat
0
${iter}
1
0 0 $depth
0.0
1
0
EOF

    makeCutawayForGMT param_H_${depth}km.dat

    mv \
        resistivity_GMT_iter${iter}.dat \
        resistivity_GMT_iter${iter}_${depth}km.dat

    in_file="resistivity_GMT_iter${iter}_${depth}km.dat"


    # ========================================================
    # SPREAD BELOW
    #
    # m4 / minimum
    # ========================================================

    cat << EOF > param_H_${depth}km_spreadbelow.dat
0
${iter3}
1
0 0 $depth
0.0
1
0
EOF

    makeCutawayForGMT \
        param_H_${depth}km_spreadbelow.dat

    mv \
        resistivity_GMT_iter${iter3}.dat \
        resistivity_GMT_iter${iter3}_${depth}km.dat

    in_file_below="resistivity_GMT_iter${iter3}_${depth}km.dat"


    # ========================================================
    # SPREAD ABOVE
    #
    # maximum / m4
    # ========================================================

    cat << EOF > param_H_${depth}km_spreadabove.dat
0
${iter2}
1
0 0 $depth
0.0
1
0
EOF

    makeCutawayForGMT \
        param_H_${depth}km_spreadabove.dat

    mv \
        resistivity_GMT_iter${iter2}.dat \
        resistivity_GMT_iter${iter2}_${depth}km.dat

    in_file_above="resistivity_GMT_iter${iter2}_${depth}km.dat"



    # ========================================================
    # ROW 1: REFERENCE RESISTIVITY
    # ========================================================

    grep Z "$in_file" |
        awk '{
            z=$3;
            if (z < 0) z=0;
            else if (z > 4) z=4;
            print z
        }' > value.txt

    awk '{
        if ($1 == ">")
            print $0;
        else
            printf "%15.6e %15.6e\n", $1, $2
    }' "$in_file" > tmp.txt


    gmt subplot set $ref_panel


    gmt plot tmp.txt \
        -CRes.cpt \
        -G+z \
        -Zvalue.txt \
        -Bxa1f1+l"y (km)" \
        -Bya1f1+l"x (km)" \
        -BWSne \
        -L \
        -V0


    # Stations
    awk '{print $2, $1}' loc.txt |
        gmt plot -Gblack -Ss5p

    # Vuache fault
    gmt plot vuache_trace_modelspace.dat \
        -W2p,red

    # Profiles
    gmt plot profiles_modelspace.dat


    echo "-3.8 5.6 ${ref_label}) Res. map at ${depth} km bsl" |
        gmt text \
            -F+f14p,Helvetica-Bold,black+jTL \
            -N



    # ========================================================
    # ROW 2: SPREAD BELOW
    #
    # m4 / minimum
    # ========================================================

    grep Z "$in_file_below" |
        awk '{
            z=$3;
            if (z < 0) z=0;
            else if (z > 2) z=2;
            print z
        }' > value.txt

    awk '{
        if ($1 == ">")
            print $0;
        else
            printf "%15.6e %15.6e\n", $1, $2
    }' "$in_file_below" > tmp.txt


    gmt subplot set $below_panel


    gmt plot tmp.txt \
        -Cspread.cpt \
        -G+z \
        -Zvalue.txt \
        -Bxa1f1+l"y (km)" \
        -Bya1f1+l"x (km)" \
        -BWSne \
        -L \
        -V0


    # Stations
    awk '{print $2, $1}' loc.txt |
        gmt plot -Gblack -Ss5p

    # Vuache fault
    gmt plot vuache_trace_modelspace.dat \
        -W2p,red

    # Profiles
    gmt plot profiles_modelspace.dat


    echo "-3.8 5.6 ${below_label}) m4 / Min. res. at ${depth} km bsl" |
        gmt text \
            -F+f14p,Helvetica-Bold,white+jTL \
            -N



    # ========================================================
    # ROW 3: SPREAD ABOVE
    #
    # maximum / m4
    # ========================================================

    grep Z "$in_file_above" |
        awk '{
            z=$3;
            if (z < 0) z=0;
            else if (z > 2) z=2;
            print z
        }' > value.txt

    awk '{
        if ($1 == ">")
            print $0;
        else
            printf "%15.6e %15.6e\n", $1, $2
    }' "$in_file_above" > tmp.txt


    gmt subplot set $above_panel


    gmt plot tmp.txt \
        -Cspread.cpt \
        -G+z \
        -Zvalue.txt \
        -Bxa1f1+l"y (km)" \
        -Bya1f1+l"x (km)" \
        -BWSne \
        -L \
        -V0


    # Stations
    awk '{print $2, $1}' loc.txt |
        gmt plot -Gblack -Ss5p

    # Vuache fault
    gmt plot vuache_trace_modelspace.dat \
        -W2p,red

    # Profiles
    gmt plot profiles_modelspace.dat


    echo "-3.8 5.6 ${above_label}) Max. res. / m4 at ${depth} km bsl" |
        gmt text \
            -F+f14p,Helvetica-Bold,white+jTL \
            -N


    i=$((i + 1))

done


# ============================================================
# TWO COLOUR BARS ONLY
#
# Anchor them to the bottom-left panel, which is:
#
#     panel = 2 * ncols
#
# For 7 columns this is panel 14.
# ============================================================

bottom_left_panel=$((2 * ncols))

gmt subplot set $bottom_left_panel


# ------------------------------------------------------------
# Resistivity scale
# ------------------------------------------------------------

gmt colorbar \
    -CRes.cpt \
    -Dx0c/-1.2c+w7c/0.4c+h+jTL \
    -Bxa1f1+l"Resistivity [log(@~W@~m)]" \
    -N


# ------------------------------------------------------------
# Shared spread scale
#
# Used for both:
#     m4 / minimum
#     maximum / m4
# ------------------------------------------------------------

gmt colorbar \
    -Cspread.cpt \
    -Dx0c/-3.0c+w7c/0.4c+h+jTL \
    -Bxa0.5f0.1+l"Resistivity spread [log ratio]" \
    -N


gmt subplot end

gmt end

echo "All plots completed."
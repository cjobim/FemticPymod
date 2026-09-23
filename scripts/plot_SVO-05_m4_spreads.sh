#!/bin/bash

# ============================================================
# Models
# ============================================================

iter_ref=6      # reference resistivity model
iter_up=88       # upward spread = max / m4
iter_down=77     # downward spread = m4 / min

cut=SVO-05


# ============================================================
# Plot region (km)
# ============================================================

Hrange=-4.79/4.79
Vrange=-3/2

Xspan=9.58
Zspan=5.0


# ============================================================
# Figure geometry
# Same physical scale horizontally and vertically
# ============================================================

Wbig=14.0

Hbig=$(awk -v w="$Wbig" -v xs="$Xspan" -v zs="$Zspan" \
    'BEGIN {printf "%.4f", w*zs/xs}')

Wsmall=6.8
Xgap=0.4

Hsmall=$(awk -v w="$Wsmall" -v xs="$Xspan" -v zs="$Zspan" \
    'BEGIN {printf "%.4f", w*zs/xs}')

Xstep=$(awk -v w="$Wsmall" -v g="$Xgap" \
    'BEGIN {printf "%.4f", w+g}')

# Gap between large panel and lower panels
Vgap=0.7

Ystep=$(awk -v h="$Hsmall" -v g="$Vgap" \
    'BEGIN {printf "%.4f", h+g}')

# Initial page origin
X0=2.0
Y0=11.0

# Colour bars
CbarW=5.8
CbarH=0.35

# Vertical shift below lower plots
CbarYdrop=1.35


echo "Main panel   : ${Wbig} x ${Hbig} cm"
echo "Small panels : ${Wsmall} x ${Hsmall} cm"


# ============================================================
# Colour palettes
# ============================================================

gmt makecpt -Cjet -Z -I -T0/4/0.2 > Res.cpt
gmt makecpt -Chot -Z -T0/2/0.1 > spread.cpt
gmt makecpt -Cred2green -Z -T0/15/0.2 > days.cpt


# ============================================================
# GMT settings
# ============================================================

gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,black
gmt set FONT_ANNOT_SECONDARY 10p,Helvetica,black
gmt set FONT_HEADING 12p,Helvetica,black
gmt set FONT_LABEL 10p,Helvetica,black

gmt set MAP_FRAME_WIDTH 0.1c
gmt set MAP_TICK_LENGTH_PRIMARY 0.1c
gmt set MAP_ANNOT_OFFSET_PRIMARY 0.2c
gmt set COLOR_NAN 255/255/255
gmt set MAP_FRAME_TYPE PLAIN
gmt set MAP_FRAME_PEN thin
gmt set MAP_LABEL_OFFSET 5p


# ============================================================
# Functions
# ============================================================

make_slice () {

    local iter=$1
    local outtag=$2
    local tmp_param="param_V_${cut}_${outtag}.tmp"

    awk -v iter="$iter" '
        NR == 2 {
            print iter
            next
        }
        {
            print
        }
    ' "param_V_${cut}.dat" > "$tmp_param"

    makeCutawayForGMT "$tmp_param"

    mv \
        "resistivity_GMT_iter${iter}.dat" \
        "${outtag}.dat"

    rm -f "$tmp_param"
}


prepare_plot_files () {

    local infile=$1
    local zmax=$2
    local xyfile=$3
    local zfile=$4

    grep Z "$infile" | \
        awk -v zmax="$zmax" '
        {
            z=$3

            if (z < 0)
                z=0

            else if (z > zmax)
                z=zmax

            print z
        }
    ' > "$zfile"


    awk '
        {
            if ($1 == ">")
                print $0

            else
                printf "%15.6e %15.6e\n", $1, -$2
        }
    ' "$infile" > "$xyfile"
}


plot_common_overlays_large () {

    # Stations
    awk -F',' 'NR>1 {print $2, $3}' \
        "${cut}_stations.csv" | \
        gmt plot \
            -Sc0.10c \
            -Gblack \
            -W0.25p


    # Station names
    awk -F',' \
        'NR>1 {
            printf "%.6f %.6f %s\n",
            $2, $3+0.10, $1
        }' \
        "${cut}_stations.csv" | \
        gmt text \
            -F+f6p,Helvetica,black+jTL+a90 \
            -N


    # Profile name
    echo "0 1.75 ${cut}" | \
        gmt text \
            -F+f16p,Helvetica,black


    # NW / SE
    echo "-4.7 1.85 NW" | \
        gmt text \
            -F+f12p,Helvetica,black+jTL \
            -N

    echo "4.7 1.85 SE" | \
        gmt text \
            -F+f12p,Helvetica,black+jTR \
            -N
}


plot_common_overlays_small () {

    # Stations
    awk -F',' 'NR>1 {print $2, $3}' \
        "${cut}_stations.csv" | \
        gmt plot \
            -Sc0.08c \
            -Gblack \
            -W0.20p


    # Station names
    awk -F',' \
        'NR>1 {
            printf "%.6f %.6f %s\n",
            $2, $3+0.08, $1
        }' \
        "${cut}_stations.csv" | \
        gmt text \
            -F+f5p,Helvetica,black+jTL+a90 \
            -N


    # Profile name
    echo "0 1.75 ${cut}" | \
        gmt text \
            -F+f11p,Helvetica,black


    # NW / SE
    echo "-4.7 1.85 NW" | \
        gmt text \
            -F+f9p,Helvetica,black+jTL \
            -N

    echo "4.7 1.85 SE" | \
        gmt text \
            -F+f9p,Helvetica,black+jTR \
            -N
}


# ============================================================
# Build slices
# ============================================================

make_slice "$iter_ref"  "ref_model"
make_slice "$iter_down" "down_spread"
make_slice "$iter_up"   "up_spread"


# ============================================================
# Start GMT
# ============================================================

gmt begin Model_spread_profile_${cut}_reformatted pdf


# ============================================================
# 1) LARGE RESISTIVITY PANEL
# ============================================================

prepare_plot_files \
    "ref_model.dat" \
    4 \
    "tmp_ref_xy.txt" \
    "tmp_ref_z.txt"


gmt plot "tmp_ref_xy.txt" \
    -R$Hrange/$Vrange \
    -JX${Wbig}c/${Hbig}c \
    -X${X0}c \
    -Y${Y0}c \
    -CRes.cpt \
    -G+z \
    -Ztmp_ref_z.txt \
    -Bxa1f1+l"distance along profile (km)" \
    -Bya1f1+l"z (km)" \
    -BWSne \
    -L \
    -V0


plot_common_overlays_large


# ------------------------------------------------------------
# Earthquake hypocentres
#
# columns:
# x z dx dy magnitude days
# ------------------------------------------------------------

awk -F',' \
    'NR>1 {
        print $2,
              $3,
              $4,
              $5,
              $6*0.05+0.05,
              $7
    }' \
    "${cut}_EQ.csv" > EQ_plot.txt

dos2unix EQ_plot.txt 2>/dev/null


# Colour = days after mainshock
# Size   = magnitude
gmt plot EQ_plot.txt \
    -i0,1,5,4 \
    -Sc \
    -Cdays.cpt \
    -W0.25p,black


# Uncomment to restore error bars
# gmt plot EQ_plot.txt \
#     -i0,1,2,3 \
#     -E+cap \
#     -W0.25p,black


# Panel title
echo "-4.65 2.20 1) Resistivity (alpha = 1.5)" | \
    gmt text \
        -F+f10p,Helvetica-Bold,black+jTL \
        -N


# ============================================================
# 2) LOWER LEFT — DOWNWARD SPREAD
# ============================================================

prepare_plot_files \
    "down_spread.dat" \
    2 \
    "tmp_down_xy.txt" \
    "tmp_down_z.txt"


gmt plot "tmp_down_xy.txt" \
    -R$Hrange/$Vrange \
    -JX${Wsmall}c/${Hsmall}c \
    -X0c \
    -Y-${Ystep}c \
    -Cspread.cpt \
    -G+z \
    -Ztmp_down_z.txt \
    -Bxa1f1+l"distance along profile (km)" \
    -Bya1f1+l"z (km)" \
    -BWSne \
    -L \
    -V0


plot_common_overlays_small


echo "-4.65 2.20 2) Downward spread (m4 / min)" | \
    gmt text \
        -F+f9p,Helvetica-Bold,black+jTL \
        -N


# ============================================================
# 3) LOWER RIGHT — UPWARD SPREAD
# ============================================================

prepare_plot_files \
    "up_spread.dat" \
    2 \
    "tmp_up_xy.txt" \
    "tmp_up_z.txt"


gmt plot "tmp_up_xy.txt" \
    -R$Hrange/$Vrange \
    -JX${Wsmall}c/${Hsmall}c \
    -X${Xstep}c \
    -Y0c \
    -Cspread.cpt \
    -G+z \
    -Ztmp_up_z.txt \
    -Bxa1f1+l"distance along profile (km)" \
    -Bya1f1+l"z (km)" \
    -BWSne \
    -L \
    -V0


plot_common_overlays_small


echo "-4.65 2.20 3) Upward spread (max / m4)" | \
    gmt text \
        -F+f9p,Helvetica-Bold,black+jTL \
        -N


# ============================================================
# COLOUR BARS
#
# No surrounding basemap/frame.
# We only use origin shifts.
# ============================================================


# ------------------------------------------------------------
# Move back beneath lower-left plot
# ------------------------------------------------------------

gmt basemap \
    -R0/1/0/1 \
    -JX${Wsmall}c/0.1c \
    -X-${Xstep}c \
    -Y-${CbarYdrop}c \
    -B+n


# Resistivity colour bar
gmt colorbar \
    -CRes.cpt \
    -Dx3.4c/0c+w${CbarW}c/${CbarH}c+h+jTC \
    -Bxa1f0.2+l"Resistivity [log(@~W@~m)]" \
    --FONT_LABEL=10p,Helvetica,black \
    --FONT_ANNOT_PRIMARY=9p,Helvetica,black


# ------------------------------------------------------------
# Move below lower-right plot
# ------------------------------------------------------------

gmt basemap \
    -R0/1/0/1 \
    -JX${Wsmall}c/0.1c \
    -X${Xstep}c \
    -Y0c \
    -B+n


# Spread colour bar
gmt colorbar \
    -Cspread.cpt \
    -Dx3.4c/0c+w${CbarW}c/${CbarH}c+h+jTC \
    -Bxa0.5f0.1+l"Resistivity spread [log ratio]" \
    --FONT_LABEL=10p,Helvetica,black \
    --FONT_ANNOT_PRIMARY=9p,Helvetica,black


gmt end

echo "All plots completed."

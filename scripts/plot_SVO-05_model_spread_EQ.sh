#!/bin/bash

# Plot settings
Hrange=-4.79/4.79
Vrange=-5/3
AreaSize=7/5.84
ah=1
fh=1
av=1
fv=1
cut=SVO-05

# Prepare color palette
gmt makecpt -Cjet -Z -I -T0/4/0.2 > Res.cpt
gmt makecpt -Chot -Z -T0/2/0.1 > spread.cpt
# Create CPT for magnitude
gmt makecpt -Cred2green -Z -T0/15/0.2 > days.cpt

# --- draw color section ---
#gmt set ANOT_FONT Helvetica
gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,black
gmt set FONT_ANNOT_SECONDARY 10p,Helvetica,black
#gmt set HEADER_FONT Helvetica
gmt set FONT_HEADING 12p,Helvetica,black
#gmt set LABEL_FONT Helvetica
gmt set FONT_LABEL 10p,Helvetica,black
#gmt set ANOT_FONT_SIZE 20pt HEADER_FONT_SIZE 14pt LABEL_FONT_SIZE 20pt
gmt set MAP_FRAME_WIDTH 0.1c
gmt set MAP_TICK_LENGTH_PRIMARY 0.1c 
gmt set MAP_ANNOT_OFFSET_PRIMARY 0.2c
gmt set COLOR_NAN 255/255/255
gmt set MAP_FRAME_TYPE PLAIN
gmt set MAP_FRAME_PEN thin
gmt set MAP_LABEL_OFFSET 5p

# Start GMT session
gmt begin Model_spread_profile_${cut}_alpha1.5_EQ pdf
gmt subplot begin 2x2 -Fs$AreaSize -M0c/0c -R$Hrange/$Vrange -JX$AreaSize

iter=13		#best model, constant_start_19
iter2=9999	#9thdecile/1stdecile from constant_start ensemble

#plot resistivity

#replace value of iter in the param file:
awk -v iter="$iter" 'NR == 2 {print iter; next} {print}' param_V_${cut}.dat > tmp_param.dat
mv tmp_param.dat param_V_${cut}.dat

makeCutawayForGMT param_V_${cut}.dat

in_file="resistivity_GMT_iter${iter}.dat"
ps_data="tmp.txt"
value_file="value.txt"

grep Z $in_file | awk '{z=$3; if (z < 0) z=0; else if (z > 4) z=4; print z}' > "$value_file"
awk '{if ($1 == ">") print $0; else printf "%15.6e %15.6e\n", $1, -$2 }' $in_file > "$ps_data"

gmt subplot set 0

	gmt plot "$ps_data" -CRes.cpt -G+z -Z"$value_file" -Bxa1f1+l"distance along profile (km)" -Bya1f1+l"z (km)" -BWSne -L -V0
	
	awk -F',' 'NR>1 {print $2, $3}' ${cut}_stations.csv | \
	gmt plot -Sc0.1c -Gblack -W0.25p
	
	awk -F',' 'NR>1 {printf "%.6f %.6f %s\n", $2, $3+0.1, $1}' ${cut}_stations.csv | \
	gmt pstext -F+f6p,Helvetica,black+jTL+a90 -N
	
	# --- Earthquake hypocenters ---
	#awk -F',' 'NR>1 {print $2, $3, $4, $5, $6}' ${cut}_EQ.csv > EQ_plot.txt
	awk -F',' 'NR>1 {print $2, $3, $4, $5, $6*0.05+0.05, $7}' ${cut}_EQ.csv > EQ_plot.txt
	dos2unix EQ_plot.txt
	# --- Plot hypocenters ---
	# File order: x  z  dx  dy  mag  days
	# -E+x+y+cap   : error bars (horizontal & vertical) with caps
	# -S           : circle symbol, radius from magnitude column
	# -Cdays.cpt   : color from "days" column
	# -W0.25p,black: black outline for circles and error bars
	gmt plot EQ_plot.txt -i0,1,5,4 -Sc -Cdays.cpt -W0.25p,black
	#gmt plot EQ_plot.txt -i0,1,2,3 -E+cap -W0.25p,black 

	#plot model label
	label_num=1
	echo "-4.7 3.7 ${label_num}) Resistivity (alpha=1.5) " | gmt pstext -F+f10p,Helvetica-Bold,black+jTL -N
	#plot title
    echo "0 2.6 ${cut}" | gmt pstext -F+f16p,Helvetica,black -N
    #plot NW
	echo "-4.7 2.8 NW" | gmt pstext -F+f12p,Helvetica,black+jTL -N
    #plot SE
	echo "4.7 2.8 SE" | gmt pstext -F+f12p,Helvetica,black+jTR -N
	#plot RMS
    #echo "-4.7 2.2 RMS=${rms}" | gmt pstext -F+f9p,Helvetica,black+jTL -N


#plot spread
iter=9999	#9thdecile/1stdecile
#replace value of iter in the param file:
awk -v iter="$iter" 'NR == 2 {print iter; next} {print}' param_V_${cut}.dat > tmp_param.dat
mv tmp_param.dat param_V_${cut}.dat

makeCutawayForGMT param_V_${cut}.dat

in_file="resistivity_GMT_iter${iter}.dat"
ps_data="tmp.txt"
value_file="value.txt"

grep Z $in_file | awk '{z=$3; if (z < 0) z=0; else if (z > 4) z=4; print z}' > "$value_file"
awk '{if ($1 == ">") print $0; else printf "%15.6e %15.6e\n", $1, -$2 }' $in_file > "$ps_data"

gmt subplot set 1

	gmt plot "$ps_data" -Cspread.cpt -G+z -Z"$value_file" -Bxa1f1+l"distance along profile (km)" -Byf1 -BWSne -L -V0
	
	awk -F',' 'NR>1 {print $2, $3}' ${cut}_stations.csv | \
	gmt plot -Sc0.1c -Gblack -W0.25p
	
	awk -F',' 'NR>1 {printf "%.6f %.6f %s\n", $2, $3+0.1, $1}' ${cut}_stations.csv | \
	gmt pstext -F+f6p,Helvetica,black+jTL+a90 -N
	
	#plot model label
	label_num=2
	echo "-4.7 3.7 ${label_num}) Partial uncertainty (Res. spread)" | gmt pstext -F+f10p,Helvetica-Bold,black+jTL -N
	#plot title
    echo "0 2.6 ${cut}" | gmt pstext -F+f16p,Helvetica,black -N
    #plot NW
	echo "-4.7 2.8 NW" | gmt pstext -F+f12p,Helvetica,black+jTL -N
    #plot SE
	echo "4.7 2.8 SE" | gmt pstext -F+f12p,Helvetica,black+jTR -N
	#plot RMS
    #echo "-4.7 2.2 RMS=${rms}" | gmt pstext -F+f9p,Helvetica,black+jTL -N

gmt subplot set 2
# Horizontal colorbar below all plots, centered at x=0, y=-1.5c
gmt psscale -CRes.cpt -Dx0c/4.7c+w3c/0.4c+h -Bxa1f1+l"Resistivity [log(@~W@~m)]" -N --FONT_LABEL=16p,Helvetica,black --FONT_ANNOT_PRIMARY=16p,Helvetica,black

# Horizontal colorbar below all plots, centered at x=0, y=-1.5c
gmt psscale -Cspread.cpt -Dx4c/4.7c+w3c/0.4c+h -Bxa1f1+l"Resistivity spread [log(@~W@~m)]" -N --FONT_LABEL=16p,Helvetica,black --FONT_ANNOT_PRIMARY=16p,Helvetica,black

gmt subplot set 3
# Add magnitude colorbar
gmt psscale -Cdays.cpt -Dx0c/4.7c+w7c/0.4c+h -Bxa1f+l"Days after mainshock" -N --FONT_LABEL=10p,Helvetica,black --FONT_ANNOT_PRIMARY=10p,Helvetica,black


gmt subplot end

gmt end
echo "All plots completed."

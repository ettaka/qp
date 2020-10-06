#!/bin/env zsh 

function logit {
	echo $1
	echo $1 >> plot_tf_and_qp.log
}

#FILES=("QXF_2D_20200723-71428_s6_m250")
#QPFILES=("meas/Q22_MQXFS6b.txt")
#PKNPK="TRANSFER2_MQXFS4_PK_NPK.txt"
#PKNPK2="TRANSFER1_PK_NPK_simple.txt"
#TF_FILE="TF2_MQXFS6b.csv"

FILES=$1
QPFILES=$2
PKNPK=$3
PKNPK2=$4
TF_FILE=$5
QP_COIL_SEARCH_TERM=$6
QP_COIL_FILTER_OUT=$7
YVAR_ORDER=$8
LEGEND=$9
MEAS_LEGEND=$10
MAGNETNAME=$11

if [ "$#" -ne 11 ]; then
    echo "Illegal number of parameters"
    echo "These are the required arguments:"
    echo "---------------------------------"
    echo '$1 = 2d ansys folders'
    echo '$2 = qp training files'
    echo '$3 = PK and NPK file for TF1'
    echo '$4 = PK and NPK file for TF2'
    echo '$5 = TF measurement file'
    echo '$6 = Additional text for finding coil data'
    echo '$7 = Additional text for filtering coil data'
    echo '$8 = Order the y variable list'
    echo '$9 = legend string'
    echo '$10 = measurement legend string'
    echo '$11 = magnet name'
    exit
fi

if [ "$LEGEND" = "nolegend" ]; then
	NO_XLABEL=("--no-xlabel")
	#NO_XLABEL=("--no-xlabel" "--no-xticklabels")
else
	NO_XLABEL=""
fi

FIGWIDTH="2.7"
FIGHEIGHT="2.7"

echo "Processing TF1" > plot_tf_and_qp.log
transfer_function.py -pk $PKNPK $TF_FILE --ansys-2d-files $FILES  --operations Idling Keys Initial "Insert" --operations-ignore Warm 293 1.9 nominal --legend-location "bottom outside" -nerr --set-ylim "-120 0" --set-xlim "0 80" -uf 1 --image-name "TF1" --fig-width $FIGWIDTH --fig-height $FIGHEIGHT --legend-location $LEGEND --line-width 1 --marker-size 3 $NO_XLABEL>> plot_tf_and_qp.log

logit "Processing Shell Azimuthal Stress"
qp $QPFILES -y "Shell Azimuthal Stress" --normalize-x 320.0521 -s --ansys-2d-files $FILES --gradient-filter 100 --save-fig temp_shell.png --min-points --y-sort-labels $YVAR_ORDER --fig-width $FIGWIDTH --fig-height $FIGHEIGHT -rc --legend-location $LEGEND $MEAS_LEGEND --line-width 1 --marker-size 3 $NO_XLABEL >> plot_tf_and_qp.log
logit "Processing Coil Azimuthal Stress"
qp $QPFILES -y "Coil Azimuthal Stress $QP_COIL_SEARCH_TERM" --filter-out-any "$QP_COIL_FILTER_OUT" --normalize-x 320.0521 -s --ansys-2d-files $FILES --gradient-filter 100 --save-fig temp_pole.png -d --last-points --max-points --y-sort-labels $YVAR_ORDER --fig-width $FIGWIDTH --fig-height $FIGHEIGHT -rc --legend-location $LEGEND $MEAS_LEGEND --line-width 1 --marker-size 3 --set-ylim "0 180" --set-xlim "0 1.4" --set-xticks " 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4" --set-yticks "0 30 60 90 120 150 180" --fit $NO_XLABEL >> plot_tf_and_qp.log

logit "Processing Rod Strain"
qp $QPFILES -y "Rod" --normalize-x 320.0521 --gradient-filter 100 --save-fig temp_rod_strain.png --last-points --y-sort-labels "A B C D" --fig-width $FIGWIDTH --fig-height $FIGHEIGHT -rc --legend-location $LEGEND $MEAS_LEGEND --line-width 1 --marker-size 3 $NO_XLABEL >> plot_tf_and_qp.log
logit "Processing Rod Stress"
qp $QPFILES -y "Rod Stress" --normalize-x 320.0521 -s --gradient-filter 100 --save-fig temp_rod_stress.png --last-points --y-sort-labels "A B C D" --fig-width $FIGWIDTH --fig-height $FIGHEIGHT -rc --legend-location $LEGEND $MEAS_LEGEND --line-width 1 --marker-size 3 $NO_XLABEL >> plot_tf_and_qp.log
logit "Processing Rod Force"
qp $QPFILES -y "Rod Force" --normalize-x 320.0521 -s --rod-force --gradient-filter 100 --save-fig temp_rod_stress.png --last-points --y-sort-labels "A B C D" --fig-width $FIGWIDTH --fig-height $FIGHEIGHT -rc --legend-location $LEGEND $MEAS_LEGEND --line-width 1 --marker-size 3 $NO_XLABEL >> plot_tf_and_qp.log

logit "Processing TF2"
transfer_function.py -pk $PKNPK2 $TF_FILE --ansys-2d-files $FILES  --operations Idling Keys Initial approx "Insert" --operations-ignore Warm 293 nominal SM18 --legend-location "bottom outside" -nerr --set-ylim "-200 0" --set-xlim "0 200" --set-yticks "0 -40 -80 -120 -160 -200" -uf 1 --image-name "TF2" --TF2 --fig-width $FIGWIDTH --fig-height $FIGHEIGHT --legend-location $LEGEND --line-width 1 --marker-size 3 --plot-note $MAGNETNAME $NO_XLABEL>> plot_tf_and_qp.log

logit "Processing KP Shell"
transfer_function.py $TF_FILE  --ansys-2d-files $FILES  --operations Idling Keys Initial "Insert" --operations-ignore Warm 293 1.9 nominal --legend-location "bottom outside" -nerr --key-shell  --set-ylim "0 80" --set-xticks "13.2 13.4 13.6 13.8 14" --set-xlim "13.15 14.15" --set-ylim "0 120" -uf 1 --image-name "KP_shell" --fig-width $FIGWIDTH --fig-height $FIGHEIGHT --legend-location $LEGEND --line-width 1 --marker-size 3 $NO_XLABEL>> plot_tf_and_qp.log

logit "Processing KP Pole"
transfer_function.py $TF_FILE  --ansys-2d-files $FILES  --operations Idling Keys Initial "Insert" --operations-ignore Warm 293 1.9 nominal --legend-location "bottom outside" -nerr --key-pole --set-ylim "-150 0" --set-xticks "13.0 13.2 13.4 13.6 13.8 14" --set-xlim "13.15 14.15" --set-ylim "-120 0" -uf 1 --image-name "KP_pole" --fig-width $FIGWIDTH --fig-height $FIGHEIGHT --legend-location $LEGEND --line-width 1 --marker-size 3 $NO_XLABEL>> plot_tf_and_qp.log

awk '/max points/,/done/' plot_tf_and_qp.log > plot_tf_and_qp_data.txt
awk '/min points/,/done/' plot_tf_and_qp.log >> plot_tf_and_qp_data.txt
awk '/last points/,/done/' plot_tf_and_qp.log >> plot_tf_and_qp_data.txt
awk '/qp fit:/,/done/' plot_tf_and_qp.log >> plot_tf_and_qp_data.txt
awk '/ANSYS2D output:/' plot_tf_and_qp.log > plot_tf_and_qp_data_2d_fem.txt

logit "Processing images"
montage -tile 3x1 -mode concatenate TF2.png \
			temp_shell.png         \
			temp_pole.png    \
			shell_coil.png
#eog shell_coil.png


montage -tile 3x1 -mode concatenate TF1.png         \
			KP_shell.png    \
			KP_pole.png	\
			TF_KP_plots.png

#eog TF_KP_plots.png

montage -tile 5x1 -mode concatenate TF2.png         \
			KP_shell.png    \
			KP_pole.png	\
			temp_shell.png \
			temp_pole.png \
			TF_KP_qp_plots.png

#eog TF_KP_qp_plots_all.png

montage -tile 4x1 -mode concatenate TF2.png         \
			KP_shell.png    \
			KP_pole.png	\
			temp_pole.png \
			TF_KP_qp_plots.png

eog TF_KP_qp_plots.png

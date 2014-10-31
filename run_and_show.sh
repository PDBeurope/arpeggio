#!/bin/bash

#############################################
# NOTE: DOESN'T WORK WITH OUTPUT POSTFIXES. #
#       JUST DESIGNED FOR QUICK TESTING.    #
#                                           #
#       PASS NORMAL OPTIONS TO ARPEGGIO.    #
#                                           #
#       SHOW CONTACTS USES -bs OPTION.      #
#                                           #
#############################################

# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Running Arpeggio."
# http://stackoverflow.com/questions/3190818/pass-all-arguments-from-bash-script-to-another-command
python $SCRIPT_DIR/arpeggio.py "$@"

echo "Running PyMOL."
pymol -R &
PYMOL_PID=$!

echo "Waiting..."
echo "(Fudge factor, PyMOL has to load, increase if script fails here.)"
sleep 1

echo "Running show contacts."
python $SCRIPT_DIR/show_contacts.py $1 -xml -bs

echo "Finished running things."
echo "PyMOL PID = $PYMOL_PID"

while :
do

	echo "########################################"
	echo "#                                      #"
	echo "# Type 'exit' here to quit everything. #"
	echo "#                                      #"
	echo "# N.B.: Exiting without closing PyMOL  #"
	echo "#       (leaving more than one PyMOL   #"
	echo "#       server running) may get        #"
	echo "#       confusing!                     #"
	echo "#                                      #"
	echo "########################################"
	read OPTION

	if [ "$OPTION" == "exit" ]; then

		echo "Quitting PyMOL (PID:$PYMOL_PID)."
		pkill -TERM -P $PYMOL_PID
		kill $PYMOL_PID
		exit

	fi

done
